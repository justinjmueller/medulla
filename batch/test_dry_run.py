"""
Tests that --dry-run (utilities.DRY_RUN) changes nothing: no jobsub_submit, no
writes to project.db or the project directory, no deleted outputs, and no stray
files in the working directory. Each dry-run test has a control (same setup,
dry run off) so that the test would fail if the mutation it guards against
stopped being detectable.
"""
import hashlib
import sqlite3
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).parent))
import medulla  # noqa: E402
import utilities  # noqa: E402

BATCH_DIR = Path(__file__).parent
GOOD = b'x' * 4096     # >= 1 KB: counts as a real output
STUB = b'x' * 10       # < 1 KB: counts as a stub


# ---------------------------------------------------------------------------
# helpers / fixtures
# ---------------------------------------------------------------------------

def snapshot(directory: Path) -> dict:
    """Relative path -> (size, sha256, mtime_ns) for every file under directory."""
    snap = {}
    for p in sorted(Path(directory).rglob('*')):
        if p.is_file():
            snap[str(p.relative_to(directory))] = (
                p.stat().st_size,
                hashlib.sha256(p.read_bytes()).hexdigest(),
                p.stat().st_mtime_ns,
            )
    return snap


def make_project(root: Path, jobs, outputs=None, name='proj') -> Path:
    """
    Build a minimal project directory by hand.

    jobs    : list of (jobid, sample_tag, status)
    outputs : dict of relative filename under output/ -> bytes
    """
    project = root / name
    (project / 'output').mkdir(parents=True)
    conn = sqlite3.connect(project / 'project.db')
    conn.execute(utilities.SCHEMA_CONFIGURATION)
    conn.execute(utilities.SCHEMA_JOBS)
    for jid, tag, status in jobs:
        cfg = f'[[sample]]\nname = "s{jid}"\ntag = "{tag}"\nismc = true\n'
        conn.execute('INSERT INTO configuration VALUES (?, ?)', (jid, cfg))
        conn.execute('INSERT INTO jobs VALUES (?, ?)', (jid, status))
    conn.commit()
    conn.close()
    for rel, data in (outputs or {}).items():
        (project / 'output' / rel).write_bytes(data)
    return project


def job_statuses(project: Path) -> dict:
    conn = sqlite3.connect(project / 'project.db')
    rows = dict(conn.execute('SELECT jobid, status FROM jobs').fetchall())
    conn.close()
    return rows


@pytest.fixture(autouse=True)
def reset_flags():
    utilities.set_dry_run(False)
    utilities.set_assume_yes(False)
    yield
    utilities.set_dry_run(False)
    utilities.set_assume_yes(False)


@pytest.fixture
def work(tmp_path, monkeypatch):
    """A separate working directory (the code copies project.db into the cwd)."""
    w = tmp_path / 'work'
    w.mkdir()
    monkeypatch.chdir(w)
    return w


@pytest.fixture
def jobsub_calls(monkeypatch):
    """
    Replace utilities.subprocess.run: record jobsub_submit calls and pretend
    they succeed; run everything else (cp, mv, rm) for real.
    """
    calls = []
    real_run = subprocess.run

    def fake_run(cmd, *args, **kwargs):
        if isinstance(cmd, (list, tuple)) and cmd and cmd[0] == 'jobsub_submit':
            calls.append(list(cmd))
            return SimpleNamespace(stdout='submitted\n', stderr='', returncode=0)
        return real_run(cmd, *args, **kwargs)

    monkeypatch.setattr(utilities.subprocess, 'run', fake_run)
    return calls


@pytest.fixture
def no_os_system(monkeypatch):
    """Fail the test if anything shells out through os.system."""
    def boom(cmd):
        raise AssertionError(f'os.system called: {cmd}')
    monkeypatch.setattr(utilities.os, 'system', boom)


# ---------------------------------------------------------------------------
# check_project_status
# ---------------------------------------------------------------------------

def status_project(tmp_path):
    return make_project(
        tmp_path,
        jobs=[(0, 'nominal', 'pending'), (1, 'nominal', 'submitted'), (2, 'nominal', 'submitted')],
        outputs={
            'output_jobid0000.root': GOOD,
            'output_systematics_jobid0000.root': GOOD,
            'output_jobid0001.root': STUB,   # stub: would be deleted + reset
        },                                    # job 2 has no output: "missing"
    )


@pytest.mark.parametrize('assume_yes', [False, True])
def test_check_status_dry_run_changes_nothing(tmp_path, work, assume_yes):
    project = status_project(tmp_path)
    before_proj, before_work = snapshot(project), snapshot(work)
    utilities.set_dry_run(True)
    utilities.set_assume_yes(assume_yes)   # dry run must win over --yes
    utilities.check_project_status(project)
    assert snapshot(project) == before_proj           # db + outputs untouched, stub kept
    assert snapshot(work) == before_work              # no leftover ./project.db
    assert job_statuses(project) == {0: 'pending', 1: 'submitted', 2: 'submitted'}


def test_check_status_control_actually_mutates(tmp_path, work):
    project = status_project(tmp_path)
    utilities.set_assume_yes(True)
    utilities.check_project_status(project)
    assert not (project / 'output' / 'output_jobid0001.root').exists()   # stub deleted
    assert job_statuses(project)[0] == 'completed'
    assert job_statuses(project)[1] == 'pending'                          # reset after delete


# ---------------------------------------------------------------------------
# launch_jobsub
# ---------------------------------------------------------------------------

def launch_project(tmp_path):
    return make_project(
        tmp_path,
        jobs=[(0, 'nominal', 'pending'), (1, 'nominal', 'pending'), (2, 'nominal', 'submitted')],
        outputs={'output_jobid0000.root': GOOD},
    )


@pytest.mark.parametrize('relaunch_missing', [False, True])
def test_launch_dry_run_does_not_submit_or_write(tmp_path, work, jobsub_calls, capsys, relaunch_missing):
    project = launch_project(tmp_path)
    before_proj, before_work = snapshot(project), snapshot(work)
    utilities.set_dry_run(True)
    utilities.set_assume_yes(True)
    ok = utilities.launch_jobsub(project, exp='icarus', njobs=-1, tag='abc123',
                                 memory=10000, disk=10, relaunch_missing=relaunch_missing)
    assert ok is False
    assert jobsub_calls == []
    assert snapshot(project) == before_proj
    assert snapshot(work) == before_work
    assert job_statuses(project) == {0: 'pending', 1: 'pending', 2: 'submitted'}
    out = capsys.readouterr().out
    assert 'jobsub_submit' in out and '--tag=abc123' in out and '--memory=10000MB' in out
    assert '--disk=10GB' in out and 'Nothing submitted' in out


def test_launch_control_submits(tmp_path, work, jobsub_calls):
    project = launch_project(tmp_path)
    utilities.set_assume_yes(True)
    utilities.launch_jobsub(project, exp='icarus', njobs=-1, tag='abc123')
    assert len(jobsub_calls) == 1


# ---------------------------------------------------------------------------
# create_new_project and medulla.main
# ---------------------------------------------------------------------------

def test_create_dry_run_creates_nothing(tmp_path, work):
    tml = tmp_path / 'selection.toml'
    tml.write_text('[general]\n')
    project = tmp_path / 'newproj'
    utilities.set_dry_run(True)
    utilities.create_new_project(project, str(tml), 30, None)
    assert not project.exists()


def test_main_dry_run_sets_flag_itself(tmp_path, work):
    """medulla.main(dry_run=True) must not rely on the CLI to set the flag."""
    tml = tmp_path / 'selection.toml'
    tml.write_text('[general]\n')
    project = tmp_path / 'newproj'
    medulla.main(str(project), 'icarus', create_project=True, tml=str(tml),
                 batch_size=30, dry_run=True)
    assert not project.exists()
    assert utilities.DRY_RUN is True


def test_cli_dry_run_create(tmp_path):
    tml = tmp_path / 'selection.toml'
    tml.write_text('[general]\n')
    project = tmp_path / 'newproj'
    r = subprocess.run(
        [sys.executable, str(BATCH_DIR / 'medulla.py'), '-p', str(project), '-e', 'icarus',
         '-c', '-t', str(tml), '-b', '30', '--dry-run'],
        capture_output=True, text=True, cwd=tmp_path, stdin=subprocess.DEVNULL,
    )
    assert r.returncode == 0, r.stderr
    assert 'DRY RUN' in r.stdout
    assert not project.exists()


# ---------------------------------------------------------------------------
# Phase 1
# ---------------------------------------------------------------------------

def write_variation_toml(path: Path, splines_nue=None):
    text = '[input]\npath = "in.root"\n[output]\npath = "out.root"\n[variations]\nkeys = ["cv"]\n'
    if splines_nue is not None:
        text += f'splines_file_nue = "{splines_nue}"\n'
    text += '[[tree]]\norigin = "events/a/b"\ndestination = "events/a/"\nname = "b"\naction = "add_weights"\n'
    path.write_text(text)
    return path


def test_phase1_jobsub_dry_run_writes_nothing(tmp_path, work, jobsub_calls):
    project = make_project(tmp_path, jobs=[(0, 'detector_variation', 'completed')])
    vtoml = write_variation_toml(tmp_path / 'var.toml')
    before_proj, before_work = snapshot(project), snapshot(work)
    utilities.set_dry_run(True)
    utilities.set_assume_yes(True)
    ok = utilities.launch_variation_phase1_jobsub(project, 'root://x//merged.root', str(vtoml),
                                                  exp='icarus', tag='abc', memory=10000, disk=10)
    assert ok is False
    assert jobsub_calls == []
    assert snapshot(project) == before_proj      # variation_systematics.toml NOT copied in
    assert snapshot(work) == before_work


def test_phase1_jobsub_control_copies_and_submits(tmp_path, work, jobsub_calls):
    project = make_project(tmp_path, jobs=[(0, 'detector_variation', 'completed')])
    vtoml = write_variation_toml(tmp_path / 'var.toml')
    utilities.set_assume_yes(True)
    utilities.launch_variation_phase1_jobsub(project, 'root://x//merged.root', str(vtoml),
                                             exp='icarus', tag='abc')
    assert (project / 'variation_systematics.toml').exists()
    assert len(jobsub_calls) == 1


def test_phase1_interactive_dry_run_runs_nothing(tmp_path, work, monkeypatch, no_os_system):
    project = make_project(tmp_path, jobs=[(0, 'detector_variation', 'completed')])
    vtoml = write_variation_toml(tmp_path / 'var.toml')
    calls = []
    monkeypatch.setattr(utilities.subprocess, 'run', lambda *a, **k: calls.append(a) or None)
    # make the binary check pass so the dry-run branch itself is exercised
    monkeypatch.setattr(utilities.Path, 'exists', lambda self: True)
    before_proj = snapshot(project)
    utilities.set_dry_run(True)
    utilities.set_assume_yes(True)
    ok = utilities.run_variation_phase1_interactive(project, 'root://x//merged.root', str(vtoml))
    assert ok is False
    assert calls == []
    monkeypatch.undo()
    assert snapshot(project) == before_proj


# ---------------------------------------------------------------------------
# Phase 2
# ---------------------------------------------------------------------------

def phase2_project(tmp_path):
    project = make_project(
        tmp_path,
        jobs=[(0, 'nominal', 'completed'), (1, 'nominal', 'completed'), (2, 'detector_variation', 'completed')],
        outputs={
            'output_systematics_jobid0000.root': GOOD,
            'output_systematics_jobid0001.root': GOOD,
            'output_varsys_jobid0000.root': STUB,    # bad: would be deleted + resubmitted
        },
    )
    (project / 'variation_splines.root').write_bytes(GOOD)
    return project


def test_phase2_dry_run_changes_nothing(tmp_path, work, jobsub_calls, capsys):
    project = phase2_project(tmp_path)
    nue = tmp_path / 'splines_nue.root'
    nue.write_bytes(GOOD)
    vtoml = write_variation_toml(tmp_path / 'var.toml', splines_nue=nue)
    before_proj, before_work = snapshot(project), snapshot(work)
    utilities.set_dry_run(True)
    utilities.set_assume_yes(True)
    ok = utilities.launch_variation_phase2_jobsub(project, str(vtoml), exp='icarus', tag='abc',
                                                  memory=10000, disk=10)
    assert ok is False
    assert jobsub_calls == []
    assert snapshot(project) == before_proj      # no toml/manifest written, stub varsys kept
    assert snapshot(work) == before_work         # no leftover project_p2.db etc.
    out = capsys.readouterr().out
    assert 'Would delete 1 stub varsys' in out
    assert '-N 2' in out                          # both nominal jobs pending after the stub is dropped


def test_phase2_control_writes_deletes_submits(tmp_path, work, jobsub_calls):
    project = phase2_project(tmp_path)
    nue = tmp_path / 'splines_nue.root'
    nue.write_bytes(GOOD)
    vtoml = write_variation_toml(tmp_path / 'var.toml', splines_nue=nue)
    utilities.set_assume_yes(True)
    utilities.launch_variation_phase2_jobsub(project, str(vtoml), exp='icarus', tag='abc')
    assert not (project / 'output' / 'output_varsys_jobid0000.root').exists()
    assert (project / 'variation_phase2_manifest.txt').exists()
    assert len(jobsub_calls) == 1


# ---------------------------------------------------------------------------
# merge_list.py --dry-run
# ---------------------------------------------------------------------------

def merge_list_project(tmp_path):
    return make_project(
        tmp_path,
        jobs=[(0, 'nominal', 'completed'), (1, 'detector_variation', 'completed')],
        outputs={
            'output_systematics_jobid0000.root': GOOD,
            'output_jobid0001.root': GOOD,
        },
    )


def run_merge_list(project, work, out_file, dry_run):
    cmd = [sys.executable, str(BATCH_DIR / 'merge_list.py'), str(project),
           '--output', str(out_file)]
    if dry_run:
        cmd.append('--dry-run')
    return subprocess.run(cmd, capture_output=True, text=True, cwd=work,
                          stdin=subprocess.DEVNULL)


def test_merge_list_dry_run_writes_nothing(tmp_path, work):
    project = merge_list_project(tmp_path)
    out_file = tmp_path / 'merge_list.txt'
    before_proj = snapshot(project)
    r = run_merge_list(project, work, out_file, dry_run=True)
    assert r.returncode == 0, r.stderr
    assert 'DRY-RUN' in r.stdout
    assert not out_file.exists()              # no --output file written
    assert snapshot(work) == {}               # temp project.db copy cleaned up, cwd clean
    assert snapshot(project) == before_proj   # project untouched


def test_merge_list_control_writes_output(tmp_path, work):
    project = merge_list_project(tmp_path)
    out_file = tmp_path / 'merge_list.txt'
    r = run_merge_list(project, work, out_file, dry_run=False)
    assert r.returncode == 0, r.stderr
    assert out_file.exists()                  # the mutation the dry-run guards against
    assert out_file.read_text().strip() != ''
