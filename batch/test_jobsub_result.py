"""
jobsub_submit can exit 0 without submitting anything (2026-11-09:
"condor_vault_storer failed: exit code 256 ... BAILING OUT", then medulla.py
printed "Submitted 2471 Phase 2 jobs"). utilities._run_jobsub_submit treats a
submission as real only if the output has the "Use job id ..." line and no
"BAILING OUT"; every submit path must report failure otherwise.
"""
import subprocess
from types import SimpleNamespace

import pytest

from test_dry_run import (  # noqa: F401
    GOOD, JOBSUB_OK, launch_project, phase2_project, write_variation_toml, work, reset_flags,
)
import utilities

BAIL = ("Transferring files to web sandbox...\nSubmitting job(s)\n"
        "Failed to process job credential requests (1): 'process_job_credentials(): "
        "'/opt/jobsub_lite/bin/condor_vault_storer' failed: exit code 256\n'; BAILING OUT.\n")


def fake_jobsub(monkeypatch, stdout='', stderr='', returncode=0):
    """Replace subprocess.run for jobsub_submit only; run cp/mv/rm for real."""
    calls = []
    real_run = subprocess.run

    def fake_run(cmd, *args, **kwargs):
        if isinstance(cmd, (list, tuple)) and cmd and cmd[0] == 'jobsub_submit':
            calls.append(list(cmd))
            if returncode != 0 and kwargs.get('check'):
                raise subprocess.CalledProcessError(returncode, cmd, output=stdout, stderr=stderr)
            return SimpleNamespace(stdout=stdout, stderr=stderr, returncode=returncode)
        return real_run(cmd, *args, **kwargs)

    monkeypatch.setattr(utilities.subprocess, 'run', fake_run)
    return calls


@pytest.mark.parametrize('stdout,stderr,rc,ok', [
    (JOBSUB_OK, '', 0, True),
    (BAIL, '', 0, False),                         # the 2026-11-09 case
    ('', BAIL, 0, False),                         # same message on stderr
    ('Submitting job(s)\n', '', 0, False),        # exit 0, no job id
    ('', 'ExpiredSignatureError', 1, False),
    ('', 'something else', 2, False),
])
def test_run_jobsub_submit(monkeypatch, capsys, stdout, stderr, rc, ok):
    fake_jobsub(monkeypatch, stdout, stderr, rc)
    result = utilities._run_jobsub_submit(['jobsub_submit', '-G', 'icarus'])
    assert (result is not None) == ok
    out = capsys.readouterr().out
    if not ok:
        assert '[ERROR]' in out or 'ERROR' in out
    if 'condor_vault_storer' in stdout + stderr:
        assert 'kinit' in out


def test_launch_reports_bail_out(tmp_path, work, monkeypatch, capsys):
    project = launch_project(tmp_path)
    calls = fake_jobsub(monkeypatch, stdout=BAIL)
    utilities.set_assume_yes(True)
    ok = utilities.launch_jobsub(project, exp='icarus', njobs=-1, tag='abc123')
    out = capsys.readouterr().out
    assert len(calls) == 1
    assert ok is False
    assert 'Launched' not in out and 'did NOT submit' in out


def test_phase2_reports_bail_out(tmp_path, work, monkeypatch, capsys):
    project = phase2_project(tmp_path)
    nue = tmp_path / 'splines_nue.root'
    nue.write_bytes(GOOD)
    vtoml = write_variation_toml(tmp_path / 'var.toml', splines_nue=nue)
    calls = fake_jobsub(monkeypatch, stdout=BAIL)
    utilities.set_assume_yes(True)
    ok = utilities.launch_variation_phase2_jobsub(project, str(vtoml), exp='icarus', tag='abc')
    out = capsys.readouterr().out
    assert len(calls) == 1
    assert ok is False
    assert 'Submitted' not in out and 'did NOT submit' in out


def test_phase2_reports_success(tmp_path, work, monkeypatch, capsys):
    project = phase2_project(tmp_path)
    nue = tmp_path / 'splines_nue.root'
    nue.write_bytes(GOOD)
    vtoml = write_variation_toml(tmp_path / 'var.toml', splines_nue=nue)
    fake_jobsub(monkeypatch, stdout=JOBSUB_OK)
    utilities.set_assume_yes(True)
    ok = utilities.launch_variation_phase2_jobsub(project, str(vtoml), exp='icarus', tag='abc')
    assert ok is True
    assert 'Submitted 2 Phase 2 jobs' in capsys.readouterr().out


@pytest.mark.parametrize('stdout,flag', [(BAIL, True), (JOBSUB_OK, False)])
def test_submit_failed_flag(monkeypatch, stdout, flag):
    monkeypatch.setattr(utilities, 'SUBMIT_FAILED', False)
    fake_jobsub(monkeypatch, stdout=stdout)
    utilities._run_jobsub_submit(['jobsub_submit'])
    assert utilities.SUBMIT_FAILED is flag


@pytest.mark.parametrize('fake_out,rc', [(BAIL, 1), (JOBSUB_OK, 0)])
def test_medulla_cli_exit_code(tmp_path, work, fake_out, rc):
    """End to end: medulla.py --launch-jobs with a fake jobsub_submit first on PATH."""
    import os
    import sys
    from pathlib import Path
    project = launch_project(tmp_path)
    bindir = tmp_path / 'bin'
    bindir.mkdir()
    (bindir / 'out.txt').write_text(fake_out)
    fake = bindir / 'jobsub_submit'
    fake.write_text(f'#!/bin/bash\ncat {bindir / "out.txt"}\nexit 0\n')
    fake.chmod(0o755)
    env = dict(os.environ, PATH=f'{bindir}:{os.environ["PATH"]}')
    medulla_py = Path(utilities.__file__).parent / 'medulla.py'
    res = subprocess.run([sys.executable, str(medulla_py), '-p', str(project), '-e', 'icarus',
                          '--launch-jobs', '-y'], cwd=work, env=env, capture_output=True,
                         text=True, stdin=subprocess.DEVNULL)
    assert res.returncode == rc, res.stdout + res.stderr
    assert ('did NOT submit' in res.stdout) == (rc == 1)
