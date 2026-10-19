"""
Tests for analysis/plot_all.py (the plotting driver) and the per-family CLI
dry-runs. The memory gate and job dispatch are exercised with an injected
available-memory value and a fake runner, so no ROOT and no 30 GB input are
needed. The CLI dry-runs are run as real subprocesses (they defer the heavy
imports, so --dry-run needs no ROOT) and must write nothing.
"""
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import plot_all  # noqa: E402
import plot_cli_common  # noqa: E402

GOOD = b'x' * 4096


# ---------------------------------------------------------------------------
# pure logic: fits / load_jobs / build_cmd
# ---------------------------------------------------------------------------

def test_fits_threshold():
    assert plot_all.fits(30, 30, 1.0) is True
    assert plot_all.fits(29, 30, 1.0) is False
    assert plot_all.fits(45, 30, 1.5) is True
    assert plot_all.fits(44, 30, 1.5) is False
    assert plot_all.fits(None, 30) is False
    assert plot_all.fits(30, None) is False


def write_manifest(tmp_path, inp):
    m = tmp_path / 'jobs.toml'
    m.write_text(
        '[[job]]\n'
        'name = "selection"\ncli = "analysis/makeSelectionPlotsCLI.py"\n'
        'toml = "analysis/nue_plots.toml"\n'
        f'input = "{inp}"\noutput_dir = "out/selection"\n\n'
        '[[job]]\n'
        'name = "eff"\ncli = "analysis/plotEfficiencyCLI.py"\n'
        'toml = "analysis/eff.toml"\n'
        f'input = "{inp}"\noutput_dir = "out/eff"\nfigures = ["efficiency_vs_energy"]\n'
    )
    return m


def test_load_jobs_ok(tmp_path):
    jobs = plot_all.load_jobs(write_manifest(tmp_path, tmp_path / 'in.root'))
    assert [j['name'] for j in jobs] == ['selection', 'eff']


def test_load_jobs_rejects_missing_field(tmp_path):
    m = tmp_path / 'jobs.toml'
    m.write_text('[[job]]\nname = "x"\ncli = "a.py"\n')  # no toml/input/output_dir
    with pytest.raises(ValueError):
        plot_all.load_jobs(m)


def test_load_jobs_rejects_duplicate_names(tmp_path):
    m = tmp_path / 'jobs.toml'
    body = ('[[job]]\nname = "dup"\ncli = "a.py"\ntoml = "t"\ninput = "i"\noutput_dir = "o"\n')
    m.write_text(body + body)
    with pytest.raises(ValueError):
        plot_all.load_jobs(m)


def test_build_cmd_shapes():
    job = {'name': 'j', 'cli': 'analysis/x.py', 'toml': 't.toml',
           'input': 'i.root', 'output_dir': 'o'}
    cmd = plot_all.build_cmd(job, dry_run=False)
    assert cmd[1:] == ['analysis/x.py', '--toml', 't.toml', '--input', 'i.root', '--output-dir', 'o']
    assert '--dry-run' not in cmd
    cmd_fig = plot_all.build_cmd({**job, 'figures': ['a', 'b']}, dry_run=True)
    assert '--figures' in cmd_fig and 'a' in cmd_fig and 'b' in cmd_fig
    assert cmd_fig[-1] == '--dry-run'


def test_build_cmd_args_passthrough():
    job = {'name': 'j', 'cli': 'analysis/x.py', 'toml': 't', 'input': 'i', 'output_dir': 'o',
           'args': ['--systematics-set', 'xsec']}
    cmd = plot_all.build_cmd(job, dry_run=False)
    assert cmd[-2:] == ['--systematics-set', 'xsec']


# ---------------------------------------------------------------------------
# apply_systematics_set (pure; no ROOT)
# ---------------------------------------------------------------------------

def test_apply_systematics_set_overrides_only_spine_systematics():
    cfg = {
        'figure': [
            {'name': 'a', 'artists': [{'type': 'SpineSystematics', 'systematics': ['old']},
                                      {'type': 'SpineHistogram'}]},
            {'name': 'b', 'artists': [{'type': 'SpineSystematics', 'systematics': ['old']}]},
        ]
    }
    n = plot_cli_common.apply_systematics_set(cfg, ['x', 'y', 'z'])
    assert n == 2
    assert cfg['figure'][0]['artists'][0]['systematics'] == ['x', 'y', 'z']
    assert 'systematics' not in cfg['figure'][0]['artists'][1]   # non-systematics artist untouched
    assert cfg['figure'][1]['artists'][0]['systematics'] == ['x', 'y', 'z']


def test_apply_systematics_set_none_present():
    cfg = {'figure': [{'name': 'a', 'artists': [{'type': 'SpineHistogram'}]}]}
    assert plot_cli_common.apply_systematics_set(cfg, ['x']) == 0


# ---------------------------------------------------------------------------
# run_jobs dispatch (fake runner, injected memory)
# ---------------------------------------------------------------------------

def fake_runner_factory():
    calls = []

    def runner(cmd):
        calls.append(cmd)
        return SimpleNamespace(returncode=0)
    return calls, runner


def jobs_with_input(inp):
    return [
        {'name': 'a', 'cli': 'analysis/x.py', 'toml': 't', 'input': str(inp), 'output_dir': 'o'},
    ]


def test_run_jobs_skips_when_memory_low(tmp_path):
    inp = tmp_path / 'in.root'
    inp.write_bytes(GOOD)
    calls, runner = fake_runner_factory()
    res = plot_all.run_jobs(jobs_with_input(inp), dry_run=False,
                            available_bytes=1, runner=runner, out=lambda *_: None)
    assert res == [('a', 'skipped-memory')]
    assert calls == []                      # nothing actually run


def test_run_jobs_runs_when_memory_ok(tmp_path):
    inp = tmp_path / 'in.root'
    inp.write_bytes(GOOD)
    calls, runner = fake_runner_factory()
    res = plot_all.run_jobs(jobs_with_input(inp), dry_run=False,
                            available_bytes=10 ** 9, runner=runner, out=lambda *_: None)
    assert res == [('a', 'ok')]
    assert len(calls) == 1 and '--dry-run' not in calls[0]


def test_run_jobs_dry_run_ignores_gate_and_never_runs_real(tmp_path):
    inp = tmp_path / 'in.root'
    inp.write_bytes(GOOD)
    calls, runner = fake_runner_factory()
    res = plot_all.run_jobs(jobs_with_input(inp), dry_run=True,
                            available_bytes=1, runner=runner, out=lambda *_: None)
    assert res == [('a', 'dry-run')]         # ran despite tiny memory
    assert calls and calls[0][-1] == '--dry-run'


def test_run_jobs_skips_missing_input(tmp_path):
    calls, runner = fake_runner_factory()
    jobs = jobs_with_input(tmp_path / 'nope.root')
    res = plot_all.run_jobs(jobs, dry_run=False, available_bytes=10 ** 12,
                            runner=runner, out=lambda *_: None)
    assert res == [('a', 'missing-input')]
    assert calls == []


def test_run_jobs_only_filter(tmp_path):
    inp = tmp_path / 'in.root'
    inp.write_bytes(GOOD)
    jobs = [
        {'name': 'a', 'cli': 'x.py', 'toml': 't', 'input': str(inp), 'output_dir': 'o'},
        {'name': 'b', 'cli': 'y.py', 'toml': 't', 'input': str(inp), 'output_dir': 'o'},
    ]
    calls, runner = fake_runner_factory()
    res = plot_all.run_jobs(jobs, dry_run=True, available_bytes=1, only=['b'],
                            runner=runner, out=lambda *_: None)
    assert ('a', 'not-selected') in res and ('b', 'dry-run') in res
    assert len(calls) == 1 and 'y.py' in calls[0]


# ---------------------------------------------------------------------------
# CLI dry-runs as real subprocesses (defer heavy imports -> no ROOT needed)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('cli', [
    'makeSelectionPlotsCLI.py',
    'makeSystematicPlotsCLI.py',
    'makeNMinus1PlotsCLI.py',
])
def test_cli_dry_run_writes_nothing(tmp_path, cli):
    tml = tmp_path / 'plots.toml'
    tml.write_text('[samples]\n[variables]\n[[figure]]\nname = "reco_nu"\n')
    out_dir = tmp_path / 'out'
    r = subprocess.run(
        [sys.executable, str(HERE / cli), '--toml', str(tml),
         '--input', str(tmp_path / 'missing.root'), '--output-dir', str(out_dir),
         '--figures', 'reco_nu', '--dry-run'],
        capture_output=True, text=True, cwd=HERE.parent, stdin=subprocess.DEVNULL,
    )
    assert r.returncode == 0, r.stderr
    assert 'DRY RUN' in r.stdout
    assert not out_dir.exists()             # nothing created


@pytest.mark.parametrize('cli', [
    'makeSelectionPlotsCLI.py',
    'makeSystematicPlotsCLI.py',
    'makeNMinus1PlotsCLI.py',
])
def test_cli_help(tmp_path, cli):
    r = subprocess.run([sys.executable, str(HERE / cli), '--help'],
                       capture_output=True, text=True, stdin=subprocess.DEVNULL)
    assert r.returncode == 0
    assert '--figures' in r.stdout and '--dry-run' in r.stdout


def test_systematic_cli_set_dry_run(tmp_path):
    tml = tmp_path / 'sys.toml'
    tml.write_text(
        '[samples]\n'
        '[[figure]]\nname = "sys_x"\n'
        '[[figure.artists]]\ntype = "SpineSystematics"\nsystematics = ["old"]\n'
    )
    sets = tmp_path / 'sets.toml'
    sets.write_text('[sets]\nmyset = ["a", "b", "c"]\n')
    out_dir = tmp_path / 'out'
    r = subprocess.run(
        [sys.executable, str(HERE / 'makeSystematicPlotsCLI.py'),
         '--toml', str(tml), '--input', str(tmp_path / 'missing.root'),
         '--output-dir', str(out_dir), '--figures', 'sys_x',
         '--systematics-set', 'myset', '--sets-file', str(sets), '--dry-run'],
        capture_output=True, text=True, cwd=HERE.parent, stdin=subprocess.DEVNULL,
    )
    assert r.returncode == 0, r.stderr
    assert 'Systematics set: myset (3 names)' in r.stdout
    assert 'applied to 1' in r.stdout
    assert not out_dir.exists()


def test_systematic_cli_only_has_set_flag(tmp_path):
    """--systematics-set is on the systematic CLI but not on selection/n-1."""
    r_sys = subprocess.run([sys.executable, str(HERE / 'makeSystematicPlotsCLI.py'), '--help'],
                           capture_output=True, text=True, stdin=subprocess.DEVNULL)
    r_sel = subprocess.run([sys.executable, str(HERE / 'makeSelectionPlotsCLI.py'), '--help'],
                           capture_output=True, text=True, stdin=subprocess.DEVNULL)
    assert '--systematics-set' in r_sys.stdout
    assert '--systematics-set' not in r_sel.stdout
