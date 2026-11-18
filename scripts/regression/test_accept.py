"""
Tests for `run_chain.py accept` (reference runs + CHANGELOG.md) on fake runs.
Run: python3 -m pytest scripts/regression/ -v   (host python with uproot, toml)
"""
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import toml
import uproot

sys.path.insert(0, str(Path(__file__).parent))
import run_chain  # noqa: E402


def fake_run(root, label, x, tier='full'):
    d = root / label / tier
    (d / 'job' / 'data').mkdir(parents=True)
    with uproot.recreate(d / 'job' / 'output.root') as f:
        f['t'] = {'x': np.asarray(x, dtype=float)}
    (d / 'job' / 'data' / 'input.root').symlink_to('/nonexistent/input.root')
    with open(d / 'run_info.toml', 'w') as f:
        toml.dump({'tier': tier, 'label': label, 'sha': 'a' * 40, 'ref': 'WORKTREE',
                   'dirty_files': ['x.cc'], 'configs_from': '/repo'}, f)
    return d


@pytest.fixture
def env(tmp_path, monkeypatch):
    log = tmp_path / 'CHANGELOG.md'
    monkeypatch.setattr(run_chain, 'CHANGELOG', log)
    return tmp_path, log


def accept(root, run_dir, note='ok', dry_run=False):
    run_chain.cmd_accept(SimpleNamespace(root=root, run_dir=run_dir, note=note, dry_run=dry_run))


def test_first_then_changed(env):
    root, log = env
    r1 = fake_run(root, 'v1', [1, 2, 3])
    accept(root, r1, note='initial (owner)')
    cur = root / 'reference' / 'full' / 'current'
    assert cur.is_symlink() and (cur / 'job' / 'output.root').exists()
    assert (cur / 'job' / 'data' / 'input.root').is_symlink()          # inputs stay symlinks
    assert toml.load(cur / 'accepted.toml')['note'] == 'initial (owner)'
    text = log.read_text()
    assert text.startswith('# Regression sandbox reference log')
    assert 'full | v1' in text and 'first reference' in text and 'uncommitted changes to 1 files' in text

    r2 = fake_run(root, 'v2', [1, 2, 4])
    accept(root, r2, note='x changed on purpose (owner)')
    assert cur.resolve().name.endswith('_v2')
    assert len(list((root / 'reference' / 'full').glob('*_v*'))) == 2  # old one kept
    text = log.read_text()
    assert 'x changed on purpose' in text and 'RUNS DIFFER' in text and 't/x: 1 of 3 values differ' in text


def test_condense():
    lines = ['  DIFF events/a/all: branch w only in B', '  DIFF events/b/all: branch w only in B',
             '  DIFF h: 1 bin values differ (first bin 3: 1.0 vs 2.0)']
    assert run_chain.condense(lines) == ['2 x branch w only in B', 'DIFF h: 1 bin values differ (first bin 3: 1.0 vs 2.0)']


def test_reference_is_a_copy(env):
    root, _ = env
    r1 = fake_run(root, 'v1', [1, 2, 3])
    accept(root, r1)
    (r1 / 'job' / 'output.root').unlink()                               # e.g. a later --force rerun
    assert (root / 'reference' / 'full' / 'current' / 'job' / 'output.root').exists()


def test_dry_run_changes_nothing(env, capsys):
    root, log = env
    accept(root, fake_run(root, 'v1', [1, 2, 3]))
    before = log.read_text()
    accept(root, fake_run(root, 'v2', [9, 2, 3]), dry_run=True)
    assert log.read_text() == before
    assert (root / 'reference' / 'full' / 'current').resolve().name.endswith('_v1')
    assert 'RUNS DIFFER' in capsys.readouterr().out


def test_refuses_unfinished_or_reference(env):
    root, _ = env
    r1 = fake_run(root, 'v1', [1])
    (r1 / 'run_info.toml').unlink()
    with pytest.raises(SystemExit):
        accept(root, r1)
    r2 = fake_run(root, 'v2', [1])
    accept(root, r2)
    with pytest.raises(SystemExit):
        accept(root, root / 'reference' / 'full' / 'current')
