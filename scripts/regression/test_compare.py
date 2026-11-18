"""
Tests for scripts/regression/compare.py on small synthetic ROOT files.
Run: python3 -m pytest scripts/regression/ -v   (host python with uproot)
"""
import sys
from pathlib import Path

import numpy as np
import pytest
import uproot

sys.path.insert(0, str(Path(__file__).parent))
import compare  # noqa: E402

BASE = {'Evt': np.arange(5, dtype=np.int32),
        'energy': np.array([1.0, 2.0, np.nan, 4.0, 5.0]),
        'weight': np.array([0.9, 1.0, 1.1, 1.2, 1.3])}


def make_run(root: Path, tree=None, hist_counts=None, extra_file=False):
    d = root / 'job'
    (d / 'data').mkdir(parents=True)
    with uproot.recreate(d / 'output.root') as f:
        f['events/NuMIFull/selected'] = dict(tree if tree is not None else BASE)
        counts = np.array(hist_counts if hist_counts is not None else [1.0, 2.0, 3.0])
        f['hist'] = (counts, np.array([0.0, 1.0, 2.0, 3.0]))
    # an input symlink target under data/ must be ignored
    with uproot.recreate(d / 'data' / 'input.root') as f:
        f['t'] = {'x': np.arange(3) + (10 if extra_file else 0)}
    if extra_file:
        with uproot.recreate(d / 'extra.root') as f:
            f['t'] = {'x': np.arange(3)}
    return root


def run(a, b, *args):
    return compare.main([str(a), str(b)] + list(args))


def test_identical_including_nan(tmp_path, capsys):
    a, b = make_run(tmp_path / 'a'), make_run(tmp_path / 'b')
    assert run(a, b) == 0
    assert 'IDENTICAL' in capsys.readouterr().out


def test_changed_value(tmp_path, capsys):
    t = dict(BASE, weight=BASE['weight'] * np.array([1, 1, 1, 1, 1.001]))
    a, b = make_run(tmp_path / 'a'), make_run(tmp_path / 'b', tree=t)
    assert run(a, b) == 1
    out = capsys.readouterr().out
    assert 'events/NuMIFull/selected/weight: 1 of 5 values differ' in out
    assert '[4]' in out


def test_tolerance(tmp_path):
    t = dict(BASE, weight=BASE['weight'] * (1 + 1e-9))
    a, b = make_run(tmp_path / 'a'), make_run(tmp_path / 'b', tree=t)
    assert run(a, b) == 1
    assert run(a, b, '--rtol', '1e-6') == 0


def test_nan_vs_number(tmp_path, capsys):
    t = dict(BASE, energy=np.array([1.0, 2.0, 3.0, 4.0, 5.0]))
    a, b = make_run(tmp_path / 'a'), make_run(tmp_path / 'b', tree=t)
    assert run(a, b) == 1
    assert 'energy: 1 of 5' in capsys.readouterr().out


def test_missing_branch_and_entries(tmp_path, capsys):
    t = {k: v[:4] for k, v in BASE.items() if k != 'weight'}
    a, b = make_run(tmp_path / 'a'), make_run(tmp_path / 'b', tree=t)
    assert run(a, b) == 1
    out = capsys.readouterr().out
    assert 'branch weight only in A' in out
    assert 'entries 5 vs 4' in out


def test_histogram_and_ignore(tmp_path, capsys):
    a, b = make_run(tmp_path / 'a'), make_run(tmp_path / 'b', hist_counts=[1.0, 2.0, 4.0])
    assert run(a, b) == 1
    assert 'hist: 1 bin values differ' in capsys.readouterr().out
    assert run(a, b, '--ignore', '^hist$') == 0


def test_extra_output_file_and_inputs_ignored(tmp_path, capsys):
    a, b = make_run(tmp_path / 'a'), make_run(tmp_path / 'b', extra_file=True)
    assert run(a, b) == 1
    out = capsys.readouterr().out
    assert 'job/extra.root' in out and 'file only in B' in out
    assert 'input.root' not in out  # data/ is skipped even though it differs


class FakeGraph:
    """Stands in for uproot's TGraph model: values() -> (x, y)."""

    def __init__(self, y):
        self.y = y

    def values(self):
        return np.array([-1.0, 0.0, 1.0]), np.array(self.y)


def test_object_numbers_flattens_graph_containers():
    # to_gundam writes a TClonesArray of TGraph per entry; compare their points
    assert compare._numbers([FakeGraph([0.9, 1.0, 1.1])]) == [-1.0, 0.0, 1.0, 0.9, 1.0, 1.1]
    assert compare._numbers([[FakeGraph([1.0, 1.0, 1.0])], 2.5]) == [-1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.5]
