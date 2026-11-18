#!/usr/bin/env python3
"""
Diff two regression-sandbox runs (scripts/regression/run_chain.py), e.g.

    python3 scripts/regression/compare.py \
        /exp/icarus/data/users/micarrig/regression/<sha_a>/quick \
        /exp/icarus/data/users/micarrig/regression/<sha_b>/quick

Every *.root file under A (except the inputs) is matched by relative path in
B. For each file it reports objects added/removed; for each TTree the entry
count, branches added/removed, and per branch whether the values are equal
(exact by default, --rtol/--atol to allow rounding; NaN == NaN), with the first
differing entries; for histograms the bin contents and errors. Exits 1 if
anything differs, 0 if the runs are identical.

Host python (/usr/bin/python3, uproot + numpy + awkward) is enough.
"""
import argparse
import re
import sys
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

SKIP_DIRS = {'data', 'logs'}  # symlinked inputs and logs are not outputs


def outputs(run_dir: Path):
    return sorted(p.relative_to(run_dir) for p in run_dir.rglob('*.root')
                  if not (set(p.relative_to(run_dir).parts) & SKIP_DIRS))


def objects(f):
    """{path: classname} for every object, without cycle numbers (last cycle wins)."""
    return {k.rsplit(';', 1)[0]: c for k, c in f.classnames(recursive=True).items()}


def flat(arr):
    a = ak.flatten(arr, axis=None) if isinstance(arr, ak.Array) else np.asarray(arr).ravel()
    return np.asarray(a, dtype=np.float64)


def _numbers(obj):
    """All numbers in a ROOT object read by uproot: TGraph-like objects via
    values() (x then y), containers (TClonesArray, lists) element by element."""
    if hasattr(obj, 'values') and callable(obj.values) and not isinstance(obj, dict):
        vals = obj.values()
        parts = vals if isinstance(vals, tuple) else (vals,)
        return [float(v) for p in parts for v in np.asarray(p, dtype=float).ravel()]
    if isinstance(obj, (int, float, np.number)):
        return [float(obj)]
    return [v for item in obj for v in _numbers(item)]


def object_numbers(branch, n):
    """Per entry, the flattened numbers of an object branch, as an awkward array."""
    return ak.Array([_numbers(o) for o in branch.array(library='np', entry_stop=n)])


def equal(x, y, rtol, atol):
    if x.shape != y.shape:
        return np.zeros(max(len(x), len(y)), dtype=bool)
    with np.errstate(invalid='ignore'):
        return np.isclose(x, y, rtol=rtol, atol=atol, equal_nan=True) if (rtol or atol) else \
            ((x == y) | (np.isnan(x) & np.isnan(y)))


class Report:
    def __init__(self, show):
        self.show, self.ndiff, self.lines = show, 0, []

    def diff(self, msg):
        self.ndiff += 1
        self.lines.append('  DIFF ' + msg)


def compare_tree(rep, path, ta, tb, rtol, atol, ignore):
    if ta.num_entries != tb.num_entries:
        rep.diff(f'{path}: entries {ta.num_entries} vs {tb.num_entries}')
    ka, kb = set(ta.keys()), set(tb.keys())
    for k in sorted(ka - kb):
        rep.diff(f'{path}: branch {k} only in A')
    for k in sorted(kb - ka):
        rep.diff(f'{path}: branch {k} only in B')
    n = min(ta.num_entries, tb.num_entries)
    for k in sorted(ka & kb):
        if ignore and ignore.search(f'{path}/{k}'):
            continue
        try:
            xa, xb = ta[k].array(entry_stop=n), tb[k].array(entry_stop=n)
        except Exception:
            # Object branches (to_gundam's TClonesArray of TGraph splines):
            # compare the numbers inside each entry's objects instead.
            try:
                xa, xb = object_numbers(ta[k], n), object_numbers(tb[k], n)
            except Exception as e:  # still unreadable: report, keep going
                rep.diff(f'{path}/{k}: cannot read ({type(e).__name__}: {str(e).splitlines()[0]})')
                continue
        if isinstance(xa, ak.Array) and xa.ndim > 1:
            ca, cb = ak.num(xa, axis=1), ak.num(xb, axis=1)
            if not np.array_equal(np.asarray(ca), np.asarray(cb)):
                bad = np.nonzero(np.asarray(ca) != np.asarray(cb))[0]
                rep.diff(f'{path}/{k}: element counts differ in {len(bad)} entries (first entry {bad[0]})')
                continue
        fa, fb = flat(xa), flat(xb)
        ok = equal(fa, fb, rtol, atol)
        if not ok.all():
            bad = np.nonzero(~ok)[0]
            with np.errstate(invalid='ignore'):
                maxd = np.nanmax(np.abs(fa[bad] - fb[bad])) if fa.shape == fb.shape else float('nan')
            ex = ', '.join(f'[{i}] {fa[i]!r} vs {fb[i]!r}' for i in bad[:rep.show]) if fa.shape == fb.shape else ''
            rep.diff(f'{path}/{k}: {len(bad)} of {len(fa)} values differ, max |diff| {maxd:.3g}; {ex}')


def compare_hist(rep, path, ha, hb, rtol, atol):
    for what, fa, fb in (('values', ha.values(flow=True), hb.values(flow=True)),
                         ('errors', ha.errors(flow=True), hb.errors(flow=True))):
        fa, fb = np.asarray(fa, float).ravel(), np.asarray(fb, float).ravel()
        if fa.shape != fb.shape:
            rep.diff(f'{path}: binning differs ({fa.shape} vs {fb.shape})')
            return
        ok = equal(fa, fb, rtol, atol)
        if not ok.all():
            bad = np.nonzero(~ok)[0]
            rep.diff(f'{path}: {len(bad)} bin {what} differ (first bin {bad[0]}: {fa[bad[0]]!r} vs {fb[bad[0]]!r})')


def compare_file(rep, fa_path, fb_path, rtol, atol, ignore):
    with uproot.open(fa_path) as fa, uproot.open(fb_path) as fb:
        oa, ob = objects(fa), objects(fb)
        for k in sorted(set(oa) - set(ob)):
            rep.diff(f'{k} ({oa[k]}) only in A')
        for k in sorted(set(ob) - set(oa)):
            rep.diff(f'{k} ({ob[k]}) only in B')
        for k in sorted(set(oa) & set(ob)):
            ca, cb = oa[k], ob[k]
            if ca != cb:
                rep.diff(f'{k}: class {ca} vs {cb}')
            elif ca.startswith('TTree'):
                compare_tree(rep, k, fa[k], fb[k], rtol, atol, ignore)
            elif ca.startswith(('TH1', 'TH2', 'TH3', 'TProfile')):
                if not (ignore and ignore.search(k)):
                    compare_hist(rep, k, fa[k], fb[k], rtol, atol)
            # TDirectory, TGraph, splines, ...: presence is checked above


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('run_a', type=Path)
    p.add_argument('run_b', type=Path)
    p.add_argument('--rtol', type=float, default=0.0, help='relative tolerance (default 0: exact)')
    p.add_argument('--atol', type=float, default=0.0, help='absolute tolerance (default 0: exact)')
    p.add_argument('--ignore', help='regex of object or tree/branch paths to skip')
    p.add_argument('--show', type=int, default=3, help='differing entries shown per branch')
    a = p.parse_args(argv)
    ignore = re.compile(a.ignore) if a.ignore else None
    for d in (a.run_a, a.run_b):
        if not d.is_dir():
            print(f'not a directory: {d}')
            return 2
    fa, fb = outputs(a.run_a), outputs(a.run_b)
    total = 0
    for rel in sorted(set(fa) | set(fb)):
        rep = Report(a.show)
        if rel not in fb:
            rep.diff('file only in A')
        elif rel not in fa:
            rep.diff('file only in B')
        else:
            compare_file(rep, a.run_a / rel, a.run_b / rel, a.rtol, a.atol, ignore)
        print(f'{"OK  " if rep.ndiff == 0 else "DIFF"} {rel}' + (f' ({rep.ndiff})' if rep.ndiff else ''))
        for line in rep.lines:
            print(line)
        total += rep.ndiff
    print(f'\n{len(set(fa) | set(fb))} files, {total} differences: '
          f'{"IDENTICAL" if total == 0 else "RUNS DIFFER"}')
    return 0 if total == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
