#!/usr/bin/env python3
"""
Add one TTree from a small ROOT file to a COPY of a big merged file.

Made for the v31 nc_sideband patch: phase 2 was run locally on the merged
nc_sideband tree, giving events/NuMIFull/nc_sideband_variationTree in a small
file; to_gundam needs that tree next to nc_sideband in the stage-5 file.

    python3 scripts/graft_tree.py \
        --base   /exp/icarus/data/users/micarrig/merged_final_v31d_flux.root \
        --extra  /exp/icarus/data/users/micarrig/v31_sideband/sideband_varsys.root \
        --tree   events/NuMIFull/nc_sideband_variationTree \
        --ref    events/NuMIFull/nc_sideband \
        --output /exp/icarus/data/users/micarrig/merged_final_v31e_flux.root [--dry-run]

What it does:
  1. Checks (read-only) that --tree is in --extra and NOT in --base, that --ref
     is in --base, that both trees have the same number of entries, and that
     Run/Subrun/Evt match entry by entry (to_gundam pairs them by entry index).
  2. Stops here with --dry-run.
  3. Copies --base to <output>.partial, adds the tree in UPDATE mode, re-opens
     it and re-checks, then renames <output>.partial to --output.

It never writes to --base or --extra, refuses an existing --output or
<output>.partial, and never deletes anything (a failed run leaves
<output>.partial for the owner to look at). Needs PyROOT (SL7 container:
scripts/in_sl7.sh).
"""
import argparse
import os
import shutil
import sys
import time

import ROOT

ROOT.gROOT.SetBatch(True)
IDS = ['Run', 'Subrun', 'Evt']


def fail(msg):
    print(f'[ERROR] {msg}', file=sys.stderr)
    sys.exit(1)


def open_read(path):
    f = ROOT.TFile.Open(path, 'READ')
    if not f or f.IsZombie():
        fail(f'cannot open {path}')
    if f.TestBit(ROOT.TFile.kRecovered):
        fail(f'{path} was not closed properly (recovered)')
    return f


def get_tree(f, path, label):
    t = f.Get(path)
    if not t or not t.InheritsFrom('TTree'):
        return None
    return t


def ids(f, path):
    df = ROOT.RDataFrame(path, f)
    cols = df.AsNumpy(IDS)
    return [cols[c] for c in IDS]


def check_alignment(extra_file, tree, ref_file, ref):
    """Entry count and Run/Subrun/Evt equality, entry by entry."""
    a = ids(extra_file, tree)
    b = ids(ref_file, ref)
    if len(a[0]) != len(b[0]):
        fail(f'entry counts differ: {tree} has {len(a[0])}, {ref} has {len(b[0])}')
    for name, x, y in zip(IDS, a, b):
        bad = (x != y).nonzero()[0]
        if len(bad):
            fail(f'{name} differs at {len(bad)} entries (first at entry {bad[0]}: {x[bad[0]]} vs {y[bad[0]]})')
    print(f'[INFO] {tree} and {ref}: {len(a[0])} entries, Run/Subrun/Evt identical entry by entry.')
    return len(a[0])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--base', required=True, help='big merged file (read only)')
    p.add_argument('--extra', required=True, help='small file holding --tree (read only)')
    p.add_argument('--tree', required=True, help='tree path, the same in --extra and in the output')
    p.add_argument('--ref', required=True, help='tree in --base that --tree must line up with')
    p.add_argument('--output', required=True, help='new file; must not exist')
    p.add_argument('--dry-run', action='store_true', help='run the checks only, copy and write nothing')
    a = p.parse_args()

    base, extra, out = (os.path.abspath(x) for x in (a.base, a.extra, a.output))
    partial = out + '.partial'
    for path in (base, extra):
        if not os.path.isfile(path):
            fail(f'not found: {path}')
    if out in (base, extra):
        fail('--output must differ from --base and --extra')
    for path in (out, partial):
        if os.path.exists(path):
            fail(f'refusing to overwrite existing {path}')
    out_dir = os.path.dirname(out)
    need = os.path.getsize(base) + os.path.getsize(extra)
    free = shutil.disk_usage(out_dir).free
    print(f'[INFO] base {os.path.getsize(base)/1e9:.2f} GB; free space in {out_dir}: {free/1e9:.1f} GB '
          f'(filesystem; a user quota may be lower)')
    if free < need:
        fail(f'not enough free space for the copy ({need/1e9:.1f} GB needed)')

    fb, fe = open_read(base), open_read(extra)
    if get_tree(fe, a.tree, '--extra') is None:
        fail(f'{a.tree} not found in {extra}')
    if get_tree(fb, a.tree, '--base') is not None:
        fail(f'{a.tree} already exists in {base}; nothing to graft')
    if get_tree(fb, a.ref, '--base') is None:
        fail(f'{a.ref} not found in {base}')
    dest_dir_path = os.path.dirname(a.tree)
    if dest_dir_path and not fb.GetDirectory(dest_dir_path):
        fail(f'directory {dest_dir_path} not found in {base}')
    n = check_alignment(fe, a.tree, fb, a.ref)
    fb.Close()

    if a.dry_run:
        print(f'[DRY RUN] would copy {base} -> {partial}, add {a.tree} ({n} entries), re-check, rename to {out}.')
        return 0

    t0 = time.time()
    print(f'[INFO] copying {base} -> {partial} ...', flush=True)
    shutil.copyfile(base, partial)
    print(f'[INFO] copy done in {time.time() - t0:.0f} s', flush=True)

    fo = ROOT.TFile.Open(partial, 'UPDATE')
    if not fo or fo.IsZombie():
        fail(f'cannot open {partial} for update (left in place)')
    d = fo.GetDirectory(dest_dir_path) if dest_dir_path else fo
    d.cd()
    src = fe.Get(a.tree)
    clone = src.CloneTree(-1, 'fast')
    clone.SetDirectory(d)
    clone.Write(os.path.basename(a.tree))
    fo.Close()
    fe.Close()

    # Re-open and re-check before giving the file its final name.
    fo = open_read(partial)
    if get_tree(fo, a.tree, 'output') is None:
        fail(f'{a.tree} missing after write; {partial} left in place')
    check_alignment(fo, a.tree, fo, a.ref)
    fo.Close()
    os.rename(partial, out)
    print(f'[INFO] wrote {out} ({os.path.getsize(out)/1e9:.2f} GB) in {time.time() - t0:.0f} s')
    return 0


if __name__ == '__main__':
    sys.exit(main())
