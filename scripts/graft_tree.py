#!/usr/bin/env python3
"""
Add TTrees from small ROOT files to a COPY of a big merged file.

Made for the v31 nc_sideband patch: phase 2 was run locally on the merged
nc_sideband tree, giving events/NuMIFull/nc_sideband_variationTree in a small
file; to_gundam needs that tree next to nc_sideband in the stage-5 file.

One tree:

    python3 scripts/graft_tree.py \
        --base   /exp/icarus/data/users/micarrig/merged_final_v31d_flux.root \
        --extra  /exp/icarus/data/users/micarrig/v31_sideband/sideband_varsys.root \
        --tree   events/NuMIFull/nc_sideband_variationTree \
        --ref    events/NuMIFull/nc_sideband \
        --output /exp/icarus/data/users/micarrig/merged_final_v31e_flux.root [--dry-run]

Several trees, one copy of --base (repeat --graft EXTRA TREE REF; it can be
combined with --extra/--tree/--ref):

    python3 scripts/graft_tree.py --base BIG.root --output NEW.root --replace \
        --graft varsys.root  events/NuMIFull/selected_variationTree   events/NuMIFull/selected \
        --graft varsys.root  events/NuMIFull/all_signal_variationTree events/NuMIFull/all_signal \
        --graft sideband.root events/NuMIFull/photon_sideband_variationTree events/NuMIFull/photon_sideband

With --replace, every tree must ALREADY be in --base (e.g. v32: phase 2 writes
events/NuMIFull/photon_sideband_variationTree with electron-energy weights
that are meaningless for the sideband) and must have the same branch names as
the one it replaces; every cycle of it is deleted in the copy and the new
version is written in its place. Without --replace, every tree must be new.

What it does:
  1. Checks (read-only), for every tree: that it is in its extra file and
     NOT in --base (with --replace: that it IS in --base, with the same
     branches), that its ref is in --base, that both have the same number of
     entries, and that Run/Subrun/Evt match entry by entry (to_gundam pairs
     them by entry index). The same tree may not be given twice.
  2. Stops here with --dry-run.
  3. Copies --base ONCE to <output>.partial, adds all trees in UPDATE mode,
     re-opens it and re-checks every tree, then renames <output>.partial to
     --output.

It never writes to --base or the extra files, refuses an existing --output or
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


def get_tree(f, path):
    t = f.Get(path)
    if not t or not t.InheritsFrom('TTree'):
        return None
    return t


def branch_names(t):
    return sorted(b.GetName() for b in t.GetListOfBranches())


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
            fail(f'{tree}: {name} differs at {len(bad)} entries (first at entry {bad[0]}: {x[bad[0]]} vs {y[bad[0]]})')
    print(f'[INFO] {tree} and {ref}: {len(a[0])} entries, Run/Subrun/Evt identical entry by entry.')
    return len(a[0])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--base', required=True, help='big merged file (read only)')
    p.add_argument('--extra', help='small file holding --tree (read only)')
    p.add_argument('--tree', help='tree path, the same in --extra and in the output')
    p.add_argument('--ref', help='tree in --base that --tree must line up with')
    p.add_argument('--graft', nargs=3, action='append', default=[], metavar=('EXTRA', 'TREE', 'REF'),
                   help='one more tree to add: file, tree path, ref tree in --base (repeatable)')
    p.add_argument('--output', required=True, help='new file; must not exist')
    p.add_argument('--dry-run', action='store_true', help='run the checks only, copy and write nothing')
    p.add_argument('--replace', action='store_true',
                   help='every tree must already exist in --base; replace it (all cycles) in the copy')
    a = p.parse_args()

    single = (a.extra, a.tree, a.ref)
    if any(single) and not all(single):
        fail('--extra, --tree and --ref must be given together')
    grafts = ([single] if all(single) else []) + [tuple(g) for g in a.graft]
    if not grafts:
        fail('nothing to graft: give --extra/--tree/--ref and/or --graft EXTRA TREE REF')
    trees = [g[1] for g in grafts]
    if len(set(trees)) != len(trees):
        fail(f'a tree is given more than once: {trees}')

    base, out = os.path.abspath(a.base), os.path.abspath(a.output)
    grafts = [(os.path.abspath(e), t, r) for e, t, r in grafts]
    extras = sorted({e for e, _, _ in grafts})
    partial = out + '.partial'
    for path in [base] + extras:
        if not os.path.isfile(path):
            fail(f'not found: {path}')
    if out == base or out in extras:
        fail('--output must differ from --base and the extra files')
    for path in (out, partial):
        if os.path.exists(path):
            fail(f'refusing to overwrite existing {path}')
    out_dir = os.path.dirname(out)
    need = os.path.getsize(base) + sum(os.path.getsize(e) for e in extras)
    free = shutil.disk_usage(out_dir).free
    print(f'[INFO] base {os.path.getsize(base)/1e9:.2f} GB; free space in {out_dir}: {free/1e9:.1f} GB '
          f'(filesystem; a user quota may be lower)')
    if free < need:
        fail(f'not enough free space for the copy ({need/1e9:.1f} GB needed)')

    fb = open_read(base)
    fext = {e: open_read(e) for e in extras}
    counts = {}
    for extra, tree, ref in grafts:
        fe = fext[extra]
        new_t = get_tree(fe, tree)
        if new_t is None:
            fail(f'{tree} not found in {extra}')
        old_t = get_tree(fb, tree)
        if old_t is not None and not a.replace:
            fail(f'{tree} already exists in {base}; use --replace to overwrite it in the copy')
        if a.replace and old_t is None:
            fail(f'--replace given but {tree} is not in {base}')
        if a.replace:
            if branch_names(old_t) != branch_names(new_t):
                only_old = sorted(set(branch_names(old_t)) - set(branch_names(new_t)))
                only_new = sorted(set(branch_names(new_t)) - set(branch_names(old_t)))
                fail(f'{tree}: branches differ from the one in --base (only old: {only_old}, only new: {only_new})')
            print(f'[INFO] --replace: existing {tree} in --base has {old_t.GetEntries()} entries; '
                  f'same branches; it will be replaced in the copy.')
        if get_tree(fb, ref) is None:
            fail(f'{ref} not found in {base}')
        dest_dir_path = os.path.dirname(tree)
        if dest_dir_path and not fb.GetDirectory(dest_dir_path):
            fail(f'directory {dest_dir_path} not found in {base}')
        counts[tree] = check_alignment(fe, tree, fb, ref)
    fb.Close()

    verb = 'replace' if a.replace else 'add'
    if a.dry_run:
        print(f'[DRY RUN] would copy {base} -> {partial}, {verb} {len(grafts)} tree(s), re-check, rename to {out}:')
        for extra, tree, _ in grafts:
            print(f'[DRY RUN]   {tree} ({counts[tree]} entries) from {extra}')
        return 0

    t0 = time.time()
    print(f'[INFO] copying {base} -> {partial} ...', flush=True)
    shutil.copyfile(base, partial)
    print(f'[INFO] copy done in {time.time() - t0:.0f} s', flush=True)

    fo = ROOT.TFile.Open(partial, 'UPDATE')
    if not fo or fo.IsZombie():
        fail(f'cannot open {partial} for update (left in place)')
    for extra, tree, _ in grafts:
        dest_dir_path = os.path.dirname(tree)
        d = fo.GetDirectory(dest_dir_path) if dest_dir_path else fo
        d.cd()
        if a.replace:
            d.Delete(os.path.basename(tree) + ';*')
            if d.Get(os.path.basename(tree)):
                fail(f'could not delete the old {tree} in {partial} (left in place)')
        clone = fext[extra].Get(tree).CloneTree(-1, 'fast')
        clone.SetDirectory(d)
        clone.Write(os.path.basename(tree))
        print(f'[INFO] {verb}: {tree}', flush=True)
    fo.Close()
    for f in fext.values():
        f.Close()

    # Re-open and re-check every tree before giving the file its final name.
    fo = open_read(partial)
    for _, tree, ref in grafts:
        if get_tree(fo, tree) is None:
            fail(f'{tree} missing after write; {partial} left in place')
        check_alignment(fo, tree, fo, ref)
        dest_dir_path = os.path.dirname(tree)
        d = fo.GetDirectory(dest_dir_path) if dest_dir_path else fo
        cycles = [k.GetCycle() for k in d.GetListOfKeys() if k.GetName() == os.path.basename(tree)]
        if len(cycles) != 1:
            fail(f'{tree} has {len(cycles)} key cycles after write (expected 1); {partial} left in place')
    fo.Close()
    os.rename(partial, out)
    print(f'[INFO] wrote {out} ({os.path.getsize(out)/1e9:.2f} GB, {len(grafts)} tree(s)) in {time.time() - t0:.0f} s')
    return 0


if __name__ == '__main__':
    sys.exit(main())
