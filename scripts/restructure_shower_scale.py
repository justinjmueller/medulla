#!/usr/bin/env python3
"""
Copy a merged stage-0 file and add the shower-scale trees under the layout
phase 1 (detsys.cc) expects.

selection_shower_scale.toml (and the nue version) has ONE sample, shower_scale,
with trees named after the scale (medulla applies shower_energy_scale per tree),
so the merge holds:
    events/shower_scale/shower_scale_up[_exposure]
    events/shower_scale/true_signal_shower_scale_up[_exposure]
    events/shower_scale/photon_sideband_shower_scale_up[_exposure]
    (and the same for _down)
detsys reads events/<key>/<tree> and events/<key>/<tree>_exposure with
key = shower_scale_up / shower_scale_down, so this writes:
    events/shower_scale_up/selected[_exposure]
    events/shower_scale_up/true_signal[_exposure]
    events/shower_scale_up/photon_sideband[_exposure]
    (and the same for _down)
plus a copy of events/shower_scale/POT into each new directory if it exists.
The original events/shower_scale/ trees are left as they are.

    python3 scripts/restructure_shower_scale.py \
        --input  <proj>/merged_stage0.root \
        --output <proj>/merged_stage0_v32.root [--dry-run]

Never writes to --input, refuses an existing --output or <output>.partial,
never deletes anything (a failed run leaves <output>.partial). Needs PyROOT
(SL7 container: scripts/in_sl7.sh).
"""
import argparse
import os
import shutil
import sys
import time

import ROOT

ROOT.gROOT.SetBatch(True)
SRC_DIR = 'events/shower_scale'
SCALES = ['shower_scale_up', 'shower_scale_down']


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


def plan(f):
    """Return [(src_name, dest_dir, dest_name, class)] for every object to copy."""
    d = f.GetDirectory(SRC_DIR)
    if not d:
        fail(f'{SRC_DIR} not found')
    names = {k.GetName(): k.GetClassName() for k in d.GetListOfKeys()}
    out = []
    for scale in SCALES:
        dest_dir = f'events/{scale}'
        if f.GetDirectory(dest_dir):
            fail(f'{dest_dir} already exists in the input; nothing to do')
        found = 0
        for name, cls in sorted(names.items()):
            base, exp = (name[:-len('_exposure')], '_exposure') if name.endswith('_exposure') else (name, '')
            if base == scale:
                dest = 'selected'
            elif base.endswith('_' + scale):
                dest = base[:-len('_' + scale)]
            else:
                continue
            out.append((name, dest_dir, dest + exp, cls))
            found += 1
        if not any(dn == 'selected' for _, dd, dn, _ in out if dd == dest_dir):
            fail(f'no {SRC_DIR}/{scale} tree found (keys: {sorted(names)})')
        if 'POT' in names:
            out.append(('POT', dest_dir, 'POT', names['POT']))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', required=True, help='merged stage-0 file (read only)')
    p.add_argument('--output', required=True, help='new file; must not exist')
    p.add_argument('--dry-run', action='store_true', help='list what would be copied; write nothing')
    a = p.parse_args()

    inp, out = os.path.abspath(a.input), os.path.abspath(a.output)
    partial = out + '.partial'
    if not os.path.isfile(inp):
        fail(f'not found: {inp}')
    if out == inp:
        fail('--output must differ from --input')
    for path in (out, partial):
        if os.path.exists(path):
            fail(f'refusing to overwrite existing {path}')
    free = shutil.disk_usage(os.path.dirname(out)).free
    print(f'[INFO] input {os.path.getsize(inp)/1e9:.2f} GB; free space: {free/1e9:.1f} GB')
    if free < 1.1 * os.path.getsize(inp):
        fail('not enough free space for the copy')

    fi = open_read(inp)
    todo = plan(fi)
    counts = {}
    for src, dd, dn, cls in todo:
        o = fi.Get(f'{SRC_DIR}/{src}')
        n = o.GetEntries() if o.InheritsFrom('TTree') else o.GetBinContent(1)
        counts[(dd, dn)] = n
        print(f'[PLAN] {SRC_DIR}/{src} ({cls}, {"entries" if o.InheritsFrom("TTree") else "bin1"}={n:g}) -> {dd}/{dn}')
    fi.Close()
    if a.dry_run:
        print(f'[DRY RUN] would copy {inp} -> {partial}, add {len(todo)} objects, re-check, rename to {out}.')
        return 0

    t0 = time.time()
    print(f'[INFO] copying {inp} -> {partial} ...', flush=True)
    shutil.copyfile(inp, partial)
    print(f'[INFO] copy done in {time.time() - t0:.0f} s', flush=True)

    fo = ROOT.TFile.Open(partial, 'UPDATE')
    if not fo or fo.IsZombie():
        fail(f'cannot open {partial} for update (left in place)')
    for scale in SCALES:
        fo.GetDirectory('events').mkdir(scale)
    for src, dd, dn, cls in todo:
        d = fo.GetDirectory(dd)
        d.cd()
        o = fo.Get(f'{SRC_DIR}/{src}')
        if o.InheritsFrom('TTree'):
            c = o.CloneTree(-1, 'fast')
            c.SetDirectory(d)
            c.Write(dn)
        else:
            c = o.Clone(dn)
            c.SetDirectory(d)
            c.Write(dn)
    fo.Close()

    fo = open_read(partial)
    for (dd, dn), n in counts.items():
        o = fo.Get(f'{dd}/{dn}')
        if not o:
            fail(f'{dd}/{dn} missing after write; {partial} left in place')
        m = o.GetEntries() if o.InheritsFrom('TTree') else o.GetBinContent(1)
        if m != n:
            fail(f'{dd}/{dn}: {m:g} != {n:g} after write; {partial} left in place')
    fo.Close()
    os.rename(partial, out)
    print(f'[INFO] wrote {out} ({os.path.getsize(out)/1e9:.2f} GB) in {time.time() - t0:.0f} s; all {len(counts)} objects re-checked')
    return 0


if __name__ == '__main__':
    sys.exit(main())
