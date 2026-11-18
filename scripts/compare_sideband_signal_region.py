#!/usr/bin/env python3
"""Photon sideband vs signal-region (SR) background, area-normalized shape comparisons (v32+).

Read-only. Reads two small reco trees of the merged file:
  SR background  events/NuMIFull/selected        with true_category in --cats (default: all != 0)
  sideband       events/NuMIFull/photon_sideband with the same true_category filter
Shapes are weighted by ppfx_cv_weight (unless --unweighted) and normalized per region.

--sideband-from-all rebuilds the sideband from events/NuMIFull/all (read in chunks) instead of the
photon_sideband tree, to try other definitions without a new production: fiducial, flash,
containment, no_muons, leading-photon primary softmax > 0.95, leading-photon energy >
--photon-threshold (default 400 MeV), and no reco electron, or with --allow-electron-below X a leading
reco electron below X MeV (X <= 300 keeps it orthogonal to the SR's > 300 MeV electron cut). `all`
has no ppfx_cv_weight, so both regions are unweighted in this mode. With the defaults it reproduces
the photon_sideband tree (4583 vs 4572 entries in v32b).

Variables: reco visible energy (photon-sideband fit bins marked), reco pi0 mass (two-photon events
only), reco leading-photon energy (events with a reco photon), true W, true Q^2 and reco - true
visible energy. Modes from true_interaction_mode: QE 0, MEC 10, RES 1, DIS 2, COH 3, Other
(CC and NC pooled within a mode).

Output in --outdir (suffix _all, or _cats<list> for --cats):
  compare_<sfx>.png       each variable, SR background vs sideband
  modes_stacked_<sfx>.png each variable stacked by mode, one row per region, plus mode fractions
  modes_shapes_<sfx>.png  per mode, SR background vs sideband shapes
  per_fitbin_<sfx>.png    per reco_visible_energy fit bin: mode fractions, W, Q^2, reco-true E_vis
  summary_<sfx>.txt       counts and mode fractions (inclusive and per fit bin)

  spineplot/.pixi/envs/default/bin/python scripts/compare_sideband_signal_region.py \\
      --outdir plots_v32/sideband_vs_sr [--cats 2,4] [--unweighted] [--dry-run]
      [--sideband-from-all [--photon-threshold 300] [--allow-electron-below 300]]
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'spineplot'))
from style import load_style_from_toml, stamp  # noqa: E402

FILE = '/exp/icarus/data/users/micarrig/merged_final_v32b_flux.root'
TREES = {'SR background': 'events/NuMIFull/selected', 'Photon sideband': 'events/NuMIFull/photon_sideband'}
BRANCHES = ['true_category', 'true_interaction_mode', 'reco_visible_energy', 'true_visible_energy',
            'reco_pi0_invariant_mass', 'reco_leading_photon_energy', 'true_hadronic_invariant_mass',
            'true_momentum_transfer', 'ppfx_cv_weight']
CATNAMES = {1: r'other $\nu_e$ CC', 2: r'$\nu_\mu$ CC', 3: 'other CC', 4: 'NC', 5: 'cosmic', 6: 'other'}
MODES = [(0, 'QE'), (10, 'MEC'), (1, 'RES'), (2, 'DIS'), (3, 'COH')]
FITBINS = np.array([0.0, 0.81, 1.1, 1.4, 5.0])  # photon-sideband reco_visible_energy fit bins [GeV]
# (key, axis label, bin edges, positive-only: values <= 0 count as undefined)
VARS = [('evis', 'Reco visible energy [GeV]', np.linspace(0, 4, 21), False),
        ('mpi0', r'Reco $\pi^0$ invariant mass [MeV]', np.linspace(0, 400, 21), True),
        ('egam', 'Reco leading-photon energy [MeV]', np.linspace(0, 2500, 26), True),
        ('W', r'True hadronic invariant mass $W$ [GeV]', np.linspace(0.9, 3.5, 27), False),
        ('Q2', r'True $Q^2$ [GeV$^2$]', np.linspace(0, 3, 31), False),
        ('dE', r'Reco $-$ true visible energy [GeV]', np.linspace(-1.5, 1.0, 26), False)]
STYLE = {'SR background': dict(histtype='stepfilled', alpha=0.6),
         'Photon sideband': dict(histtype='step', color='black', lw=2)}


def load(path, tree, cats, weighted):
    import uproot
    t = uproot.open(path)[tree]
    missing = [b for b in BRANCHES if b not in t.keys()]
    if missing:
        sys.exit(f'[ERROR] {tree}: missing branches {missing}')
    a = t.arrays(BRANCHES, library='np')
    d, keep = columns(a, cats, weighted)
    print(f'[INFO] {tree}: {len(keep):,} entries, {keep.sum():,} kept (categories {cats})')
    return {k: v[keep] for k, v in d.items()}, len(keep)


def columns(a, cats, weighted):
    cat = np.where(np.isfinite(a['true_category']), a['true_category'], -1).astype(int)
    keep = np.isin(cat, cats)
    w = a['ppfx_cv_weight'] if weighted else np.ones(len(cat))
    w = np.where(np.isfinite(w), w, 1.0)
    d = dict(evis=a['reco_visible_energy'], mpi0=a['reco_pi0_invariant_mass'],
             egam=a['reco_leading_photon_energy'], W=a['true_hadronic_invariant_mass'],
             Q2=a['true_momentum_transfer'], dE=a['reco_visible_energy'] - a['true_visible_energy'],
             mode=a['true_interaction_mode'], cat=cat, w=w)
    return d, keep


CUT_BRANCHES = ['reco_fiducial_cut', 'reco_flash_cut', 'reco_containment_cut', 'reco_no_muons',
                'reco_electron_multiplicity', 'reco_leading_electron_energy', 'reco_leading_photon_primary_softmax']


def load_sideband_from_all(path, cats, photon_thr, ele_below):
    """Rebuild a photon sideband from the `all` tree, in chunks (keeps only passing rows)."""
    import uproot
    br = [b for b in BRANCHES if b != 'ppfx_cv_weight'] + CUT_BRANCHES
    parts, n = [], 0
    for a in uproot.iterate(f'{path}:events/NuMIFull/all', br, step_size='500 MB', library='np'):
        n += len(a['true_category'])
        m = ((a['reco_fiducial_cut'] == 1) & (a['reco_flash_cut'] == 1) & (a['reco_containment_cut'] == 1)
             & (a['reco_no_muons'] == 1) & (a['reco_leading_photon_primary_softmax'] > 0.95)
             & (a['reco_leading_photon_energy'] > photon_thr))
        ele = a['reco_electron_multiplicity'] == 0
        if ele_below is not None:
            ele = ele | ((a['reco_electron_multiplicity'] > 0) & (a['reco_leading_electron_energy'] < ele_below))
        parts.append({k: v[m & ele] for k, v in a.items()})
    a = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    d, keep = columns(a, cats, False)
    print(f'[INFO] events/NuMIFull/all: {n:,} entries scanned, {len(keep):,} pass the sideband cuts, '
          f'{keep.sum():,} kept (categories {cats})')
    return {k: v[keep] for k, v in d.items()}, len(keep)


def defined(v, positive):
    ok = np.isfinite(v)
    return ok & (v > 0) if positive else ok


def mode_masks(mode):
    masks = [mode == m for m, _ in MODES]
    masks.append(~np.any(masks, axis=0))
    return masks, [m for _, m in MODES] + ['Other']


def shape(ax, v, w, edges, positive, label, total_w, **kw):
    """Area-normalized histogram of the defined values; legend shows count and defined fraction."""
    ok = defined(v, positive)
    if not ok.any():
        return
    x = np.clip(v[ok], edges[0], edges[-1] - 1e-9)
    frac = w[ok].sum() / total_w if total_w > 0 else 0
    lab = f'{label} ({ok.sum():,}' + (f', {100 * frac:.0f}% defined)' if positive else ')')
    ax.hist(x, bins=edges, weights=w[ok] / w[ok].sum(), label=lab, **kw)


def header(fig, title, left, right, y):
    fig.suptitle(title, y=0.985)
    hdr = fig.add_axes([left, y, right - left, 0.0001])
    hdr.axis('off')
    stamp(hdr)


def fractions(d, sel=None):
    sel = np.ones(len(d['w']), bool) if sel is None else sel
    masks, labels = mode_masks(d['mode'])
    tot = d['w'][sel].sum()
    return [100 * d['w'][sel & m].sum() / tot if tot > 0 else 0 for m in masks], labels, sel.sum()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input', default=FILE)
    ap.add_argument('--outdir', default=os.path.join(REPO, 'plots_v32', 'sideband_vs_sr'))
    ap.add_argument('--cats', default=None, help='true_category values, e.g. 2,4 (default: 1-6, all backgrounds)')
    ap.add_argument('--unweighted', action='store_true', help='do not weight by ppfx_cv_weight')
    ap.add_argument('--dry-run', action='store_true', help='print what would be read and written, open nothing')
    ap.add_argument('--sideband-from-all', action='store_true',
                    help='rebuild the sideband from the `all` tree (unweighted) instead of photon_sideband')
    ap.add_argument('--photon-threshold', type=float, default=400.0,
                    help='with --sideband-from-all: leading-photon energy cut [MeV] (default 400)')
    ap.add_argument('--allow-electron-below', type=float, default=None,
                    help='with --sideband-from-all: also allow a leading reco electron below this energy '
                         '[MeV]; keep it <= 300 for orthogonality with the SR')
    args = ap.parse_args()
    if not args.sideband_from_all and (args.photon_threshold != 400.0 or args.allow_electron_below is not None):
        ap.error('--photon-threshold / --allow-electron-below need --sideband-from-all')
    if args.allow_electron_below is not None and args.allow_electron_below > 300:
        print('[WARN] --allow-electron-below > 300 MeV overlaps the SR (leading electron > 300 MeV)')
    if args.sideband_from_all:
        args.unweighted = True
    cats = [int(c) for c in args.cats.split(',')] if args.cats else [1, 2, 3, 4, 5, 6]
    sfx = 'all' if args.cats is None else 'cats' + ''.join(map(str, cats))
    sblab = ''
    if args.sideband_from_all:
        sfx += f'_sbAll_g{args.photon_threshold:.0f}'
        sblab = f'$E_\\gamma$ > {args.photon_threshold:.0f} MeV, ' + (
            'no reco e' if args.allow_electron_below is None else f'reco e < {args.allow_electron_below:.0f} MeV')
        if args.allow_electron_below is not None:
            sfx += f'_e{args.allow_electron_below:.0f}'
    elif args.unweighted:
        sfx += '_unweighted'
    catlab = 'all backgrounds' if args.cats is None else ' + '.join(CATNAMES.get(c, str(c)) for c in cats)
    outs = [f'{n}_{sfx}.{e}' for n, e in [('compare', 'png'), ('modes_stacked', 'png'),
                                          ('modes_shapes', 'png'), ('per_fitbin', 'png'), ('summary', 'txt')]]
    if args.dry_run:
        print(f'[DRY RUN] input {args.input} ({"exists" if os.path.exists(args.input) else "MISSING"})')
        for r, t in TREES.items():
            if r != 'SR background' and args.sideband_from_all:
                t = f'events/NuMIFull/all, rebuilt in chunks ({sblab})'
            print(f'[DRY RUN] {r}: {t}, true_category in {cats}')
        print(f'[DRY RUN] branches: {BRANCHES}; weights: {"none" if args.unweighted else "ppfx_cv_weight"}')
        for o in outs:
            p = os.path.join(args.outdir, o)
            print(f'[DRY RUN] would write {p}' + (' (OVERWRITE)' if os.path.exists(p) else ''))
        return

    os.makedirs(args.outdir, exist_ok=True)
    style = load_style_from_toml('nue_icarus')
    plt.style.use(style.get_style())
    regs = {}
    for r, t in TREES.items():
        if r != 'SR background' and args.sideband_from_all:
            regs[r], _ = load_sideband_from_all(args.input, cats, args.photon_threshold, args.allow_electron_below)
        else:
            regs[r], _ = load(args.input, t, cats, not args.unweighted)
    STYLE['SR background']['color'] = style.get_color(0)
    wlab = (', unweighted' if args.unweighted else ', PPFX-weighted') + (f'; sideband: {sblab}' if sblab else '')
    mcol = [style.get_color(i) for i in range(len(MODES))] + ['lightgray']

    # 1) inclusive comparison
    fig, axes = plt.subplots(2, 3, figsize=(19, 11), layout='none')
    fig.subplots_adjust(left=0.06, right=0.98, top=0.87, bottom=0.07, wspace=0.25, hspace=0.3)
    for ax, (k, xlab, edges, pos) in zip(axes.flat, VARS):
        for r, d in regs.items():
            shape(ax, d[k], d['w'], edges, pos, r, d['w'].sum(), **STYLE[r])
        if k == 'evis':
            for x in FITBINS[1:-1]:
                ax.axvline(x, color='gray', ls=':', lw=1)
        ax.set_xlabel(xlab)
        ax.set_ylabel('Fraction of events')
        ax.legend(fontsize=9)
    header(fig, f'Signal-Region Background vs Photon Sideband ({catlab}{wlab})', 0.06, 0.98, 0.915)
    fig.savefig(os.path.join(args.outdir, outs[0]), dpi=300)
    plt.close(fig)

    # 2) stacked by mode, one row per region, last column the mode fractions
    fig, axes = plt.subplots(2, len(VARS) + 1, figsize=(34, 11), layout='none')
    fig.subplots_adjust(left=0.035, right=0.99, top=0.84, bottom=0.08, wspace=0.3, hspace=0.38)
    for row, (r, d) in zip(axes, regs.items()):
        masks, mlab = mode_masks(d['mode'])
        for ax, (k, xlab, edges, pos) in zip(row, VARS):
            ok = defined(d[k], pos)
            tot = d['w'][ok].sum()
            vals = [np.clip(d[k][ok & m], edges[0], edges[-1] - 1e-9) for m in masks]
            wts = [d['w'][ok & m] / tot for m in masks]
            ax.hist(vals, bins=edges, weights=wts, stacked=True, histtype='stepfilled', color=mcol, label=mlab)
            ax.set_xlabel(xlab)
            ax.set_ylabel('Fraction of events')
            ax.set_title(r, loc='center', fontsize=12)
        fr, labels, n = fractions(d)
        ax = row[-1]
        ax.bar(labels, fr, color=mcol)
        for i, f in enumerate(fr):
            ax.text(i, f + 1, f'{f:.0f}%', ha='center', fontsize=10)
        ax.set_ylim(0, 100)
        ax.set_ylabel('Fraction of events [%]')
        ax.set_title(f'{r} ({n:,})', loc='center', fontsize=12)
    axes[0][0].legend(fontsize=10)
    header(fig, f'True Interaction Mode: Signal-Region Background vs Photon Sideband ({catlab}{wlab})',
           0.035, 0.99, 0.9)
    fig.savefig(os.path.join(args.outdir, outs[1]), dpi=300)
    plt.close(fig)

    # 3) per-mode shapes
    fig, axes = plt.subplots(len(MODES), len(VARS), figsize=(30, 23), layout='none')
    fig.subplots_adjust(left=0.04, right=0.99, top=0.9, bottom=0.04, wspace=0.27, hspace=0.45)
    for row, (m, ml) in zip(axes, MODES):
        for ax, (k, xlab, edges, pos) in zip(row, VARS):
            for r, d in regs.items():
                sel = d['mode'] == m
                shape(ax, d[k][sel], d['w'][sel], edges, pos, r, d['w'][sel].sum(), **STYLE[r])
            ax.set_xlabel(xlab)
            ax.set_ylabel('Fraction of events')
            ax.set_title(ml, loc='center', fontsize=12)
            ax.legend(fontsize=8)
    header(fig, f'Per-Mode Shapes: Signal-Region Background vs Photon Sideband ({catlab}{wlab})',
           0.04, 0.99, 0.94)
    fig.savefig(os.path.join(args.outdir, outs[2]), dpi=300)
    plt.close(fig)

    # 4) per reco_visible_energy fit bin
    pv = [v for v in VARS if v[0] in ('W', 'Q2', 'dE')]
    nb = len(FITBINS) - 1
    fig, axes = plt.subplots(nb, 1 + len(pv), figsize=(24, 5.2 * nb), layout='none')
    fig.subplots_adjust(left=0.05, right=0.98, top=0.91, bottom=0.04, wspace=0.27, hspace=0.45)
    lines = []
    for i, row in enumerate(axes):
        lo, hi = FITBINS[i], FITBINS[i + 1]
        blab = f'{lo:.2f} < $E_{{vis}}^{{reco}}$ < {hi:.2f} GeV'
        ax = row[0]
        x = np.arange(len(MODES) + 1)
        for j, (r, d) in enumerate(regs.items()):
            sel = (d['evis'] >= lo) & (d['evis'] < hi)
            fr, labels, n = fractions(d, sel)
            ax.bar(x + (j - 0.5) * 0.4, fr, width=0.4, label=f'{r} ({n:,})',
                   color=STYLE[r].get('color', 'white'), edgecolor='black',
                   alpha=0.6 if j == 0 else 1.0, hatch=None if j == 0 else '//')
            lines.append(f'{blab.replace("$E_{vis}^{reco}$", "Evis")} | {r:16s} n={n:6d} | ' +
                         '  '.join(f'{lb} {f:5.1f}%' for lb, f in zip(labels, fr)))
        ax.set_xticks(x, labels)
        ax.set_ylim(0, 100)
        ax.set_ylabel('Fraction of events [%]')
        ax.set_title(blab, loc='center', fontsize=12)
        ax.legend(fontsize=9)
        for ax, (k, xlab, edges, pos) in zip(row[1:], pv):
            for r, d in regs.items():
                sel = (d['evis'] >= lo) & (d['evis'] < hi)
                shape(ax, d[k][sel], d['w'][sel], edges, pos, r, d['w'][sel].sum(), **STYLE[r])
            ax.set_xlabel(xlab)
            ax.set_ylabel('Fraction of events')
            ax.set_title(blab, loc='center', fontsize=12)
            ax.legend(fontsize=8)
    header(fig, f'Per Fit Bin: Signal-Region Background vs Photon Sideband ({catlab}{wlab})', 0.05, 0.98, 0.945)
    fig.savefig(os.path.join(args.outdir, outs[3]), dpi=300)
    plt.close(fig)

    # 5) summary
    with open(os.path.join(args.outdir, outs[4]), 'w') as f:
        f.write(f'# {args.input}\n# categories {cats} ({catlab}); weights: '
                f'{"none" if args.unweighted else "ppfx_cv_weight"}\n')
        f.write(f'# sideband: {"rebuilt from all: " + sblab if sblab else "photon_sideband tree"}\n\n')
        for r, d in regs.items():
            cs = ', '.join(f'{CATNAMES.get(c, c)}: {(d["cat"] == c).sum()}' for c in cats)
            fr, labels, n = fractions(d)
            f.write(f'{r}: {n} events (weighted {d["w"].sum():.1f}); by category: {cs}\n')
            f.write('  modes: ' + '  '.join(f'{lb} {v:.1f}%' for lb, v in zip(labels, fr)) + '\n')
        sr, sb = (fractions(d)[0] for d in regs.values())
        f.write('  SR/sideband mode-fraction ratio: ' + '  '.join(
            f'{lb} {a / b:.2f}' if b > 0 else f'{lb} -' for lb, a, b in zip(fractions(regs['SR background'])[1], sr, sb)) + '\n')
        f.write('\nPer reco_visible_energy fit bin:\n' + '\n'.join(lines) + '\n')
    for o in outs:
        print('[INFO] wrote', os.path.join(args.outdir, o))


if __name__ == '__main__':
    main()
