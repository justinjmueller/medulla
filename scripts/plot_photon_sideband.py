#!/usr/bin/env python3
"""Photon sideband vs the signal-region background (NuMIFull nominal MC), area normalized.

Read-only. Reads scan_events.npz from photon_sideband_scan.py (the loose preselection of
events/NuMIFull/all) and events/NuMIFull/nc_sideband for reference. Regions:
  SR background   signal-region cuts
  photon sideband no_muons, no reco electron, leading photon E > --threshold MeV and primary
                  softmax > 0.95
  NC pi0 sideband the current events/NuMIFull/nc_sideband tree (reference)
All three keep only the true_category values in --cats (default 2,4: numu CC FV and NC FV;
1,2,3,4 adds other nue CC and other CC).
Variables: reco visible energy, reco leading-shower energy (the electron in the SR, the photon in
the sideband), true W and true Q^2.

Output in --outdir:
  photon_sb_compare.png       the four variables, SR background vs the two sidebands
  photon_sb_modes_stacked.png each variable stacked by true interaction mode, one row per region,
                              plus the mode fractions
  photon_sb_modes_shapes.png  per mode (QE / MEC / RES / DIS), SR background vs photon sideband shapes
  photon_sb_evis_reweighted.png  W, Q^2 and leading-shower energy with the photon sideband
                              reweighted to the SR background's reco visible energy distribution
                              (what an E_vis-binned fit sees within each bin)
Each file name gets a _cats<list> suffix when --cats is not the default.

  spineplot/.pixi/envs/default/bin/python scripts/plot_photon_sideband.py \
      --outdir plots_v31/photon_sideband
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import uproot

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'spineplot'))
sys.path.insert(0, os.path.join(REPO, 'scripts'))
from style import load_style_from_toml, stamp  # noqa: E402
from photon_sideband_scan import regions, gt  # noqa: E402

FILE = '/exp/icarus/data/users/micarrig/merged_final_v31e_flux.root'
CATNAMES = {1: r'other $\nu_e$ CC', 2: r'$\nu_\mu$ CC', 3: 'other CC', 4: 'NC'}
MODES = [(0, 'QE'), (10, 'MEC'), (1, 'RES'), (2, 'DIS')]
# (key, axis label, bin edges); the key is resolved per region by columns()
VARS = [('evis', 'Reco visible energy [GeV]', np.linspace(0, 4, 21)),
        ('shower', 'Reco leading-shower energy [MeV]', np.linspace(0, 2500, 26)),
        ('W', r'True hadronic invariant mass $W$ [GeV]', np.linspace(0.9, 3.5, 27)),
        ('Q2', r'True $Q^2$ [GeV$^2$]', np.linspace(0, 3, 31))]


def columns(d, shower):
    return dict(evis=d['reco_visible_energy'], shower=d[shower], W=d['true_hadronic_invariant_mass'],
                Q2=d['true_momentum_transfer'], mode=d['true_interaction_mode'])


def norm_hist(ax, v, edges, label, **kw):
    v = v[np.isfinite(v)]
    v = np.clip(v, edges[0], edges[-1] - 1e-9)
    ax.hist(v, bins=edges, weights=np.full(len(v), 1 / max(len(v), 1)), label=f'{label} ({len(v):,})', **kw)


def header(fig, title, left, right):
    fig.suptitle(title, y=0.975)
    hdr = fig.add_axes([left, 0.905, right - left, 0.0001])
    hdr.axis('off')
    stamp(hdr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--threshold', type=float, default=400.0, help='leading photon energy [MeV]')
    ap.add_argument('--cats', default='2,4', help='true_category values kept in every region')
    args = ap.parse_args()
    cats = [int(c) for c in args.cats.split(',')]
    sfx = '' if cats == [2, 4] else '_cats' + ''.join(map(str, cats))
    catlab = ' + '.join(CATNAMES.get(c, str(c)) for c in cats)
    style = load_style_from_toml('nue_icarus')
    plt.style.use(style.get_style())

    d = dict(np.load(os.path.join(args.outdir, 'scan_events.npz')))
    cat = np.where(np.isfinite(d['true_category']), d['true_category'], -1).astype(int)
    sr = regions(d)['SR'] & np.isin(cat, cats)
    sb = (d['reco_electron_multiplicity'] == 0) & gt(d['reco_leading_photon_energy'], args.threshold) \
        & gt(d['reco_leading_photon_primary_softmax'], 0.95) & np.isin(cat, cats)
    a = uproot.open(FILE)['events/NuMIFull/nc_sideband'].arrays(
        ['true_category', 'reco_visible_energy', 'reco_leading_photon_energy', 'true_hadronic_invariant_mass',
         'true_momentum_transfer', 'true_interaction_mode'], library='np')
    a = {k: v[np.isin(a['true_category'], cats)] for k, v in a.items()}
    sbl = f'Photon sideband ($E_\\gamma$ > {args.threshold:.0f} MeV)'
    regs = [('SR background', {k: v[sr] for k, v in columns(d, 'reco_leading_electron_energy').items()},
             dict(histtype='stepfilled', color=style.get_color(0), alpha=0.6)),
            (sbl, {k: v[sb] for k, v in columns(d, 'reco_leading_photon_energy').items()},
             dict(histtype='step', color='black', lw=2)),
            (r'NC $\pi^0$ sideband', columns(a, 'reco_leading_photon_energy'),
             dict(histtype='step', color=style.get_color(1), lw=1.5, ls='--'))]

    # 1) the four variables
    fig, axes = plt.subplots(2, 2, figsize=(13, 10.5), layout='none')
    fig.subplots_adjust(left=0.08, right=0.97, top=0.88, bottom=0.07, wspace=0.25, hspace=0.3)
    for ax, (k, xlab, edges) in zip(axes.flat, VARS):
        for lab, c, kw in regs:
            norm_hist(ax, c[k], edges, lab, **kw)
        ax.set_xlabel(xlab)
        ax.set_ylabel('Fraction of events')
        ax.legend(fontsize=9)
    header(fig, f'Signal-Region Background vs Sidebands: {catlab}', 0.08, 0.97)
    out = os.path.join(args.outdir, f'photon_sb_compare{sfx}.png')
    fig.savefig(out)
    plt.close(fig)
    print('[INFO] wrote', out)

    # 2) stacked by mode, one row per region, last column the mode fractions
    mcol = [style.get_color(i) for i in range(len(MODES))] + ['lightgray']
    mlab = [m for _, m in MODES] + ['Other']
    fig, axes = plt.subplots(3, 5, figsize=(24, 14), layout='none')
    fig.subplots_adjust(left=0.05, right=0.98, top=0.84, bottom=0.06, wspace=0.28, hspace=0.38)
    for row, (lab, c, _) in zip(axes, regs):
        n = len(c['mode'])
        masks = [c['mode'] == m for m, _ in MODES]
        masks.append(~np.any(masks, axis=0))
        for ax, (k, xlab, edges) in zip(row, VARS):
            vals, wts = [], []
            for msk in masks:
                v = c[k][msk]
                v = np.clip(v[np.isfinite(v)], edges[0], edges[-1] - 1e-9)
                vals.append(v)
                wts.append(np.full(len(v), 1 / n))
            ax.hist(vals, bins=edges, weights=wts, stacked=True, histtype='stepfilled', color=mcol, label=mlab)
            ax.set_xlabel(xlab)
            ax.set_ylabel('Fraction of events')
            ax.set_title(lab, loc='center', fontsize=12)
        ax = row[-1]
        fr = [100 * np.mean(m) for m in masks]
        ax.bar(mlab, fr, color=mcol)
        for i, f in enumerate(fr):
            ax.text(i, f + 1, f'{f:.0f}%', ha='center', fontsize=10)
        ax.set_ylim(0, 85)
        ax.set_ylabel('Fraction of events [%]')
        ax.set_title(f'{lab} ({n:,})', loc='center', fontsize=12)
    axes[0][0].legend(fontsize=10)
    header(fig, f'True Interaction Mode: Signal-Region Background vs Sidebands ({catlab})', 0.05, 0.98)
    out = os.path.join(args.outdir, f'photon_sb_modes_stacked{sfx}.png')
    fig.savefig(out)
    plt.close(fig)
    print('[INFO] wrote', out)

    # 3) per-mode shapes, SR background vs photon sideband
    fig, axes = plt.subplots(len(MODES), len(VARS), figsize=(22, 19), layout='none')
    fig.subplots_adjust(left=0.05, right=0.98, top=0.88, bottom=0.05, wspace=0.25, hspace=0.42)
    for row, (m, ml) in zip(axes, MODES):
        for ax, (k, xlab, edges) in zip(row, VARS):
            for lab, c, kw in regs[:2]:
                norm_hist(ax, c[k][c['mode'] == m], edges, lab, **kw)
            ax.set_xlabel(xlab)
            ax.set_ylabel('Fraction of events')
            ax.set_title(ml, loc='center', fontsize=12)
            ax.legend(fontsize=8)
    fig.suptitle(f'Per-Mode Shapes: Signal-Region Background vs Photon Sideband ({catlab})', y=0.985)
    hdr = fig.add_axes([0.05, 0.93, 0.93, 0.0001])
    hdr.axis('off')
    stamp(hdr)
    out = os.path.join(args.outdir, f'photon_sb_modes_shapes{sfx}.png')
    fig.savefig(out)
    plt.close(fig)
    print('[INFO] wrote', out)

    # 4) photon sideband reweighted to the SR background's reco visible energy
    sr_c, sb_c = regs[0][1], regs[1][1]
    ev = np.array([0, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 1e9])
    h_sr = np.histogram(sr_c['evis'], ev)[0] / len(sr_c['evis'])
    h_sb = np.histogram(sb_c['evis'], ev)[0] / len(sb_c['evis'])
    ratio = np.divide(h_sr, h_sb, out=np.zeros_like(h_sr), where=h_sb > 0)
    w = ratio[np.clip(np.digitize(sb_c['evis'], ev) - 1, 0, len(ratio) - 1)]
    fig, axes = plt.subplots(1, 3, figsize=(19, 6.5), layout='none')
    fig.subplots_adjust(left=0.05, right=0.98, top=0.8, bottom=0.12, wspace=0.25)
    for ax, (k, xlab, edges) in zip(axes, [VARS[2], VARS[3], VARS[1]]):
        norm_hist(ax, sr_c[k], edges, 'SR background', **regs[0][2])
        norm_hist(ax, sb_c[k], edges, sbl, histtype='step', color='gray', lw=1.5, ls=':')
        v = sb_c[k]
        ok = np.isfinite(v)
        ax.hist(np.clip(v[ok], edges[0], edges[-1] - 1e-9), bins=edges, weights=w[ok] / w[ok].sum(),
                histtype='step', color='black', lw=2, label=r'Photon sideband, reweighted to SR $E_{vis}$')
        ax.set_xlabel(xlab)
        ax.set_ylabel('Fraction of events')
        ax.legend(fontsize=9)
    fig.suptitle(f'Photon Sideband Reweighted in Reco Visible Energy ({catlab})', y=0.975)
    hdr = fig.add_axes([0.05, 0.84, 0.93, 0.0001])
    hdr.axis('off')
    stamp(hdr)
    out = os.path.join(args.outdir, f'photon_sb_evis_reweighted{sfx}.png')
    fig.savefig(out)
    plt.close(fig)
    print('[INFO] wrote', out)
    print('[INFO] E_vis bins', ev[:-1], 'SR/SB shape ratio', np.round(ratio, 2))


if __name__ == '__main__':
    main()
