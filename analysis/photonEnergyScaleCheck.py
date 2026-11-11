#!/usr/bin/env python3
"""
Photon energy scale: det-var CV vs each main detector variation.

For each sample, reco_leading_photon_energy / true_leading_photon_energy per
event, in bins of the TRUE leading photon energy. The scale is the median of
that ratio within a window around the peak (default 0.5-1.5, to drop
mismatched photons); its error is 1.2533 * sigma / sqrt(n).

Use a tree WITHOUT a reco photon-energy cut (nc_sideband) to avoid threshold
bias; photon_sideband (reco photon > 400 MeV) is biased near the threshold.

    spineplot/.pixi/envs/default/bin/python analysis/photonEnergyScaleCheck.py \
        [--stage0 .../merged_stage0_v32.root] [--tree nc_sideband] [--outdir plots_v32/photon_sideband_detsys]
"""
import argparse
import os
import sys

import numpy as np
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'spineplot'))
os.environ.setdefault('MEDULLA_PLOT_DIR', os.path.join(REPO, 'spineplot') + '/')
from style import stamp  # noqa: E402

VARS = ['induction_high', 'induction_low', 'coherent_high', 'coherent_low', 'intrinsic_high', 'intrinsic_low',
        'lifetime_high', 'lifetime_low', 'scintillation', 'wire_gap', 'yz_uniformity', 'cathode_bend',
        'recombination_residual', 'shower_scale_up', 'shower_scale_down', 'null']
STD = VARS[:12] + ['null']
TBINS = np.array([0, 150, 300, 500, 800, 1200, 3000])


def scale(r):
    if len(r) < 10:
        return np.nan, np.nan, len(r)
    return np.median(r), 1.2533 * np.std(r) / np.sqrt(len(r)), len(r)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--stage0', default='/exp/icarus/data/users/micarrig/merged_stage0_v32.root')
    p.add_argument('--tree', default='nc_sideband')
    p.add_argument('--outdir', default=os.path.join(REPO, 'plots_v32/photon_sideband_detsys'))
    p.add_argument('--window', type=float, nargs=2, default=[0.5, 1.5])
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    plt.style.use(os.path.join(REPO, 'spineplot/styles/nue2024.mplstyle'))
    f = uproot.open(a.stage0)

    data = {}
    missing = [s for s in ['cv'] + VARS if f'events/{s}/{a.tree}' not in f]
    if missing:
        print(f'[INFO] no {a.tree} tree for {missing}; skipped')
    VARS[:] = [s for s in VARS if s not in missing]
    STD[:] = [s for s in STD if s not in missing]
    for s in ['cv'] + VARS:
        d = f[f'events/{s}/{a.tree}'].arrays(['reco_leading_photon_energy', 'true_leading_photon_energy'], library='np')
        reco, true = d['reco_leading_photon_energy'], d['true_leading_photon_energy']
        ok = np.isfinite(reco) & np.isfinite(true) & (true > 0) & (reco > 0)
        r = reco[ok] / true[ok]
        t = true[ok]
        w = (r > a.window[0]) & (r < a.window[1])
        data[s] = (r[w], t[w], ok.sum(), w.sum())

    print(f'{a.tree}: median(reco/true leading photon E) in [{a.window[0]}, {a.window[1]}]')
    print(f"{'sample':24s}{'n':>7s}{'in win':>8s}{'scale':>9s}{'+-':>7s}{'vs cv':>9s}" +
          ''.join(f'{f"{int(lo)}-{int(hi)}":>11s}' for lo, hi in zip(TBINS[:-1], TBINS[1:])))
    cv_s, _, _ = scale(data['cv'][0])
    per_bin = {}
    for s in ['cv'] + VARS:
        r, t, n, nw = data[s]
        sc, er, _ = scale(r)
        bins = []
        for lo, hi in zip(TBINS[:-1], TBINS[1:]):
            m = (t >= lo) & (t < hi)
            bins.append(scale(r[m]))
        per_bin[s] = bins
        print(f'{s:24s}{n:7d}{nw:8d}{sc:9.4f}{er:7.4f}{sc / cv_s:9.4f}' +
              ''.join(f'{b[0]:11.4f}' for b in bins))

    # pooled standard variations vs cv, per true-energy bin
    print('\nper true-E bin: median scale, cv vs the 13 standard variations pooled (incl. null)')
    pooled_r = np.concatenate([data[s][0] for s in STD])
    pooled_t = np.concatenate([data[s][1] for s in STD])
    rows = []
    for lo, hi in zip(TBINS[:-1], TBINS[1:]):
        c = scale(data['cv'][0][(data['cv'][1] >= lo) & (data['cv'][1] < hi)])
        v = scale(pooled_r[(pooled_t >= lo) & (pooled_t < hi)])
        rows.append((lo, hi, c, v))
        print(f'  true {int(lo):5d}-{int(hi):5d} MeV: cv {c[0]:.4f}+-{c[1]:.4f} (n={c[2]})   '
              f'variations {v[0]:.4f}+-{v[1]:.4f} (n={v[2]})   var/cv {v[0] / c[0]:.4f}')

    # Plot 1: ratio distributions, cv vs pooled variations (area-normalized)
    fig, ax = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={'height_ratios': [2.2, 1]}, sharex=True)
    rb = np.linspace(a.window[0], a.window[1], 51)
    hc, _ = np.histogram(data['cv'][0], bins=rb, density=True)
    hv, _ = np.histogram(pooled_r, bins=rb, density=True)
    ax[0].stairs(hc, rb, color='black', lw=2, label=f'det-var CV (n={len(data["cv"][0])}, median {cv_s:.4f})')
    ax[0].stairs(hv, rb, color='tab:red', lw=2,
                 label=f'13 variations pooled (n={len(pooled_r)}, median {np.median(pooled_r):.4f})')
    ax[0].set_ylabel('Normalized events')
    ax[0].legend(fontsize=12)
    stamp(ax[0], work_in_progress=True)
    with np.errstate(divide='ignore', invalid='ignore'):
        ax[1].stairs(hv / hc, rb, color='tab:red', lw=1.5)
    ax[1].axhline(1, color='black', lw=1)
    ax[1].set_ylim(0.5, 1.5)
    ax[1].set_ylabel('variations / CV')
    ax[1].set_xlabel('reco / true leading photon energy')
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, f'photon_energy_scale_{a.tree}_distribution.png'), dpi=110)
    plt.close(fig)

    # Plot 2: median scale vs true photon energy, cv and every variation
    fig, ax = plt.subplots(figsize=(12, 8))
    cen = 0.5 * (TBINS[:-1] + TBINS[1:])
    cmap = plt.get_cmap('tab20')
    for i, s in enumerate(VARS):
        y = np.array([b[0] for b in per_bin[s]])
        ax.plot(cen, y, 'o-', color=cmap(i % 20), ms=4, lw=1, alpha=0.8, label=s)
    y = np.array([b[0] for b in per_bin['cv']])
    e = np.array([b[1] for b in per_bin['cv']])
    ax.errorbar(cen, y, yerr=e, fmt='s-', color='black', lw=2.5, ms=7, capsize=3, label='det-var CV', zorder=10)
    ax.set_xlabel('True leading photon energy [MeV]')
    ax.set_ylabel('median(reco / true leading photon E)')
    ax.legend(fontsize=9, ncol=2)
    stamp(ax, work_in_progress=True)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, f'photon_energy_scale_{a.tree}_vs_trueE.png'), dpi=110)
    plt.close(fig)
    print(f'wrote {a.outdir}/photon_energy_scale_{a.tree}_{{distribution,vs_trueE}}.png')


if __name__ == '__main__':
    main()
