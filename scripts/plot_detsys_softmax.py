#!/usr/bin/env python3
"""Softmax scores of the leading particles in selected events: CV vs every detector variation.

Read-only. Uses events/<sample>/selected from the merged files and the corrected POT from the
POT CSV; every sample is scaled to 2.5e20 POT. CV = cv (nominal) / NuMI_nue (nue-enhanced).

Output, per set, in --outdir:
  softmax_<set>.png          particle-ID softmax of the leading electron, proton, muon, photon
                             and pion (each particle's own class), CV as a black line, each
                             variation as a coloured line; variation/CV ratio underneath.
  primary_softmax_<set>.png  the same for the primary softmax of each leading particle.
Events without that particle (NaN) are left out of its panel.

  spineplot/.pixi/envs/default/bin/python scripts/plot_detsys_softmax.py \
      --csv plots_v31/pot/detsys_pot_and_counts_norm2p5e20_pctcv.csv --outdir plots_v31/pot/softmax
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'spineplot'))
from style import load_style_from_toml, stamp  # noqa: E402

POT = 2.5e20
FILES = {'main': '/exp/icarus/data/users/micarrig/merged_final_v31e_flux.root',
         'nue': '/exp/icarus/data/users/micarrig/merged_nue_v31.root'}
SET_TITLE = {'main': 'Nominal', 'nue': r'$\nu_e$-Enhanced'}
PARTICLES = ['electron', 'proton', 'muon', 'photon', 'pion']
# (branch, panel title, x range) per page
PAGES = {
    'softmax': ('Leading-Particle Softmax Score',
                [(f'reco_leading_{p}_{p}_softmax', f'Leading {p}: {p} softmax',
                  (0.6, 1.0) if p == 'electron' else (0.5, 1.0) if p == 'photon' else (0.3, 1.0))
                 for p in PARTICLES]),
    'primary_softmax': ('Leading-Particle Primary Softmax Score',
                        [(f'reco_leading_{p}_primary_softmax', f'Leading {p}: primary softmax',
                          (0.95, 1.0) if p == 'electron' else (0.0, 1.0)) for p in PARTICLES]),
}
NBINS = 20
# 14 distinct line colours for the variations (tab20 without its greys); the CV is black.
_T20 = plt.get_cmap('tab20').colors
VAR_COLORS = [_T20[i] for i in (0, 2, 4, 6, 8, 10, 12, 16, 18, 1, 3, 5, 7, 9, 11, 13, 17, 19)]


def load(setname, rows, branches):
    f = uproot.open(FILES[setname])
    out = {}
    for _, r in rows.iterrows():
        a = f[f"events/{r['sample']}/selected"].arrays(branches, library='np')
        a['scale'] = POT / float(r['corrected_pot'])
        out[r['sample']] = a
    return out


def draw_page(data, cv, setname, page_key):
    title, panels = PAGES[page_key]
    variations = [s for s in data if s != cv]
    fig = plt.figure(figsize=(15, 10.5), layout='none')
    outer = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.28, left=0.06, right=0.98, top=0.83, bottom=0.06)
    for i, (branch, ptitle, (lo, hi)) in enumerate(panels):
        r, c = divmod(i, 3)
        inner = outer[r, c].subgridspec(2, 1, height_ratios=[3, 1.3], hspace=0.08)
        ax = fig.add_subplot(inner[0])
        rax = fig.add_subplot(inner[1], sharex=ax)
        edges = np.linspace(lo, hi, NBINS + 1)

        def hist(a):
            v = a[branch]
            v = np.clip(v[np.isfinite(v)], lo, hi - 1e-9)  # out-of-range values go to the end bins
            return np.histogram(v, bins=edges)[0].astype(float)

        n_cv = hist(data[cv])
        cv_scaled = n_cv * data[cv]['scale']
        for k, var in enumerate(variations):
            n_v = hist(data[var])
            col = VAR_COLORS[k]
            ax.stairs(n_v * data[var]['scale'], edges, color=col, lw=1.2, label=var)
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = np.where(n_cv >= 5, n_v * data[var]['scale'] / cv_scaled, np.nan)
            rax.stairs(ratio, edges, color=col, lw=1.0)
        ax.stairs(cv_scaled, edges, color='black', lw=2.2, label='CV', zorder=5)
        rax.axhline(1, color='black', lw=1.5, ls='--', zorder=5)
        ax.set_yscale('log')
        ax.set_title(f'{ptitle} ({int(n_cv.sum())} CV events)', fontsize=12, pad=6)
        ax.set_xlim(lo, hi)
        ax.tick_params(axis='x', labelbottom=False)
        ax.tick_params(axis='y', labelsize=10)
        ax.set_ylabel('Events', fontsize=11)
        rax.set_ylim(0.5, 1.5)
        rax.set_xlabel('Softmax score', fontsize=11)
        rax.tick_params(labelsize=10)
        rax.set_ylabel('Var./CV', fontsize=10)
        if i == 0:
            handles, labs = ax.get_legend_handles_labels()
    # The sixth cell holds the legend (CV first).
    lax = fig.add_subplot(outer[1, 2])
    lax.axis('off')
    order = [labs.index('CV')] + [j for j, l in enumerate(labs) if l != 'CV']
    lax.legend([handles[j] for j in order], [labs[j] for j in order], loc='center', ncol=2,
               fontsize=11, frameon=False, handlelength=2.0)
    fig.suptitle(f'{SET_TITLE[setname]} Detector Variations: {title}', y=0.975, fontsize=15)
    hdr = fig.add_axes([0.06, 0.875, 0.92, 0.0001])
    hdr.axis('off')
    stamp(hdr, pot=POT)
    return fig


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--csv', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--sets', nargs='+', default=['main', 'nue'])
    args = p.parse_args()
    style = load_style_from_toml('nue_icarus')
    plt.style.use(style.get_style())
    d = pd.read_csv(args.csv, keep_default_na=False)
    branches = sorted({b for _, panels in PAGES.values() for b, _, _ in panels})
    os.makedirs(args.outdir, exist_ok=True)
    for s in args.sets:
        rows = d[d.set == s]
        cv = rows.cv_reference.iloc[0]
        data = load(s, rows, branches)
        for key in PAGES:
            fig = draw_page(data, cv, s, key)
            out = os.path.join(args.outdir, f'{key}_{s}.png')
            fig.savefig(out)
            plt.close(fig)
            print('[INFO] wrote', out, flush=True)


if __name__ == '__main__':
    main()
