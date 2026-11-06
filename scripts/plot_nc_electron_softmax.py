#!/usr/bin/env python3
"""Softmax scores of the leading electron in selected NC events (NuMIFull nominal MC).

Read-only. Three panels: electron softmax, photon softmax and primary softmax of the leading
reco electron. NC (true_category 4) is the filled histogram; signal (category 0) is overlaid as
a line for comparison. Both are normalized to unit area (shape comparison).

  spineplot/.pixi/envs/default/bin/python scripts/plot_nc_electron_softmax.py \
      --outdir plots_v31/background_topologies
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
from style import load_style_from_toml, stamp  # noqa: E402

FILE = '/exp/icarus/data/users/micarrig/merged_final_v31e_flux.root'
TREE = 'events/NuMIFull/selected'
PANELS = [('reco_leading_electron_electron_softmax', 'Electron softmax', (0.6, 1.0)),
          ('reco_leading_electron_photon_softmax', 'Photon softmax', (0.0, 0.4)),
          ('reco_leading_electron_primary_softmax', 'Primary softmax', (0.95, 1.0))]
NBINS = 20


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    style = load_style_from_toml('nue_icarus')
    plt.style.use(style.get_style())

    a = uproot.open(FILE)[TREE].arrays(['true_category'] + [p[0] for p in PANELS], library='np')
    cat = a['true_category']
    groups = [('NC', cat == 4, dict(histtype='stepfilled', color=style.get_color(0), alpha=0.6)),
              (r'Signal $\nu_e$ CC', cat == 0, dict(histtype='step', color='black', lw=2))]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6), layout='none')
    fig.subplots_adjust(left=0.06, right=0.98, top=0.84, bottom=0.13, wspace=0.25)
    for ax, (br, title, (lo, hi)) in zip(axes, PANELS):
        edges = np.linspace(lo, hi, NBINS + 1)
        for name, sel, kw in groups:
            v = a[br][sel]
            v = np.clip(v[np.isfinite(v)], lo, hi - 1e-9)  # out-of-range values go to the end bins
            ax.hist(v, bins=edges, weights=np.full(len(v), 1 / len(v)), label=f'{name} ({len(v):,})', **kw)
        ax.set_yscale('log')
        ax.set_xlim(lo, hi)
        ax.set_xlabel(f'Leading electron: {title.lower()}')
        ax.set_ylabel('Fraction of events')
        ax.legend(fontsize=10, loc='upper right' if 'Photon' in title else 'upper left')
    fig.suptitle(r'Leading-Electron Softmax Scores: Selected NC vs Signal (area normalized)', y=0.98)
    hdr = fig.add_axes([0.06, 0.855, 0.92, 0.0001])
    hdr.axis('off')
    stamp(hdr)
    out = os.path.join(args.outdir, 'nc_leading_electron_softmax.png')
    fig.savefig(out)
    print('[INFO] wrote', out)


if __name__ == '__main__':
    main()
