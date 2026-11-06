#!/usr/bin/env python3
"""Uncorrected vs corrected POT per detector-variation sample, with file and event counts.

Read-only. Input: plots_v31/pot/detsys_pot_and_counts.csv (from run_corrected_pot.py).
  uncorrected POT = pot_th1d (CAF POT histograms summed over every file)
  failed POT      = POT of subruns where SPINE produced no truth interactions despite GENIE neutrinos
  corrected POT   = corrected_pot (the value used for normalization)
One table per set (main, nue): pot_table_<set>.png, plus pot_table.csv.

  spineplot/.pixi/envs/default/bin/python scripts/plot_pot_table.py --outdir plots_v31/pot/png
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'spineplot'))
from style import load_style_from_toml, stamp  # noqa: E402

CSV = os.path.join(REPO, 'plots_v31/pot/detsys_pot_and_counts.csv')
TITLES = {'main': 'Nominal Detector-Variation Samples: POT', 'nue': r'$\nu_e$-Enhanced Detector-Variation Samples: POT'}
CVS = {'main': 'cv', 'nue': 'NuMI_nue'}


def table_page(d, setname):
    hdr = ['Sample', 'Files', 'Subruns', 'Failed\nsubruns', 'Events', 'Failed\nevents',
           'Uncorrected\nPOT', 'Failed\nPOT', 'Corrected\nPOT', 'Change\n[%]']
    rows = [hdr]
    for r in d.itertuples():
        rows.append([f'CV ({r.sample})' if r.sample == CVS[setname] else r.sample,
                     f'{r.files_ok:,}', f'{r.subruns:,}', f'{r.failed_subruns:,}',
                     f'{r.events:,}', f'{r.failed_events:,}',
                     f'{r.pot_th1d:.4e}', f'{r.failed_pot:.3e}' if r.failed_pot > 0 else '0',
                     f'{r.corrected_pot:.4e}', f'{r.change_pct:+.2f}' if abs(r.change_pct) >= 0.005 else '0.00'])
    ncol = len(hdr)
    fig, ax = plt.subplots(figsize=(19, 0.42 * len(rows) + 2.0))
    ax.axis('off')
    tab = ax.table(cellText=rows, cellLoc='center', edges='T', bbox=[0, 0, 1, 1],
                   colWidths=[0.15] + [0.85 / (ncol - 1)] * (ncol - 1))
    tab.auto_set_font_size(False)
    tab.set_fontsize(11)
    for i in range(2, len(rows)):
        for j in range(ncol):
            tab[i, j].visible_edges = 'B' if i == len(rows) - 1 else 'open'
    for j in range(ncol):
        tab[0, j].get_text().set_weight('bold')
    ax.set_title(TITLES[setname], loc='center')
    stamp(ax)
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--csv', default=CSV)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    style = load_style_from_toml('nue_icarus')
    plt.style.use(style.get_style())

    p = pd.read_csv(args.csv, keep_default_na=False)
    for c in ('pot_th1d', 'failed_pot', 'corrected_pot'):
        p[c] = p[c].astype(float)
    p['change_pct'] = 100 * (p.corrected_pot / p.pot_th1d - 1)
    cols = ['set', 'sample', 'files_listed', 'files_ok', 'files_fail', 'subruns', 'failed_subruns', 'events',
            'failed_events', 'pot_th1d', 'failed_pot', 'corrected_pot', 'change_pct']
    p[cols].to_csv(os.path.join(args.outdir, 'pot_table.csv'), index=False)
    print(p[cols].to_string(index=False))
    for s in ('main', 'nue'):
        d = p[p.set == s]
        d = pd.concat([d[d['sample'] == CVS[s]], d[d['sample'] != CVS[s]]])  # CV first
        fig = table_page(d, s)
        out = os.path.join(args.outdir, f'pot_table_{s}.png')
        fig.savefig(out)
        plt.close(fig)
        print('[INFO] wrote', out)


if __name__ == '__main__':
    main()
