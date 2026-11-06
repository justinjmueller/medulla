#!/usr/bin/env python3
"""Plots and table images of the detector-variation event counts vs the CV.

Input: the CSV from the POT/count step (detsys_pot_and_counts_norm2p5e20_pctcv.csv).
Output (in --outdir):
  detsys_pct_vs_cv_<set>.png   percent change from the CV of the 2.5e20-POT-scaled
                               selected / signal / background / non-nue counts, one row per
                               sample; error bars are the statistical uncertainty of the
                               ratio, sqrt(1/N_var + 1/N_cv) (raw counts).
  detsys_table_<set>.png       table: corrected POT, all-tree entries, and per category
                               raw count, (scaled to 2.5e20), and % change from the CV.
Styled with spineplot's nue_icarus style and style.stamp (ICARUS label + POT); the tables
copy the layout of spineplot's selection efficiency table (efficiency.py).

  spineplot/.pixi/envs/default/bin/python scripts/plot_detsys_counts.py \
      --csv plots_v31/pot/detsys_pot_and_counts_norm2p5e20_pctcv.csv --outdir plots_v31/pot/png
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'spineplot'))
from style import load_style_from_toml, stamp  # noqa: E402

POT = 2.5e20
CATS = [('n_selected', 'Selected'), ('n_sel_signal', 'Signal (cat = 0)'),
        ('n_sel_bkgd', 'Background (cat > 0)'), ('n_sel_non_nue', r'Non-$\nu_e$ (cat > 1)')]
TITLES = {'main': 'Nominal Detector Variations',
          'nue': r'$\nu_e$-Enhanced Detector Variations'}


def pct_plot(d, setname, out, style):
    cv = d.cv_reference.iloc[0]
    ref = d[d['sample'] == cv].iloc[0]
    v = d[d['sample'] != cv].reset_index(drop=True)
    y = np.arange(len(v))[::-1]
    fig, ax = plt.subplots(figsize=(9, 0.45 * len(v) + 2.2))
    off = np.linspace(-0.27, 0.27, len(CATS))
    for i, ((c, lab), o) in enumerate(zip(CATS, off)):
        pct = v[c + '_pct_vs_cv'].astype(float)
        err = 100 * (1 + pct / 100) * np.sqrt(1 / v[c].astype(float) + 1 / float(ref[c]))
        ax.errorbar(pct, y + o, xerr=err, fmt=style.get_marker(i), ms=6, color=style.get_color(i),
                    elinewidth=1.5, capsize=0, label=lab)
    ax.axvline(0, color='gray', lw=1, ls='--', zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(v['sample'])
    ax.tick_params(axis='y', which='minor', left=False, right=False)
    ax.set_ylim(-0.7, len(v) - 0.3)
    ax.set_xlabel('Change from CV [%]')
    ax.grid(axis='x')
    ax.set_title(TITLES[setname], loc='center')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.09), ncol=4, fontsize=10,
              handletextpad=0.3, columnspacing=1.0)
    stamp(ax, pot=POT)
    fig.savefig(out)
    plt.close(fig)


def table_png(d, setname, out):
    cv = d.cv_reference.iloc[0]
    header = ['Sample', 'Corrected POT', 'All Tree'] + [lab for _, lab in CATS]
    rows = [header]
    for _, r in d.iterrows():
        is_cv = r['sample'] == cv
        cells = ['CV' if is_cv else r['sample'], f"{float(r.corrected_pot):.3e}", f"{int(r.n_all):,}"]
        for c, _ in CATS:
            pct = '' if is_cv else f"  {float(r[c + '_pct_vs_cv']):+.1f}%"
            cells.append(f"{int(r[c]):,} ({float(r[c + '_2p5e20']):,.1f}){pct}")
        rows.append(cells)

    fig, ax = plt.subplots(figsize=(16, 0.34 * len(rows) + 1.1))
    ax.axis('off')
    # Same layout as spineplot's selection efficiency table: header as the first data row,
    # rules above and below it and at the bottom, no other cell edges.
    # bbox = the whole axes, so the table sits below the stamp and fills the figure.
    table = ax.table(cellText=rows, colLabels=None, cellLoc='center', edges='T',
                     colWidths=[0.11, 0.09, 0.09] + [0.1775] * 4, bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    for i in range(2, len(rows)):
        for j in range(len(rows[i])):
            table[i, j].visible_edges = 'B' if i == len(rows) - 1 else 'open'
    for j in range(len(header)):
        table[0, j].get_text().set_weight('bold')
    ax.set_title(TITLES[setname] + '\nraw count (scaled to 2.5e20 POT) and change from CV',
                 loc='center')
    stamp(ax, pot=POT)
    fig.savefig(out)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--csv', required=True)
    p.add_argument('--outdir', required=True)
    args = p.parse_args()
    style = load_style_from_toml('nue_icarus')
    plt.style.use(style.get_style())
    d = pd.read_csv(args.csv, keep_default_na=False)
    os.makedirs(args.outdir, exist_ok=True)
    for s in ('main', 'nue'):
        ds = d[d.set == s]
        out = os.path.join(args.outdir, f'detsys_pct_vs_cv_{s}.png')
        pct_plot(ds, s, out, style)
        print('[INFO] wrote', out)
        out = os.path.join(args.outdir, f'detsys_table_{s}.png')
        table_png(ds, s, out)
        print('[INFO] wrote', out)


if __name__ == '__main__':
    main()
