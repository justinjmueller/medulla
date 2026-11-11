#!/usr/bin/env python3
"""
Every main detector variation's photon-sideband rate vs the det-var CV, per POT,
in the sideband's reco_visible_energy bins. POT per sample is in the legends.

    spineplot/.pixi/envs/default/bin/python analysis/photonSidebandVariationsVsCV.py \
        [--stage0 /exp/icarus/data/users/micarrig/merged_stage0_v32.root] \
        [--outdir plots_v32/photon_sideband_detsys] [--tree photon_sideband]

POT is the CSV's corrected_pot (set main), as used by phase 1: recombination_residual
uses the recombination row, shower_scale_* the cv row. Writes
variations_vs_cv_grid.png, variations_vs_cv_overlay.png, variations_vs_cv_totals.png.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'spineplot'))
os.environ.setdefault('MEDULLA_PLOT_DIR', os.path.join(REPO, 'spineplot') + '/')
from style import stamp  # noqa: E402

CSV = os.path.join(REPO, 'plots_v31/pot/detsys_corrected_pot.csv')
ALIAS = {'recombination_residual': 'recombination', 'shower_scale_up': 'cv', 'shower_scale_down': 'cv'}
VARS = ['induction_high', 'induction_low', 'coherent_high', 'coherent_low', 'intrinsic_high', 'intrinsic_low',
        'lifetime_high', 'lifetime_low', 'scintillation', 'wire_gap', 'yz_uniformity', 'cathode_bend',
        'recombination_residual', 'shower_scale_up', 'shower_scale_down', 'null']
BINS = np.array([0.0, 0.81, 1.1, 1.4, 5.0])
SCALE = 2.5e20


def label(s):
    return f'{s} (no detector effects)' if s == 'null' else s


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--stage0', default='/exp/icarus/data/users/micarrig/merged_stage0_v32.root')
    p.add_argument('--outdir', default=os.path.join(REPO, 'plots_v32/photon_sideband_detsys'))
    p.add_argument('--tree', default='photon_sideband')
    p.add_argument('--var', default='reco_visible_energy', help='branch to histogram')
    p.add_argument('--bins', type=float, nargs='+', default=None, help='bin edges (default: the sideband visE bins)')
    p.add_argument('--xlabel', default=None)
    a = p.parse_args()
    bins = np.array(a.bins) if a.bins else BINS
    xlabel = a.xlabel or a.var
    tag = '' if a.var == 'reco_visible_energy' else f'_{a.var}'
    os.makedirs(a.outdir, exist_ok=True)
    plt.style.use(os.path.join(REPO, 'spineplot/styles/nue2024.mplstyle'))

    pot_tab = pd.read_csv(CSV, keep_default_na=False)
    pot_tab = pot_tab[pot_tab['set'] == 'main'].set_index('sample')['corrected_pot'].astype(float)
    f = uproot.open(a.stage0)

    n, pot = {}, {}
    for s in ['cv'] + VARS:
        x = f[f'events/{s}/{a.tree}'].arrays([a.var], library='np')[a.var]
        n[s] = np.histogram(x[np.isfinite(x)], bins=bins)[0].astype(float)
        pot[s] = float(pot_tab[ALIAS.get(s, s)])
    h = {s: n[s] * SCALE / pot[s] for s in n}
    eh = {s: np.sqrt(n[s]) * SCALE / pot[s] for s in n}
    tot = {s: n[s].sum() / (pot[s] / 1e18) for s in n}
    etot = {s: np.sqrt(n[s].sum()) / (pot[s] / 1e18) for s in n}

    print(f"{'sample':24s}{'events':>8s}{'POT':>11s}{'per1e18':>9s}{'vs cv':>14s}")
    for s in ['cv'] + VARS:
        r = tot[s] / tot['cv']
        er = r * np.sqrt(1 / n[s].sum() + 1 / n['cv'].sum())
        print(f'{s:24s}{int(n[s].sum()):8d}{pot[s]:11.3e}{tot[s]:9.3f}{r:8.3f}+-{er:.3f}')

    centers = 0.5 * (bins[:-1] + bins[1:])
    print(f'{a.var}: totals above count only events inside the binning')

    # 1. One panel per variation: CV vs variation, with a ratio sub-panel
    fig = plt.figure(figsize=(24, 26))
    outer = fig.add_gridspec(4, 4, hspace=0.32, wspace=0.25, top=0.94)
    for i, s in enumerate(VARS):
        sub = outer[i // 4, i % 4].subgridspec(2, 1, height_ratios=[3, 1.3], hspace=0.05)
        ax0 = fig.add_subplot(sub[0])
        ax1 = fig.add_subplot(sub[1], sharex=ax0)
        r_tot = tot[s] / tot['cv']
        e_tot = r_tot * np.sqrt(1 / n[s].sum() + 1 / n['cv'].sum())
        col = 'tab:purple' if s == 'null' else 'tab:red'
        ax0.stairs(h['cv'], bins, color='black', lw=1.8, label=f'cv: {pot["cv"]:.3g} POT')
        ax0.errorbar(centers, h['cv'], yerr=eh['cv'], fmt='none', ecolor='black', capsize=2)
        ax0.stairs(h[s], bins, color=col, lw=1.8, label=f'{label(s)}: {pot[s]:.3g} POT')
        ax0.errorbar(centers, h[s], yerr=eh[s], fmt='none', ecolor=col, capsize=2)
        ax0.set_ylim(0, 1.35 * max(h['cv'].max(), h[s].max()))
        ax0.set_title(f'{label(s)}   total/cv = {r_tot:.3f} $\\pm$ {e_tot:.3f}', fontsize=13)
        ax0.legend(fontsize=10, loc='lower left')
        ax0.tick_params(labelbottom=False)
        if i % 4 == 0:
            ax0.set_ylabel(f'Events / {SCALE:.1e} POT', fontsize=12)
        with np.errstate(divide='ignore', invalid='ignore'):
            r = h[s] / h['cv']
            er = r * np.sqrt(1 / n[s] + 1 / n['cv'])
        ax1.stairs(r, bins, color=col, lw=1.8)
        ax1.errorbar(centers, r, yerr=er, fmt='none', ecolor=col, capsize=2)
        ax1.axhline(1, color='black', lw=1)
        ax1.set_ylim(0.7, 1.4)
        ax1.set_xlabel(xlabel, fontsize=11)
        if i % 4 == 0:
            ax1.set_ylabel('var / cv', fontsize=11)
    # (no suptitle: the style's layout puts it on top of the first row; each panel is titled)
    fig.savefig(os.path.join(a.outdir, f'variations_vs_cv{tag}_grid.png'), dpi=90, bbox_inches='tight')
    plt.close(fig)

    # 2. All variations' ratio to the CV in one panel, CV stat band around 1
    fig, ax = plt.subplots(figsize=(14, 9))
    cv_err = 1 / np.sqrt(n['cv'])
    ax.stairs(1 + cv_err, bins, baseline=1 - cv_err, fill=True, color='gray', alpha=0.3,
              label=f'cv stat. ($\\pm 1\\sigma$): {pot["cv"]:.3g} POT')
    cmap = plt.get_cmap('tab20')
    for i, s in enumerate(VARS):
        ls = ':' if s == 'null' else ('--' if s.startswith('shower') or s.startswith('recomb') else '-')
        ax.stairs(h[s] / h['cv'], bins, color=cmap(i % 20), lw=1.8, ls=ls,
                  label=f'{label(s)}: {pot[s]:.3g} POT')
    ax.axhline(1, color='black', lw=1.2)
    ax.set_ylim(0.75, 1.35)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Variation / det-var CV (per POT)')
    ax.legend(fontsize=9, ncol=2, loc='upper left', framealpha=0.9)
    stamp(ax, pot=SCALE, work_in_progress=True)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, f'variations_vs_cv{tag}_overlay.png'), dpi=110)
    plt.close(fig)

    # 3. Total rate per 1e18 POT for every sample
    order = sorted(['cv'] + VARS, key=lambda s: tot[s])
    fig, ax = plt.subplots(figsize=(11, 10))
    y = np.arange(len(order))
    for i, s in enumerate(order):
        col = 'black' if s == 'cv' else ('tab:purple' if s == 'null' else 'tab:red')
        ax.errorbar(tot[s], i, xerr=etot[s], fmt='o', color=col, capsize=3, ms=7)
    ax.axvline(tot['cv'], color='black', lw=1, ls='--')
    ax.axvspan(tot['cv'] - etot['cv'], tot['cv'] + etot['cv'], color='gray', alpha=0.25)
    ax.set_yticks(y)
    ax.set_yticklabels([f'{label(s)}  ({pot[s]:.3g} POT)' for s in order], fontsize=11)
    ax.set_xlabel(f'{a.tree} events per 1e18 POT (in the {a.var} binning)')
    ax.set_title('Total photon-sideband rate per sample (det-var CV in black, its stat. band in gray)', fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, f'variations_vs_cv{tag}_totals.png'), dpi=110)
    plt.close(fig)
    print(f'wrote {a.outdir}/variations_vs_cv{tag}_{{grid,overlay,totals}}.png')


if __name__ == '__main__':
    main()
