#!/usr/bin/env python3
"""Reco leading-electron energy and cos(theta_beam) of the same true neutrinos in the CV and a variation (nue-enhanced set).

Read-only. Neutrinos in the all_signal tree are matched between NuMI_nue and the variation on truth
(E_nu, Q2, true electron energy, interaction mode; Run/Subrun/Evt are not unique in these samples).
Left: matched neutrinos with a reco electron energy in both samples. Right: matched neutrinos that
pass the full selection in both samples. Each CV bin (column) is normalized to 1. Phase-1 energy bins plus underflow (< 300 MeV) and overflow (> 5000 MeV), equal width per bin.

  spineplot/.pixi/envs/default/bin/python scripts/plot_detsys_matched_energy.py \
      --variation all --outdir plots_v31/detsys_matched

Also writes shared_events_{nue,main}.png and shared_events.csv: per variation, how many signal neutrinos
(all_signal) are shared with the CV or found in only one sample (nominal set: cv vs variation).
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

FILE = '/exp/icarus/data/users/micarrig/merged_nue_v31.root'
CV = 'NuMI_nue'
MAIN_FILE = '/exp/icarus/data/users/micarrig/merged_final_v31e_flux.root'
MAIN_CV = 'cv'
VARIATIONS = ['induction_high', 'induction_low', 'coherent_high', 'coherent_low', 'intrinsic_high',
              'intrinsic_low', 'lifetime_high', 'lifetime_low', 'scintillation', 'wire_gap',
              'yz_uniformity', 'recombination', 'cathode_bend', 'null']
E_BINS = np.array([300, 390, 490, 600, 730, 870, 1030, 1240, 1500, 2000, 5000], dtype=float)
C_BINS = np.array([-1.0, -0.4, 0.0, 0.2, 0.35, 0.5, 0.625, 0.725, 0.8, 0.85, 0.9, 0.94, 0.97, 1.0])
# branch, bins, tick format, axis label, positive-only, output name
VARS = {'energy': ('reco_leading_electron_energy', E_BINS, '{:.0f}', 'reco electron energy [MeV]', True, 'energy'),
        'angle': ('reco_ele_beam_open_angle', C_BINS, '{:g}', r'reco electron cos$\theta_{beam}$', False, 'angle')}
CUTS = ['reco_fiducial_cut', 'reco_flash_cut', 'reco_containment_cut', 'reco_has_particle_multiplicity',
        'reco_electron_softmax_cut', 'reco_primary_softmax_cut', 'reco_no_muons', 'reco_vertex_distance_cut',
        'reco_leading_ele_energy_cut']
TRUTH = ['true_neutrino_energy', 'true_momentum_transfer', 'true_leading_electron_energy', 'true_interaction_mode']
KEY = ['E', 'Q2', 'Ee', 'true_interaction_mode', 'k']


RECO = [v[0] for v in VARS.values()]


def load(f, s):
    d = f[f'events/{s}/all_signal'].arrays(TRUTH + CUTS + RECO, library='pd')
    d['E'] = d.true_neutrino_energy.round(4)
    d['Q2'] = d.true_momentum_transfer.round(4)
    d['Ee'] = d.true_leading_electron_energy.round(2)
    d['k'] = d.groupby(KEY[:-1]).cumcount()
    d['sel'] = np.all([np.nan_to_num(d[c].astype(float)) > 0.5 for c in CUTS], axis=0)
    return d[KEY + RECO + ['sel']]


def panel(fig, ax, x, y, title, variation, bins, fmt, xlabel):
    # bin index space: 0 = underflow, 1..N = phase-1 bins, N+1 = overflow; equal width per bin
    n = len(bins) + 1
    ix = np.digitize(x, bins)
    iy = np.digitize(y, bins)
    h = np.zeros((n, n))
    np.add.at(h, (ix, iy), 1)
    # column-normalized: fraction of the CV bin's neutrinos that land in each variation bin
    col = h.sum(axis=1, keepdims=True)
    frac = np.divide(h, col, out=np.zeros_like(h), where=col > 0)
    m = ax.pcolormesh(np.arange(n + 1), np.arange(n + 1), np.where(h > 0, frac, np.nan).T,
                      cmap='Blues', vmin=0, vmax=1)
    fig.colorbar(m, ax=ax, label='Fraction of CV bin')
    for i in range(n):
        for j in range(n):
            if h[i, j] > 0:
                txt = f'{frac[i, j]:.2f}' if frac[i, j] >= 0.005 else '<.01'
                ax.text(i + 0.5, j + 0.5, txt, ha='center', va='center', fontsize=7,
                        color='white' if frac[i, j] > 0.5 else 'black')
    for e in (1, n - 1):
        ax.axvline(e, color='black', lw=1)
        ax.axhline(e, color='black', lw=1)
    ax.plot([0, n], [0, n], color='gray', ls='--', lw=1)
    ticks = np.arange(1, n)
    labels = [fmt.format(e) for e in bins]
    ax.set_xticks(ticks, labels, rotation=45, fontsize=9)
    ax.set_yticks(ticks, labels, fontsize=9)
    flow = [0.5, n - 0.5]
    flab = ['<' + fmt.format(bins[0]), '>' + fmt.format(bins[-1])]
    ax.set_xticks(flow, flab, minor=True, rotation=45, fontsize=9)
    ax.set_yticks(flow, flab, minor=True, fontsize=9)
    ax.tick_params(which='minor', length=0)
    ax.set_xlim(0, n)
    ax.set_ylim(0, n)
    ax.set_xlabel(r'CV ($\nu_e$-Enhanced) ' + xlabel)
    ax.set_ylabel(f'{variation} ' + xlabel)
    ax.set_title(title, loc='center', fontsize=13)
    stamp(ax)
    diag = np.trace(h) / max(h.sum(), 1)
    print(f'[INFO] {title}: {len(x)} neutrinos, underflow cv {int((ix == 0).sum())} var {int((iy == 0).sum())}, '
          f'overflow cv {int((ix == n - 1).sum())} var {int((iy == n - 1).sum())}; same bin {100 * diag:.1f}%; '
          f'var bin lower {100 * np.tril(h, -1).sum() / h.sum():.1f}%, higher {100 * np.triu(h, 1).sum() / h.sum():.1f}%')
    print(f'[INFO]   median var - cv = {np.median(y - x):+.4g}, median var/cv - 1 = {100 * np.median(y / x - 1):+.2f}%')
    return h


def match_counts(a, b):
    m = a[KEY].merge(b[KEY], on=KEY, how='outer', indicator=True)['_merge'].value_counts()
    return dict(n_cv=len(a), n_var=len(b), shared=int(m.get('both', 0)),
                cv_only=int(m.get('left_only', 0)), var_only=int(m.get('right_only', 0)))


def table_page(df, title):
    hdr = ['Variation', 'CV signal $\\nu$', 'Variation signal $\\nu$', 'Shared', 'CV only', 'Variation only',
           'Shared / CV [%]', 'Shared / Var. [%]']
    rows = [hdr] + [[r.variation, f'{r.n_cv:,}', f'{r.n_var:,}', f'{r.shared:,}', f'{r.cv_only:,}', f'{r.var_only:,}',
                     f'{100 * r.shared / r.n_cv:.1f}', f'{100 * r.shared / r.n_var:.1f}'] for r in df.itertuples()]
    fig, ax = plt.subplots(figsize=(15, 0.42 * len(rows) + 2.0))
    ax.axis('off')
    ncol = len(hdr)
    tab = ax.table(cellText=rows, cellLoc='center', edges='T', bbox=[0, 0, 1, 1],
                   colWidths=[0.18] + [0.82 / (ncol - 1)] * (ncol - 1))
    tab.auto_set_font_size(False)
    tab.set_fontsize(11)
    for i in range(2, len(rows)):
        for j in range(ncol):
            tab[i, j].visible_edges = 'B' if i == len(rows) - 1 else 'open'
    for j in range(ncol):
        tab[0, j].get_text().set_weight('bold')
    ax.set_title(title, loc='center')
    stamp(ax)
    return fig


def matrix_page(m, var, variation):
    br, bins, fmt, xlabel, positive, name = VARS[var]
    xc, xv = m[br + '_cv'].to_numpy(), m[br + '_var'].to_numpy()
    reco = np.isfinite(xc) & np.isfinite(xv)
    if positive:
        reco &= (xc > 0) & (xv > 0)
    sel = reco & m.sel_cv.to_numpy() & m.sel_var.to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(17, 7.6), layout='none')
    fig.subplots_adjust(left=0.07, right=0.97, top=0.83, bottom=0.14, wspace=0.22)
    panel(fig, axes[0], xc[reco], xv[reco], 'Matched, Reconstructed in Both', variation, bins, fmt, xlabel)
    panel(fig, axes[1], xc[sel], xv[sel], 'Matched, Selected in Both', variation, bins, fmt, xlabel)
    what = 'Electron Energy' if var == 'energy' else r'Electron cos$\theta_{beam}$'
    fig.suptitle(rf'Reco {what} of Matched $\nu_e$ CC Signal: CV vs {variation}', y=0.985)
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--variation', nargs='+', default=['recombination'], help="variation name(s), or 'all'")
    ap.add_argument('--var', nargs='+', default=['energy', 'angle'], choices=list(VARS))
    ap.add_argument('--no-table', action='store_true', help='skip the shared/separate count tables')
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    style = load_style_from_toml('nue_icarus')
    plt.style.use(style.get_style())
    variations = VARIATIONS if args.variation == ['all'] else args.variation

    f = uproot.open(FILE)
    a = load(f, CV)
    rows = []
    for variation in variations:
        b = load(f, variation)
        rows.append(dict(set='nue', variation=variation, **match_counts(a, b)))
        m = a.merge(b, on=KEY, suffixes=('_cv', '_var'))
        print(f'[INFO] {variation}: {CV} {len(a)}, {variation} {len(b)}, matched {len(m)}', flush=True)
        for var in args.var:
            fig = matrix_page(m, var, variation)
            out = os.path.join(args.outdir, f'matched_{VARS[var][5]}_{variation}.png')
            fig.savefig(out)
            plt.close(fig)
            print('[INFO] wrote', out, flush=True)

    if args.no_table:
        return
    fm = uproot.open(MAIN_FILE)
    am = load(fm, MAIN_CV)
    for variation in variations:
        rows.append(dict(set='main', variation=variation, **match_counts(am, load(fm, variation))))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.outdir, 'shared_events.csv'), index=False)
    print(df.to_string(index=False))
    for sname, title in (('nue', r'$\nu_e$-Enhanced Signal Neutrinos Shared with the CV (NuMI_nue)'),
                         ('main', 'Nominal Signal Neutrinos Shared with the CV (cv)')):
        fig = table_page(df[df.set == sname], title)
        out = os.path.join(args.outdir, f'shared_events_{sname}.png')
        fig.savefig(out)
        plt.close(fig)
        print('[INFO] wrote', out)


if __name__ == '__main__':
    main()
