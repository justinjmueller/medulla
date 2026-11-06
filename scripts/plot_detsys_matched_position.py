#!/usr/bin/env python3
"""Position dependence of a detector variation, using the neutrinos shared by the CV and the variation.

Read-only. Signal neutrinos in the all_signal tree of the nue-enhanced set are matched between NuMI_nue and
each variation (same key as plot_detsys_matched_energy.py). For the shared neutrinos, the CV population is
those selected in the CV (at their CV reco position) and the variation population those selected in the
variation (at their variation reco position). Position = reco leading-electron start (within 3.5 cm of the
vertex by the selection; the v31 trees have no vertex branch).

Per variation, in --outdir/<variation>/:
  pos2d_{xy,xz,yz}_<variation>.png   CV | variation | CV/var
  pos1d_<variation>.png              x, y, z: CV and variation overlaid, CV/var below
  pos1d_<variation>.csv

Ratio error: the two populations share events, so only neutrinos that land in a bin in one sample but not
the other fluctuate: sigma(CV/var) = sqrt(n_cv + n_var - 2 n_common) / n_var.

  spineplot/.pixi/envs/default/bin/python scripts/plot_detsys_matched_position.py \
      --variation all --outdir plots_v31/detsys_matched_position
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
sys.path.insert(0, os.path.join(REPO, 'scripts'))
from style import load_style_from_toml, stamp  # noqa: E402
from plot_detsys_matched_energy import FILE, CV, VARIATIONS, CUTS, TRUTH, KEY  # noqa: E402

AX = ['x', 'y', 'z']
POS = [f'reco_leading_electron_start_{a}' for a in AX]
BINS1 = {'x': np.arange(-370, 371, 10.0), 'y': np.arange(-190, 141, 10.0), 'z': np.arange(-900, 901, 50.0)}
# 2D bins follow the selected-event volume: 6 bins per cryostat in x (the central gap is one empty bin),
# 6 in y and 12 in z
BINS2 = {'x': np.concatenate([np.linspace(-337, -86, 7), np.linspace(86, 337, 7)]),
         'y': np.linspace(-157, 110, 7), 'z': np.linspace(-865, 845, 13)}
MIN_CV = 10      # 1D ratio points need this many CV events
MIN_CV_2D = 30   # 2D ratio cells need this many CV events


def load(f, s):
    d = f[f'events/{s}/all_signal'].arrays(TRUTH + CUTS + POS, library='pd')
    d['E'] = d.true_neutrino_energy.round(4)
    d['Q2'] = d.true_momentum_transfer.round(4)
    d['Ee'] = d.true_leading_electron_energy.round(2)
    d['k'] = d.groupby(KEY[:-1]).cumcount()
    d['sel'] = np.all([np.nan_to_num(d[c].astype(float)) > 0.5 for c in CUTS], axis=0)
    return d[KEY + POS + ['sel']]


def counts(m, axes, bins):
    """Histograms of the CV- and var-selected populations, and of events in the same bin in both."""
    def idx(sfx):
        ii = [np.digitize(m[f'reco_leading_electron_start_{a}_{sfx}'].to_numpy(), bins[a]) - 1 for a in axes]
        ok = np.all([(i >= 0) & (i < len(bins[a]) - 1) for i, a in zip(ii, axes)], axis=0)
        return ii, ok
    shape = tuple(len(bins[a]) - 1 for a in axes)
    (ic, okc), (iv, okv) = idx('cv'), idx('var')
    sc, sv = m.sel_cv.to_numpy() & okc, m.sel_var.to_numpy() & okv
    hc, hv, hb = np.zeros(shape), np.zeros(shape), np.zeros(shape)
    np.add.at(hc, tuple(i[sc] for i in ic), 1)
    np.add.at(hv, tuple(i[sv] for i in iv), 1)
    same = sc & sv & np.all([a == b for a, b in zip(ic, iv)], axis=0)
    np.add.at(hb, tuple(i[same] for i in ic), 1)
    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.where(hv > 0, hc / hv, np.nan)
        e = np.where(hv > 0, np.sqrt(np.clip(hc + hv - 2 * hb, 0, None)) / hv, np.nan)
    return hc, hv, hb, r, e


def fig2d(m, a1, a2, variation):
    hc, hv, _, r, _ = counts(m, [a1, a2], BINS2)
    b1, b2 = BINS2[a1], BINS2[a2]
    fig, axes = plt.subplots(1, 3, figsize=(20, 6.6), layout='none')
    fig.subplots_adjust(left=0.06, right=0.97, top=0.78, bottom=0.12, wspace=0.38)
    vmax = max(hc.max(), hv.max(), 1)
    for ax, h, lab in ((axes[0], hc, r'CV ($\nu_e$-Enhanced)'), (axes[1], hv, variation)):
        mm = ax.pcolormesh(b1, b2, np.where(h > 0, h, np.nan).T, cmap='Blues', vmin=0, vmax=vmax)
        fig.colorbar(mm, ax=ax, label='Selected neutrinos', fraction=0.05, pad=0.03)
        ax.set_title(lab, loc='center', fontsize=13)
    rr = np.where(hc >= MIN_CV_2D, r, np.nan)
    mm = axes[2].pcolormesh(b1, b2, rr.T, cmap='RdBu_r', vmin=0.85, vmax=1.15)
    fig.colorbar(mm, ax=axes[2], label=f'CV / {variation}', fraction=0.05, pad=0.03)
    axes[2].set_title(f'CV / {variation}', loc='center', fontsize=13)
    for ax in axes:
        ax.set_xlabel(f'Leading-electron start {a1} [cm]')
        ax.set_ylabel(f'Leading-electron start {a2} [cm]')
        ax.set_facecolor('#f2f2f2')
    fig.suptitle(rf'Selected $\nu_e$ CC Signal (Shared Neutrinos): {a1}-{a2}, CV vs {variation}', y=0.975)
    hdr = fig.add_axes([0.06, 0.885, 0.91, 0.0001])
    hdr.axis('off')
    stamp(hdr)
    return fig


def draw_points(ax, edges, y, err, **kw):
    c = 0.5 * (edges[1:] + edges[:-1])
    ok = np.isfinite(y)
    ax.errorbar(c[ok], y[ok], yerr=err[ok], fmt='o', ms=3.5, capsize=0, **kw)


def fig1d(m, variation, style):
    fig, axes = plt.subplots(2, 3, figsize=(19, 7.8), layout='none', sharex='col',
                             gridspec_kw=dict(height_ratios=[3, 1.3]))
    fig.subplots_adjust(left=0.06, right=0.98, top=0.83, bottom=0.10, wspace=0.25, hspace=0.06)
    rows = []
    for j, a in enumerate(AX):
        edges = BINS1[a]
        hc, hv, hb, r, e = counts(m, [a], BINS1)
        top, bot = axes[0, j], axes[1, j]
        top.stairs(hc, edges, color='black', lw=1.5, label=rf'CV ($\nu_e$-Enhanced) ({int(hc.sum()):,})')
        top.stairs(hv, edges, color=style.get_color(0), lw=1.5, label=f'{variation} ({int(hv.sum()):,})')
        draw_points(top, edges, np.where(hc > 0, hc, np.nan), np.sqrt(hc), color='black')
        draw_points(top, edges, np.where(hv > 0, hv, np.nan), np.sqrt(hv), color=style.get_color(0))
        top.set_ylabel('Selected neutrinos')
        top.set_ylim(0, 1.25 * max(hc.max(), hv.max(), 1))
        top.legend(fontsize=9, loc='upper right')
        rr = np.where(hc >= MIN_CV, r, np.nan)
        draw_points(bot, edges, rr, e, color='black')
        bot.axhline(1, color='gray', ls='--', lw=1)
        ok = np.isfinite(rr)
        if ok.any():
            lo, hi = np.nanmin(rr[ok] - e[ok]), np.nanmax(rr[ok] + e[ok])
            pad = 0.1 * max(hi - lo, 0.02)
            bot.set_ylim(min(lo, 1) - pad, max(hi, 1) + pad)
        bot.set_ylabel(f'CV / {variation}', fontsize=10)
        bot.set_xlabel(f'Leading-electron start {a} [cm]')
        bot.set_xlim(edges[0], edges[-1])
        for k in range(len(hc)):
            rows.append(dict(axis=a, lo=edges[k], hi=edges[k + 1], n_cv=int(hc[k]), n_var=int(hv[k]),
                             n_same_bin=int(hb[k]), ratio=r[k], ratio_err=e[k]))
    fig.suptitle(rf'Selected $\nu_e$ CC Signal (Shared Neutrinos): Position, CV vs {variation}', y=0.975)
    hdr = fig.add_axes([0.06, 0.885, 0.92, 0.0001])
    hdr.axis('off')
    stamp(hdr)
    return fig, pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--variation', nargs='+', default=['recombination'], help="variation name(s), or 'all'")
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    style = load_style_from_toml('nue_icarus')
    plt.style.use(style.get_style())
    variations = VARIATIONS if args.variation == ['all'] else args.variation

    f = uproot.open(FILE)
    a = load(f, CV)
    for variation in variations:
        out = os.path.join(args.outdir, variation)
        os.makedirs(out, exist_ok=True)
        m = a.merge(load(f, variation), on=KEY, suffixes=('_cv', '_var'))
        print(f'[INFO] {variation}: shared {len(m)}, selected cv {int(m.sel_cv.sum())}, var {int(m.sel_var.sum())}, '
              f'both {int((m.sel_cv & m.sel_var).sum())}', flush=True)
        for a1, a2 in (('x', 'y'), ('x', 'z'), ('z', 'y')):
            fig = fig2d(m, a1, a2, variation)
            fig.savefig(os.path.join(out, f'pos2d_{a1}{a2}_{variation}.png'))
            plt.close(fig)
        fig, df = fig1d(m, variation, style)
        fig.savefig(os.path.join(out, f'pos1d_{variation}.png'))
        plt.close(fig)
        df.to_csv(os.path.join(out, f'pos1d_{variation}.csv'), index=False, float_format='%.5g')
        for ax_ in AX:
            d = df[(df.axis == ax_) & (df.n_cv >= MIN_CV)]
            tot = d.n_cv.sum() / d.n_var.sum()
            chi2 = np.nansum(((d.ratio - 1) / d.ratio_err) ** 2)
            shape = np.nansum(((d.ratio - tot) / d.ratio_err) ** 2)
            print(f'[INFO]   {ax_}: total CV/var {tot:.4f}, chi2 vs 1 = {chi2:.1f}, chi2 vs total (shape) = '
                  f'{shape:.1f}, nbins {len(d)}', flush=True)
        print(f'[INFO] wrote {out}', flush=True)


if __name__ == '__main__':
    main()
