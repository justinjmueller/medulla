#!/usr/bin/env python3
"""Reco particle multiplicities of the selected events: CV vs each detector variation.

Read-only. Uses events/<sample>/selected from the merged files and the corrected POT from
the POT CSV (scripts/run_corrected_pot.py, scripts/detsys_event_counts.py); every sample is
scaled to 2.5e20 POT. CV = cv (nominal) and NuMI_nue (nue-enhanced), as in the CSV.

Output, per set (main / nue), in --outdir/<set>/:
  interaction_type_pct_vs_cv.png
      % change from the CV in the number of selected events of each reco interaction type
      (1e1p0pi, 1e2p0pi, 1e>=3p0pi, 1e1pi, Other), one row per variation (stat. errors).
  all/<variation>.png
      6 panels: electron, proton, pion, photon, muon multiplicity and the truth category
      (signal / other nue CC / non-nue); CV histogram vs variation, ratio underneath.
  <type>/<variation>.png
      the same, for selected events of one reco interaction type only.
  multiplicity_<set>.pdf   all of the above in one file.

  spineplot/.pixi/envs/default/bin/python scripts/plot_detsys_multiplicity.py \
      --csv plots_v31/pot/detsys_pot_and_counts_norm2p5e20_pctcv.csv \
      --outdir plots_v31/pot/multiplicity
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
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
# Reco interaction type (selection/include/variables.h reco_interaction_type), grouped as in
# analysis/nue_plots_reco.toml, with its labels and colours.
TYPES = [('1e1p0pi', '1e1p0$\\pi$', [1]), ('1e2p0pi', '1e2p0$\\pi$', [2]),
         ('1e3p0pi', '1e$\\geq$3p0$\\pi$', [3]), ('1e1pi', '1e1$\\pi$', [4]), ('other', 'Other', [5, 0])]
TYPE_COLORS = ['#2C194C', '#49749F', '#83A1BE', '#8EB38D', '#4F894E']
# (branch, panel title, first bin, last bin; the last bin is an overflow ">= last")
MULTS = [('reco_electron_multiplicity', 'Electron multiplicity', 1, 5),
         ('reco_proton_multiplicity', 'Proton multiplicity', 0, 5),
         ('reco_pion_multiplicity', 'Pion multiplicity', 0, 4),
         ('reco_photon_multiplicity', 'Photon multiplicity', 0, 4),
         ('reco_muon_multiplicity', 'Muon multiplicity', 0, 1)]
TRUTH = [('Signal', lambda c: c == 0), (r'Other $\nu_e$ CC', lambda c: c == 1),
         (r'Non-$\nu_e$', lambda c: c > 1)]
BRANCHES = [m[0] for m in MULTS] + ['reco_reco_interaction_type', 'true_category']


def load(setname, rows):
    f = uproot.open(FILES[setname])
    out = {}
    for _, r in rows.iterrows():
        a = f[f"events/{r['sample']}/selected"].arrays(BRANCHES, library='np')
        a['scale'] = POT / float(r['corrected_pot'])
        out[r['sample']] = a
    return out


def hist_mult(a, branch, lo, hi, mask):
    x = np.clip(a[branch][mask], lo, hi)
    return np.array([np.sum(x == k) for k in range(lo, hi + 1)], dtype=float)


def hist_truth(a, mask):
    c = a['true_category'][mask]
    return np.array([np.sum(fn(c)) for _, fn in TRUTH], dtype=float)


def type_mask(a, codes):
    return np.isin(a['reco_reco_interaction_type'], codes)


def panel(ax, rax, n_cv, n_var, s_cv, s_var, labels, title, var_name, style):
    """CV (filled step) vs variation (points), both scaled; ratio var/CV with stat errors."""
    x = np.arange(len(labels))
    cv, var = n_cv * s_cv, n_var * s_var
    edges = np.append(x - 0.5, x[-1] + 0.5)
    ax.stairs(cv, edges, fill=True, color=style.get_color(1), alpha=0.35, label='CV')
    ax.stairs(cv, edges, color=style.get_color(1), lw=1.5)
    ax.errorbar(x, var, yerr=np.sqrt(n_var) * s_var, fmt='o', ms=4, color=style.get_color(0),
                elinewidth=1.2, capsize=0, label=var_name)
    ax.set_title(title, fontsize=12, pad=6)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0, max(cv.max(), var.max(), 1e-9) * 1.25)
    ax.tick_params(axis='x', labelbottom=False, which='both', bottom=False, top=False)
    ax.tick_params(axis='y', labelsize=10)
    ax.set_ylabel('Events', fontsize=11)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = var / cv
        err = ratio * np.sqrt(1 / n_var + 1 / n_cv)
    ok = (n_cv > 0) & (n_var > 0)
    rax.axhline(1, color='gray', ls='--', lw=1)
    rax.errorbar(x[ok], ratio[ok], yerr=err[ok], fmt='o', ms=4, color=style.get_color(0),
                 elinewidth=1.2, capsize=0)
    rax.set_ylim(0.6, 1.4)
    rax.set_xlim(edges[0], edges[-1])
    rax.set_xticks(x)
    rax.set_xticklabels(labels, fontsize=10)
    rax.tick_params(axis='x', which='minor', bottom=False, top=False)
    rax.tick_params(axis='y', labelsize=9)
    rax.set_ylabel('Var./CV', fontsize=10)


def page(data, cv, var, setname, codes, type_label, style):
    fig = plt.figure(figsize=(15, 9.5), layout='none')  # manual gridspec; the style turns on constrained
    outer = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.28, left=0.06, right=0.98, top=0.83, bottom=0.06)
    a_cv, a_var = data[cv], data[var]
    m_cv = type_mask(a_cv, codes) if codes else np.ones(len(a_cv['true_category']), bool)
    m_var = type_mask(a_var, codes) if codes else np.ones(len(a_var['true_category']), bool)
    for i in range(6):
        r, c = divmod(i, 3)
        inner = outer[r, c].subgridspec(2, 1, height_ratios=[3, 1.1], hspace=0.08)
        ax = fig.add_subplot(inner[0])
        rax = fig.add_subplot(inner[1], sharex=ax)
        if i < len(MULTS):
            b, title, lo, hi = MULTS[i]
            labels = [str(k) for k in range(lo, hi)] + [f'$\\geq${hi}']
            n_cv, n_var = hist_mult(a_cv, b, lo, hi, m_cv), hist_mult(a_var, b, lo, hi, m_var)
        else:
            title, labels = 'Truth category', [t for t, _ in TRUTH]
            n_cv, n_var = hist_truth(a_cv, m_cv), hist_truth(a_var, m_var)
        panel(ax, rax, n_cv, n_var, a_cv['scale'], a_var['scale'], labels, title, var, style)
        if i == 0:
            ax.legend(loc='upper right', fontsize=10)
    sel = 'all selected events' if not codes else f'selected {type_label} events'
    fig.suptitle(f'{SET_TITLE[setname]} Detector Variations: {var} vs CV ({sel})', y=0.975, fontsize=15)
    hdr = fig.add_axes([0.06, 0.875, 0.92, 0.0001])
    hdr.axis('off')
    stamp(hdr, pot=POT)
    return fig


def combined_panel(ax, rax, hists_cv, hists_var, s_cv, s_var, labels, title):
    """CV as unstacked step lines (each topology + total); variation as points (each topology +
    total) spread across the bin; ratio per topology + total with the same x offsets."""
    x = np.arange(len(labels))
    edges = np.append(x - 0.5, x[-1] + 0.5)
    n_var_tot = np.sum(hists_var, axis=0)
    n_cv_tot = np.sum(hists_cv, axis=0)
    series = list(zip(hists_cv, hists_var, TYPE_COLORS))
    series.append((n_cv_tot, n_var_tot, 'black'))
    off = np.linspace(-0.3, 0.3, len(series))
    ymax = 0
    for (n_cv, n_var, col), o in zip(series, off):
        lw = 2.0 if col == 'black' else 1.5
        ax.stairs(n_cv * s_cv, edges, color=col, lw=lw)
        ax.errorbar(x + o, n_var * s_var, yerr=np.sqrt(n_var) * s_var, fmt='o', ms=4, color=col,
                    elinewidth=1, capsize=0)
        ymax = max(ymax, (n_cv * s_cv).max(), (n_var * s_var).max())
    ax.set_title(title, fontsize=12, pad=6)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0, max(ymax, 1e-9) * 1.15)
    ax.tick_params(axis='x', labelbottom=False, which='both', bottom=False, top=False)
    ax.tick_params(axis='y', labelsize=10)
    ax.set_ylabel('Events', fontsize=11)
    rax.axhline(1, color='gray', ls='--', lw=1)
    for (n_cv, n_var, col), o in zip(series, off):
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = (n_var * s_var) / (n_cv * s_cv)
            err = ratio * np.sqrt(1 / n_var + 1 / n_cv)
        ok = (n_cv >= 5) & (n_var > 0)  # skip bins with almost no CV events
        rax.errorbar(x[ok] + o, ratio[ok], yerr=err[ok], fmt='o', ms=3.5, color=col,
                     elinewidth=1, capsize=0)
    rax.set_ylim(0.5, 1.5)
    rax.set_xlim(edges[0], edges[-1])
    rax.set_xticks(x)
    rax.set_xticklabels(labels, fontsize=10)
    rax.tick_params(axis='x', which='minor', bottom=False, top=False)
    rax.tick_params(axis='y', labelsize=9)
    rax.set_ylabel('Var./CV', fontsize=10)


def combined_page(data, cv, var, setname, style):
    """One page per variation: every panel split by reco topology."""
    fig = plt.figure(figsize=(15, 10.2), layout='none')
    outer = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.28, left=0.06, right=0.98, top=0.83, bottom=0.11)
    a_cv, a_var = data[cv], data[var]
    for i in range(6):
        r, c = divmod(i, 3)
        inner = outer[r, c].subgridspec(2, 1, height_ratios=[3, 1.3], hspace=0.08)
        ax = fig.add_subplot(inner[0])
        rax = fig.add_subplot(inner[1], sharex=ax)
        if i < len(MULTS):
            b, title, lo, hi = MULTS[i]
            labels = [str(k) for k in range(lo, hi)] + [f'$\\geq${hi}']
            h_cv = [hist_mult(a_cv, b, lo, hi, type_mask(a_cv, codes)) for _, _, codes in TYPES]
            h_var = [hist_mult(a_var, b, lo, hi, type_mask(a_var, codes)) for _, _, codes in TYPES]
        else:
            title, labels = 'Truth category', [t for t, _ in TRUTH]
            h_cv = [hist_truth(a_cv, type_mask(a_cv, codes)) for _, _, codes in TYPES]
            h_var = [hist_truth(a_var, type_mask(a_var, codes)) for _, _, codes in TYPES]
        combined_panel(ax, rax, np.array(h_cv), np.array(h_var), a_cv['scale'], a_var['scale'], labels, title)
    # One legend for the page: colour = topology; line = CV, point = variation.
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=col, lw=1.5, marker='o', ms=5) for col in TYPE_COLORS]
    handles.append(Line2D([], [], color='black', lw=2, marker='o', ms=5))
    labs = [lab for _, lab, _ in TYPES] + ['Total']
    handles += [Line2D([], [], color='gray', lw=1.5), Line2D([], [], color='gray', ls='none', marker='o', ms=5)]
    labs += ['CV (line)', f'{var} (points)']
    fig.legend(handles, labs, loc='lower center', ncol=8, fontsize=11, bbox_to_anchor=(0.5, 0.005),
               handlelength=1.8, columnspacing=1.2)
    fig.suptitle(f'{SET_TITLE[setname]} Detector Variations: {var} vs CV, by Reco Topology',
                 y=0.975, fontsize=15)
    hdr = fig.add_axes([0.06, 0.875, 0.92, 0.0001])
    hdr.axis('off')
    stamp(hdr, pot=POT)
    return fig


def type_pct_plot(data, cv, setname, style):
    vars_ = [s for s in data if s != cv]
    y = np.arange(len(vars_))[::-1]
    fig, ax = plt.subplots(figsize=(9, 0.45 * len(vars_) + 2.2))
    off = np.linspace(-0.28, 0.28, len(TYPES))
    a_cv = data[cv]
    n_cv = np.array([np.sum(type_mask(a_cv, codes)) for _, _, codes in TYPES], float)
    for k, ((_, lab, codes), o, col) in enumerate(zip(TYPES, off, TYPE_COLORS)):
        pct, err = [], []
        for v in vars_:
            n_v = np.sum(type_mask(data[v], codes))
            r = (n_v * data[v]['scale']) / (n_cv[k] * a_cv['scale'])
            pct.append(100 * (r - 1))
            err.append(100 * r * np.sqrt(1 / n_v + 1 / n_cv[k]))
        ax.errorbar(pct, y + o, xerr=err, fmt=style.get_marker(k), ms=6, color=col,
                    elinewidth=1.5, capsize=0, label=lab)
    ax.axvline(0, color='gray', lw=1, ls='--', zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(vars_)
    ax.tick_params(axis='y', which='minor', left=False, right=False)
    ax.set_ylim(-0.7, len(vars_) - 0.3)
    ax.set_xlabel('Change from CV [%]')
    ax.grid(axis='x')
    ax.set_title(f'{SET_TITLE[setname]} Detector Variations by Reco Interaction Type', loc='center')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.09), ncol=5, fontsize=10,
              handletextpad=0.3, columnspacing=1.0)
    stamp(ax, pot=POT)
    return fig


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--csv', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--sets', nargs='+', default=['main', 'nue'])
    p.add_argument('--pages', choices=['all', 'combined'], default='all',
                   help='combined: only the one-page-per-variation topology plots')
    args = p.parse_args()
    style = load_style_from_toml('nue_icarus')
    plt.style.use(style.get_style())
    d = pd.read_csv(args.csv, keep_default_na=False)
    for s in args.sets:
        rows = d[d.set == s]
        cv = rows.cv_reference.iloc[0]
        data = load(s, rows)
        outdir = os.path.join(args.outdir, s)
        os.makedirs(outdir, exist_ok=True)
        variations = [x for x in data if x != cv]
        os.makedirs(os.path.join(outdir, 'combined'), exist_ok=True)
        with PdfPages(os.path.join(outdir, f'multiplicity_combined_{s}.pdf')) as pdf:
            for var in variations:
                fig = combined_page(data, cv, var, s, style)
                fig.savefig(os.path.join(outdir, 'combined', f'{var}.png'))
                pdf.savefig(fig)
                plt.close(fig)
        print(f'[INFO] {s}: {len(variations)} combined pages in {outdir}/combined', flush=True)
        if args.pages == 'combined':
            continue
        n = 0
        with PdfPages(os.path.join(outdir, f'multiplicity_{s}.pdf')) as pdf:
            fig = type_pct_plot(data, cv, s, style)
            fig.savefig(os.path.join(outdir, 'interaction_type_pct_vs_cv.png'))
            pdf.savefig(fig)
            plt.close(fig)
            n += 1
            for slug, label, codes in [('all', None, None)] + TYPES:
                os.makedirs(os.path.join(outdir, slug), exist_ok=True)
                for var in [x for x in data if x != cv]:
                    fig = page(data, cv, var, s, codes, label, style)
                    fig.savefig(os.path.join(outdir, slug, f'{var}.png'))
                    pdf.savefig(fig)
                    plt.close(fig)
                    n += 1
        print(f'[INFO] {s}: {n} pages in {outdir}', flush=True)


if __name__ == '__main__':
    main()
