#!/usr/bin/env python3
"""Electron energy and cos(theta) resolution: CV vs each nue-enhanced detector variation.

Read-only on the ROOT files. Two modes (--mode, default both), each in its own subdirectory:

splines  1. For every variation v and each variable (reco leading-electron energy, reco cos(theta)
            to the beam): POT-normalized histogram of SIGNAL events (true_category == 0) in the
            phase-1 bins, divided by the same for the CV (NuMI_nue). This ratio is the phase-1
            spline value at z = +1 (systematics/src/detsys.cc: 1 + 1*(ratio - 1)). Empty CV bins
            and values outside the bins get weight 1, as in DetsysCalculator::get_weight.
         2. The +1 sigma weights are applied to the signal events of the v31 MC
            (events/NuMIFull/selected): the energy resolution uses the energy-binned weight, the
            angle resolution the cos(theta)-binned one. CV = NuMIFull with weight 1.
direct   The resolution of the signal events in the variation sample itself vs the CV sample
         (NuMI_nue), no splines.

POT: det-var samples use corrected_pot from plots_v31/pot/detsys_pot_and_counts.csv; NuMIFull uses
its POT histogram. Everything is scaled to 2.5e20 POT.
Resolutions: (E_reco - E_true) / E_true and dcos = cos_reco - cos_true, each fitted with a
double-sided crystal ball (Gaussian core, power-law tails on both sides).

Output in --outdir/<mode>/:
  <variation>/{energy_total,angle_total,energy_bins,angle_bins}.png
  fit_params_{energy,angle}.png      tables of the total-resolution fits (CV + every variation)
  fit_params_total.csv, fit_params_bins.csv (+ spline_ratios.csv for splines)
  detsys_resolution_<mode>.pdf       all pages

  spineplot/.pixi/envs/default/bin/python scripts/detsys_resolution.py --outdir plots_v31/detsys_resolution
"""
import argparse
import os
import sys
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import uproot
from scipy.optimize import curve_fit

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'spineplot'))
from style import load_style_from_toml, stamp  # noqa: E402

POT = 2.5e20
NUE_FILE = '/exp/icarus/data/users/micarrig/merged_nue_v31.root'
MC_FILE = '/exp/icarus/data/users/micarrig/merged_final_v31e_flux.root'
POT_CSV = os.path.join(REPO, 'plots_v31/pot/detsys_pot_and_counts.csv')
CV = 'NuMI_nue'
VARIATIONS = ['induction_high', 'induction_low', 'coherent_high', 'coherent_low', 'intrinsic_high',
              'intrinsic_low', 'lifetime_high', 'lifetime_low', 'scintillation', 'wire_gap',
              'yz_uniformity', 'recombination', 'cathode_bend', 'null']
# Phase-1 spline bins (systematics/toml/NuMI_nue_gundam.toml [variations]).
E_BINS = np.array([300, 390, 490, 600, 730, 870, 1030, 1240, 1500, 2000, 5000], dtype=float)
C_BINS = np.array([-1.0, -0.4, 0.0, 0.2, 0.35, 0.5, 0.625, 0.725, 0.8, 0.85, 0.9, 0.94, 0.97, 1.0])
# kind: (reco branch, true branch, spline bins, resolution label, range, nbins, title, true-bin label)
KINDS = {
    'energy': ('reco_leading_electron_energy', 'true_leading_electron_energy', E_BINS,
               r'$(E^{reco}_e - E^{true}_e)\,/\,E^{true}_e$', (-0.6, 0.3), 90,
               'Electron Energy Resolution', 'True $E_e$ [MeV]'),
    'angle': ('reco_ele_beam_open_angle', 'true_ele_beam_open_angle', C_BINS,
              r'$\cos\theta^{reco} - \cos\theta^{true}$', (-0.15, 0.06), 105,
              r'Electron $\cos\theta$ Resolution', r'True $\cos\theta$'),
}
BRANCHES = ['true_category'] + [b for k in KINDS.values() for b in k[:2]]
MODE_TITLE = {'splines': 'signal-only splines (+1$\\sigma$) on the v31 MC',
              'direct': 'detector-variation samples directly'}
CV_LABEL = {'splines': 'CV (NuMIFull)', 'direct': 'CV (NuMI_nue)'}
VAR_COLOR = '#C0392B'
PARS = ['mu', 'sigma', 'aL', 'nL', 'aR', 'nR']
# Minimum tail-start alpha per kind: the summed cos(theta) resolution has no Gaussian core, so a
# free alpha runs to ~0.2 and sigma becomes meaningless; alpha >= 1 keeps sigma a core width.
ALPHA_MIN = {'energy': 0.1, 'angle': 1.0}


# ---------------------------------------------------------------- inputs
def load_pot():
    d = pd.read_csv(POT_CSV, keep_default_na=False)
    return d[d.set == 'nue'].set_index('sample').corrected_pot.astype(float)


def load_signal(fname, sample):
    a = uproot.open(fname)[f'events/{sample}/selected'].arrays(BRANCHES, library='np')
    sig = a['true_category'] == 0
    a = {k: v[sig] for k, v in a.items()}
    a['res_energy'] = (a['reco_leading_electron_energy'] - a['true_leading_electron_energy']) \
        / a['true_leading_electron_energy']
    a['res_angle'] = a['reco_ele_beam_open_angle'] - a['true_ele_beam_open_angle']
    return a


def spline_ratios(nue, pot):
    rows, ratios = [], {}
    for v in VARIATIONS:
        for kind, (rb, _, bins, *_r) in KINDS.items():
            nv = np.histogram(nue[v][rb], bins=bins)[0].astype(float)
            nc = np.histogram(nue[CV][rb], bins=bins)[0].astype(float)
            with np.errstate(divide='ignore', invalid='ignore'):
                r = np.where(nc > 0, (nv / pot[v]) / (nc / pot[CV]), 1.0)
                err = np.where((nc > 0) & (nv > 0), r * np.sqrt(1 / np.maximum(nv, 1) + 1 / np.maximum(nc, 1)), 0.0)
            ratios[(v, kind)] = r
            for i in range(len(r)):
                rows.append(dict(variation=v, variable=kind, bin_lo=bins[i], bin_hi=bins[i + 1],
                                 n_var=int(nv[i]), n_cv=int(nc[i]), ratio=r[i], ratio_err=err[i]))
    return ratios, pd.DataFrame(rows)


def weights(x, bins, ratio):
    i = np.digitize(x, bins) - 1
    ok = (i >= 0) & (i < len(ratio)) & np.isfinite(x)
    w = np.ones_like(x, dtype=float)
    w[ok] = ratio[i[ok]]
    return w


# ---------------------------------------------------------------- double-sided crystal ball
def dscb(x, amp, mu, sigma, aL, nL, aR, nR):
    """Double-sided crystal ball, peak height amp (counts per bin)."""
    t = (x - mu) / sigma
    out = np.exp(-0.5 * t ** 2)
    lo, hi = t < -aL, t > aR
    if lo.any():
        A, B = (nL / aL) ** nL * np.exp(-0.5 * aL ** 2), nL / aL - aL
        out = np.where(lo, A * np.abs(B - t) ** (-nL), out)
    if hi.any():
        A, B = (nR / aR) ** nR * np.exp(-0.5 * aR ** 2), nR / aR - aR
        out = np.where(hi, A * np.abs(B + t) ** (-nR), out)
    return amp * out


def wquantile(x, w, q):
    o = np.argsort(x)
    c = np.cumsum(w[o])
    return np.interp(np.asarray(q) * c[-1], c, x[o])


def fit_dscb(x, w, rng, nbins, amin=0.1):
    """Binned chi2 fit of a double-sided crystal ball to weighted values."""
    edges = np.linspace(*rng, nbins + 1)
    centers = 0.5 * (edges[1:] + edges[:-1])
    inr = (x >= rng[0]) & (x < rng[1])
    h = np.histogram(x, bins=edges, weights=w)[0]
    e2 = np.histogram(x, bins=edges, weights=w ** 2)[0]
    out = dict(n_events=int(inr.sum()), ok=False)
    if inr.sum() > 0:
        q16, q50, q84 = wquantile(x[inr], w[inr], [0.16, 0.5, 0.84])
        out.update(median=q50, sigma68=0.5 * (q84 - q16))
    keep = e2 > 0
    if keep.sum() < 10 or inr.sum() < 30:
        return out
    med = np.median(x[inr])
    q16, q84 = np.percentile(x[inr], [16, 84])
    s0 = max(0.5 * (q84 - q16) * 0.7, 1e-3)
    p0 = [h.max(), med, s0, 1.0, 3.0, 1.0, 3.0]
    p0[3] = p0[5] = max(1.0, amin + 0.2)
    bounds = ([0, rng[0], 1e-4, amin, 1.01, amin, 1.01], [np.inf, rng[1], rng[1] - rng[0], 10, 150, 10, 150])
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            p, cov = curve_fit(dscb, centers[keep], h[keep], p0=p0, sigma=np.sqrt(e2[keep]),
                               absolute_sigma=True, bounds=bounds, maxfev=40000)
        perr = np.sqrt(np.abs(np.diag(cov)))
        chi2 = np.sum(((h[keep] - dscb(centers[keep], *p)) ** 2) / e2[keep])
        ndf = int(keep.sum() - len(p))
        out.update(ok=True, chi2=chi2, ndf=ndf, chi2_ndf=chi2 / max(ndf, 1), _p=p)
        for i, name in enumerate(PARS):
            out[name], out[name + '_err'] = p[i + 1], perr[i + 1]
    except (RuntimeError, ValueError) as err:
        out['error'] = str(err)[:80]
    return out


def fit_text(name, r, short=False):
    if not r.get('ok'):
        return f'{name}: fit failed'
    if short:
        return f"{name}: $\\mu$={r['mu']:+.4f}  $\\sigma$={r['sigma']:.4f}"
    return (f"{name}\n$\\mu$ = {r['mu']:+.4f} $\\pm$ {r['mu_err']:.4f}\n"
            f"$\\sigma$ = {r['sigma']:.4f} $\\pm$ {r['sigma_err']:.4f}\n"
            f"$\\alpha_L$ = {r['aL']:.2f} $\\pm$ {r['aL_err']:.2f},  $n_L$ = {r['nL']:.1f} $\\pm$ {r['nL_err']:.1f}\n"
            f"$\\alpha_R$ = {r['aR']:.2f} $\\pm$ {r['aR_err']:.2f},  $n_R$ = {r['nR']:.1f} $\\pm$ {r['nR_err']:.1f}\n"
            f"$\\chi^2$/ndf = {r['chi2']:.1f}/{r['ndf']}\n"
            f"$\\sigma_{{68}}$ (16-84%) = {r['sigma68']:.4f}")


def draw_points(ax, x, w, edges, color, marker, ms=4, label=None):
    """Weighted histogram as points with sqrt(sum w^2) error bars; returns the bin contents."""
    h = np.histogram(x, bins=edges, weights=w)[0]
    e = np.sqrt(np.histogram(x, bins=edges, weights=w ** 2)[0])
    c = 0.5 * (edges[1:] + edges[:-1])
    nz = h > 0
    ax.errorbar(c[nz], h[nz], yerr=e[nz], fmt=marker, color=color, ms=ms,
                mfc=color if marker == 'o' else 'white', elinewidth=0.8, capsize=0, label=label, zorder=3)
    return h


def draw_fit(ax, r, edges, **kw):
    if r.get('ok'):
        xs = np.linspace(edges[0], edges[-1], 800)
        ax.plot(xs, dscb(xs, *r['_p']), **kw)


# ---------------------------------------------------------------- plots
def header(fig, left, right, y):
    hdr = fig.add_axes([left, y, right - left, 0.0001])
    hdr.axis('off')
    stamp(hdr, pot=POT)


def page_total(var, kind, mode, cv, vv, r_cv, r_v):
    """cv / vv = (resolution values, true values, weights)."""
    _, _, _, xlabel, rng, nbins, title, _ = KINDS[kind]
    fig = plt.figure(figsize=(9, 8.6), layout='none')
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.05, left=0.12, right=0.97, top=0.85, bottom=0.08)
    ax, rax = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    edges = np.linspace(*rng, nbins + 1)
    hc = draw_points(ax, cv[0], cv[2], edges, 'black', 'o', label='CV')
    hv = draw_points(ax, vv[0], vv[2], edges, VAR_COLOR, 's', label=var)
    draw_fit(ax, r_cv, edges, color='black', lw=2.2, label='CV fit', zorder=4)
    draw_fit(ax, r_v, edges, color=VAR_COLOR, lw=2.2, label=f'{var} fit', zorder=4)
    ax.text(0.02, 0.97, fit_text('CV', r_cv), transform=ax.transAxes, va='top', fontsize=9)
    ax.text(0.02, 0.62, fit_text(var, r_v), transform=ax.transAxes, va='top', fontsize=9, color=VAR_COLOR)
    ax.set_xlim(*rng)
    ax.set_ylim(0, max(hc.max(), hv.max()) * 1.08)
    ax.set_ylabel(r'Signal events / 2.5$\times10^{20}$ POT')
    ax.tick_params(labelbottom=False)
    ax.legend(loc='upper right', fontsize=9.5)
    # Ratio as points, only where the CV bin has >= 10 raw events; the y range follows the
    # ratio values (min to max, padded) so every plotted point is visible. In the
    # splines mode CV and variation are the same events, so independent errors would be
    # wrong: no error bars there.
    nraw = np.histogram(cv[0], bins=edges)[0]
    ok = (nraw >= 10) & (hc > 0) & (hv > 0)
    ratio = hv[ok] / hc[ok]
    if mode == 'direct':
        ec = np.sqrt(np.histogram(cv[0], bins=edges, weights=cv[2] ** 2)[0])[ok]
        ev = np.sqrt(np.histogram(vv[0], bins=edges, weights=vv[2] ** 2)[0])[ok]
        rerr = ratio * np.sqrt((ev / hv[ok]) ** 2 + (ec / hc[ok]) ** 2)
    else:
        rerr = None
    centers = 0.5 * (edges[1:] + edges[:-1])
    rax.errorbar(centers[ok], ratio, yerr=rerr, fmt='s', color=VAR_COLOR, ms=3, mfc='white',
                 elinewidth=0.8, capsize=0)
    rax.axhline(1, color='black', ls='--', lw=1)
    if ok.any():
        lo, hi = min(ratio.min(), 1.0), max(ratio.max(), 1.0)
        pad = max(0.12 * (hi - lo), 0.01)
        rax.set_ylim(max(lo - pad, 0), hi + pad)
    rax.set_xlim(*rng)
    rax.set_xlabel(xlabel)
    rax.set_ylabel('Var./CV')
    fig.suptitle(f'{title}: {var} vs CV\n{MODE_TITLE[mode]}', y=0.985, fontsize=14)
    header(fig, 0.12, 0.97, 0.862)
    return fig


def page_bins(var, kind, mode, cv, vv, fits):
    _, _, bins, xlabel, rng, nbins, title, tlabel = KINDS[kind]
    nb = len(bins) - 1
    ncol = 5 if nb == 10 else 4
    nrow = int(np.ceil(nb / ncol))
    H = 3.4 * nrow + 2.6
    fig = plt.figure(figsize=(4.2 * ncol, H), layout='none')
    top = 1 - 2.5 / H
    gs = fig.add_gridspec(nrow, ncol, hspace=0.45, wspace=0.28, left=0.05, right=0.98, top=top, bottom=0.06)
    edges = np.linspace(*rng, max(nbins // 2, 25) + 1)
    for i in range(nb):
        ax = fig.add_subplot(gs[i // ncol, i % ncol])
        mc = (cv[1] >= bins[i]) & (cv[1] < bins[i + 1])
        mv = (vv[1] >= bins[i]) & (vv[1] < bins[i + 1])
        hc = draw_points(ax, cv[0][mc], cv[2][mc], edges, 'black', 'o', ms=3)
        hv = draw_points(ax, vv[0][mv], vv[2][mv], edges, VAR_COLOR, 's', ms=3)
        rc, rv = fits[('CV', i)], fits[(var, i)]
        draw_fit(ax, rc, edges, color='black', lw=1.8, zorder=4)
        draw_fit(ax, rv, edges, color=VAR_COLOR, lw=1.8, zorder=4)
        ax.text(0.03, 0.96, fit_text('CV', rc, short=True), transform=ax.transAxes, va='top', fontsize=7.5)
        ax.text(0.03, 0.86, fit_text('var', rv, short=True), transform=ax.transAxes, va='top', fontsize=7.5,
                color=VAR_COLOR)
        fmt = '{:.0f}' if kind == 'energy' else '{:.3g}'
        ax.set_title(f'{tlabel.split(" [")[0]}: {fmt.format(bins[i])}–{fmt.format(bins[i + 1])}', fontsize=9.5,
                     loc='center')
        ax.set_xlim(*rng)
        ax.set_ylim(0, max(hc.max(), hv.max(), 1e-9) * 1.45)
        ax.tick_params(labelsize=8)
        if i + ncol >= nb:
            ax.set_xlabel(xlabel, fontsize=9)
    fig.suptitle(f'{title} per {tlabel.split(" [")[0]} Bin: {var} (red squares) vs CV (black dots), lines = fits\n'
                 f'{MODE_TITLE[mode]}', y=0.99, fontsize=14)
    header(fig, 0.05, 0.98, 1 - 1.3 / H)
    return fig


def table_page(df, kind, mode):
    title = KINDS[kind][6]
    d = df[df.variable == kind].set_index('sample')
    cv = d.loc['CV']
    # Comparison columns first. For energy the resolution is fractional, so 100*dmu is an
    # energy-scale shift in %.
    dmu_lab = r'$\Delta\mu$ (scale)' if kind == 'energy' else r'$\Delta\mu$'
    hdr = ['Sample', dmu_lab, r'$\Delta\mu$ / $\sigma_{CV}$', r'$\Delta\sigma$ / $\sigma_{CV}$',
           r'$\Delta\sigma_{68}$ / $\sigma_{68,CV}$', r'$\mu$', 'median', r'$\sigma$', r'$\sigma_{68}$',
           r'$\alpha_L$', r'$n_L$', r'$\alpha_R$', r'$n_R$', r'$\chi^2$/ndf']
    rows = [hdr]
    for s in ['CV'] + VARIATIONS:
        r = d.loc[s]
        label = CV_LABEL[mode] if s == 'CV' else s
        if not bool(r.ok):
            rows.append([label, 'fit failed'] + [''] * (len(hdr) - 2))
            continue
        if s == 'CV':
            comp = ['', '', '', '']
        else:
            dmu, dsig = r.mu - cv.mu, r.sigma - cv.sigma
            comp = [f'{100 * dmu:+.2f}%' if kind == 'energy' else f'{dmu:+.4f}', f'{dmu / cv.sigma:+.2f}',
                    f'{100 * dsig / cv.sigma:+.1f}%', f'{100 * (r.sigma68 - cv.sigma68) / cv.sigma68:+.1f}%']
        rows.append([label] + comp + [f'{r.mu:+.4f} ± {r.mu_err:.4f}', f"{r['median']:+.4f}",
                                      f'{r.sigma:.4f} ± {r.sigma_err:.4f}', f'{r.sigma68:.4f}',
                                      f'{r.aL:.2f}', f'{r.nL:.1f}', f'{r.aR:.2f}', f'{r.nR:.1f}',
                                      f'{r.chi2:.0f}/{int(r.ndf)}'])
    fig, ax = plt.subplots(figsize=(21, 0.36 * len(rows) + 1.8))
    ax.axis('off')
    tab = ax.table(cellText=rows, cellLoc='center', edges='T', bbox=[0, 0, 1, 1],
                   colWidths=[0.1, 0.075, 0.07, 0.07, 0.08, 0.1, 0.06, 0.1, 0.055, 0.045, 0.045, 0.045, 0.045, 0.06])
    tab.auto_set_font_size(False)
    tab.set_fontsize(10.5)
    for i in range(1, len(rows)):
        for j in range(len(hdr)):
            tab[i, j].visible_edges = 'TB' if i == 1 else 'B' if i == len(rows) - 1 else 'open'
    for j in range(len(hdr)):
        tab[0, j].get_text().set_weight('bold')
    ax.set_title(f'{title}: Double-Sided Crystal-Ball Fits\n{MODE_TITLE[mode]}', loc='center')
    stamp(ax, pot=POT)
    return fig


# ---------------------------------------------------------------- per-mode driver
def run_mode(mode, outdir, nue, pot, mc, mc_pot):
    os.makedirs(outdir, exist_ok=True)
    # samples[name][kind] = (resolution, true value, weight scaled to 2.5e20)
    samples = {}
    if mode == 'splines':
        ratios, dr = spline_ratios(nue, pot)
        dr.to_csv(os.path.join(outdir, 'spline_ratios.csv'), index=False, float_format='%.5f')
        scale = POT / mc_pot
        one = np.full(len(mc['res_energy']), scale)
        samples['CV'] = {k: (mc[f'res_{k}'], mc[KINDS[k][1]], one) for k in KINDS}
        for v in VARIATIONS:
            samples[v] = {k: (mc[f'res_{k}'], mc[KINDS[k][1]],
                              scale * weights(mc[KINDS[k][0]], KINDS[k][2], ratios[(v, k)])) for k in KINDS}
    else:
        for s in [CV] + VARIATIONS:
            a = nue[s]
            w = np.full(len(a['res_energy']), POT / pot[s])
            samples['CV' if s == CV else s] = {k: (a[f'res_{k}'], a[KINDS[k][1]], w) for k in KINDS}

    tot_rows, bin_rows, tot, binfits = [], [], {}, {}
    cols = ['ok', 'n_events'] + [p for q in PARS for p in (q, q + '_err')] + ['chi2', 'ndf', 'chi2_ndf', 'median', 'sigma68']
    for k, (_, _, bins, _, rng, nbins, *_r) in KINDS.items():
        for s, d in samples.items():
            x, tv, w = d[k]
            r = fit_dscb(x, w, rng, nbins, ALPHA_MIN[k])
            tot[(s, k)] = r
            tot_rows.append(dict(variable=k, sample=s, sum_weights=w.sum(), **{q: r.get(q, np.nan) for q in cols}))
            if not r['ok'] or r['chi2_ndf'] > 5:
                print(f"[WARN] {mode} total {k} {s}: ok={r['ok']} chi2/ndf={r.get('chi2_ndf', np.nan):.2f} "
                      f"{r.get('error', '')}", flush=True)
            for i in range(len(bins) - 1):
                m = (tv >= bins[i]) & (tv < bins[i + 1])
                rb = fit_dscb(x[m], w[m], rng, max(nbins // 2, 25), ALPHA_MIN[k])
                binfits[(k, s, i)] = rb
                bin_rows.append(dict(variable=k, sample=s, bin_lo=bins[i], bin_hi=bins[i + 1],
                                     **{q: rb.get(q, np.nan) for q in cols}))
    dt, db = pd.DataFrame(tot_rows), pd.DataFrame(bin_rows)
    dt.to_csv(os.path.join(outdir, 'fit_params_total.csv'), index=False, float_format='%.6g')
    db.to_csv(os.path.join(outdir, 'fit_params_bins.csv'), index=False, float_format='%.6g')
    okb = db.ok.astype(bool)
    print(f'[INFO] {mode}: bin fits {len(db)}, failed {(~okb).sum()}, chi2/ndf > 5: {(db.chi2_ndf > 5).sum()}',
          flush=True)

    with PdfPages(os.path.join(outdir, f'detsys_resolution_{mode}.pdf')) as pdf:
        for k in KINDS:
            fig = table_page(dt, k, mode)
            fig.savefig(os.path.join(outdir, f'fit_params_{k}.png'))
            pdf.savefig(fig)
            plt.close(fig)
        for v in VARIATIONS:
            os.makedirs(os.path.join(outdir, v), exist_ok=True)
            for k in KINDS:
                fig = page_total(v, k, mode, samples['CV'][k], samples[v][k], tot[('CV', k)], tot[(v, k)])
                fig.savefig(os.path.join(outdir, v, f'{k}_total.png'))
                pdf.savefig(fig)
                plt.close(fig)
            for k in KINDS:
                fits = {(s, i): binfits[(k, s, i)] for s in ('CV', v) for i in range(len(KINDS[k][2]) - 1)}
                fig = page_bins(v, k, mode, samples['CV'][k], samples[v][k], fits)
                fig.savefig(os.path.join(outdir, v, f'{k}_bins.png'))
                pdf.savefig(fig)
                plt.close(fig)
    print(f'[INFO] {mode}: wrote {outdir}', flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--mode', nargs='+', default=['splines', 'direct'], choices=['splines', 'direct'])
    args = ap.parse_args()
    style = load_style_from_toml('nue_icarus')
    plt.style.use(style.get_style())

    pot = load_pot()
    nue = {s: load_signal(NUE_FILE, s) for s in [CV] + VARIATIONS}
    print('[INFO] nue samples (signal): ' + ', '.join(f'{s}={len(a["res_energy"])}' for s, a in nue.items()), flush=True)
    mc = load_signal(MC_FILE, 'NuMIFull')
    mc_pot = float(uproot.open(MC_FILE)['events/NuMIFull/POT'].values()[0])
    print(f'[INFO] NuMIFull signal: {len(mc["res_energy"])}, POT {mc_pot:.4e}', flush=True)
    for mode in args.mode:
        run_mode(mode, os.path.join(args.outdir, mode), nue, pot, mc, mc_pot)


if __name__ == '__main__':
    main()
