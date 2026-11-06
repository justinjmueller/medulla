#!/usr/bin/env python3
"""Cutflow tables for one detector variation compared with the CV.

Read-only. Uses the `all` tree and its per-cut flag branches, with the cut order and labels of the
CV selection-efficiency table (analysis/efficiency_analysis_allTree.toml). Signal = true_category 0,
Background = categories 1-6. Raw (unweighted) counts, like the CV table.

For each set (nue: NuMI_nue vs variation; main: cv vs variation), group and cut i:
  step efficiency        N(cuts 1..i) / N(cuts 1..i-1)
  cumulative efficiency  N(cuts 1..i) / N(total)
  differential           N(cut i alone) / N(total)          (CSV only)
  delta = variation - CV in percentage points, binomial errors added in quadrature
  N at cut i scaled to 2.5e20 POT (corrected_pot) and its % change
  purity after cut i (main set, signal table)

Output in --outdir/<variation>/: cutflow_<set>_{signal,background}.png, cutflow_<set>.csv,
cutflow_<variation>.pdf

  spineplot/.pixi/envs/default/bin/python scripts/detsys_cutflow.py --variation recombination \
      --outdir plots_v31/detsys_cutflow
"""
import argparse
import copy
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
POT_CSV = os.path.join(REPO, 'plots_v31/pot/detsys_pot_and_counts.csv')
SETS = {'nue': ('/exp/icarus/data/users/micarrig/merged_nue_v31.root', 'NuMI_nue', r'$\nu_e$-Enhanced'),
        'main': ('/exp/icarus/data/users/micarrig/merged_final_v31e_flux.root', 'cv', 'Nominal')}
# Cut order and labels from the CV table (analysis/efficiency_analysis_allTree.toml).
CUTS = [('reco_no_cut', 'Total Events'),
        ('reco_fiducial_cut', 'Fiducial Volume'),
        ('reco_flash_cut', 'Flash Cut'),
        ('reco_containment_cut', 'Interaction Containment Cut'),
        ('reco_has_particle_multiplicity', r'$\geq$1 Electron-like Particle'),
        ('reco_electron_softmax_cut', 'Leading Electron Electron-Softmax > 0.6'),
        ('reco_primary_softmax_cut', 'Leading Electron Primary-Softmax > 0.95'),
        ('reco_no_muons', 'No Muon-like Particles'),
        ('reco_vertex_distance_cut', 'Vertex Distance < 3.5 cm'),
        ('reco_leading_ele_energy_cut', 'Leading Electron Energy > 300 MeV')]
VARIATIONS = ['induction_high', 'induction_low', 'coherent_high', 'coherent_low', 'intrinsic_high',
              'intrinsic_low', 'lifetime_high', 'lifetime_low', 'scintillation', 'wire_gap',
              'yz_uniformity', 'recombination', 'cathode_bend', 'null']
GROUPS = [('signal', 'Signal', lambda c: c == 0), ('background', 'Background', lambda c: c > 0)]


_CACHE = {}


def counts(fname, sample):
    """Cached wrapper: the CV trees are read once per run, and callers get a fresh copy."""
    if (fname, sample) not in _CACHE:
        _CACHE[(fname, sample)] = _counts(fname, sample)
    return copy.deepcopy(_CACHE[(fname, sample)])


def _counts(fname, sample):
    """Per group: (total, cumulative N after each cut, N passing each cut alone); plus totals for purity."""
    a = uproot.open(fname)[f'events/{sample}/all'].arrays(['true_category'] + [c for c, _ in CUTS], library='np')
    cat = a['true_category']
    flags = [np.nan_to_num(a[c].astype(float)) > 0.5 for c, _ in CUTS]
    out = {}
    for key, _, sel in GROUPS + [('all', 'All', lambda c: np.isfinite(c) | ~np.isfinite(c))]:
        g = sel(cat)
        cum, alone, m = [], [], g.copy()
        for f in flags:
            m &= f
            cum.append(int(m.sum()))
            alone.append(int((g & f).sum()))
        out[key] = dict(total=int(g.sum()), cum=np.array(cum), alone=np.array(alone))
    return out


def eff(k, n):
    k, n = np.asarray(k, float), np.asarray(n, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        e = np.where(n > 0, k / n, np.nan)
        err = np.where(n > 0, np.sqrt(np.clip(e * (1 - e), 0, None) / n), np.nan)
    return e, err


def build(setname, variation, pot):
    fname, cv, _ = SETS[setname]
    res = {s: counts(fname, s) for s in (cv, variation)}
    rows = []
    for key, glabel, _ in GROUPS:
        for name, smp in (('cv', cv), ('var', variation)):
            d = res[smp][key]
            prev = np.concatenate([[d['total']], d['cum'][:-1]])
            d['step'], d['step_err'] = eff(d['cum'], prev)
            d['cumeff'], d['cumeff_err'] = eff(d['cum'], d['total'])
            d['diff'], d['diff_err'] = eff(d['alone'], d['total'])
            d['scaled'] = d['cum'] * POT / pot[smp]
            d['purity'], d['purity_err'] = eff(d['cum'], res[smp]['all']['cum']) if key == 'signal' else (None, None)
        c, v = res[cv][key], res[variation][key]
        for i, (branch, label) in enumerate(CUTS):
            r = dict(set=setname, group=key, cut=branch, label=label, n_cv=int(c['cum'][i]), n_var=int(v['cum'][i]))
            for q in ('step', 'cumeff', 'diff'):
                r[f'{q}_cv'], r[f'{q}_var'] = 100 * c[q][i], 100 * v[q][i]
                r[f'{q}_delta'] = 100 * (v[q][i] - c[q][i])
                r[f'{q}_delta_err'] = 100 * np.hypot(c[q + '_err'][i], v[q + '_err'][i])
            r['scaled_cv'], r['scaled_var'] = c['scaled'][i], v['scaled'][i]
            r['scaled_pct'] = 100 * (v['scaled'][i] / c['scaled'][i] - 1) if c['scaled'][i] > 0 else np.nan
            if key == 'signal':
                r['purity_cv'], r['purity_var'] = 100 * c['purity'][i], 100 * v['purity'][i]
                r['purity_delta'] = 100 * (v['purity'][i] - c['purity'][i])
                r['purity_delta_err'] = 100 * np.hypot(c['purity_err'][i], v['purity_err'][i])
            rows.append(r)
    return pd.DataFrame(rows)


def table_page(df, setname, key, variation, with_purity):
    glabel = dict((k, l) for k, l, _ in GROUPS)[key]
    d = df[df.group == key]
    # Background cumulative efficiencies are tiny (the total includes every cosmic), so they are
    # printed with significant figures and their change as a relative %.
    rel = key == 'background'
    hdr = ['Cut', 'Step eff.\nCV [%]', f'Step eff.\n{variation} [%]', '$\\Delta$ step\n[pp]',
           'Cum. eff.\nCV [%]', f'Cum. eff.\n{variation} [%]', '$\\Delta$ cum.\n' + ('[rel. %]' if rel else '[pp]'),
           'N CV\n(2.5e20)', f'N {variation}\n(2.5e20)', '$\\Delta$N\n[%]']
    if with_purity:
        hdr += ['Purity\nCV [%]', f'Purity\n{variation} [%]', '$\\Delta$ purity\n[pp]']
    rows = [hdr]
    cfmt = (lambda x: f'{x:.3g}') if rel else (lambda x: f'{x:.2f}')
    for _, r in d.iterrows():
        first = r.cut == CUTS[0][0]
        if first:
            dcum = ''
        elif rel:
            rc = r.cumeff_delta / r.cumeff_cv * 100
            rerr = r.cumeff_delta_err / r.cumeff_cv * 100
            dcum = f'{rc:+.1f} ± {rerr:.1f}'
        else:
            dcum = f'{r.cumeff_delta:+.2f} ± {r.cumeff_delta_err:.2f}'
        row = [r.label,
               f'{r.step_cv:.2f}', f'{r.step_var:.2f}',
               '' if first else f'{r.step_delta:+.2f} ± {r.step_delta_err:.2f}',
               cfmt(r.cumeff_cv), cfmt(r.cumeff_var), dcum,
               f'{r.scaled_cv:,.0f}', f'{r.scaled_var:,.0f}', f'{r.scaled_pct:+.1f}']
        if with_purity:
            row += [f'{r.purity_cv:.2f}', f'{r.purity_var:.2f}', f'{r.purity_delta:+.2f} ± {r.purity_delta_err:.2f}']
        rows.append(row)
    ncol = len(hdr)
    first_w = 0.2 if with_purity else 0.22
    fig, ax = plt.subplots(figsize=(22 if with_purity else 18, 0.42 * len(rows) + 2.0))
    ax.axis('off')
    tab = ax.table(cellText=rows, cellLoc='center', edges='T', bbox=[0, 0, 1, 1],
                   colWidths=[first_w] + [(1 - first_w) / (ncol - 1)] * (ncol - 1))
    tab.auto_set_font_size(False)
    tab.set_fontsize(10.5)
    for i in range(2, len(rows)):
        for j in range(ncol):
            tab[i, j].visible_edges = 'B' if i == len(rows) - 1 else 'open'
    for j in range(ncol):
        tab[0, j].get_text().set_weight('bold')
    ax.set_title(f'{SETS[setname][2]} {glabel} Cutflow: {variation} vs CV', loc='center')
    stamp(ax, pot=POT)
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--variation', nargs='+', required=True, help="variation name(s), or 'all' for all 14")
    ap.add_argument('--sets', nargs='+', default=['nue', 'main'], choices=list(SETS))
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--check-numifull', action='store_true', help='print the NuMIFull signal cutflow (sanity check)')
    args = ap.parse_args()
    style = load_style_from_toml('nue_icarus')
    plt.style.use(style.get_style())
    variations = VARIATIONS if args.variation == ['all'] else args.variation

    p = pd.read_csv(POT_CSV, keep_default_na=False)
    if args.check_numifull:
        c = counts(SETS['main'][0], 'NuMIFull')['signal']
        print('[CHECK] NuMIFull signal cumulative:', ', '.join(f'{n}/{c["total"]}' for n in c['cum']), flush=True)

    for variation in variations:
        out = os.path.join(args.outdir, variation)
        os.makedirs(out, exist_ok=True)
        with PdfPages(os.path.join(out, f'cutflow_{variation}.pdf')) as pdf:
            for s in args.sets:
                pot = p[p.set == s].set_index('sample').corrected_pot.astype(float)
                df = build(s, variation, pot)
                df.to_csv(os.path.join(out, f'cutflow_{s}.csv'), index=False, float_format='%.5g')
                last = df.groupby('group').tail(1)
                print(f'[INFO] {variation} {s}: final ' + '; '.join(
                    f'{r.group} cv={r.n_cv} var={r.n_var}' for _, r in last.iterrows()), flush=True)
                for key, _, _ in GROUPS:
                    fig = table_page(df, s, key, variation, with_purity=(key == 'signal' and s == 'main'))
                    fig.savefig(os.path.join(out, f'cutflow_{s}_{key}.png'))
                    pdf.savefig(fig)
                    plt.close(fig)
        print(f'[INFO] wrote {out}', flush=True)


if __name__ == '__main__':
    main()
