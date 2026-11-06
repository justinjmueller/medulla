#!/usr/bin/env python3
"""Kinematics of each background in the signal region vs the NC pi0 sideband (NuMIFull nominal MC).

Read-only. For each background true_category, one figure comparing events/NuMIFull/selected
(signal region) with events/NuMIFull/nc_sideband, area normalized:
  true neutrino energy, true hadronic invariant mass W, true Q^2, reco and true visible energy,
  and the true interaction mode (QE / RES / DIS / COH / MEC / other).
Raw selected-event counts (no POT scaling, no flux weights).

Output in --outdir:
  kinematics_cat<c>.png   one figure per category
  kinematics_summary.csv  per category and region: n, medians, interaction-mode fractions

  spineplot/.pixi/envs/default/bin/python scripts/compare_region_kinematics.py \
      --outdir plots_v31/background_topologies_sideband
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

FILE = '/exp/icarus/data/users/micarrig/merged_final_v31e_flux.root'
REGIONS = [('selected', r'Signal region'), ('nc_sideband', r'NC $\pi^0$ sideband')]
CATS = [(2, r'$\nu_\mu$ CC (FV)'), (3, 'Other CC (OOFV / uncontained)'), (4, 'NC (FV)')]
# (branch, axis label, bin edges)
VARS = [('true_neutrino_energy', r'True $E_\nu$ [GeV]', np.linspace(0, 8, 33)),
        ('true_hadronic_invariant_mass', r'True hadronic invariant mass $W$ [GeV]', np.linspace(0.9, 3.5, 27)),
        ('true_momentum_transfer', r'True $Q^2$ [GeV$^2$]', np.linspace(0, 4, 33)),
        ('reco_visible_energy', 'Reco visible energy [GeV]', np.linspace(0, 5, 26)),
        ('true_visible_energy', 'True visible energy [GeV]', np.linspace(0, 5, 26))]
# caf::genie interaction mode
MODES = [(0, 'QE'), (1, 'RES'), (2, 'DIS'), (3, 'COH'), (10, 'MEC')]


def norm_hist(ax, v, edges, label, **kw):
    v = v[np.isfinite(v)]
    v = np.clip(v, edges[0], edges[-1] - 1e-9)
    ax.hist(v, bins=edges, weights=np.full(len(v), 1 / max(len(v), 1)), label=f'{label} ({len(v):,})', **kw)


def mode_fractions(mode):
    f = [100 * np.mean(mode == k) for k, _ in MODES]
    return f + [100 - sum(f)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    style = load_style_from_toml('nue_icarus')
    plt.style.use(style.get_style())

    f = uproot.open(FILE)
    br = ['true_category', 'true_interaction_mode'] + [v for v, _, _ in VARS]
    data = {t: f[f'events/NuMIFull/{t}'].arrays(br, library='np') for t, _ in REGIONS}
    kws = [dict(histtype='stepfilled', color=style.get_color(0), alpha=0.6),
           dict(histtype='step', color='black', lw=2)]
    mode_labels = [m for _, m in MODES] + ['Other']

    recs = []
    for c, name in CATS:
        fig, axes = plt.subplots(2, 3, figsize=(17, 10.5), layout='none')
        fig.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.07, wspace=0.25, hspace=0.3)
        for (t, rlab), kw in zip(REGIONS, kws):
            a = data[t]
            sel = a['true_category'] == c
            rec = dict(category=c, name=name, region=t, n=int(sel.sum()))
            for ax, (v, xlab, edges) in zip(axes.flat, VARS):
                norm_hist(ax, a[v][sel], edges, rlab, **kw)
                ax.set_xlabel(xlab)
                ax.set_ylabel('Fraction of events')
                x = a[v][sel]
                rec[f'median_{v}'] = np.median(x[np.isfinite(x)]) if np.isfinite(x).any() else np.nan
            for lab, fr in zip(mode_labels, mode_fractions(a['true_interaction_mode'][sel])):
                rec[f'mode_{lab}_pct'] = fr
            recs.append(rec)

        # interaction mode, side-by-side bars
        ax = axes.flat[len(VARS)]
        x = np.arange(len(mode_labels))
        for i, ((t, rlab), kw) in enumerate(zip(REGIONS, kws)):
            r = recs[-2 + i]
            fr = [r[f'mode_{lab}_pct'] for lab in mode_labels]
            ax.bar(x + (i - 0.5) * 0.4, fr, width=0.4, label=f"{rlab} ({r['n']:,})",
                   color=kw['color'], alpha=kw.get('alpha', 1) if i == 0 else 1,
                   fill=(i == 0), edgecolor=kw['color'], lw=kw.get('lw', 0))
        ax.set_xticks(x, mode_labels)
        ax.set_ylabel('Fraction of events [%]')
        ax.set_xlabel('True interaction mode')
        for a_ in axes.flat:
            a_.legend(fontsize=10)

        fig.suptitle(f'{name}: Signal Region vs NC $\\pi^0$ Sideband', y=0.975)
        hdr = fig.add_axes([0.06, 0.905, 0.92, 0.0001])
        hdr.axis('off')
        stamp(hdr)
        out = os.path.join(args.outdir, f'kinematics_cat{c}.png')
        fig.savefig(out)
        plt.close(fig)
        print('[INFO] wrote', out, flush=True)

    df = pd.DataFrame(recs)
    df.to_csv(os.path.join(args.outdir, 'kinematics_summary.csv'), index=False, float_format='%.3f')
    with pd.option_context('display.width', 250, 'display.max_columns', 30):
        print(df.drop(columns='name').round(2).to_string(index=False))


if __name__ == '__main__':
    main()
