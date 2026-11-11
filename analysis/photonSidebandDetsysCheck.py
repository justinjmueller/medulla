#!/usr/bin/env python3
"""
Photon-sideband version of the last cells of detectorSystematics.ipynb: compare
each detector-variation sample's photon_sideband rate with the detector-variation
CV (and, if given, the nominal NuMIFull MC), per POT, in reco_visible_energy.

Why: in v32 the main photon-sideband detector splines sit above 1 for nearly all
variations, i.e. the det-var CV's photon_sideband rate looks low. If the CV is low
against BOTH the variations and NuMIFull, the CV sample is the odd one out.

    /usr/bin/python3 analysis/photonSidebandDetsysCheck.py \
        --stage0 /exp/icarus/data/users/micarrig/merged_stage0_v32.root \
        [--final /exp/icarus/data/users/micarrig/merged_final_v32.root] \
        --outdir plots_v32/photon_sideband_detsys [--tree photon_sideband]

Reads only a few branches of small trees. POT for the det-var samples comes from
the CSV (corrected_pot, set main; recombination_residual -> recombination,
shower_scale_* -> cv); NuMIFull POT is the sum of its <tree>_exposure 'pot'.
Writes rates.csv, summary.txt and one PNG per variation.
"""
import argparse
import os

import numpy as np
import pandas as pd
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CSV = '/exp/icarus/app/users/micarrig/nuESpine/medulla/plots_v31/pot/detsys_corrected_pot.csv'
ALIAS = {'recombination_residual': 'recombination', 'shower_scale_up': 'cv', 'shower_scale_down': 'cv'}
VARS = ['induction_high', 'induction_low', 'coherent_high', 'coherent_low', 'intrinsic_high',
        'intrinsic_low', 'lifetime_high', 'lifetime_low', 'scintillation', 'wire_gap', 'yz_uniformity',
        'recombination_residual', 'cathode_bend', 'null', 'shower_scale_up', 'shower_scale_down']
BINS = np.array([0.0, 0.81, 1.1, 1.4, 5.0])
SCALE = 2.5e20


def hist(f, path, var):
    a = f[path].arrays([var, 'true_category'], library='np')
    x = np.asarray(a[var], dtype=float)
    x = x[np.isfinite(x)]
    n, _ = np.histogram(x, bins=BINS)
    return n, len(a[var])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--stage0', required=True)
    p.add_argument('--final', default=None, help='merged_final with events/NuMIFull/<tree> (optional)')
    p.add_argument('--outdir', required=True)
    p.add_argument('--tree', default='photon_sideband')
    p.add_argument('--var', default='reco_visible_energy')
    p.add_argument('--numi-pot', type=float, default=7.095834e20,
                   help='NuMIFull POT when --final has no <tree>_exposure (the value the plotting TOMLs use)')
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    pot = pd.read_csv(CSV, keep_default_na=False)
    pot = pot[pot['set'] == 'main'].set_index('sample')['corrected_pot'].astype(float)

    f0 = uproot.open(a.stage0)
    rows = []
    counts = {}
    for s in ['cv'] + VARS:
        n, tot = hist(f0, f'events/{s}/{a.tree}', a.var)
        counts[s] = (n, pot[ALIAS.get(s, s)])
        rows.append((s, tot, pot[ALIAS.get(s, s)], *n))

    if a.final:
        ff = uproot.open(a.final)
        n, tot = hist(ff, f'events/NuMIFull/{a.tree}', a.var)
        ek = f'events/NuMIFull/{a.tree}_exposure'
        if ek in ff:
            epot = float(np.sum(ff[ek].arrays(['pot'], library='np')['pot']))
        else:
            epot = a.numi_pot
            print(f'[INFO] no {ek}; NuMIFull POT = {epot:.6g} (--numi-pot; POT histogram says {ff["events/NuMIFull/POT"].values()[0]:.6g})')
        counts['NuMIFull'] = (n, epot)
        rows.append(('NuMIFull', tot, epot, *n))

    df = pd.DataFrame(rows, columns=['sample', 'entries', 'pot'] + [f'bin{i}' for i in range(len(BINS) - 1)])
    df['per1e18'] = df['entries'] / (df['pot'] / 1e18)
    df.to_csv(os.path.join(a.outdir, 'rates.csv'), index=False)

    ref = 'NuMIFull' if 'NuMIFull' in counts else None
    ncv, pcv = counts['cv']
    lines = [f'{a.tree} in {a.var}, rates per 1e18 POT (CSV corrected POT for det-var samples)',
             f"{'sample':<24}{'per1e18':>9}{'vs cv':>8}" + (f"{'vs NuMIFull':>13}" if ref else '')]
    rate = lambda s: counts[s][0].sum() / (counts[s][1] / 1e18)
    for s in (['NuMIFull'] if ref else []) + ['cv'] + VARS:
        r = rate(s)
        line = f'{s:<24}{r:9.3f}{r / rate("cv"):8.3f}'
        if ref:
            line += f'{r / rate(ref):13.3f}'
        lines.append(line)
    var_over_cv = np.array([rate(s) / rate('cv') for s in VARS])
    lines.append(f'mean(variation/cv) over {len(VARS)} variations: {var_over_cv.mean():.3f} '
                 f'(stat error on cv alone: {1 / np.sqrt(ncv.sum()):.3f})')
    open(os.path.join(a.outdir, 'summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))

    centers = 0.5 * (BINS[:-1] + BINS[1:])
    for s in VARS:
        nv, pv = counts[s]
        hv, hc = nv * SCALE / pv, ncv * SCALE / pcv
        fig, ax = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
        ax[0].stairs(hc, BINS, color='tab:blue', lw=2, label='det-var CV')
        ax[0].stairs(hv, BINS, color='tab:red', lw=2, label=s)
        if ref:
            nr, pr = counts[ref]
            ax[0].stairs(nr * SCALE / pr, BINS, color='black', lw=2, ls='--', label='NuMIFull (nominal)')
        with np.errstate(divide='ignore', invalid='ignore'):
            r = hv / hc
            e = r * np.sqrt(1 / np.maximum(nv, 1) + 1 / np.maximum(ncv, 1))
        ax[1].stairs(r, BINS, color='tab:red', lw=2, label=f'{s} / det-var CV')
        ax[1].errorbar(centers, r, yerr=e, fmt='none', ecolor='tab:red')
        if ref:
            ax[1].stairs((nr * SCALE / pr) / hc, BINS, color='black', lw=2, ls='--', label='NuMIFull / det-var CV')
        ax[1].axhline(1, color='gray', ls='--', lw=1)
        ax[1].set_ylim(0.7, 1.3)
        ax[0].set_ylabel(f'Events / {SCALE:.1e} POT')
        ax[1].set_xlabel('Reco visible energy [GeV]')
        ax[0].legend()
        ax[1].legend(fontsize=9)
        ax[0].set_title(f'{a.tree}: {s}')
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, f'{a.tree}_{s}.png'))
        plt.close(fig)
    print(f'wrote {a.outdir}')


if __name__ == '__main__':
    main()
