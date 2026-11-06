#!/usr/bin/env python3
"""Emulate an inverted-primary-score sideband from events/NuMIFull/all (NuMIFull nominal MC).

Read-only. The proposed sideband is the signal selection with the leading-electron energy cut,
the muon veto and the electron softmax cut removed, and the primary softmax cut inverted:
  fiducial, flash, containment, vertex distance, >= 1 reco electron,
  leading-electron primary softmax < @primary_softmax.
It is orthogonal to the signal region by construction (same leading electron, opposite side of
the same threshold). Variants that put back the energy cut and/or the muon veto are also made.
The signal region (SR) is rebuilt from the same flags and its count is checked against
events/NuMIFull/selected.

For each region it reports the true_category composition, the true identity of the reco leading
"electron" (true_pdg_from_reco_ele), and per-category medians and interaction modes.
Raw counts (no POT scaling, no flux weights). The tree is read in chunks.

Output in --outdir:
  sideband_composition.csv  per region and category: n, fraction
  sideband_leading_pdg.csv  per region and category: true particle behind the reco leading electron
  sideband_kinematics.csv   per region and category: medians and interaction-mode fractions
  sideband_events.npz       per-event arrays of every event in the loosest region (for plotting)

  spineplot/.pixi/envs/default/bin/python scripts/inverted_primary_sideband.py \
      --outdir plots_v31/inverted_primary_sideband
"""
import argparse
import os

import numpy as np
import pandas as pd
import uproot

FILE = '/exp/icarus/data/users/micarrig/merged_final_v31e_flux.root'
TREE = 'events/NuMIFull/all'
PRIMARY_SOFTMAX = 0.95
ELECTRON_SOFTMAX = 0.6
ELECTRON_THRESHOLD = 300.0  # MeV
CATS = {0: r'Signal (nue CC, FV)', 1: 'Other nue CC', 2: 'numu CC (FV)', 3: 'Other CC (OOFV / uncontained)',
        4: 'NC (FV)', 5: 'Cosmic', 6: 'Other', -1: 'No category'}
PDG = {11: 'e', 13: 'mu', 22: 'gamma', 211: 'pi', 2212: 'p'}
MODES = [(0, 'QE'), (1, 'RES'), (2, 'DIS'), (3, 'COH'), (10, 'MEC')]
FLAGS = ['reco_fiducial_cut', 'reco_flash_cut', 'reco_containment_cut', 'reco_vertex_distance_cut',
         'reco_has_particle_multiplicity', 'reco_no_muons']
KEEP = ['true_category', 'true_interaction_mode', 'true_pdg_from_reco_ele', 'true_neutrino_energy',
        'true_hadronic_invariant_mass', 'true_momentum_transfer', 'reco_visible_energy',
        'reco_leading_electron_energy', 'reco_leading_electron_electron_softmax',
        'reco_leading_electron_primary_softmax', 'reco_muon_multiplicity', 'reco_photon_multiplicity',
        'reco_no_muons']
MEDIANS = ['true_neutrino_energy', 'true_hadronic_invariant_mass', 'true_momentum_transfer',
           'reco_visible_energy', 'reco_leading_electron_energy']


def regions(d):
    """Boolean masks for each region, from one dict of per-event arrays."""
    ps, es, e = (d['reco_leading_electron_primary_softmax'], d['reco_leading_electron_electron_softmax'],
                 d['reco_leading_electron_energy'])
    with np.errstate(invalid='ignore'):
        inv = np.isfinite(ps) & (ps < PRIMARY_SOFTMAX)
        hi_e = np.isfinite(e) & (e > ELECTRON_THRESHOLD)
        nomu = d['reco_no_muons'] == 1
        sr = (d['base'] & nomu & np.isfinite(es) & (es > ELECTRON_SOFTMAX)
              & np.isfinite(ps) & (ps > PRIMARY_SOFTMAX) & hi_e)
    sb = d['base'] & inv
    return {'SR': sr, 'SB proposed': sb, 'SB + E>300': sb & hi_e, 'SB + mu veto': sb & nomu,
            'SB + E>300 + mu veto': sb & hi_e & nomu}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--step', default='2000000', help='entries per chunk')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    f = uproot.open(FILE)
    n_sel = f['events/NuMIFull/selected'].num_entries
    parts = []
    for i, a in enumerate(uproot.iterate(f'{FILE}:{TREE}', FLAGS + KEEP, step_size=int(args.step), library='np')):
        base = np.ones(len(a['true_category']), bool)
        for k in FLAGS[:-1]:
            base &= a[k] == 1
        a['base'] = base
        m = regions(a)
        loose = m['SR'] | m['SB proposed']
        parts.append({k: a[k][loose] for k in KEEP + ['base']})
        print(f'[INFO] chunk {i}: {len(base):,} entries, kept {loose.sum():,}', flush=True)
    d = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    np.savez_compressed(os.path.join(args.outdir, 'sideband_events.npz'), **d)

    cat = np.where(np.isfinite(d['true_category']), d['true_category'], -1).astype(int)
    pdg = np.abs(np.nan_to_num(d['true_pdg_from_reco_ele'], nan=0)).astype(int)
    comp, lead, kin = [], [], []
    for rname, r in regions(d).items():
        n = int(r.sum())
        print(f'\n[REGION] {rname}: {n:,} events')
        for c in sorted(set(cat[r])):
            s = r & (cat == c)
            ns = int(s.sum())
            comp.append(dict(region=rname, category=c, name=CATS[c], n=ns, pct=100 * ns / n))
            row = dict(region=rname, category=c, name=CATS[c], n=ns)
            for p, lab in PDG.items():
                row[f'lead_{lab}_pct'] = 100 * np.mean(pdg[s] == p)
            row['lead_other_pct'] = 100 - sum(row[f'lead_{lab}_pct'] for lab in PDG.values())
            lead.append(row)
            k = dict(region=rname, category=c, name=CATS[c], n=ns)
            for v in MEDIANS:
                x = d[v][s]
                x = x[np.isfinite(x)]
                k[f'median_{v}'] = np.median(x) if len(x) else np.nan
            for mc, ml in MODES:
                k[f'mode_{ml}_pct'] = 100 * np.mean(d['true_interaction_mode'][s] == mc)
            kin.append(k)
            print(f"  {c:2d} {CATS[c]:32s} {ns:8,d} ({100 * ns / n:5.1f}%)  lead e/mu/gamma/pi/p/other: "
                  + '/'.join(f"{row[f'lead_{lab}_pct']:.0f}" for lab in list(PDG.values()) + ['other']) + '%')
    n_sr = int(regions(d)['SR'].sum())
    print(f'\n[CHECK] SR rebuilt from flags: {n_sr:,}   events/NuMIFull/selected: {n_sel:,}')
    pd.DataFrame(comp).to_csv(os.path.join(args.outdir, 'sideband_composition.csv'), index=False, float_format='%.3f')
    pd.DataFrame(lead).to_csv(os.path.join(args.outdir, 'sideband_leading_pdg.csv'), index=False, float_format='%.3f')
    pd.DataFrame(kin).to_csv(os.path.join(args.outdir, 'sideband_kinematics.csv'), index=False, float_format='%.3f')
    print(f'[INFO] wrote CSVs and sideband_events.npz to {args.outdir}')


if __name__ == '__main__':
    main()
