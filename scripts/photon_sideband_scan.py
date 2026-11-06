#!/usr/bin/env python3
"""Scan candidate sidebands that mimic the signal-region background (NuMIFull nominal MC).

Read-only. The selected-sample background is mostly an energetic primary photon reconstructed
as the electron. The candidates keep the signal-region preselection (fiducial, flash,
containment, no muons) and the energetic primary shower, but require it to be identified
differently:
  photon        no reco electron; leading photon E > 300 MeV and primary softmax > 0.95
  photon 1g     same, exactly one reco photon (>= @photon_threshold)
  photon >=2g   same, at least two reco photons
  low e-score   signal region with the electron softmax cut inverted (< 0.6)
All are orthogonal to the signal region (SR). For each region, and for the SR background
(true_category 1-4) as the reference, it reports the composition, the true particle behind the
leading shower, and the kinematics of the neutrino background (true_category 1-4); plus KS
distances to the SR background in true W, true E_nu and reco visible energy.

Output in --outdir:
  scan_summary.csv     one row per region
  scan_categories.csv  per region and category: n, %
  scan_events.npz      per-event arrays of the loose preselection (for further scans)

  spineplot/.pixi/envs/default/bin/python scripts/photon_sideband_scan.py \
      --outdir plots_v31/photon_sideband
  (add --from-npz to rescan the saved events without reading the tree)
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
SHOWER_THRESHOLD = 300.0  # MeV, the signal-region electron threshold
CATS = {0: 'Signal', 1: 'Other nue CC', 2: 'numu CC (FV)', 3: 'Other CC', 4: 'NC (FV)', 5: 'Cosmic',
        6: 'Other', -1: 'No category'}
BKG = [1, 2, 3, 4]
MODES = [(0, 'QE'), (1, 'RES'), (2, 'DIS'), (3, 'COH'), (10, 'MEC')]
FLAGS = ['reco_fiducial_cut', 'reco_flash_cut', 'reco_containment_cut', 'reco_no_muons']
KEEP = ['true_category', 'true_interaction_mode', 'true_neutrino_energy', 'true_hadronic_invariant_mass',
        'true_momentum_transfer', 'true_visible_energy', 'reco_visible_energy', 'reco_vertex_distance_cut', 'reco_has_particle_multiplicity',
        'reco_electron_multiplicity', 'reco_photon_multiplicity',
        'reco_leading_electron_energy', 'reco_leading_electron_electron_softmax',
        'reco_leading_electron_primary_softmax', 'reco_leading_photon_energy',
        'reco_leading_photon_primary_softmax', 'reco_leading_photon_photon_softmax',
        'true_pdg_from_reco_ele', 'true_pdg_from_reco_photon']


def gt(x, t):
    return np.isfinite(x) & (x > t)


def regions(d):
    ele = (d['reco_vertex_distance_cut'] == 1) & (d['reco_has_particle_multiplicity'] == 1) \
        & gt(d['reco_leading_electron_primary_softmax'], PRIMARY_SOFTMAX) \
        & gt(d['reco_leading_electron_energy'], SHOWER_THRESHOLD)
    es = d['reco_leading_electron_electron_softmax']
    pho = (d['reco_electron_multiplicity'] == 0) & gt(d['reco_leading_photon_energy'], SHOWER_THRESHOLD) \
        & gt(d['reco_leading_photon_primary_softmax'], PRIMARY_SOFTMAX)
    with np.errstate(invalid='ignore'):
        return {'SR': ele & gt(es, ELECTRON_SOFTMAX),
                'low e-score': ele & np.isfinite(es) & (es < ELECTRON_SOFTMAX),
                'photon': pho,
                'photon 1g': pho & (d['reco_photon_multiplicity'] == 1),
                'photon >=2g': pho & (d['reco_photon_multiplicity'] >= 2)}


def ks(a, b):
    a, b = np.sort(a[np.isfinite(a)]), np.sort(b[np.isfinite(b)])
    if not len(a) or not len(b):
        return np.nan
    x = np.concatenate([a, b])
    return np.max(np.abs(np.searchsorted(a, x, 'right') / len(a) - np.searchsorted(b, x, 'right') / len(b)))


def load(outdir):
    parts = []
    for i, a in enumerate(uproot.iterate(f'{FILE}:{TREE}', FLAGS + KEEP, step_size=2000000, library='np')):
        base = np.ones(len(a['true_category']), bool)
        for k in FLAGS:
            base &= a[k] == 1
        loose = base & (gt(a['reco_leading_electron_energy'], 100) | gt(a['reco_leading_photon_energy'], 100))
        parts.append({k: a[k][loose] for k in KEEP})
        print(f'[INFO] chunk {i}: kept {loose.sum():,} of {len(base):,}', flush=True)
    d = {k: np.concatenate([p[k] for p in parts]) for k in KEEP}
    np.savez_compressed(os.path.join(outdir, 'scan_events.npz'), **d)
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--from-npz', action='store_true')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    if args.from_npz:
        d = dict(np.load(os.path.join(args.outdir, 'scan_events.npz')))
    else:
        d = load(args.outdir)

    cat = np.where(np.isfinite(d['true_category']), d['true_category'], -1).astype(int)
    mode, W, Enu, Evis = (d['true_interaction_mode'], d['true_hadronic_invariant_mass'],
                          d['true_neutrino_energy'], d['reco_visible_energy'])
    regs = regions(d)
    ref = regs['SR'] & np.isin(cat, BKG)
    summ, cats = [], []
    for name, r in regs.items():
        n = int(r.sum())
        b = r & np.isin(cat, BKG)
        # the leading shower: the electron in the electron regions, the photon in the photon ones
        pdg = d['true_pdg_from_reco_photon'] if name.startswith('photon') else d['true_pdg_from_reco_ele']
        pdg = np.abs(np.nan_to_num(pdg[b], nan=0))
        row = dict(region=name, n=n, signal_pct=100 * np.mean(cat[r] == 0), n_bkg=int(b.sum()),
                   cosmic_nocat_pct=100 * np.mean(np.isin(cat[r], [-1, 5])),
                   shower_true_gamma_pct=100 * np.mean(pdg == 22), shower_true_e_pct=100 * np.mean(pdg == 11),
                   median_W=np.nanmedian(W[b]), median_Enu=np.nanmedian(Enu[b]), median_Evis=np.nanmedian(Evis[b]),
                   ks_W=ks(W[b], W[ref]), ks_Enu=ks(Enu[b], Enu[ref]), ks_Evis=ks(Evis[b], Evis[ref]))
        for mc, ml in MODES:
            row[f'{ml}_pct'] = 100 * np.mean(mode[b] == mc)
        for c in BKG:
            row[f'cat{c}_of_bkg_pct'] = 100 * np.mean(cat[b] == c)
        summ.append(row)
        for c in sorted(set(cat[r])):
            cats.append(dict(region=name, category=c, name=CATS[c], n=int(np.sum(r & (cat == c))),
                             pct=100 * np.mean(cat[r] == c)))
    s = pd.DataFrame(summ)
    s.to_csv(os.path.join(args.outdir, 'scan_summary.csv'), index=False, float_format='%.3f')
    pd.DataFrame(cats).to_csv(os.path.join(args.outdir, 'scan_categories.csv'), index=False, float_format='%.3f')
    with pd.option_context('display.width', 250, 'display.max_columns', 40):
        print(s.round(2).T.to_string())


if __name__ == '__main__':
    main()
