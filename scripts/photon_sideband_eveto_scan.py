#!/usr/bin/env python3
"""Does the photon sideband gain from allowing soft primary reco electrons? (NuMIFull nominal MC)

Read-only. The current photon_sideband requires no primary reco electron at all
(no_electrons, threshold 0) plus a leading photon > 400 MeV with primary softmax > 0.95, after
fiducial, flash, containment and no_muons. This scans the electron-veto threshold T: the event
passes if it has no primary reco electron with KE >= T. Any T <= 300 MeV stays orthogonal to the
signal region, which needs a leading electron > 300 MeV (checked below as overlap_SR).

The flat `all` tree only has the leading electron (by KE, primary or not), so for T > 0:
  sure pass   no primary electron at all, or the leading electron is below T
  sure fail   the leading electron is >= T and primary (primary softmax > 0.5)
  ambiguous   the leading electron is >= T but not primary, and a primary electron exists
              somewhere below it with unknown KE. Not counted in the region; reported as n_ambig.
Regions:
  SR                      signal region (reference; its background is true_category 1-4)
  veto T                  photon sideband with electron veto threshold T (T = 0 is the current one)
  added T                 only the events veto T adds over veto 0 (they have a soft primary electron)

Output in --outdir:
  eveto_summary.csv   one row per region: n, signal %, background makeup, modes, W, KS to SR background
  eveto_bins.csv      per region and reco_visible_energy fit bin: n, signal %, DIS % of background
  eveto_events.npz    the loose preselection (rerun with --from-npz to skip reading the tree)

  spineplot/.pixi/envs/default/bin/python scripts/photon_sideband_eveto_scan.py \\
      --outdir plots_v31/photon_sideband_eveto
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
SR_THRESHOLD = 300.0      # MeV, the signal-region electron threshold
PHOTON_SB = 400.0         # MeV, @photon_sb_threshold
VETOES = [0.0, 50.0, 100.0, 200.0, 300.0]
EVIS_BINS = [0.0, 0.81, 1.1, 1.4, 5.0]  # GeV, the photon sideband fit binning
BKG = [1, 2, 3, 4]
MODES = [(0, 'QE'), (1, 'RES'), (2, 'DIS'), (3, 'COH'), (10, 'MEC')]
FLAGS = ['reco_fiducial_cut', 'reco_flash_cut', 'reco_containment_cut', 'reco_no_muons']
KEEP = ['true_category', 'true_interaction_mode', 'true_neutrino_energy', 'true_hadronic_invariant_mass',
        'reco_visible_energy', 'reco_vertex_distance_cut', 'reco_has_particle_multiplicity',
        'reco_electron_multiplicity', 'reco_photon_multiplicity',
        'reco_leading_electron_energy', 'reco_leading_electron_electron_softmax',
        'reco_leading_electron_primary_softmax', 'reco_leading_photon_energy',
        'reco_leading_photon_primary_softmax', 'true_pdg_from_reco_ele', 'true_pdg_from_reco_photon']


def gt(x, t):
    return np.isfinite(x) & (x > t)


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
    np.savez_compressed(os.path.join(outdir, 'eveto_events.npz'), **d)
    return d


def regions(d):
    """Returns {name: mask} and {name: number of ambiguous events}."""
    leadE = d['reco_leading_electron_energy']
    lead_primary = gt(d['reco_leading_electron_primary_softmax'], 0.5)
    no_primary_e = d['reco_electron_multiplicity'] == 0
    sr = (d['reco_vertex_distance_cut'] == 1) & (d['reco_has_particle_multiplicity'] == 1) \
        & gt(d['reco_leading_electron_primary_softmax'], PRIMARY_SOFTMAX) \
        & gt(d['reco_leading_electron_energy'], SR_THRESHOLD) \
        & gt(d['reco_leading_electron_electron_softmax'], ELECTRON_SOFTMAX)
    pho = gt(d['reco_leading_photon_energy'], PHOTON_SB) & gt(d['reco_leading_photon_primary_softmax'], PRIMARY_SOFTMAX)
    regs, ambig = {'SR': sr}, {'SR': 0}
    for t in VETOES:
        if t == 0:
            ok, amb = no_primary_e, np.zeros_like(no_primary_e)
        else:
            below = ~gt(leadE, t - 1e-9)  # leading electron KE < t (or no electron)
            ok = no_primary_e | below
            amb = ~ok & ~lead_primary
        regs[f'veto {t:g}'] = pho & ok
        ambig[f'veto {t:g}'] = int((pho & amb).sum())
        if t > 0:
            regs[f'added {t:g}'] = pho & ok & ~no_primary_e
            ambig[f'added {t:g}'] = ambig[f'veto {t:g}']
    return regs, ambig


def ks(a, b):
    a, b = np.sort(a[np.isfinite(a)]), np.sort(b[np.isfinite(b)])
    if not len(a) or not len(b):
        return np.nan
    x = np.concatenate([a, b])
    return np.max(np.abs(np.searchsorted(a, x, 'right') / len(a) - np.searchsorted(b, x, 'right') / len(b)))


def pct(x):
    return 100 * np.mean(x) if len(x) else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--from-npz', action='store_true')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    d = dict(np.load(os.path.join(args.outdir, 'eveto_events.npz'))) if args.from_npz else load(args.outdir)

    cat = np.where(np.isfinite(d['true_category']), d['true_category'], -1).astype(int)
    mode, W, Enu, Evis = (d['true_interaction_mode'], d['true_hadronic_invariant_mass'],
                          d['true_neutrino_energy'], d['reco_visible_energy'])
    regs, ambig = regions(d)
    ref = regs['SR'] & np.isin(cat, BKG)
    ebin = np.digitize(Evis, EVIS_BINS) - 1
    summ, bins = [], []
    for name, r in regs.items():
        b = r & np.isin(cat, BKG)
        pho_pdg = np.abs(np.nan_to_num(d['true_pdg_from_reco_photon'][b], nan=0))
        ele_pdg = np.abs(np.nan_to_num(d['true_pdg_from_reco_ele'][b], nan=0))
        row = dict(region=name, n=int(r.sum()), n_ambig=ambig[name], overlap_SR=int((r & regs['SR']).sum()) if name != 'SR' else 0,
                   signal_pct=pct(cat[r] == 0), n_bkg=int(b.sum()),
                   cosmic_nocat_pct=pct(np.isin(cat[r], [-1, 5])),
                   photon_true_gamma_pct=pct(pho_pdg == 22), lead_ele_true_gamma_pct=pct(ele_pdg == 22),
                   lead_ele_true_e_pct=pct(ele_pdg == 11),
                   median_W=np.nanmedian(W[b]) if b.any() else np.nan,
                   ks_W=ks(W[b], W[ref]), ks_Enu=ks(Enu[b], Enu[ref]), ks_Evis=ks(Evis[b], Evis[ref]))
        for mc, ml in MODES:
            row[f'{ml}_pct'] = pct(mode[b] == mc)
        for c in BKG:
            row[f'cat{c}_of_bkg_pct'] = pct(cat[b] == c)
        summ.append(row)
        for i in range(len(EVIS_BINS) - 1):
            rb, bb = r & (ebin == i), b & (ebin == i)
            bins.append(dict(region=name, bin=f'{EVIS_BINS[i]:g}-{EVIS_BINS[i + 1]:g}', n=int(rb.sum()),
                             stat_err_pct=100 / np.sqrt(rb.sum()) if rb.any() else np.nan,
                             signal_pct=pct(cat[rb] == 0), DIS_of_bkg_pct=pct(mode[bb] == 2),
                             RES_of_bkg_pct=pct(mode[bb] == 1),
                             median_W=np.nanmedian(W[bb]) if bb.any() else np.nan))
    s, bt = pd.DataFrame(summ), pd.DataFrame(bins)
    s.to_csv(os.path.join(args.outdir, 'eveto_summary.csv'), index=False, float_format='%.3f')
    bt.to_csv(os.path.join(args.outdir, 'eveto_bins.csv'), index=False, float_format='%.3f')
    with pd.option_context('display.width', 250, 'display.max_columns', 40, 'display.max_rows', 200):
        print(s.round(2).T.to_string())
        print()
        print(bt.round(2).to_string(index=False))


if __name__ == '__main__':
    main()
