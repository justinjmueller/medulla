#!/usr/bin/env python3
"""
Event-matched comparison of the photon sideband between the det-var CV and each
detector variation (same generated neutrino interactions, different detector).

Interactions are keyed by (Evt, true_neutrino_energy) like cell 38 of
detectorSystematics.ipynb (Run/Subrun are ~constant in this MC). Only
interactions present in the `all` tree of BOTH samples are used, so the
CV's larger exposure and interaction-level statistics drop out. For those
shared interactions we count photon_sideband membership:
  both     in the sideband in the CV and in the variation
  cv_only  in the CV's sideband only (migrated out in the variation)
  var_only in the variation's sideband only (migrated in)
and the matched ratio var/cv = (both + var_only) / (both + cv_only), per
reco_visible_energy bin (bin taken from each sample's own reco value).

    /usr/bin/python3 analysis/photonSidebandMatchedCheck.py \
        --stage0 /exp/icarus/data/users/micarrig/merged_stage0_v32.root \
        [--variations null induction_high ...] [--tree photon_sideband]
"""
import argparse

import numpy as np
import pandas as pd
import uproot

BINS = np.array([0.0, 0.81, 1.1, 1.4, 5.0])
STD = ['induction_high', 'induction_low', 'coherent_high', 'coherent_low', 'intrinsic_high', 'intrinsic_low',
       'lifetime_high', 'lifetime_low', 'scintillation', 'wire_gap', 'yz_uniformity', 'cathode_bend', 'null']
KEY = ['Evt', 'true_neutrino_energy']


def interactions(f, sample):
    """Set of (Evt, true E) of all interactions with a true neutrino in the `all` tree."""
    keys = set()
    for chunk in f[f'events/{sample}/all'].iterate(KEY, library='np', step_size='200 MB'):
        e = chunk['true_neutrino_energy']
        ok = np.isfinite(e) & (e > 0)
        keys.update(zip(chunk['Evt'][ok].tolist(), np.round(e[ok], 6).tolist()))
    return keys


def sideband(f, sample, tree):
    a = f[f'events/{sample}/{tree}'].arrays(KEY + ['reco_visible_energy', 'true_category'], library='np')
    df = pd.DataFrame({'Evt': a['Evt'], 'E': np.round(a['true_neutrino_energy'], 6),
                       'visE': a['reco_visible_energy'], 'cat': a['true_category']})
    df = df[np.isfinite(df['E']) & (df['E'] > 0)]
    return df.drop_duplicates(['Evt', 'E'])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--stage0', required=True)
    p.add_argument('--variations', nargs='*', default=STD)
    p.add_argument('--tree', default='photon_sideband')
    a = p.parse_args()

    f = uproot.open(a.stage0)
    cv_int = interactions(f, 'cv')
    cv_sb = sideband(f, 'cv', a.tree)
    print(f'cv: {len(cv_int)} interactions with a true neutrino, {len(cv_sb)} in {a.tree}')
    tot = {'both': np.zeros(len(BINS) - 1), 'cv_only': np.zeros(len(BINS) - 1), 'var_only': np.zeros(len(BINS) - 1)}
    for s in a.variations:
        v_int = interactions(f, s)
        common = cv_int & v_int
        v_sb = sideband(f, s, a.tree)
        c = cv_sb[[k in common for k in zip(cv_sb['Evt'], cv_sb['E'])]]
        v = v_sb[[k in common for k in zip(v_sb['Evt'], v_sb['E'])]]
        m = c.merge(v, on=['Evt', 'E'], how='outer', suffixes=('_cv', '_var'), indicator=True)
        both = m[m['_merge'] == 'both']
        cvo = m[m['_merge'] == 'left_only']
        vo = m[m['_merge'] == 'right_only']
        h = lambda x: np.histogram(x[np.isfinite(x)], bins=BINS)[0].astype(float)
        hb_cv, hb_v = h(both['visE_cv'].to_numpy()), h(both['visE_var'].to_numpy())
        hc, hv = h(cvo['visE_cv'].to_numpy()), h(vo['visE_var'].to_numpy())
        num, den = hb_v + hv, hb_cv + hc
        r = num / den
        print(f'{s:16s} shared interactions {len(common):7d} ({len(common)/len(v_int):.0%} of the variation\'s) | '
              f'sideband both {len(both):4d}  cv-only {len(cvo):4d}  var-only {len(vo):4d} | '
              f'matched var/cv per bin ' + ' '.join(f'{x:.3f}' for x in r))
        tot['both'] += hb_cv
        tot['cv_only'] += hc
        tot['var_only'] += hv
        tot.setdefault('num', np.zeros(len(BINS) - 1))
        tot.setdefault('den', np.zeros(len(BINS) - 1))
        tot['num'] += num
        tot['den'] += den
    r = tot['num'] / tot['den']
    # Binomial-like error on the migrations only (the shared 'both' events cancel)
    err = np.sqrt(tot['cv_only'] + tot['var_only']) / tot['den']
    print('pooled matched var/cv per bin: ' + '  '.join(f'{x:.3f}+-{e:.3f}' for x, e in zip(r, err)))


if __name__ == '__main__':
    main()
