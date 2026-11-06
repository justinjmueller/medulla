#!/usr/bin/env python3
"""Leading-electron energy and vertex distance in selected NC events vs signal (NuMIFull nominal MC).

Read-only. Four panels:
  (a) reco leading-electron energy, NC vs signal (area normalized);
  (b) NC with >= 2 true photons: reco electron energy / true leading-photon energy and
      / (leading + secondary photon energy). A single photon reconstructed as the electron
      peaks near 1 in the first ratio; two merged photons peak near 1 in the second;
  (c) NC with >= 2 true photons: reco E_e / E_gamma1 vs the true photon opening angle;
  (d) leading-electron vertex distance, NC vs signal (area normalized).

  spineplot/.pixi/envs/default/bin/python scripts/plot_nc_electron_energy.py \
      --outdir plots_v31/background_topologies
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import uproot

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'spineplot'))
from style import load_style_from_toml, stamp  # noqa: E402

FILE = '/exp/icarus/data/users/micarrig/merged_final_v31e_flux.root'
TREE = 'events/NuMIFull/selected'
BR = ['true_category', 'reco_leading_electron_energy', 'reco_leading_electron_vertex_distance',
      'true_leading_photon_energy', 'true_secondary_photon_energy', 'true_photon_opening_angle',
      'true_photon_multiplicity']


def norm_hist(ax, v, edges, label, **kw):
    v = v[np.isfinite(v)]
    v = np.clip(v, edges[0], edges[-1] - 1e-9)
    ax.hist(v, bins=edges, weights=np.full(len(v), 1 / max(len(v), 1)), label=f'{label} ({len(v):,})', **kw)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    style = load_style_from_toml('nue_icarus')
    plt.style.use(style.get_style())

    a = uproot.open(FILE)[TREE].arrays(BR, library='np')
    cat = a['true_category']
    nc, sig = cat == 4, cat == 0
    fill = dict(histtype='stepfilled', color=style.get_color(0), alpha=0.6)
    line = dict(histtype='step', color='black', lw=2)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10.5), layout='none')
    fig.subplots_adjust(left=0.08, right=0.97, top=0.88, bottom=0.07, wspace=0.25, hspace=0.3)
    (ax_e, ax_r), (ax_2d, ax_v) = axes

    # (a) electron energy
    edges = np.linspace(0, 3000, 31)
    norm_hist(ax_e, a['reco_leading_electron_energy'][nc], edges, 'NC', **fill)
    norm_hist(ax_e, a['reco_leading_electron_energy'][sig], edges, r'Signal $\nu_e$ CC', **line)
    ax_e.set_xlabel('Reco leading-electron energy [MeV]')
    ax_e.set_ylabel('Fraction of events')
    ax_e.legend(fontsize=10)

    # (b) energy ratios, NC with two true photons
    e = a['reco_leading_electron_energy']
    g1, g2 = a['true_leading_photon_energy'], a['true_secondary_photon_energy']
    two = nc & np.isfinite(g1) & np.isfinite(g2) & (g1 > 0) & (g2 > 0) & (a['true_photon_multiplicity'] >= 2)
    edges = np.linspace(0, 2.5, 51)
    norm_hist(ax_r, (e / g1)[two], edges, r'$E_e^{reco}$ / $E_{\gamma 1}^{true}$', **fill)
    norm_hist(ax_r, (e / (g1 + g2))[two], edges, r'$E_e^{reco}$ / ($E_{\gamma 1}^{true}$ + $E_{\gamma 2}^{true}$)',
              histtype='step', color=style.get_color(1), lw=2)
    ax_r.axvline(1, color='gray', ls='--', lw=1)
    ax_r.set_xlabel('Energy ratio')
    ax_r.set_ylabel('Fraction of events')
    ax_r.set_title(r'NC with $\geq$2 true photons', loc='center', fontsize=12)
    ax_r.legend(fontsize=10)

    # (c) ratio vs opening angle
    ang = np.degrees(a['true_photon_opening_angle'][two])
    h = ax_2d.hist2d(np.clip(ang, 0, 179.9), np.clip((e / g1)[two], 0, 2.49), bins=[36, 50],
                     range=[[0, 180], [0, 2.5]], cmap='Blues', cmin=1)
    fig.colorbar(h[3], ax=ax_2d, label='Events')
    ax_2d.axhline(1, color='gray', ls='--', lw=1)
    ax_2d.set_xlabel(r'True photon opening angle [$^\circ$]')
    ax_2d.set_ylabel(r'$E_e^{reco}$ / $E_{\gamma 1}^{true}$')
    ax_2d.set_title(r'NC with $\geq$2 true photons', loc='center', fontsize=12)

    # (d) vertex distance
    edges = np.linspace(0, 3.5, 36)
    norm_hist(ax_v, a['reco_leading_electron_vertex_distance'][nc], edges, 'NC', **fill)
    norm_hist(ax_v, a['reco_leading_electron_vertex_distance'][sig], edges, r'Signal $\nu_e$ CC', **line)
    ax_v.set_yscale('log')
    ax_v.set_xlabel('Reco leading-electron vertex distance [cm]')
    ax_v.set_ylabel('Fraction of events')
    ax_v.legend(fontsize=10)

    fig.suptitle('Leading Electron in Selected NC Events', y=0.975)
    hdr = fig.add_axes([0.08, 0.905, 0.89, 0.0001])
    hdr.axis('off')
    stamp(hdr)
    out = os.path.join(args.outdir, 'nc_leading_electron_energy.png')
    fig.savefig(out)
    print('[INFO] wrote', out)

    # numbers for the chat
    r1, rs = (e / g1)[two], (e / (g1 + g2))[two]
    print(f'[INFO] NC with 2 true photons: {two.sum()} of {nc.sum()}')
    print(f'[INFO] |E/Eg1-1|<0.2: {100*np.mean(abs(r1-1)<0.2):.1f}%   |E/(Eg1+Eg2)-1|<0.2: {100*np.mean(abs(rs-1)<0.2):.1f}%'
          f'   E/Eg1>1.2: {100*np.mean(r1>1.2):.1f}%')
    small = ang < 20
    print(f'[INFO] opening angle <20deg: {small.sum()} events, E/Eg1>1.2 there: {100*np.mean(r1[small]>1.2):.1f}%;'
          f' >=20deg: {100*np.mean(r1[~small]>1.2):.1f}%')
    for lab, s in (('NC', nc), ('signal', sig)):
        print(f'[INFO] {lab}: median E_e {np.median(e[s]):.0f} MeV, median vertex dist '
              f'{np.median(a["reco_leading_electron_vertex_distance"][s]):.2f} cm, '
              f'dist>1cm {100*np.mean(a["reco_leading_electron_vertex_distance"][s]>1):.1f}%')


if __name__ == '__main__':
    main()
