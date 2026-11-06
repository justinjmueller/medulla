#!/usr/bin/env python3
"""Particle content of each background in the selected sample (NuMIFull nominal MC).

Read-only. For every true_category of events/NuMIFull/<--tree> (default: selected, the signal
region; nc_sideband for the NC pi0 sideband) it counts, with the reco and
the true multiplicities side by side:
  - events with >= 1 electron / proton / pion / photon / muon, and the mean multiplicity,
  - the most common (Ne, Np, Npi, Ngamma) combinations,
  - reco vs true electron count (how often the reco electron has no true electron behind it).
Raw selected-event counts (no POT scaling).

Output in --outdir:
  counts_reco_true.csv   per category: n, and per particle n/% with >= 1 and mean (reco, true)
  top_combos.csv         per category and reco/true: the --top most common combinations
  electron_match.csv     per category: reco electron count vs true electron count
  table_<cat>.png        one table per category (efficiency-table style)

  spineplot/.pixi/envs/default/bin/python scripts/background_topologies.py \
      --outdir plots_v31/background_topologies
  spineplot/.pixi/envs/default/bin/python scripts/background_topologies.py \
      --tree nc_sideband --outdir plots_v31/background_topologies_sideband
"""
import argparse
import os
import sys
from collections import Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'spineplot'))
from style import load_style_from_toml, stamp  # noqa: E402

POT = 2.5e20
FILE = '/exp/icarus/data/users/micarrig/merged_final_v31e_flux.root'
REGION = {'selected': r'Selected $\nu_e$', 'nc_sideband': r'NC $\pi^0$ Sideband'}
PARTS = ['electron', 'proton', 'pion', 'photon', 'muon']
SYM = {'electron': 'e', 'proton': 'p', 'pion': r'$\pi$', 'photon': r'$\gamma$', 'muon': r'$\mu$'}
COMBO = ['electron', 'proton', 'pion', 'photon']  # particles that define a combination
CATS = [(0, r'Signal ($\nu_e$ CC, FV)'), (1, r'Other $\nu_e$ CC'), (2, r'$\nu_\mu$ CC (FV)'),
        (3, 'Other CC (OOFV / uncontained)'), (4, 'NC (FV)'), (5, 'Cosmic'), (6, 'Other'),
        (-1, 'No category')]


def mult(a, kind, p):
    """Multiplicity array; NaN/placeholder values (no truth match) become -1."""
    x = a[f'{kind}_{p}_multiplicity'].astype(float)
    bad = ~np.isfinite(x) | (x < 0) | (x > 1000)
    return np.where(bad, -1, x).astype(int)


def combo_label(row):
    return ''.join(f"{n}{SYM[p]}" for n, p in zip(row, COMBO))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--top', type=int, default=10)
    ap.add_argument('--tree', default='selected', choices=sorted(REGION))
    args = ap.parse_args()
    tree, region = f'events/NuMIFull/{args.tree}', REGION[args.tree]
    os.makedirs(args.outdir, exist_ok=True)
    style = load_style_from_toml('nue_icarus')
    plt.style.use(style.get_style())

    t = uproot.open(FILE)[tree]
    br = ['true_category'] + [f'{k}_{p}_multiplicity' for k in ('reco', 'true') for p in PARTS]
    a = t.arrays(br, library='np')
    cat = a['true_category'].astype(float)
    cat = np.where(np.isfinite(cat), cat, -1).astype(int)
    m = {(k, p): mult(a, k, p) for k in ('reco', 'true') for p in PARTS}
    print(f'[INFO] {tree}: {len(cat)} selected events', flush=True)

    counts, combos, ematch = [], [], []
    for c, name in CATS:
        sel = cat == c
        n = int(sel.sum())
        if n == 0:
            continue
        row = dict(category=c, name=name, n=n)
        for k in ('reco', 'true'):
            has_truth = m[(k, 'electron')][sel] >= 0
            row[f'{k}_n_valid'] = int(has_truth.sum())
            for p in PARTS:
                x = m[(k, p)][sel]
                x = x[x >= 0]
                row[f'{k}_{p}_ge1'] = int((x >= 1).sum())
                row[f'{k}_{p}_ge1_pct'] = 100 * (x >= 1).mean() if len(x) else np.nan
                row[f'{k}_{p}_mean'] = x.mean() if len(x) else np.nan
            x = m[(k, 'electron')][sel]
            x = x[x >= 0]
            row[f'{k}_electron_ge2_pct'] = 100 * (x >= 2).mean() if len(x) else np.nan
            mat = np.stack([m[(k, p)][sel] for p in COMBO], axis=1)
            mat = mat[(mat >= 0).all(axis=1)]
            for (combo, cnt) in Counter(map(tuple, mat)).most_common(args.top):
                combos.append(dict(category=c, name=name, kind=k, combo=combo_label(combo),
                                   ne=combo[0], np=combo[1], npi=combo[2], ngamma=combo[3],
                                   count=cnt, pct=100 * cnt / n))
        counts.append(row)
        re, te = np.clip(m[('reco', 'electron')][sel], -1, 2), np.clip(m[('true', 'electron')][sel], -1, 2)
        for (r_, t_), cnt in sorted(Counter(zip(re, te)).items()):
            ematch.append(dict(category=c, name=name, reco_e=('>=2' if r_ == 2 else r_),
                               true_e=('none' if t_ < 0 else '>=2' if t_ == 2 else t_),
                               count=cnt, pct=100 * cnt / n))
        print(f"[CAT] {c:2d} {name:32s} n={n:6d}  " + '  '.join(
            f"{p[:2]} r{row[f'reco_{p}_ge1_pct']:5.1f}%/t{row[f'true_{p}_ge1_pct']:5.1f}%" for p in PARTS), flush=True)

    dc, dk, de = pd.DataFrame(counts), pd.DataFrame(combos), pd.DataFrame(ematch)
    assert dc.n.sum() == len(cat), 'category counts do not add up'
    dc.to_csv(os.path.join(args.outdir, 'counts_reco_true.csv'), index=False, float_format='%.3f')
    dk.to_csv(os.path.join(args.outdir, 'top_combos.csv'), index=False, float_format='%.2f')
    de.to_csv(os.path.join(args.outdir, 'electron_match.csv'), index=False, float_format='%.2f')

    # One table per category: particle content (reco | true), then the top combinations.
    for _, r in dc.iterrows():
        c = int(r.category)
        rows = [['', r'Reco: $\geq$1', 'Reco: mean', r'True: $\geq$1', 'True: mean']]
        for p in PARTS:
            cell = lambda k: (f"{int(r[f'{k}_{p}_ge1']):,} ({r[f'{k}_{p}_ge1_pct']:.1f}%)"
                              if np.isfinite(r[f'{k}_{p}_ge1_pct']) else '-')
            mean = lambda k: f"{r[f'{k}_{p}_mean']:.2f}" if np.isfinite(r[f'{k}_{p}_mean']) else '-'
            rows.append([p.capitalize(), cell('reco'), mean('reco'), cell('true'), mean('true')])
        rows.append(['Top combinations', 'Reco', 'Events', 'True', 'Events'])
        rk = dk[(dk.category == c) & (dk.kind == 'reco')].reset_index(drop=True)
        tk = dk[(dk.category == c) & (dk.kind == 'true')].reset_index(drop=True)
        for i in range(max(len(rk), len(tk))):
            rc = (rk.combo[i], f"{rk['count'][i]:,} ({rk.pct[i]:.1f}%)") if i < len(rk) else ('', '')
            tc = (tk.combo[i], f"{tk['count'][i]:,} ({tk.pct[i]:.1f}%)") if i < len(tk) else ('', '')
            rows.append([f'#{i + 1}', rc[0], rc[1], tc[0], tc[1]])
        fig, ax = plt.subplots(figsize=(11, 0.34 * len(rows) + 1.3))
        ax.axis('off')
        tab = ax.table(cellText=rows, cellLoc='center', edges='T', colWidths=[0.2, 0.2, 0.2, 0.2, 0.2],
                       bbox=[0, 0, 1, 1])
        tab.auto_set_font_size(False)
        tab.set_fontsize(11)
        sub = len(PARTS) + 1  # the second header row
        for i in range(2, len(rows)):
            for j in range(5):
                tab[i, j].visible_edges = 'B' if i in (len(rows) - 1, sub - 1, sub) else 'open'
        for j in range(5):
            tab[0, j].get_text().set_weight('bold')
            tab[sub, j].get_text().set_weight('bold')
        ax.set_title(f"{r['name']}: {int(r.n):,} {region} events", loc='center')
        stamp(ax, pot=POT)
        out = os.path.join(args.outdir, f'table_cat{c}.png')
        fig.savefig(out)
        plt.close(fig)
        print('[INFO] wrote', out, flush=True)
    summary_table(dc, dk, os.path.join(args.outdir, 'summary_backgrounds.png'), region)
    topology_tables(cat, m, dc, args.outdir, region)
    print(f'[INFO] wrote CSVs to {args.outdir}')


BINS = ['0', '1', r'$\geq$2']


def draw_table(rows, title, out, colw, sep_rows=()):
    """Efficiency-table layout: bold header, rules under the header, after sep_rows and at the bottom."""
    fig, ax = plt.subplots(figsize=(14, 0.36 * len(rows) + 1.3))
    ax.axis('off')
    tab = ax.table(cellText=rows, cellLoc='center', edges='T', colWidths=colw, bbox=[0, 0, 1, 1])
    tab.auto_set_font_size(False)
    tab.set_fontsize(11)
    for i in range(2, len(rows)):
        for j in range(len(rows[0])):
            tab[i, j].visible_edges = 'B' if (i == len(rows) - 1 or i in sep_rows) else 'open'
    for j in range(len(rows[0])):
        tab[0, j].get_text().set_weight('bold')
    ax.set_title(title, loc='center')
    stamp(ax, pot=POT)
    fig.savefig(out)
    plt.close(fig)
    print('[INFO] wrote', out, flush=True)


def topology_tables(cat, m, dc, outdir, region, min_pct=1.0):
    """Reco topology fractions per category (columns: signal + each background).
    1) per particle, the fraction with 0 / 1 / >=2 (12 rows);
    2) every (Ne, Np, Npi, Ngamma) combination with 0/1/>=2 bins that is >= min_pct of any category."""
    cats = [(int(r.category), r['name'], int(r.n)) for _, r in dc.iterrows()]
    header = ['Reco topology'] + [f'{name}\n({n:,})' for _, name, n in cats]
    binned = {p: np.clip(m[('reco', p)], 0, 2) for p in COMBO}

    # 1) per particle
    rows, recs, seps = [header], [], []
    for p in COMBO:
        for b in range(3):
            row = [f'{BINS[b]} {SYM[p]}']
            for c, name, n in cats:
                f = 100 * np.sum(binned[p][cat == c] == b) / n
                row.append(f'{f:.1f}%')
                recs.append(dict(particle=p, bin=['0', '1', '>=2'][b], category=c, name=name, pct=f))
            rows.append(row)
        seps.append(len(rows) - 1)
    pd.DataFrame(recs).to_csv(os.path.join(outdir, 'topology_per_particle.csv'), index=False, float_format='%.3f')
    draw_table(rows, f'{region} Events: Reco Multiplicity per Particle',
               os.path.join(outdir, 'topology_per_particle.png'),
               [0.16] + [0.84 / len(cats)] * len(cats), sep_rows=seps[:-1])

    # 2) joint topologies
    key = binned['electron'] * 27 + binned['proton'] * 9 + binned['pion'] * 3 + binned['photon']
    frac = {}
    for k in np.unique(key):
        frac[k] = [100 * np.sum((key == k) & (cat == c)) / n for c, _, n in cats]
    lab = lambda k: ''.join(f"{BINS[(k // d) % 3].replace(r'$\geq$', '≥')}{s}"
                            for d, s in zip((27, 9, 3, 1), ('e', 'p', 'π', 'γ')))
    pd.DataFrame([dict(topology=lab(k), **{name: f for (_, name, _), f in zip(cats, v)})
                  for k, v in frac.items()]).to_csv(os.path.join(outdir, 'topology_joint.csv'),
                                                    index=False, float_format='%.3f')
    keep = sorted([k for k, v in frac.items() if max(v[1:] + [v[0]]) >= min_pct], key=lambda k: -sum(frac[k][1:]))
    rows = [header]
    for k in keep:
        rows.append([lab(k)] + [f'{f:.1f}%' for f in frac[k]])
    other = [100 - sum(frac[k][i] for k in keep) for i in range(len(cats))]
    rows.append([f'Other ({len(frac) - len(keep)} topologies)'] + [f'{f:.1f}%' for f in other])
    draw_table(rows, f'{region} Events: Reco Topology Fractions (0 / 1 / ' + r'$\geq$2 of each particle)',
               os.path.join(outdir, 'topology_joint.png'),
               [0.16] + [0.84 / len(cats)] * len(cats), sep_rows=(len(rows) - 2,))


def summary_table(dc, dk, out, region):
    """One row per background (category > 0): reco/true particle fractions and top reco topologies."""
    pct = lambda v: f'{v:.0f}%' if v >= 10 or v == 0 else f'{v:.1f}%'
    header = ['Background', 'Events', r'Reco: $\geq$1 e', r'Reco: $\geq$2 e', r'Reco: $\geq$1 $\gamma$',
              r'Reco: $\geq$1 $\pi$', r'Reco: $\geq$1 p', r'True: $\geq$1 e', 'Most common reco topologies']
    rows = [header]
    for _, r in dc[dc.category > 0].iterrows():
        top = dk[(dk.category == r.category) & (dk.kind == 'reco')].head(2)
        topo = ', '.join(f'{c} ({p:.0f}%)' for c, p in zip(top.combo, top.pct))
        rows.append([r['name'], f'{int(r.n):,}', pct(r.reco_electron_ge1_pct), pct(r.reco_electron_ge2_pct),
                     pct(r.reco_photon_ge1_pct), pct(r.reco_pion_ge1_pct), pct(r.reco_proton_ge1_pct),
                     pct(r.true_electron_ge1_pct), topo])
    fig, ax = plt.subplots(figsize=(17, 0.45 * len(rows) + 1.3))
    ax.axis('off')
    tab = ax.table(cellText=rows, cellLoc='center', edges='T',
                   colWidths=[0.17, 0.06, 0.075, 0.075, 0.075, 0.075, 0.075, 0.075, 0.30], bbox=[0, 0, 1, 1])
    tab.auto_set_font_size(False)
    tab.set_fontsize(11)
    for i in range(2, len(rows)):
        for j in range(len(header)):
            tab[i, j].visible_edges = 'B' if i == len(rows) - 1 else 'open'
    for j in range(len(header)):
        tab[0, j].get_text().set_weight('bold')
    ax.set_title(f'{region} Backgrounds: Particle Content', loc='center')
    stamp(ax, pot=POT)
    fig.savefig(out)
    plt.close(fig)
    print('[INFO] wrote', out, flush=True)


if __name__ == '__main__':
    main()
