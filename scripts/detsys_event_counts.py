#!/usr/bin/env python3
"""Add event counts per detector-variation sample to the corrected-POT table.

Read-only on the ROOT files. For each row of the POT CSV (scripts/run_corrected_pot.py)
it opens events/<sample>/ in the merged file for that set and adds:
  n_all           entries in the `all` tree
  n_selected      entries in the `selected` tree
  n_sel_signal    selected with true_category == 0
  n_sel_bkgd      selected with true_category > 0
  n_sel_non_nue   selected with true_category > 1 (not an electron-neutrino CC event)
  n_sel_nocat     selected with true_category = NaN (no category; in neither sum above)
Only the true_category branch is read.

  spineplot/.pixi/envs/default/bin/python scripts/detsys_event_counts.py \
      --pot plots_v31/pot/detsys_corrected_pot.csv \
      --output plots_v31/pot/detsys_pot_and_counts.csv
"""
import argparse
import csv
import os
import sys

import numpy as np
import uproot

DEFAULT_FILES = {
    'main': '/exp/icarus/data/users/micarrig/merged_final_v31e_flux.root',
    'nue': '/exp/icarus/data/users/micarrig/merged_nue_v31.root',
}
NEW = ['input_file', 'n_all', 'n_selected', 'n_sel_signal', 'n_sel_bkgd', 'n_sel_non_nue', 'n_sel_nocat']


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--pot', required=True, help='CSV from run_corrected_pot.py')
    p.add_argument('--output', required=True, help='output CSV (must not exist)')
    p.add_argument('--main-file', default=DEFAULT_FILES['main'])
    p.add_argument('--nue-file', default=DEFAULT_FILES['nue'])
    args = p.parse_args()
    if os.path.exists(args.output):
        sys.exit(f'[ERROR] {args.output} exists; pick a new name')
    files = {'main': args.main_file, 'nue': args.nue_file}

    with open(args.pot) as fh:
        rows = list(csv.DictReader(fh))
    opened = {k: uproot.open(v) for k, v in files.items()}
    for r in rows:
        f = opened[r['set']]
        d = f'events/{r["sample"]}'
        cat = f[f'{d}/selected']['true_category'].array(library='np')
        r.update(input_file=files[r['set']],
                 n_all=f[f'{d}/all'].num_entries,
                 n_selected=len(cat),
                 n_sel_signal=int(np.sum(cat == 0)),
                 n_sel_bkgd=int(np.sum(cat > 0)),
                 n_sel_non_nue=int(np.sum(cat > 1)),
                 n_sel_nocat=int(np.sum(np.isnan(cat))))
        print(f"{r['set']:4s} {r['sample']:15s} all={r['n_all']:>9} sel={r['n_selected']:>6} "
              f"sig={r['n_sel_signal']:>6} bkgd={r['n_sel_bkgd']:>5} non_nue={r['n_sel_non_nue']:>5} "
              f"nocat={r['n_sel_nocat']}", flush=True)
    fields = list(rows[0].keys())
    fields = [k for k in fields if k not in NEW] + NEW
    with open(args.output, 'x', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f'[INFO] wrote {args.output}')


if __name__ == '__main__':
    main()
