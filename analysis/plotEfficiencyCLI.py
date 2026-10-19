#!/usr/bin/env python3
"""
Command-line version of plotEfficiency.py: make the main efficiency figures for
one TOML and one input file.

Usage:
    python3 analysis/plotEfficiencyCLI.py --toml analysis/efficiency_analysis_allTree.toml \
        --input /exp/icarus/data/users/micarrig/merged_v30.root --output-dir plots_v30/efficiency_allTree

Each figure is saved as <name>.png (300 dpi) and <name>.pkl in the output dir.
Output file names match plotEfficiency.py (e.g. efficiency_vs_trueInvMass).
--figures accepts any [[figure]] name defined in the TOML; the default is the
8 main ones below.
"""
import argparse
import os
import pickle
import sys
from pathlib import Path

os.environ.setdefault('MPLBACKEND', 'Agg')
# Style sheets in the TOMLs are relative to spineplot/ (same as plot_cli_common).
os.environ.setdefault('MEDULLA_PLOT_DIR', str(Path(__file__).resolve().parent.parent / 'spineplot') + '/')
sys.path.append(str(Path(__file__).resolve().parent.parent / 'spineplot'))

import matplotlib.pyplot as plt
import ROOT as r

from efficiency import SpineEfficiency  # noqa: F401  (registers the artist type)
from analysis import Analysis
from utilities import relayout_figure

r.gStyle.SetOptStat(0)

# Figure name in the TOML -> output file stem. Same list and names as the
# uncommented part of plotEfficiency.py.
FIGURES = {
    'selection_efficiency_table': 'selection_efficiency_table',
    'efficiency_vs_energy': 'efficiency_vs_energy',
    'efficiency_vs_electronEnergy': 'efficiency_vs_electronEnergy',
    'efficiency_vs_electronAngle': 'efficiency_vs_electronAngle',
    'efficiency_vs_trueMomentumTransfer': 'efficiency_vs_trueMomentumTransfer',
    'efficiency_vs_trueInvariantMass': 'efficiency_vs_trueInvMass',
    'efficiency_vs_trueLeadingElectronEnergy': 'efficiency_vs_trueLeadingElectronEnergy',
    'efficiency_vs_trueLeadingElectronCosTheta': 'efficiency_vs_trueLeadingElectronAngle',
}


def stem_for(name):
    """Output file stem: the historical name for the main figures, else the figure name."""
    return FIGURES.get(name, name)


def dry_run(args, out):
    """Report the plan without opening the input file (no data is loaded)."""
    import toml
    cfg = toml.load(args.toml)
    inp = Path(args.input)
    size = f'{inp.stat().st_size / 1e9:.1f} GB' if inp.exists() else 'MISSING'
    print(f'[DRY RUN] TOML:       {args.toml}')
    print(f'[DRY RUN] Input:      {inp} ({size}) -- would be read with uproot')
    for name, smp in cfg.get('samples', {}).items():
        print(f"[DRY RUN] Sample {name}: key={smp.get('key')} trees={smp.get('trees')} "
              f"{len(smp.get('branches', []))} branches read, {len(smp.get('precompute', {}))} precomputed")
    print(f'[DRY RUN] Output dir: {out} ({"exists" if out.exists() else "would be created"})')
    existing = [f for f in args.figures if (out / f'{stem_for(f)}.png').exists()]
    for f in args.figures:
        flag = ' (would OVERWRITE)' if f in existing else ''
        print(f'[DRY RUN] Figure {f} -> {stem_for(f)}.png/.pkl{flag}')
    print('[DRY RUN] Memory: each sample loads its tree(s) with the listed branches into pandas '
          '(all rows). Check the free memory on this node before a real run.')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--toml', required=True, help='Analysis TOML with the efficiency figures.')
    p.add_argument('--input', required=True, help='Input ROOT file.')
    p.add_argument('--output-dir', required=True, help='Directory for the png/pkl files (created if missing).')
    p.add_argument('--figures', nargs='+', default=list(FIGURES),
                   help='TOML [[figure]] names to make (default: the 8 main ones).')
    p.add_argument('--dry-run', action='store_true',
                   help='Print what would be loaded and written, then exit. Opens no ROOT '
                        'file and creates nothing.')
    args = p.parse_args()

    import toml
    defined = [f.get('name') for f in toml.load(args.toml).get('figure', [])]
    unknown = [f for f in args.figures if f not in defined]
    if unknown:
        p.error(f'figure(s) {unknown} not defined in {args.toml}; it defines {defined}')

    out = Path(args.output_dir)

    if args.dry_run:
        dry_run(args, out)
        return
    out.mkdir(parents=True, exist_ok=True)

    ana = Analysis(args.toml, args.input)
    failed = []
    for name in args.figures:
        stem = stem_for(name)
        print(f'[INFO] -- Making {name} -> {out / stem}.png', flush=True)
        try:
            fig = ana.run_interactively(name)
            # Same size and trimming as the other plot CLIs (plot_cli_common);
            # the legend and headroom are redone for the new size.
            fig.set_size_inches(10, 8)
            relayout_figure(fig)
            plt.savefig(out / f'{stem}.png', dpi=300, bbox_inches='tight')
            with open(out / f'{stem}.pkl', 'wb') as f:
                pickle.dump(fig, f)
        except Exception as e:  # keep going so one bad figure does not lose the rest
            print(f'[ERROR] -- {name} failed: {type(e).__name__}: {e}', flush=True)
            failed.append(name)
        finally:
            plt.close('all')

    print(f'[INFO] -- Done. {len(args.figures) - len(failed)}/{len(args.figures)} figures written to {out}')
    if failed:
        print(f'[ERROR] -- Failed: {failed}')
        sys.exit(1)


if __name__ == '__main__':
    main()
