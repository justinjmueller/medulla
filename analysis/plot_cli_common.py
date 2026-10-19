#!/usr/bin/env python3
"""
Shared pieces for the per-family plotting CLIs (selection, systematic, n-1).

Each family CLI supplies only a default list of figure/variable names and a one
line description; everything else - argument parsing, the --dry-run report, and
the save loop - lives here. These CLIs are the command-line equivalents of the
hardcoded make*Plots.py scripts: same TOML, same `run_interactively`, same
plot_<name>.png/.pkl output, but TOML/input/output/figures are arguments so the
same TOML can be run many ways (different figure subsets, different inputs)
without editing the script. The originals are left untouched (CLAUDE.md rule 9).

Heavy imports (matplotlib, ROOT, spineplot) are deferred into run_figures so
that --dry-run, --help and the unit tests work without ROOT. A plot variant is
selected by (TOML, figure-name subset): run_interactively takes a figure name
and the TOML defines each figure's tree and category `groups`, so there is no
separate --tree/--groups flag.
"""
import argparse
import os
import sys
from pathlib import Path

SPINEPLOT = str(Path(__file__).resolve().parent.parent / 'spineplot')


def build_parser(description, default_figures, with_sys_set=False):
    p = argparse.ArgumentParser(description=description,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--toml', required=True, help='Analysis TOML with the figures/variables.')
    p.add_argument('--input', required=True, help='Input ROOT file.')
    p.add_argument('--output-dir', required=True,
                   help='Directory for the png/pkl files (created if missing).')
    p.add_argument('--figures', nargs='+', default=list(default_figures),
                   help='Figure/variable names to make (default: the full family list).')
    p.add_argument('--dry-run', action='store_true',
                   help='Print what would be loaded and written, then exit. Opens no '
                        'ROOT file and creates nothing.')
    if with_sys_set:
        p.add_argument('--systematics-set', default=None,
                       help='Name of a set in --sets-file; overrides the `systematics` '
                            'list of every SpineSystematics artist, so the one TOML can '
                            'produce the xsec / flux / detector / ... breakdowns.')
        p.add_argument('--sets-file', default=str(Path(__file__).resolve().parent / 'systematic_sets.toml'),
                       help='TOML file with a [sets] table of named systematic lists.')
    return p


def stem_for(name):
    """Output file stem for a figure/variable name (matches make*Plots.py)."""
    return f'plot_{name}'


def apply_systematics_set(cfg, names):
    """
    Set the `systematics` list of every SpineSystematics artist in a loaded
    analysis-TOML dict to `names`. Returns the number of artists changed. Pure:
    mutates and returns a count, no I/O, so it is unit-testable without ROOT.
    """
    count = 0
    for fig in cfg.get('figure', []):
        for art in fig.get('artists', []):
            if art.get('type') == 'SpineSystematics':
                art['systematics'] = list(names)
                count += 1
    return count


def load_set(sets_file, name):
    """Return the list of systematic names for `name` from a [sets] TOML file."""
    import toml
    sets = toml.load(sets_file).get('sets', {})
    if name not in sets:
        raise SystemExit(f"[ERROR] systematics set {name!r} not in {sets_file}; "
                         f"available: {sorted(sets)}")
    return sets[name]


def dry_run(args):
    """Report the plan without opening the input file (no data is loaded)."""
    import toml
    cfg = toml.load(args.toml)
    inp = Path(args.input)
    size = f'{inp.stat().st_size / 1e9:.1f} GB' if inp.exists() else 'MISSING'
    out = Path(args.output_dir)
    known = set(cfg.get('variables', {})) | {f.get('name') for f in cfg.get('figure', [])}
    print(f'[DRY RUN] TOML:       {args.toml}')
    print(f'[DRY RUN] Input:      {inp} ({size}) -- would be read with ROOT')
    for name, smp in cfg.get('samples', {}).items():
        print(f"[DRY RUN] Sample {name}: key={smp.get('key')} trees={smp.get('trees')} "
              f"{len(smp.get('branches', []))} branches")
    print(f'[DRY RUN] Output dir: {out} ({"exists" if out.exists() else "would be created"})')
    if getattr(args, 'systematics_set', None):
        knobs = load_set(args.sets_file, args.systematics_set)
        n_artists = apply_systematics_set(cfg, knobs)
        print(f'[DRY RUN] Systematics set: {args.systematics_set} ({len(knobs)} names) '
              f'-> applied to {n_artists} SpineSystematics artist(s)')
    for name in args.figures:
        miss = '' if name in known else ' (NOT in TOML variables/figures - would be skipped)'
        over = ' (would OVERWRITE)' if (out / f'{stem_for(name)}.png').exists() else ''
        print(f'[DRY RUN] Figure {name} -> {stem_for(name)}.png/.pkl{miss}{over}')
    print('[DRY RUN] Memory: each sample loads its tree(s) with the listed branches into '
          'pandas (all rows). Check free memory before a real run.')


def run_figures(args):
    """Real run. Heavy imports are here so --dry-run / tests need no ROOT."""
    import pickle
    import tempfile
    os.environ.setdefault('MPLBACKEND', 'Agg')
    os.environ.setdefault('MEDULLA_PLOT_DIR', SPINEPLOT + '/')
    if SPINEPLOT not in sys.path:
        sys.path.append(SPINEPLOT)
    import matplotlib.pyplot as plt
    from efficiency import SpineEfficiency  # noqa: F401  (registers the artist type)
    from analysis import Analysis
    from utilities import relayout_figure

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # If a systematics set was chosen, write a temp TOML with the overridden
    # `systematics` lists and hand that to Analysis. The original is untouched.
    toml_for_analysis = args.toml
    tmp_toml = None
    if getattr(args, 'systematics_set', None):
        import toml as _toml
        cfg = _toml.load(args.toml)
        knobs = load_set(args.sets_file, args.systematics_set)
        n = apply_systematics_set(cfg, knobs)
        print(f'[INFO] -- systematics set {args.systematics_set}: {len(knobs)} names '
              f'applied to {n} artist(s)', flush=True)
        fd = tempfile.NamedTemporaryFile('w', suffix='.toml', delete=False, dir=out)
        _toml.dump(cfg, fd)
        fd.close()
        tmp_toml = fd.name
        toml_for_analysis = tmp_toml

    ana = Analysis(toml_for_analysis, args.input)
    failed = []
    for name in args.figures:
        stem = stem_for(name)
        print(f'[INFO] -- Making {name} -> {out / stem}.png', flush=True)
        try:
            fig = ana.run_interactively(name)
            fig.set_size_inches(10, 8)
            # Legend size and headroom depend on the figure size.
            relayout_figure(fig)
            for ax in fig.axes:
                leg = ax.get_legend()
                if leg is not None:
                    leg.get_frame().set_facecolor('white')
                    leg.get_frame().set_alpha(0.8)  # hides cut lines behind it; data is kept below by the headroom
            plt.savefig(out / f'{stem}.png', dpi=300, bbox_inches='tight')
            with open(out / f'{stem}.pkl', 'wb') as f:
                pickle.dump(fig, f)
        except Exception as e:  # keep going so one bad figure does not lose the rest
            print(f'[ERROR] -- {name} failed: {type(e).__name__}: {e}', flush=True)
            failed.append(name)
        finally:
            plt.close('all')
    if tmp_toml is not None:
        try:
            os.unlink(tmp_toml)
        except OSError:
            pass
    print(f'[INFO] -- Done. {len(args.figures) - len(failed)}/{len(args.figures)} '
          f'figures written to {out}')
    return failed


def main(description, default_figures, with_sys_set=False):
    args = build_parser(description, default_figures, with_sys_set=with_sys_set).parse_args()
    if args.dry_run:
        dry_run(args)
        return
    failed = run_figures(args)
    if failed:
        print(f'[ERROR] -- Failed: {failed}')
        sys.exit(1)
