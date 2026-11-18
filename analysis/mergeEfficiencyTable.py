#!/usr/bin/env python3
"""
Merge two efficiency tables made by plotEfficiencyCLI.py: copy the table of
--base and replace one column (default 'Purity [%]') with the same column of
--column-from. Typical use: efficiencies from the true_signal tree, purity from
the all tree (the true_signal tree has almost no background, so its purity is
meaningless).

Reads only the two pickled figures (pkl/<name>.pkl), no ROOT input, so it is
cheap to run anywhere. The rows are matched by the 'Cut' column and must agree
in order. The ICARUS/POT labels and the figure size are copied from --base.

Usage:
    python3 analysis/mergeEfficiencyTable.py \
        --base plots_v32/efficiency_trueSignalTree/pkl/selection_efficiency_table_signalOnly.pkl \
        --column-from plots_v32/efficiency_allTree/pkl/selection_efficiency_table_signalOnly.pkl \
        --output-dir plots_v32/efficiency_merged \
        --name selection_efficiency_table_signalOnly_trueSignalEff_allTreePurity

Writes <output-dir>/png/<name>.png (300 dpi) and <output-dir>/pkl/<name>.pkl.
"""
import argparse
import os
import pickle
import sys
from pathlib import Path

os.environ.setdefault('MPLBACKEND', 'Agg')
HERE = Path(__file__).resolve().parent
sys.path.append(str(HERE.parent / 'spineplot'))

import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox

from efficiency import fit_table

STYLE = HERE.parent / 'spineplot' / 'styles' / 'muon2024.mplstyle'


def read_table(path):
    """Return (figure, axis, table, rows as lists of str) from a pickled table figure."""
    with open(path, 'rb') as f:
        fig = pickle.load(f)
    tables = [(ax, t) for ax in fig.axes for t in ax.tables]
    if len(tables) != 1:
        sys.exit(f'[ERROR] {path}: expected exactly one table, found {len(tables)}')
    ax, table = tables[0]
    cells = table.get_celld()
    nr = 1 + max(r for r, _ in cells)
    nc = 1 + max(c for _, c in cells)
    rows = [[cells[(r, c)].get_text().get_text() for c in range(nc)] for r in range(nr)]
    return fig, ax, table, rows


def crop_box(fig, pad=0.1):
    """Full width; vertically the table and labels (same as plotEfficiencyCLI.table_bbox)."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = [a.get_window_extent(renderer=renderer)
             for ax in fig.axes for a in (*ax.tables, *ax.texts) if a.get_visible()]
    b = Bbox.union(boxes)
    return Bbox.from_extents(0, b.y0 / fig.dpi - pad, fig.get_size_inches()[0], b.y1 / fig.dpi + pad)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--base', required=True, help='Pickled table figure to copy.')
    p.add_argument('--column-from', required=True, help='Pickled table figure to take the column from.')
    p.add_argument('--column', default='Purity [%]', help="Header of the column to replace (default 'Purity [%%]').")
    p.add_argument('--output-dir', required=True, help='Directory for png/ and pkl/ (created if missing).')
    p.add_argument('--name', required=True, help='Output file stem.')
    args = p.parse_args()

    base_fig, base_ax, base_table, base = read_table(args.base)
    _, _, _, other = read_table(args.column_from)

    hdr = base[0]
    if args.column not in hdr or args.column not in other[0]:
        sys.exit(f"[ERROR] column '{args.column}' missing: base has {hdr}, other has {other[0]}")
    if 'Cut' not in hdr or 'Cut' not in other[0]:
        sys.exit("[ERROR] both tables need a 'Cut' column to match rows")
    bcut, ocut = hdr.index('Cut'), other[0].index('Cut')
    bcol, ocol = hdr.index(args.column), other[0].index(args.column)
    base_cuts = [r[bcut] for r in base[1:]]
    other_cuts = [r[ocut] for r in other[1:]]
    if base_cuts != other_cuts:
        sys.exit(f'[ERROR] cut rows differ:\n  base:  {base_cuts}\n  other: {other_cuts}')

    merged = [list(r) for r in base]
    for i in range(1, len(merged)):
        merged[i][bcol] = other[i][ocol]
    for b, m in zip(base[1:], merged[1:]):
        print(f'[INFO] {b[bcut]:45s} {args.column}: {b[bcol]}  ->  {m[bcol]}')

    # Same style sheet as the plots, but the efficiency tables come out in
    # DejaVu Sans (the sheet's Helvetica-like fonts are not used there).
    plt.style.use(str(STYLE))
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
    fig, ax = plt.subplots(figsize=tuple(base_fig.get_size_inches()))
    ax.axis('off')
    table = ax.table(cellText=merged, colLabels=None, loc='upper center', cellLoc='center', edges='T')
    # Same edges as the base table (header lines, group separators, bottom line).
    bcells = base_table.get_celld()
    for key, cell in table.get_celld().items():
        cell.visible_edges = bcells[key].visible_edges
    fit_table(ax, table)
    for t in base_ax.texts:  # ICARUS label and POT
        ax.text(*t.get_position(), t.get_text(), transform=ax.transAxes, fontsize=t.get_fontsize(),
                color=t.get_color(), ha=t.get_ha(), va=t.get_va(), clip_on=False)

    out = Path(args.output_dir)
    (out / 'png').mkdir(parents=True, exist_ok=True)
    (out / 'pkl').mkdir(parents=True, exist_ok=True)
    png, pkl = out / 'png' / f'{args.name}.png', out / 'pkl' / f'{args.name}.pkl'
    fig.savefig(png, dpi=300, bbox_inches=crop_box(fig))
    with open(pkl, 'wb') as f:
        pickle.dump(fig, f)
    print(f'[INFO] wrote {png} and {pkl}')


if __name__ == '__main__':
    main()
