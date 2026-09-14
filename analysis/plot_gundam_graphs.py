#!/usr/bin/env python3
"""
Plot per-event weight-vs-sigma TGraph splines from a to_gundam output file.

For each systematic branch (TClonesArray<TGraph>), overlays the response
curves for all events and draws a median + 16-84th-percentile band.
Saves one PDF page per systematic.

Usage:
    python plot_gundam_graphs.py output.root
    python plot_gundam_graphs.py output.root -o my_graphs.pdf
    python plot_gundam_graphs.py output.root -t events/NuMIFull/selected -n 1000
"""
import argparse
import sys
from pathlib import Path

import ROOT
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT.gErrorIgnoreLevel = ROOT.kFatal


def iter_trees(directory, path=''):
    """Recursively yield (path, TTree) pairs from a TFile/TDirectory."""
    for key in directory.GetListOfKeys():
        classname = key.GetClassName()
        name = key.GetName()
        full_path = f'{path}/{name}' if path else name
        if classname in ('TTree', 'TNtuple'):
            yield full_path, key.ReadObj()
        elif 'Directory' in classname:
            yield from iter_trees(key.ReadObj(), full_path)


def find_tgraph_branches(tree):
    """Return names of branches that are TClonesArray<TGraph>."""
    names = []
    for branch in tree.GetListOfBranches():
        if 'TClonesArray' in branch.GetClassName():
            names.append(branch.GetName())
    return names


def collect_graphs(tree, branch_name, max_events=None):
    """
    Read TGraphs from a TClonesArray branch.
    Returns a list of (x_array, y_array) tuples, one per event.
    """
    arr = ROOT.TClonesArray('TGraph')
    # ROOT.AddressOf gives a TClonesArray** equivalent, which is what
    # TBranchObject::SetAddress expects (not the object pointer itself).
    rc = tree.SetBranchAddress(branch_name, ROOT.AddressOf(arr))
    if rc < 0:
        print(f'  WARNING: SetBranchAddress failed (rc={rc})', file=sys.stderr)
        return []

    n = tree.GetEntries()
    if max_events is not None:
        n = min(n, max_events)

    graphs = []
    for i in range(n):
        tree.GetEntry(i)
        if arr.GetEntries() < 1:
            continue
        g = arr.At(0)
        npts = g.GetN()
        if npts == 0:
            continue
        xs = np.array([g.GetX()[j] for j in range(npts)])
        ys = np.array([g.GetY()[j] for j in range(npts)])
        graphs.append((xs, ys))

    tree.ResetBranchAddresses()
    return graphs


def plot_systematic(ax, graphs, title, max_overlay=300):
    """
    Overlay individual event TGraphs and draw median + percentile band.
    """
    if not graphs:
        return

    # Build a common x-grid from the union of all x-values.
    # For a given systematic all events share the same nsigma points,
    # so in practice this is just one set of values.
    all_x = sorted({float(x) for xs, _ in graphs for x in xs})
    x_grid = np.array(all_x)

    # Interpolate every event onto the common grid.
    interp_ys = []
    for xs, ys in graphs:
        if len(xs) >= 2:
            interp_ys.append(np.interp(x_grid, xs, ys))
    if not interp_ys:
        return
    ys_arr = np.array(interp_ys)

    # Mask out sentinel non-nu values (-5) before computing statistics.
    valid = np.all(ys_arr > -4, axis=1)
    ys_valid = ys_arr[valid]

    step = max(1, len(graphs) // max_overlay)
    for xs, ys in graphs[::step]:
        color = 'steelblue' if np.all(ys > -4) else 'grey'
        ax.plot(xs, ys, color=color, alpha=0.06, linewidth=0.5)

    if ys_valid.size > 0:
        median = np.median(ys_valid, axis=0)
        p16    = np.percentile(ys_valid, 16, axis=0)
        p84    = np.percentile(ys_valid, 84, axis=0)
        ax.fill_between(x_grid, p16, p84, color='orange', alpha=0.35, label='16–84th pct')
        ax.plot(x_grid, median, color='darkorange', linewidth=1.5, label='median')

    ax.axhline(1.0, color='black', linewidth=0.8, linestyle='--', alpha=0.6)
    ax.axvline(0.0, color='black', linewidth=0.8, linestyle='--', alpha=0.6)
    ax.set_xlabel('σ')
    ax.set_ylabel('weight')
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=8, loc='best')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='to_gundam output ROOT file')
    parser.add_argument('-o', '--output', default=None,
                        help='Output PDF path (default: <input>_graphs.pdf)')
    parser.add_argument('-t', '--tree', default=None,
                        help='Specific tree path to read (default: all trees)')
    parser.add_argument('-n', '--max-events', type=int, default=None,
                        help='Max events to read per systematic (default: all)')
    parser.add_argument('--max-overlay', type=int, default=300,
                        help='Max individual curves drawn per plot (default: 300)')
    parser.add_argument('-s', '--sys', default=None,
                        help='Only plot systematics whose branch name contains this string')
    args = parser.parse_args()

    out_path = args.output or (str(Path(args.input).with_suffix('')) + '_graphs.pdf')

    f = ROOT.TFile.Open(args.input)
    if f is None or f.IsZombie():
        print(f'ERROR: cannot open {args.input}', file=sys.stderr)
        sys.exit(1)

    trees = list(iter_trees(f))
    if not trees:
        print('ERROR: no TTrees found in file', file=sys.stderr)
        sys.exit(1)

    if args.tree:
        trees = [(p, t) for p, t in trees if p == args.tree]
        if not trees:
            print(f'ERROR: tree {args.tree!r} not found. Available trees:', file=sys.stderr)
            for p, _ in list(iter_trees(f)):
                print(f'  {p}', file=sys.stderr)
            sys.exit(1)

    n_plots = 0
    with PdfPages(out_path) as pdf:
        for tree_path, tree in trees:
            branch_names = find_tgraph_branches(tree)
            if args.sys:
                branch_names = [b for b in branch_names if args.sys in b]
            if not branch_names:
                continue

            print(f'Tree: {tree_path}  ({len(branch_names)} systematic branches, '
                  f'{tree.GetEntries()} events)')

            for bname in branch_names:
                print(f'  {bname} ...', end='', flush=True)
                graphs = collect_graphs(tree, bname, max_events=args.max_events)
                print(f' {len(graphs)} events read')
                if not graphs:
                    continue

                fig, ax = plt.subplots(figsize=(7, 5))
                plot_systematic(ax, graphs,
                                title=f'{bname}  [{tree_path}]',
                                max_overlay=args.max_overlay)
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
                n_plots += 1

    f.Close()
    print(f'\nSaved {n_plots} plot(s) to: {out_path}')


if __name__ == '__main__':
    main()
