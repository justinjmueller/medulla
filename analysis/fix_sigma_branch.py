#!/usr/bin/env python3
"""
Fix _sigma branches in an intermediate systematics ROOT file.

Replaces the sigma vector stored in matching _sigma branches with a corrected
set of values for every event, then writes a new output file. All other
branches and non-tree objects (POT, Livetime histograms, etc.) are copied
unchanged.

Typical use: the ZExpPCA weights were computed at [-3,-2,-1,0,1,2,3] but
the TOML had the wrong ordering [-1,1,-2,2,-3,3], producing scrambled TGraphs
in the to_gundam output. Fix the intermediate file with this script, then
re-run to_gundam.

Usage:
    python fix_sigma_branch.py input.root output.root
    python fix_sigma_branch.py input.root output.root --filter ZExpPCA
    python fix_sigma_branch.py input.root output.root --sigmas -3 -2 -1 0 1 2 3
"""
import argparse
import sys
import ROOT

ROOT.gErrorIgnoreLevel = ROOT.kFatal


def make_sigma_vec(values):
    vec = ROOT.std.vector('double')()
    for v in values:
        vec.push_back(v)
    return vec


def latest_keys(tdir):
    """Return one TKey per unique name — the highest cycle (latest write)."""
    return { k.GetName(): k for k in tdir.GetListOfKeys() }.values()


def copy_non_tree_objects(src_dir, dst_dir):
    """Copy histograms and other non-tree, non-directory objects."""
    for key in latest_keys(src_dir):
        cls = key.GetClassName()
        if 'Directory' in cls or cls in ('TTree', 'TNtuple'):
            continue
        obj = key.ReadObj()
        dst_dir.cd()
        obj.Write(key.GetName())


def log(msg):
    print(msg, flush=True)


def process_directory(src_dir, dst_dir, correct_sigmas, name_filter):
    copy_non_tree_objects(src_dir, dst_dir)

    for key in latest_keys(src_dir):
        cls  = key.GetClassName()
        name = key.GetName()

        if 'Directory' in cls:
            log(f'[DIR ] {src_dir.GetPath()}/{name}')
            sub = dst_dir.mkdir(name)
            process_directory(key.ReadObj(), sub, correct_sigmas, name_filter)

        elif cls in ('TTree', 'TNtuple'):
            tree_in = key.ReadObj()
            n_entries = tree_in.GetEntries()

            sigma_branches = [
                b.GetName() for b in tree_in.GetListOfBranches()
                if name_filter in b.GetName() and b.GetName().endswith('_sigma')
            ]

            dst_dir.cd()

            if not sigma_branches:
                log(f'[COPY] {name}  ({n_entries} entries) ...')
                tree_in.CloneTree(-1).Write()
                log(f'[COPY] {name}  done')
                continue

            log(f'[FIX ] {name}  ({n_entries} entries, {len(sigma_branches)} sigma branch(es)) ...')
            for bname in sigma_branches:
                log(f'       {bname}')

            correct_vec = make_sigma_vec(correct_sigmas)

            # Phase 1: bulk-copy all branches EXCEPT sigma at C++ speed.
            # Disabling a branch in the source causes CloneTree to skip it
            # entirely, so the output tree won't have sigma branches yet.
            tree_in.SetBranchStatus('*', 1)
            for bname in sigma_branches:
                tree_in.SetBranchStatus(bname, 0)
            log(f'[FIX ] {name}  cloning non-sigma branches ...')
            tree_out = tree_in.CloneTree(-1)

            # Phase 2: add the sigma branches back with the corrected constant
            # vector. correct_vec is the same for every event so we never read
            # from tree_in — just create each branch and fill it N times.
            # This is fast: only a small std::vector<double> is serialized per entry.
            log(f'[FIX ] {name}  filling sigma branches ...')
            n = tree_out.GetEntries()
            for bname in sigma_branches:
                br_out = tree_out.Branch(bname, correct_vec)
                for _ in range(n):
                    br_out.Fill()

            tree_out.Write()
            log(f'[FIX ] {name}  done')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input',  help='Intermediate ROOT file to fix')
    parser.add_argument('output', help='Output path for the fixed file')
    parser.add_argument('--filter', '-f', default='ZExpPCA',
                        help='Fix _sigma branches whose name contains this string '
                             '(default: ZExpPCA)')
    parser.add_argument('--sigmas', type=float, nargs='+',
                        default=[-3., -2., -1., 0., 1., 2., 3.],
                        help='Correct sigma values to write '
                             '(default: -3 -2 -1 0 1 2 3)')
    args = parser.parse_args()

    print(f'Filter  : {args.filter!r}')
    print(f'Sigmas  : {args.sigmas}')
    print(f'Input   : {args.input}')
    print(f'Output  : {args.output}')

    f_in = ROOT.TFile.Open(args.input, 'READ')
    if f_in is None or f_in.IsZombie():
        print(f'ERROR: cannot open {args.input}', file=sys.stderr)
        sys.exit(1)

    f_out = ROOT.TFile.Open(args.output, 'RECREATE')

    process_directory(f_in, f_out, args.sigmas, args.filter)

    f_out.Close()
    f_in.Close()
    print(f'Done → {args.output}')


if __name__ == '__main__':
    main()
