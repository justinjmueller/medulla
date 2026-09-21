#!/usr/bin/env python3
"""
Account for every job ID and input file in one or more medulla projects.

The same reconciliation as `campaign.py reconcile`, for projects that are not
part of a campaign -- for example ones created and launched with medulla.py.
Read-only and needs no ROOT, so it is safe to run while jobs are still
landing.

Usage:
    python3 reconcile.py PROJECT_DIR [PROJECT_DIR ...] [--out DIR] [--max-ids N]

Prints a report per project; with --out, also writes each to
DIR/<project name>.txt. Exits 1 if any project has a gap that nothing else in
the workflow would report (output lost after success, a stale record, a job ID
processed twice within one launch, events read != events written, or input
files that do not add up), so it can gate a merge step.
"""
import sys
from argparse import ArgumentParser
from pathlib import Path

from utilities import reconcile_project, format_reconcile_report


def gaps(r : dict) -> list:
    """The findings that mean part of the dataset is silently wrong or missing."""
    c = r['categories']
    out = []
    if c['output_lost']:
        out.append(f"{len(c['output_lost'])} job ID(s) succeeded but their output is gone")
    if c['complete_stale_record']:
        out.append(f"{len(c['complete_stale_record'])} complete job ID(s) shadowed by a failed attempt's record")
    if r['duplicate_claims']:
        out.append(f"{len(r['duplicate_claims'])} job ID(s) processed twice within one launch")
    if r['events']['mismatched']:
        out.append(f"{len(r['events']['mismatched'])} job ID(s) wrote a different number of events than they read")
    if not r['files']['balanced']:
        out.append("input files do not add up")
    return out


def main(argv=None) -> int:
    p = ArgumentParser(description=__doc__.strip().splitlines()[0])
    p.add_argument('project_dirs', nargs='+', type=Path, help='Project directories to reconcile')
    p.add_argument('--out', type=Path, help='Also write each report to OUT/<project>.txt')
    p.add_argument('--max-ids', type=int, default=50, help='Job IDs to list per category (default 50)')
    args = p.parse_args(argv)

    any_gap = False
    for proj in args.project_dirs:
        r = reconcile_project(proj)
        report = format_reconcile_report(proj.name, r, max_ids=args.max_ids)
        print(report)
        if args.out:
            args.out.mkdir(parents=True, exist_ok=True)
            (args.out / f"{proj.name}.txt").write_text(report)
        found = gaps(r)
        if found:
            any_gap = True
            print(f"GAPS in {proj.name}:")
            for g in found:
                print(f"  - {g}")
            print()
    return 1 if any_gap else 0


if __name__ == '__main__':
    sys.exit(main())
