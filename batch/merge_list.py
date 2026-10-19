#!/usr/bin/env python3
"""
Generate a list of output ROOT files to merge, grouped by sample name.

For each job in project.db:
  - tag "data", "detector_variation", or "rock" -> output_jobid{N:04d}.root
  - tag "nominal", "nue_enhanced", or "g4"       -> output_systematics_jobid{N:04d}.root

Also includes any output_varsys*.root files found in output/.

Usage:
    python merge_list.py <project_dir> [--output merge_list.txt]
                         [--tags TAG [TAG ...]] [--no-varsys] [--dry-run]

--tags restricts the list to jobs whose sample tag is one of TAG (e.g.
`--tags detector_variation` for the phase 1 input, which covers the variation
samples and the CV sample). --no-varsys leaves out output_varsys*.root files.

--dry-run reads project.db and stats the outputs exactly as a real run does, but
writes no --output file and prints only a per-sample count summary. It leaves no
file behind (the temporary project.db copy is always removed).
"""
import argparse
import sqlite3
import toml
from glob import glob
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('project_dir', type=Path)
    parser.add_argument('--output', type=Path, default=None,
                        help='Write file list to this path (default: print to stdout)')
    parser.add_argument('--tags', nargs='+', default=None,
                        help='Only include jobs whose sample tag is in this list '
                             '(e.g. detector_variation). Default: all known tags.')
    parser.add_argument('--no-varsys', action='store_true',
                        help='Do not add output_varsys*.root files to the list.')
    parser.add_argument('--dry-run', action='store_true',
                        help='Preview the file list (per-sample counts) without '
                             'writing the --output file.')
    args = parser.parse_args()

    project_dir = args.project_dir.resolve()
    output_dir = project_dir / 'output'

    db_path = project_dir / 'project.db'
    if not db_path.exists():
        raise FileNotFoundError(f"project.db not found in {project_dir}")

    # Work on a local copy to avoid dcache issues
    local_db = Path('./project_merge_tmp.db')
    subprocess.run(['cp', str(db_path), str(local_db)], check=True)
    conn = sqlite3.connect(local_db)
    curs = conn.cursor()
    curs.execute("SELECT jobid, cfg FROM configuration")
    rows = curs.fetchall()
    conn.close()
    local_db.unlink()

    files_by_sample = {}

    for jobid, cfg_str in rows:
        cfg = toml.loads(cfg_str)
        samples = cfg.get('sample', [])
        if not samples:
            continue

        sample = samples[0]
        name = sample.get('name', f'job{jobid:04d}')
        tag = sample.get('tag', '')

        if args.tags is not None and tag not in args.tags:
            continue

        if tag in ('data', 'detector_variation', 'rock'):
            candidate = output_dir / f'output_jobid{jobid:04d}.root'
        elif tag in ('nominal', 'nue_enhanced', 'g4'):
            candidate = output_dir / f'output_systematics_jobid{jobid:04d}.root'
        else:
            print(f'[WARN] Unknown tag "{tag}" for jobid {jobid} (sample "{name}"), skipping.')
            continue

        if not candidate.exists():
            print(f'[WARN] Missing file for jobid {jobid}: {candidate}')
            continue

        files_by_sample.setdefault(name, []).append(str(candidate))

    # varsys files are not per-jobid — collect them all under a special key
    varsys_files = [] if args.no_varsys else sorted(glob(str(output_dir / 'output_varsys*.root')))
    if varsys_files:
        files_by_sample['_varsys'] = varsys_files

    lines = []
    for paths in files_by_sample.values():
        for p in sorted(paths):
            lines.append(p)

    text = '\n'.join(lines)
    total = sum(len(v) for v in files_by_sample.values())

    if args.dry_run:
        print('[DRY-RUN] merge_list.py: no file will be written.')
        for name, paths in files_by_sample.items():
            print(f'[DRY-RUN]   {name}: {len(paths)} file(s)')
        dest = str(args.output) if args.output else 'stdout'
        print(f'[DRY-RUN] Would write {total} file(s) to {dest}.')
    elif args.output:
        args.output.write_text(text)
        print(f'Wrote {total} file(s) to {args.output}')
    else:
        print(text)


if __name__ == '__main__':
    main()
