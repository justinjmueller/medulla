#!/usr/bin/env python3
"""
Generate a list of output ROOT files to merge, grouped by sample name.

For each job in project.db:
  - tag "data" or "detector_variation" -> output_jobid{N:04d}.root
  - tag "nominal"                       -> output_systematics_jobid{N:04d}.root

Also includes any output_varsys*.root files found in output/.

Usage:
    python merge_list.py <project_dir> [--output merge_list.txt]
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

        if tag in ('data', 'detector_variation', 'rock'):
            candidate = output_dir / f'output_jobid{jobid:04d}.root'
        elif tag == 'nominal' or tag == "nue_enhanced":
            candidate = output_dir / f'output_systematics_jobid{jobid:04d}.root'
        else:
            print(f'[WARN] Unknown tag "{tag}" for jobid {jobid} (sample "{name}"), skipping.')
            continue

        if not candidate.exists():
            print(f'[WARN] Missing file for jobid {jobid}: {candidate}')
            continue

        files_by_sample.setdefault(name, []).append(str(candidate))

    # varsys files are not per-jobid — collect them all under a special key
    varsys_files = sorted(glob(str(output_dir / 'output_varsys*.root')))
    if varsys_files:
        files_by_sample['_varsys'] = varsys_files

    lines = []
    for paths in files_by_sample.values():
        for p in sorted(paths):
            lines.append(p)

    text = '\n'.join(lines)

    if args.output:
        args.output.write_text(text)
        print(f'Wrote {sum(len(v) for v in files_by_sample.values())} file(s) to {args.output}')
    else:
        print(text)


if __name__ == '__main__':
    main()
