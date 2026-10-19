#!/usr/bin/env python3
"""
Run every plotting job in a manifest, one family at a time.

Reads a manifest TOML (default analysis/plot_jobs.toml) and runs each job's CLI
as a separate subprocess, so memory is released between families. A memory gate
refuses a real run when this node's available memory does not comfortably exceed
the input file size -- so the SAME command is safe on the ~11 GB ICARUS GPVM (it
refuses and points you at the big-memory server) and runs for real on a big node.

--dry-run passes --dry-run to every job (no ROOT opened, nothing written) and
ignores the memory gate. A job whose input file is missing is skipped with a
warning, not an error.

Usage:
    python3 analysis/plot_all.py --dry-run
    python3 analysis/plot_all.py --only selection efficiency_allTree
    python3 analysis/plot_all.py                  # real run (needs a big node)
    python3 analysis/plot_all.py --mem-factor 1.5 # require 1.5x the input size free

Each job is one invocation of a plotting CLI (see analysis/plot_jobs.toml); the
CLIs share the --toml/--input/--output-dir/--figures/--dry-run interface.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import toml

DEFAULT_MANIFEST = Path(__file__).resolve().parent / 'plot_jobs.toml'


def load_jobs(manifest_path):
    """Return the list of job dicts from a manifest TOML."""
    cfg = toml.load(manifest_path)
    jobs = cfg.get('job', [])
    required = ('name', 'cli', 'toml', 'input', 'output_dir')
    for j in jobs:
        missing = [k for k in required if k not in j]
        if missing:
            raise ValueError(f"job {j.get('name', '?')!r} is missing fields: {missing}")
    names = [j['name'] for j in jobs]
    dupes = {n for n in names if names.count(n) > 1}
    if dupes:
        raise ValueError(f'duplicate job name(s): {sorted(dupes)}')
    return jobs


def mem_available_bytes():
    """Available memory in bytes from /proc/meminfo, or None if unreadable."""
    try:
        with open('/proc/meminfo') as f:
            for line in f:
                if line.startswith('MemAvailable:'):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    return None


def fits(available_bytes, input_bytes, factor=1.0):
    """True if there is comfortably enough memory to load input_bytes."""
    if available_bytes is None or input_bytes is None:
        return False
    return available_bytes >= input_bytes * factor


def build_cmd(job, dry_run):
    cmd = [sys.executable, job['cli'],
           '--toml', job['toml'],
           '--input', job['input'],
           '--output-dir', job['output_dir']]
    if job.get('figures'):
        cmd += ['--figures', *job['figures']]
    if job.get('args'):
        cmd += [str(a) for a in job['args']]
    if dry_run:
        cmd.append('--dry-run')
    return cmd


def run_jobs(jobs, dry_run=False, mem_factor=1.0, available_bytes=None,
             only=None, runner=subprocess.run, out=print):
    """
    Run (or dry-run) each job. Returns [(name, status)] where status is one of
    'dry-run', 'ok', 'failed', 'missing-input', 'skipped-memory', 'not-selected'.
    Pure enough to unit-test: inject available_bytes and a fake runner.
    """
    results = []
    for job in jobs:
        name = job['name']
        if only and name not in only:
            results.append((name, 'not-selected'))
            continue
        inp = Path(job['input'])
        size = inp.stat().st_size if inp.exists() else None
        if size is None:
            out(f'[SKIP] {name}: input missing: {inp}')
            results.append((name, 'missing-input'))
            continue
        if not dry_run and not fits(available_bytes, size, mem_factor):
            have = 0 if available_bytes is None else available_bytes / 1e9
            out(f'[SKIP] {name}: needs ~{size / 1e9:.0f} GB free (x{mem_factor}), '
                f'only {have:.0f} GB available. Run on the big-memory server.')
            results.append((name, 'skipped-memory'))
            continue
        cmd = build_cmd(job, dry_run)
        out(f'[{"DRY" if dry_run else "RUN"}] {name}: {" ".join(cmd)}')
        r = runner(cmd)
        rc = getattr(r, 'returncode', 0)
        if dry_run:
            results.append((name, 'dry-run'))
        else:
            results.append((name, 'ok' if rc == 0 else 'failed'))
    return results


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest', default=str(DEFAULT_MANIFEST),
                   help=f'Manifest TOML of jobs (default: {DEFAULT_MANIFEST}).')
    p.add_argument('--only', nargs='+', default=None,
                   help='Only run jobs with these names (default: all).')
    p.add_argument('--mem-factor', type=float, default=1.0,
                   help='Require available memory >= this x input size (default: 1.0).')
    p.add_argument('--dry-run', action='store_true',
                   help='Pass --dry-run to every job; ignore the memory gate.')
    args = p.parse_args()

    jobs = load_jobs(args.manifest)
    if args.only:
        known = {j['name'] for j in jobs}
        unknown = [n for n in args.only if n not in known]
        if unknown:
            p.error(f'unknown job name(s) {unknown}; choose from {sorted(known)}')

    avail = mem_available_bytes()
    results = run_jobs(jobs, dry_run=args.dry_run, mem_factor=args.mem_factor,
                       available_bytes=avail, only=args.only)

    ran = [n for n, s in results if s in ('ok', 'dry-run')]
    skipped = [(n, s) for n, s in results if s in ('missing-input', 'skipped-memory')]
    failed = [n for n, s in results if s == 'failed']
    print('\n[SUMMARY] plot_all.py')
    for name, status in results:
        if status != 'not-selected':
            print(f'  {name}: {status}')
    print(f'  -> {len(ran)} ran, {len(skipped)} skipped, {len(failed)} failed')
    if failed:
        sys.exit(1)


if __name__ == '__main__':
    main()
