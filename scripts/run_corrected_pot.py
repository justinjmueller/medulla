#!/usr/bin/env python3
"""Run analysis/compute_corrected_pot_list.C on every detector-variation sample.

Read-only on /pnfs. For each sample in the selection TOMLs with tag
detector_variation (plus the nue det-sys CV, sample NuMI_nue in selection_nue.toml):
  1. list its directory over xrootd (xrdfs ls) and keep files matching the TOML glob,
  2. split the root:// URLs into chunks of --chunk files (<out>/lists/<set>_<sample>_partNNN.txt),
  3. run the macro on each chunk in parallel (logs in <out>/logs/<set>_<sample>_partNNN.log),
  4. sum the chunks' RESULT lines per sample into <out>/detsys_corrected_pot.csv
     (rewritten as samples finish). Every macro output is a sum over files and subruns are
     handled per file, so summing chunks equals one macro run over the whole sample.
Files the macro could not open are written to <out>/failed_files.txt.

  spineplot/.pixi/envs/default/bin/python scripts/run_corrected_pot.py [--dry-run]

Needs a valid bearer token ($BEARER_TOKEN_FILE) and the pixi ROOT (with xrootd) on PATH.
"""
import argparse
import csv
import fnmatch
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import toml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOMLDIR = os.path.join(REPO, 'selection/toml/icarus_nueCC_inclusive_2026')
MACRO = os.path.join(REPO, 'analysis/compute_corrected_pot_list.C')
XRD_HOST = 'fndcadoor.fnal.gov:1094'
SETS = {'main': 'selection.toml', 'nue': 'selection_nue.toml'}
EXTRA = {('nue', 'NuMI_nue')}  # tagged nominal, but it is the nue det-sys CV
SUM_INT = ['files_listed', 'files_ok', 'files_fail', 'events', 'failed_events', 'subruns', 'failed_subruns']
SUM_FLOAT = ['pot_event', 'pot_th1d', 'failed_pot', 'corrected_pot']
FIELDS = (['set', 'sample', 'tag', 'path', 'n_files_listed', 'n_chunks'] + SUM_INT + SUM_FLOAT
          + ['failed_pot_frac', 'pot_event_over_th1d', 'status'])


def samples():
    out = []
    for setname, fname in SETS.items():
        cfg = toml.load(os.path.join(TOMLDIR, fname))
        for s in cfg['sample']:
            if s.get('tag') == 'detector_variation' or (setname, s['name']) in EXTRA:
                out.append(dict(set=setname, sample=s['name'], tag=s.get('tag'), path=s['path']))
    return out


def list_files(path):
    d, glob = os.path.dirname(path), os.path.basename(path)
    if not d.startswith('/pnfs/'):
        raise RuntimeError(f'not a /pnfs path: {path}')
    xd = '/pnfs/fnal.gov/usr/' + d[len('/pnfs/'):]
    r = subprocess.run(['xrdfs', XRD_HOST, 'ls', xd], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f'xrdfs ls {xd} failed: {r.stderr.strip() or r.stdout.strip()}')
    names = sorted(os.path.basename(l.strip()) for l in r.stdout.splitlines() if l.strip())
    return [f'root://{XRD_HOST}/{xd}/{n}' for n in names if fnmatch.fnmatch(n, glob)]


def run_chunk(listfile, logfile):
    """Run the macro on one list; return (RESULT dict or None, failed files, exit code)."""
    cmd = ['root', '-l', '-b', '-q', f'{MACRO}("{listfile}")']
    with open(logfile, 'w') as log:
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL).returncode
    res, failed = None, []
    with open(logfile) as log:
        for l in log:
            if l.startswith('RESULT '):
                res = dict(kv.split('=', 1) for kv in l.split()[1:])
            elif l.startswith('FAILED_FILE '):
                failed.append(l.split(None, 1)[1].strip())
    return res, failed, rc


def summarize(s, parts):
    """Sum the chunk results of one sample into a CSV row."""
    res = dict(s)
    bad = [i for i, (r, _, rc) in enumerate(parts) if r is None or rc != 0]
    for k in SUM_INT:
        res[k] = sum(int(r[k]) for r, _, _ in parts if r)
    for k in SUM_FLOAT:
        res[k] = sum(float(r[k]) for r, _, _ in parts if r)
    pe, th = res['pot_event'], res['pot_th1d']
    res['failed_pot_frac'] = f"{res['failed_pot'] / pe:.6f}" if pe else ''
    res['pot_event_over_th1d'] = f'{pe / th:.6f}' if th else ''
    for k in SUM_FLOAT:
        res[k] = f'{res[k]:.15e}'
    res['n_chunks'] = len(parts)
    res['status'] = 'ok' if not bad else f'{len(bad)} chunks without RESULT or nonzero exit'
    return res


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--out', default=os.path.join(REPO, 'plots_v31/pot'))
    p.add_argument('--workers', type=int, default=24, help='macro processes in parallel (default 24)')
    p.add_argument('--chunk', type=int, default=500, help='files per macro run (default 500)')
    p.add_argument('--only', nargs='+', help='set_sample names to run (e.g. main_null nue_null)')
    p.add_argument('--dry-run', action='store_true', help='list files and counts; run nothing, write nothing')
    args = p.parse_args()

    csvfile = os.path.join(args.out, 'detsys_corrected_pot.csv')
    if not args.dry_run and os.path.exists(csvfile):
        sys.exit(f'[ERROR] {csvfile} exists; pick a new --out')

    todo = [s for s in samples() if not args.only or f"{s['set']}_{s['sample']}" in args.only]
    print(f'[INFO] {len(todo)} samples', flush=True)
    listed = []
    for s in todo:
        try:
            files = list_files(s['path'])
        except Exception as err:
            print(f"[ERROR] {s['set']}/{s['sample']}: {err}", flush=True)
            files = None
        s['n_files_listed'] = len(files) if files is not None else ''
        print(f"[LIST] {s['set']:4s} {s['sample']:15s} {s['tag']:18s} {s['n_files_listed']!s:>5} files  {s['path']}",
              flush=True)
        listed.append((s, files))
    nfiles = sum(len(f) for _, f in listed if f)
    print(f'[INFO] {nfiles} files in total', flush=True)
    if args.dry_run:
        print('[DRY RUN] nothing written')
        return
    bad = [s for s, f in listed if not f]
    if bad:
        sys.exit(f"[ERROR] no files for: {', '.join(s['set'] + '/' + s['sample'] for s in bad)}; stopping")

    os.makedirs(os.path.join(args.out, 'lists'), exist_ok=True)
    os.makedirs(os.path.join(args.out, 'logs'), exist_ok=True)
    jobs = []  # (sample index, list file, log file)
    for si, (s, files) in enumerate(listed):
        tag = f"{s['set']}_{s['sample']}"
        for ci in range(0, len(files), args.chunk):
            stem = f'{tag}_part{ci // args.chunk:03d}'
            lf = os.path.join(args.out, 'lists', stem + '.txt')
            with open(lf, 'w') as fh:
                fh.write('\n'.join(files[ci:ci + args.chunk]) + '\n')
            jobs.append((si, lf, os.path.join(args.out, 'logs', stem + '.log')))
    nchunks = {si: sum(1 for j in jobs if j[0] == si) for si in range(len(listed))}
    print(f'[INFO] {len(jobs)} chunks of up to {args.chunk} files, {args.workers} workers', flush=True)

    parts = {si: [] for si in range(len(listed))}
    done = {}
    failed_all = []

    def write():
        tmp = csvfile + '.tmp'
        with open(tmp, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction='ignore')
            w.writeheader()
            w.writerows(done[i] for i in sorted(done))
        os.replace(tmp, csvfile)
        with open(os.path.join(args.out, 'failed_files.txt'), 'w') as fh:
            fh.write(''.join(f + '\n' for f in failed_all))

    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(run_chunk, lf, log): (si, log) for si, lf, log in jobs}
        for n, fut in enumerate(as_completed(futs), 1):
            si, log = futs[fut]
            r = fut.result()
            parts[si].append(r)
            failed_all.extend(r[1])
            if r[0] is None or r[2] != 0:
                print(f'[WARN] {os.path.basename(log)}: no RESULT line or exit {r[2]}', flush=True)
            if len(parts[si]) == nchunks[si]:
                row = done[si] = summarize(listed[si][0], parts[si])
                print(f"[DONE] {row['set']}/{row['sample']}: {row['status']}  files_ok={row['files_ok']}/"
                      f"{row['n_files_listed']} fail={row['files_fail']}  corrected_pot={row['corrected_pot']}  "
                      f"failed_frac={row['failed_pot_frac']}", flush=True)
            write()
            if n % 20 == 0:
                print(f'[INFO] {n}/{len(jobs)} chunks done, {len(done)}/{len(listed)} samples', flush=True)
    nok = sum(1 for r in done.values() if r['status'] == 'ok' and r['files_fail'] == 0)
    print(f'[SUMMARY] {nok}/{len(listed)} samples ok with no failed files; {len(failed_all)} unreadable files; '
          f'output {csvfile}')


if __name__ == '__main__':
    main()
