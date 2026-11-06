#!/usr/bin/env python3
"""Find which flat CAF files hold a list of (Run, Evt) events.

Read-only. Opens each CAF with uproot and reads only rec.hdr.run/subrun/evt from
recTree, matches on (Run, Evt) (the merged medulla output has Subrun = 1 for every
data event, so Subrun is not used for matching; the CAF's own subrun is reported),
and stops once every requested event is found.

  python3 scripts/find_caf_events.py \
      --events plots_v31/event_displays/v31_data_display_candidates.csv \
      --files '/pnfs/icarus/persistent/users/dcarber/spine/NuMI_data_flat_cafs_2/*.flat.root' \
      --output plots_v31/event_displays/v31_data_display_candidates_cafs.csv

--files takes globs, plain paths, root:// URLs, or @list.txt (one path/URL per line).
--dry-run prints the number of files and target events and opens nothing.
The output CSV is the events CSV plus caf_file, caf_entry, caf_subrun; it is
rewritten after every match, so an interrupted scan keeps what it found.
"""
import argparse
import csv
import glob
import os
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed


def expand_files(specs):
    files = []
    for s in specs:
        if s.startswith('@'):
            with open(s[1:]) as fh:
                files += [l.strip() for l in fh if l.strip() and not l.startswith('#')]
        elif '://' in s:
            files.append(s)
        elif any(c in s for c in '*?['):
            files += sorted(glob.glob(s))
        else:
            files.append(s)
    return files


HDR = ['rec.hdr.run', 'rec.hdr.subrun', 'rec.hdr.evt']


def read_headers(path, backend):
    if backend == 'root':
        # PyROOT reads root:// URLs through its own xrootd plugin (needs a bearer
        # token in $BEARER_TOKEN_FILE or /tmp/bt_u<uid>).
        import ROOT
        a = ROOT.RDataFrame('recTree', path).AsNumpy(HDR)
    else:
        import uproot
        with uproot.open(path) as f:
            a = f['recTree'].arrays(HDR, library='np')
    return path, a[HDR[0]], a[HDR[1]], a[HDR[2]]


def write_out(path, fields, rows):
    tmp = path + '.tmp'
    with open(tmp, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--events', required=True, help='CSV with Run and Evt columns.')
    p.add_argument('--files', required=True, nargs='+', help='CAF globs/paths/URLs or @list.txt.')
    p.add_argument('--output', required=True, help='Output CSV (must not exist).')
    p.add_argument('--workers', type=int, default=8, help='Files opened in parallel (default 8).')
    p.add_argument('--backend', choices=['auto', 'uproot', 'root'], default='auto',
                   help='auto: PyROOT for root:// URLs, uproot for local paths.')
    p.add_argument('--dry-run', action='store_true')
    args = p.parse_args()

    with open(args.events) as fh:
        events = list(csv.DictReader(fh))
    targets = {(int(e['Run']), int(e['Evt'])): e for e in events}
    files = expand_files(args.files)
    print(f'[INFO] {len(targets)} target events, {len(files)} CAF files')
    if not files:
        sys.exit('[ERROR] no CAF files matched --files')
    if os.path.exists(args.output):
        sys.exit(f'[ERROR] {args.output} exists; pick a new name')
    if args.dry_run:
        print('[DRY RUN] first files:', *files[:3], sep='\n  ')
        return

    fields = list(events[0].keys()) + ['caf_file', 'caf_entry', 'caf_subrun']
    found = {}
    bad = []
    backend = args.backend
    if backend == 'auto':
        backend = 'root' if any('://' in f for f in files) else 'uproot'
    # PyROOT is not safe to drive from several threads: use processes for it.
    Pool = ProcessPoolExecutor if backend == 'root' else ThreadPoolExecutor
    print(f'[INFO] backend {backend}, {args.workers} workers')
    with Pool(max_workers=args.workers) as ex:
        futs = {ex.submit(read_headers, f, backend): f for f in files}
        for n, fut in enumerate(as_completed(futs), 1):
            try:
                path, run, sub, evt = fut.result()
            except Exception as err:  # unreadable file: note it, keep going
                bad.append((futs[fut], str(err)))
                continue
            for i, key in enumerate(zip(run.tolist(), evt.tolist())):
                if key in targets and key not in found:
                    found[key] = dict(targets[key], caf_file=path, caf_entry=i, caf_subrun=int(sub[i]))
                    print(f'[MATCH] Run {key[0]} Evt {key[1]} ({targets[key].get("topology", "")}) '
                          f'-> {os.path.basename(path)} entry {i}', flush=True)
            if found:
                write_out(args.output, fields, [found[k] for k in targets if k in found])
            if n % 50 == 0:
                print(f'[INFO] {n}/{len(files)} files read, {len(found)}/{len(targets)} found', flush=True)
            if len(found) == len(targets):
                print('[INFO] all events found; stopping early')
                for f in futs:
                    f.cancel()
                break

    missing = [k for k in targets if k not in found]
    print(f'[SUMMARY] found {len(found)}/{len(targets)}; unreadable files: {len(bad)}')
    for f, err in bad[:10]:
        print(f'[WARN] could not read {f}: {err}')
    for k in missing:
        print(f'[MISSING] Run {k[0]} Evt {k[1]}')


if __name__ == '__main__':
    main()
