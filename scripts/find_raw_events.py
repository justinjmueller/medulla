#!/usr/bin/env python3
"""Find the raw (compressed) art file, and its stage0/stage1 children, for a list of events.

Read-only: SAM queries over HTTPS and ROOT reads of EventAuxiliary over xrootd.

For each (Run, Evt) in --events:
  1. SAM: raw files in --defname with that run and first_event <= Evt <= last_event
     (several DAQ loggers interleave events, so there are usually ~10 candidates).
  2. Open each candidate's Events tree over xrootd and read EventAuxiliary until
     the event is found (about 2 s per file; only the aux branch is read).
  3. SAM: children of the matched raw file with data_tier reconstructed (stage0),
     and their children (stage1), with locations.

  python3 scripts/find_raw_events.py \
      --events plots_v31/event_displays/v31_data_display_candidates_cafs.csv \
      --output plots_v31/event_displays/v31_data_display_candidates_raw.csv

Needs a bearer token for xrootd ($BEARER_TOKEN_FILE). The SAM server is reached as
samicarus.fnal.gov because samweb.fnal.gov does not resolve on some nodes.
Output: the events CSV plus raw_file, raw_url, raw_entry, raw_subrun, stage0_files,
stage0_urls, stage1_files, stage1_urls (lists joined with ';', newest first).
"""
import argparse
import csv
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

SAM = 'https://samicarus.fnal.gov:8483/sam/icarus/api'
XRD = 'root://fndcadoor.fnal.gov:1094'


def sam(path, **params):
    cmd = ['curl', '-sS', '--fail', '--max-time', '180', '-G', f'{SAM}/{path}']
    for k, v in params.items():
        cmd += ['--data-urlencode', f'{k}={v}']
    # SAM answers 503 when busy: retry with backoff before giving up.
    for wait in (5, 15, 30, 60, 120, None):
        out = subprocess.run(cmd, capture_output=True, text=True)
        if not out.returncode:
            return out.stdout
        if wait is None:
            raise RuntimeError(f'SAM query {path} {params} failed: {out.stderr.strip()}')
        time.sleep(wait)


def sam_list(dims):
    return [l.strip() for l in sam('files/list/', dims=dims).splitlines() if l.strip()]


def xrootd_url(name):
    """First disk location of a SAM file as an xrootd URL (None if not found)."""
    for loc in sam(f'files/name/{name}/locations/').split():
        loc = loc.split('(')[0]
        if ':' in loc:
            loc = loc.split(':', 1)[1]
        if loc.startswith('/pnfs/'):
            return f'{XRD}/pnfs/fnal.gov/usr/{loc[len("/pnfs/"):]}/{name}'
    return None


def read_ids(url):
    import ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kError
    f = ROOT.TFile.Open(url)
    if not f or f.IsZombie():
        raise RuntimeError(f'cannot open {url}')
    t = f.Get('Events')
    t.SetBranchStatus('*', 0)
    t.SetBranchStatus('EventAuxiliary*', 1)
    t.SetEstimate(t.GetEntries() + 1)
    n = t.Draw('EventAuxiliary.id_.subRun_.run_.run_:EventAuxiliary.id_.subRun_.subRun_:EventAuxiliary.id_.event_',
               '', 'goff')
    ids = [(int(t.GetV1()[i]), int(t.GetV2()[i]), int(t.GetV3()[i])) for i in range(n)]
    f.Close()
    return ids


def children(names, tier_word):
    """SAM children of the given files whose file name contains tier_word."""
    out = []
    for n in names:
        out += [c for c in sam_list(f'ischildof: (file_name {n})') if tier_word in c and not c.startswith('hist_')]
    return sorted(set(out), reverse=True)  # names carry a timestamp: newest first


def find_one(ev, defname):
    try:
        return _find_one(ev, defname)
    except Exception as err:  # one bad event must not kill the whole scan
        return dict(ev, errors=f'{type(err).__name__}: {err}')


def _find_one(ev, defname):
    run, evt = int(ev['Run']), int(ev['Evt'])
    res = dict(ev)
    cands = sam_list(f'defname: {defname} and run_number {run} and first_event <= {evt} and last_event >= {evt}')
    res['n_candidates'] = len(cands)
    for name in cands:
        url = xrootd_url(name)
        if not url:
            continue
        try:
            ids = read_ids(url)
        except Exception as err:
            res.setdefault('errors', '')
            res['errors'] += f'{name}: {err}; '
            continue
        for i, (r, s, e) in enumerate(ids):
            if r == run and e == evt:
                res.update(raw_file=name, raw_url=url, raw_entry=i, raw_subrun=s)
                s0 = children([name], 'stage0')
                # stage1 art files only (their children in SAM also include the CAFs)
                s1 = [c for c in children(s0, 'stage1') if c.endswith('-stage1.root')]
                res['stage0_files'] = ';'.join(s0)
                res['stage1_files'] = ';'.join(s1)
                res['stage0_urls'] = ';'.join(filter(None, (xrootd_url(x) for x in s0[:1])))
                res['stage1_urls'] = ';'.join(filter(None, (xrootd_url(x) for x in s1[:1])))
                return res
    return res


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--events', required=True, help='CSV with Run and Evt columns.')
    p.add_argument('--output', required=True, help='Output CSV (must not exist).')
    p.add_argument('--defname', default='run2_compression_production_v09_82_02_01_numimajority_compressed_data')
    p.add_argument('--workers', type=int, default=4)
    args = p.parse_args()
    if os.path.exists(args.output):
        sys.exit(f'[ERROR] {args.output} exists; pick a new name')
    with open(args.events) as fh:
        events = list(csv.DictReader(fh))
    print(f'[INFO] {len(events)} events, defname {args.defname}', flush=True)

    extra = ['n_candidates', 'raw_file', 'raw_url', 'raw_entry', 'raw_subrun',
             'stage0_files', 'stage0_urls', 'stage1_files', 'stage1_urls', 'errors']
    fields = list(events[0].keys()) + extra
    done = {}

    def write():  # rewritten after every event (input order), so a crash keeps what was found
        tmp = args.output + '.tmp'
        with open(tmp, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(done[i] for i in sorted(done))
        os.replace(tmp, args.output)

    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(find_one, ev, args.defname): i for i, ev in enumerate(events)}
        for fut in as_completed(futs):
            res = done[futs[fut]] = fut.result()
            write()
            n0 = len(res['stage0_files'].split(';')) if res.get('stage0_files') else 0
            n1 = len(res['stage1_files'].split(';')) if res.get('stage1_files') else 0
            print(f"[{len(done)}/{len(events)}] Run {res['Run']} Evt {res['Evt']} ({res.get('topology', '')}) "
                  f"-> {res.get('raw_file', 'NOT FOUND')}; stage0: {n0}, stage1: {n1}"
                  + (f"; ERRORS: {res['errors']}" if res.get('errors') else ''), flush=True)
    rows = [done[i] for i in sorted(done)]
    found = sum(1 for r in rows if r.get('raw_file'))
    print(f'[SUMMARY] raw file found for {found}/{len(rows)}; '
          f'with stage1: {sum(1 for r in rows if r.get("stage1_files"))}')


if __name__ == '__main__':
    main()
