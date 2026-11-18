#!/usr/bin/env python3
"""
Regression sandbox: run the medulla processing chain locally on a few pinned
input files, so the outputs of two code versions can be diffed with compare.py.

Run INSIDE the SL7 container (needs ROOT, hadd and the cvmfs build setup):

    scripts/in_sl7.sh 'python3 scripts/regression/run_chain.py prepare'
    scripts/in_sl7.sh 'python3 scripts/regression/run_chain.py run --ref WORKTREE --tier quick'
    scripts/in_sl7.sh 'python3 scripts/regression/run_chain.py run --ref 4b32a0a --tier full'

prepare  copies the pinned inputs (scripts/regression/inputs.toml) and the
         production spline files to <root>/inputs/ once and writes
         <root>/inputs/manifest.toml. Runs never read /pnfs.
run      --ref WORKTREE uses this checkout and its build/ (uncommitted changes
         included; the output dir gets a -dirty suffix if there are any).
         With --build DIR it compiles this checkout into DIR instead (make is
         rerun every time, incrementally) and never touches build/.
         --ref <git ref> checks the ref out as a detached git worktree under
         <root>/src_<sha> and builds it in <root>/build_<sha> (once; reused).
         The C++ binaries and add_numi_flux_sys.py come from that tree. The
         selection/systematics TOMLs come from THIS checkout by default, so two
         refs are compared on the same configuration (--configs-from-ref to
         use the ref's own TOMLs instead).

Tiers:
  quick  selection (medulla) + run_systematics for the quick inputs.
  full   every input, then: phase 1 spline check (one variation + cv, cut-down
         [variations]), phase 2 on the NuMIFull systematics output with the
         production splines, merge B, stage 5, stage 6, the g4 chain and stage 7.

Output: <root>/<sha>[-dirty]/<tier>/ (refuses to overwrite; --force to replace
a previous run of the same ref and tier). One log per step in logs/, and
timings in timings.toml.

accept   RUN_DIR --note "..." makes a finished run the reference of its tier,
         ONLY after the owner has approved its differences: copies it to
         <root>/reference/<tier>/<time>_<label>/, points reference/<tier>/current
         at it (old references are kept) and appends an entry with the note and
         the diff vs the previous reference to scripts/regression/CHANGELOG.md.
         --dry-run shows that diff and changes nothing. Check a new run with
         compare.py <root>/reference/<tier>/current <run>.

It never submits jobs and never writes to /pnfs; the only repo file it writes
is CHANGELOG.md (accept).
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import toml

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'batch'))
import utilities  # noqa: E402  (create_systematics_cfg, make_phase1_cfg, make_phase2_cfg)

DEFAULT_ROOT = Path('/exp/icarus/data/users/micarrig/regression')
NO_SYS_TAGS = utilities.NO_SYSTEMATICS_TAGS


def info(msg):
    print(f'[regression] {msg}', flush=True)


def fail(msg):
    print(f'[regression] ERROR: {msg}', flush=True)
    sys.exit(1)


def git(*args, cwd=REPO):
    return subprocess.run(['git', *args], cwd=cwd, check=True, capture_output=True, text=True).stdout.strip()


# --------------------------------------------------------------------------- prepare

def cmd_prepare(a):
    spec = toml.load(a.inputs)
    dest = a.root / 'inputs'
    manifest_path = dest / 'manifest.toml'
    if manifest_path.exists() and not a.add:
        fail(f'{manifest_path} exists; inputs are pinned. Use --add to pin only NEW labels, '
             'or a new --root to pin different ones.')
    if a.add and not manifest_path.exists():
        fail(f'--add needs an existing {manifest_path}')
    dest.mkdir(parents=True, exist_ok=True)
    manifest = toml.load(manifest_path) if a.add else {'input': [], 'chain': dict(spec['chain'])}
    pinned = {i['label'] for i in manifest['input']}
    for inp in spec['input']:
        if inp['label'] in pinned:
            continue  # never re-pin: existing files stay as they are
        cfg = toml.load(REPO / inp['toml'])
        sample = [s for s in cfg['sample'] if s['name'] == inp['sample']]
        if len(sample) != 1:
            fail(f"{inp['toml']}: sample {inp['sample']} not found exactly once")
        if 'files' in inp:  # explicit files (e.g. a known edge case)
            pattern, chosen = inp['files'], list(inp['files'])
        else:
            pattern = sample[0]['path']
            paths = sorted(p for pat in (pattern if isinstance(pattern, list) else [pattern]) for p in glob.glob(pat))
            if len(paths) < inp['nfiles']:
                fail(f"{inp['label']}: only {len(paths)} files match {pattern}")
            chosen = paths[:inp['nfiles']]
        out = dest / inp['label']
        out.mkdir(exist_ok=True)
        files = []
        for p in chosen:
            target = out / Path(p).name
            info(f"copy {p} -> {target}")
            # timeout: a single hung /pnfs read must not block forever
            subprocess.run(['timeout', '600', 'cp', p, str(target)], check=True)
            files.append(str(target))
        manifest['input'].append({**inp, 'source': pattern, 'files': files})
    if not a.add:
        for key in ('splines_main', 'splines_nue'):
            src = spec['chain'][key]
            target = dest / f'{key}.root'
            info(f'copy {src} -> {target}')
            subprocess.run(['timeout', '600', 'cp', src, str(target)], check=True)
            manifest['chain'][key] = str(target)
            manifest['chain'][key + '_source'] = src
    with open(manifest_path, 'w') as f:
        toml.dump(manifest, f)
    info(f'wrote {manifest_path}')


# --------------------------------------------------------------------------- build

class Tree:
    """Where the binaries and python scripts of one code version live."""

    def __init__(self, src: Path, build: Path):
        self.src, self.build = src, build
        self.medulla = build / 'selection' / 'medulla'
        self.run_systematics = build / 'systematics' / 'run_systematics'
        self.to_gundam = build / 'systematics' / 'to_gundam'
        self.flux_script = src / 'systematics' / 'add_numi_flux_sys.py'
        self.env = dict(os.environ)
        libs = ':'.join(str(build / d) for d in ('selection', 'systematics', 'shared'))
        self.env['LD_LIBRARY_PATH'] = libs + ':' + self.env.get('LD_LIBRARY_PATH', '')


def resolve_tree(a):
    if a.ref == 'WORKTREE':
        sha = git('rev-parse', '--short=12', 'HEAD')
        dirty = bool(git('status', '--porcelain', '--untracked-files=no'))
        label = sha + ('-dirty' if dirty else '')
        if a.build is None:
            tree = Tree(REPO, REPO / 'build')
            info(f'using this checkout ({label}) and {tree.build}; NOT rebuilding it')
            return label, tree
        # --build: compile this checkout (uncommitted changes included) into a
        # separate dir, leaving build/ alone. make is incremental, so it always runs.
        build = a.build.resolve()
        build.mkdir(parents=True, exist_ok=True)
        info(f'building this checkout ({label}) in {build} (log {build}/build.log)')
        with open(build / 'build.log', 'w') as log:
            if not (build / 'CMakeCache.txt').exists():
                subprocess.run(['cmake', '-S', str(REPO), '-B', str(build)], check=True, stdout=log, stderr=subprocess.STDOUT)
            subprocess.run(['make', '-C', str(build), f'-j{a.jobs}'], check=True, stdout=log, stderr=subprocess.STDOUT)
        return label, Tree(REPO, build)
    sha = git('rev-parse', '--short=12', a.ref + '^{commit}')
    src, build = a.root / f'src_{sha}', a.root / f'build_{sha}'
    if not src.exists():
        info(f'git worktree add --detach {src} {sha}')
        git('worktree', 'add', '--detach', str(src), sha)
    tree = Tree(src, build)
    if not (tree.medulla.exists() and tree.run_systematics.exists() and tree.to_gundam.exists()):
        build.mkdir(exist_ok=True)
        info(f'building {sha} in {build} (log {build}/build.log)')
        with open(build / 'build.log', 'w') as log:
            subprocess.run(['cmake', '-S', str(src), '-B', str(build)], check=True, stdout=log, stderr=subprocess.STDOUT)
            subprocess.run(['make', '-C', str(build), f'-j{a.jobs}'], check=True, stdout=log, stderr=subprocess.STDOUT)
    return sha, tree


# --------------------------------------------------------------------------- run

class Chain:
    def __init__(self, a, label, tree, manifest):
        self.a, self.tree, self.m, self.label = a, tree, manifest, label
        self.chain = manifest['chain']
        self.out = a.root / label / a.tier
        self.cfg_root = tree.src if a.configs_from_ref else REPO
        self.timings = {}

    def log_path(self, step):
        return self.out / 'logs' / f'{step}.log'

    def run(self, step, cmd, cwd, env=None):
        """Run one step, its output to logs/<step>.log; stop the chain on failure."""
        t0 = time.time()
        with open(self.log_path(step), 'w') as log:
            log.write(f'# cwd: {cwd}\n# cmd: {" ".join(map(str, cmd))}\n')
            log.flush()
            rc = subprocess.run([str(c) for c in cmd], cwd=cwd, env=env or self.tree.env,
                                stdout=log, stderr=subprocess.STDOUT).returncode
        dt = time.time() - t0
        self.timings[step] = round(dt, 1)
        info(f'{step}: rc={rc} ({dt:.0f} s)')
        if rc != 0:
            fail(f'{step} failed; see {self.log_path(step)}')

    # -- selection + systematics, as batch/submit.sh does for one job --------
    def job(self, inp):
        d = self.out / inp['label']
        (d / 'data').mkdir(parents=True)
        for f in inp['files']:
            os.symlink(f, d / 'data' / Path(f).name)
        cfg = toml.load(self.cfg_root / inp['toml'])
        cfg = utilities.resolve_samples(cfg)
        sample = dict([s for s in cfg['sample'] if s['name'] == inp['sample']][0])
        sample['path'] = [f'data/{Path(f).name}' for f in inp['files']]
        cfg['general']['output'] = 'output'
        cfg['sample'] = [sample]
        with open(d / 'job_config.toml', 'w') as f:
            toml.dump(cfg, f)
        self.run(f'{inp["label"]}_selection', [self.tree.medulla, 'job_config.toml'], d)
        if sample.get('tag') in NO_SYS_TAGS:
            return d
        # The project's per-job systematics config (create_new_project).
        sys_toml = {'main': 'sys_main', 'nue': 'sys_nue', 'g4': 'sys_g4'}.get(inp['role'])
        if sys_toml is None:
            return d
        scfg = utilities.create_systematics_cfg(toml.load(self.cfg_root / self.chain[sys_toml]),
                                                cfg.get('tree', []), [sample])
        with open(d / 'systematics.toml', 'w') as f:
            toml.dump(scfg, f)
        self.run(f'{inp["label"]}_systematics', [self.tree.run_systematics, 'systematics.toml'], d)
        return d

    # -- phase 1 spline check ------------------------------------------------
    def phase1(self, dirs):
        d = self.out / 'phase1'
        d.mkdir()
        keys = self.chain['spline_keys']
        self.run('phase1_merge', ['hadd', '-f', d / 'merged_stage0.root'] + [dirs[k] / 'output.root' for k in keys], d)
        cfg = utilities.make_phase1_cfg(self.cfg_root / self.chain['sys_main'], 'merged_stage0.root')
        cfg['variations']['keys'] = keys
        cfg['variations']['seed'] = self.chain['seed']
        cfg['sys'] = [s for s in cfg.get('sys', []) if s.get('type') != 'variation' or s['name'] == self.chain['spline_sys']]
        with open(d / 'phase1.toml', 'w') as f:
            toml.dump(cfg, f)
        self.run('phase1', [self.tree.run_systematics, 'phase1.toml'], d)

    # -- phase 2 with the production splines --------------------------------
    def phase2(self, numifull_dir):
        import uproot
        d = self.out / 'phase2'
        d.mkdir()
        cfg, _ = utilities.make_phase2_cfg(self.cfg_root / self.chain['sys_main'],
                                           input_path=str(numifull_dir / 'output_sys.root'),
                                           splines=self.chain['splines_main'], splines_nue=self.chain['splines_nue'])
        cfg['variations']['seed'] = self.chain['seed']
        # Drop [[sys]] variation blocks the production splines do not have
        # (e.g. a variation added after that production); phase 2 would stop.
        result = cfg['variations'].get('result_destination', 'detsys_results/').strip('/')
        var0 = cfg['variations'].get('variables', cfg['variations'].get('variable'))
        var0 = var0[0] if isinstance(var0, list) else var0
        have = None  # spline directories present in BOTH production files
        for key in ('splines_main', 'splines_nue'):
            with uproot.open(self.chain[key]) as f:
                names = {k.split(';')[0].split('/')[-1] for k in f[result].keys(recursive=False)}
            have = names if have is None else (have & names)
        kept = []
        for s in cfg.get('sys', []):
            if s.get('type') == 'variation' and f'{s["name"]}__{var0}_splines' not in have:
                info(f'phase 2: dropping [[sys]] {s["name"]} (not in the production splines)')
                continue
            kept.append(s)
        cfg['sys'] = kept
        with open(d / 'phase2.toml', 'w') as f:
            toml.dump(cfg, f)
        self.run('phase2', [self.tree.run_systematics, 'phase2.toml'], d)
        return d / 'output_varsys.root'

    def hadd(self, step, out, files):
        self.run(step, ['hadd', '-f', out] + list(files), out.parent)

    def stage5(self, step, merged, directory, trees):
        out = merged.with_name(merged.stem + '_flux.root')
        self.run(step, [sys.executable, self.tree.flux_script, merged, '--output', out, directory] + trees,
                 merged.parent)
        return out

    def stage6(self, step, sys_key, flux_file):
        d = self.out / f'{step}_tmp'
        d.mkdir()
        self.run(step, [self.tree.to_gundam, self.cfg_root / self.chain[sys_key], flux_file], d)
        return d / 'output.root'

    def go(self):
        if self.out.exists():
            if not self.a.force:
                fail(f'{self.out} exists; use --force to replace this ref+tier run')
            shutil.rmtree(self.out)
        (self.out / 'logs').mkdir(parents=True)
        info(f'output: {self.out}')
        t0 = time.time()
        todo = [i for i in self.m['input'] if self.a.tier in i['tiers']]
        dirs = {inp['label']: self.job(inp) for inp in todo}
        if self.a.tier == 'full':
            roles = {i['label']: i['role'] for i in todo}
            self.phase1(dirs)
            main = [k for k, r in roles.items() if r == 'main'][0]
            varsys = self.phase2(dirs[main])
            mb = self.out / 'mergeB'
            mb.mkdir()
            # merge_list.py's rule: systematics file for nominal, output file for
            # data/rock/detector_variation, plus the phase 2 varsys file.
            files = [dirs[main] / 'output_sys.root', varsys]
            files += [dirs[k] / 'output.root' for k, r in roles.items() if r in ('merge', 'spline')]
            self.hadd('mergeB', mb / 'merged_final.root', files)
            flux = self.stage5('stage5', mb / 'merged_final.root', 'events/NuMIFull', self.chain['flux_trees_main'])
            g_main = self.stage6('stage6', 'sys_main', flux)
            g4 = [k for k, r in roles.items() if r == 'g4'][0]
            mg = self.out / 'mergeB_g4'
            mg.mkdir()
            self.hadd('mergeB_g4', mg / 'merged_final.root', [dirs[g4] / 'output_sys.root'])
            flux_g4 = self.stage5('stage5_g4', mg / 'merged_final.root', 'events/g4', self.chain['flux_trees_g4'])
            g_g4 = self.stage6('stage6_g4', 'sys_g4', flux_g4)
            st7 = self.out / 'stage7'
            st7.mkdir()
            self.hadd('stage7', st7 / 'merged_gundam_final.root', [g_main, g_g4])
        self.timings['total'] = round(time.time() - t0, 1)
        with open(self.out / 'timings.toml', 'w') as f:
            toml.dump(self.timings, f)
        with open(self.out / 'run_info.toml', 'w') as f:
            toml.dump(dict(self.code_state(), ref=self.a.ref, label=self.label, src=str(self.tree.src),
                           build=str(self.tree.build), configs_from=str(self.cfg_root), tier=self.a.tier,
                           finished=time.strftime('%Y-%m-%d %H:%M:%S')), f)
        info(f'done in {self.timings["total"]:.0f} s: {self.out}')

    def code_state(self):
        """Commit of the code tree; for WORKTREE also its uncommitted changes,
        saved as code.diff (tracked files) and code_untracked.txt in the run
        dir, so an accepted reference says exactly which code made it."""
        src = self.tree.src
        state = {'sha': git('rev-parse', 'HEAD', cwd=src)}
        if self.a.ref == 'WORKTREE':
            diff = git('diff', 'HEAD', cwd=src)
            untracked = git('ls-files', '--others', '--exclude-standard', cwd=src)
            (self.out / 'code.diff').write_text(diff + '\n' if diff else '')
            (self.out / 'code_untracked.txt').write_text(untracked + '\n' if untracked else '')
            state['dirty_files'] = git('diff', 'HEAD', '--name-only', cwd=src).splitlines()
        return state


def cmd_run(a):
    manifest_path = a.root / 'inputs' / 'manifest.toml'
    if not manifest_path.exists():
        fail(f'{manifest_path} not found; run `run_chain.py prepare` first')
    if shutil.which('hadd') is None:
        fail('hadd not on PATH: run inside the SL7 container (scripts/in_sl7.sh)')
    manifest = toml.load(manifest_path)
    # The manifest pins the FILES; tiers and roles come from inputs.toml so
    # they can be changed without re-pinning.
    spec = {i['label']: i for i in toml.load(Path(__file__).with_name('inputs.toml'))['input']}
    for inp in manifest['input']:
        if inp['label'] in spec:
            inp['tiers'], inp['role'] = spec[inp['label']]['tiers'], spec[inp['label']]['role']
    label, tree = resolve_tree(a)
    Chain(a, label + (f'-{a.suffix}' if a.suffix else ''), tree, manifest).go()


# --------------------------------------------------------------------------- accept

CHANGELOG = Path(__file__).with_name('CHANGELOG.md')
MAX_LOGGED_DIFF_LINES = 40


def condense(diff_lines):
    """compare.py's '  DIFF <path>: <message>' lines, with a message that recurs
    across paths (e.g. 'branch X only in B' in every tree) folded into one
    '<n> x <message>' line; one-off differences are kept whole."""
    groups = {}
    for line in diff_lines:
        path, _, msg = line[len('  DIFF '):].partition(': ')
        groups.setdefault(msg or path, []).append(line.strip())
    out = []
    for msg, lines in groups.items():
        out += [f'{len(lines)} x {msg}'] if len(lines) > 1 else lines
    return out


def reference_dir(root: Path, tier: str) -> Path:
    return root / 'reference' / tier


def cmd_accept(a):
    """Make a finished run the reference of its tier, once the owner has
    approved its differences. Copies the run (so later --force reruns cannot
    change it), repoints reference/<tier>/current, and appends an entry with
    the note and the diff against the previous reference to CHANGELOG.md.
    Old references are kept."""
    import contextlib
    import io
    sys.path.insert(0, str(Path(__file__).parent))
    import compare

    run_dir = a.run_dir.resolve()
    info_path = run_dir / 'run_info.toml'
    if not info_path.exists():
        fail(f'{info_path} not found: not a finished run_chain.py run')
    ref_root = (a.root / 'reference').resolve()
    if ref_root in run_dir.parents:
        fail(f'{run_dir} is already under {ref_root}')
    rinfo = toml.load(info_path)
    tier = rinfo['tier']
    label = rinfo.get('label', run_dir.parent.name)
    tdir = reference_dir(a.root, tier)
    current = tdir / 'current'
    prev = current.resolve() if current.exists() else None

    if prev is None:
        summary, diff_lines = 'none (first reference of this tier)', []
    else:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = compare.main([str(prev), str(run_dir)])
        out = buf.getvalue().splitlines()
        summary = f'vs {prev.name}: ' + (out[-1].strip() if out else f'compare exit {rc}')
        diff_lines = condense([l for l in out if l.startswith('  DIFF ')])

    stamp = time.strftime('%Y%m%d-%H%M%S')
    dest = tdir / f'{stamp}_{label}'
    if a.dry_run:
        info(f'[dry-run] would copy {run_dir} -> {dest}, point {current} at it, and log to {CHANGELOG}')
        info(f'[dry-run] {summary}')
        for l in diff_lines[:MAX_LOGGED_DIFF_LINES]:
            print(l)
        return
    tdir.mkdir(parents=True, exist_ok=True)
    shutil.copytree(run_dir, dest, symlinks=True)
    with open(dest / 'accepted.toml', 'w') as f:
        toml.dump({'accepted': time.strftime('%Y-%m-%d %H:%M:%S'), 'from': str(run_dir),
                   'note': a.note, 'previous': str(prev) if prev else '', 'summary': summary}, f)
    tmp = tdir / f'.current.{os.getpid()}'
    os.symlink(dest.name, tmp)
    os.replace(tmp, current)

    # runs from before 2026-11-18 have no sha; a git-ref run is still identified by its ref
    sha = rinfo.get('sha') or (rinfo['ref'] if rinfo.get('ref') not in (None, 'WORKTREE') else '')
    code = f'{sha[:12]}' if sha else '(not recorded)'
    if rinfo.get('dirty_files'):
        code += f' + uncommitted changes to {len(rinfo["dirty_files"])} files (code.diff, code_untracked.txt in the reference)'
    entry = [f'\n## {time.strftime("%Y-%m-%d %H:%M")} | {tier} | {label}\n\n',
             f'- Reference: `{dest}`\n',
             f'- Code: {code}; configs from `{rinfo.get("configs_from", "?")}`\n',
             f'- Approved: {a.note}\n',
             f'- Diff {summary}\n']
    if diff_lines:
        shown = diff_lines[:MAX_LOGGED_DIFF_LINES]
        more = len(diff_lines) - len(shown)
        entry.append('\n```\n' + '\n'.join(shown) + (f'\n... {more} more lines' if more else '') + '\n```\n')
    if not CHANGELOG.exists():
        CHANGELOG.write_text(CHANGELOG_HEADER)
    with open(CHANGELOG, 'a') as f:
        f.writelines(entry)
    info(f'accepted {run_dir} as the {tier} reference: {dest}')
    info(f'logged to {CHANGELOG}')


CHANGELOG_HEADER = """# Regression sandbox reference log

One entry per accepted reference run (`run_chain.py accept`), newest last. A run
is accepted only after the owner has looked at its differences from the previous
reference and said they are intended. The reference itself is
`/exp/icarus/data/users/micarrig/regression/reference/<tier>/current`; compare new
runs against it with `compare.py`.
"""


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--root', type=Path, default=DEFAULT_ROOT, help=f'sandbox root (default {DEFAULT_ROOT})')
    sub = p.add_subparsers(dest='cmd', required=True)
    pp = sub.add_parser('prepare', help='copy the pinned inputs and write the manifest')
    pp.add_argument('--inputs', type=Path, default=Path(__file__).with_name('inputs.toml'))
    pp.add_argument('--add', action='store_true', help='pin only labels not yet in the manifest')
    pr = sub.add_parser('run', help='run the chain for one code version')
    pr.add_argument('--ref', default='WORKTREE', help='git ref to build and run, or WORKTREE (default)')
    pr.add_argument('--tier', choices=['quick', 'full'], default='quick')
    pr.add_argument('--jobs', type=int, default=4, help='make -j for a ref build')
    pr.add_argument('--build', type=Path,
                    help='with --ref WORKTREE: compile this checkout into DIR (incremental) instead of using build/')
    pr.add_argument('--configs-from-ref', action='store_true',
                    help="use the ref's own selection/systematics TOMLs instead of this checkout's")
    pr.add_argument('--force', action='store_true', help='replace an existing run of the same ref and tier')
    pr.add_argument('--suffix', help='appended to the output label, e.g. "rep2" for a repeat run of the same ref')
    pa = sub.add_parser('accept', help="make a run its tier's reference (only after the owner approves its diffs)")
    pa.add_argument('run_dir', type=Path, help='finished run, e.g. <root>/<label>/full')
    pa.add_argument('--note', required=True, help='what changed and who approved it (goes into CHANGELOG.md)')
    pa.add_argument('--dry-run', action='store_true', help='show the diff vs the current reference, change nothing')
    a = p.parse_args()
    a.root = a.root.resolve()
    {'prepare': cmd_prepare, 'run': cmd_run, 'accept': cmd_accept}[a.cmd](a)


if __name__ == '__main__':
    main()
