"""Campaign management CLI for medulla batch processing."""
import argparse
import contextlib
import re
import shutil
import sqlite3
import sys
import tempfile
import toml
from dataclasses import dataclass, field
from datetime import datetime, timezone
from glob import glob
from pathlib import Path

from auth import authenticate
from catalog import resolve_samples
from utilities import create_new_project, check_project_status, launch_jobsub

# Repo root is two levels above this script (batch/ -> repo root).
REPO_ROOT = Path(__file__).resolve().parent.parent
TOML_DIR = REPO_ROOT / 'selection' / 'toml'


# ---------------------------------------------------------------------------
# Terminal formatting helpers
# ---------------------------------------------------------------------------

class _A:
    """ANSI escape codes (no third-party dependency needed on Linux)."""
    RESET   = '\033[0m'
    BOLD    = '\033[1m'
    DIM     = '\033[2m'
    CYAN    = '\033[96m'   # headers / created
    GREEN   = '\033[92m'   # completed
    YELLOW  = '\033[93m'   # pending
    BLUE    = '\033[94m'   # submitted
    RED     = '\033[91m'   # failed
    MAGENTA = '\033[95m'   # analysis names
    GREY    = '\033[90m'   # rules / dim chrome

_STATUS_COLOR = {
    'created':   _A.CYAN,
    'pending':   _A.YELLOW,
    'partial':   _A.YELLOW,
    'submitted': _A.BLUE,
    'completed': _A.GREEN,
    'failed':    _A.RED,
}

_ANSI_RE = re.compile(r'\033\[[0-9;]*m')


def _trunc(s, width):
    """Truncate *s* to *width* visible characters, appending '…' if needed."""
    s = str(s)
    return s if len(s) <= width else s[:width - 1] + '…'


def _cell(s, width, align='<', color=None):
    """Return a fixed-width, truncated cell, then optionally ANSI-colored."""
    padded = f"{_trunc(str(s), width):{align}{width}}"
    return f"{color}{padded}{_A.RESET}" if color else padded


def _print_table(headers, widths, rows, sep='  '):
    """
    Print a fixed-width table with a bold header and a rule beneath it.

    Parameters
    ----------
    headers : list[str]
    widths  : list[int]
    rows    : list[list]
        Each cell may be a plain string or a pre-colored string produced by
        _cell().  Pre-colored strings must already be padded to the correct
        visible width.
    sep : str
        Column separator (default: two spaces).
    """
    rule = sep.join(_A.GREY + '─' * w + _A.RESET for w in widths)
    header = sep.join(
        f"{_A.BOLD}{_trunc(h, w):<{w}}{_A.RESET}" for h, w in zip(headers, widths)
    )
    print(header)
    print(rule)
    for row in rows:
        print(sep.join(str(c) for c in row))


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class TomlEntry:
    """One [[toml]] block from a meta.toml: role, file path, experiments, and
    per-experiment enable-key lists."""
    role: str
    file: str                    # absolute path to the selection TOML
    experiments: list
    enable: dict = field(default_factory=dict)  # {experiment: [key, ...]}


@dataclass
class AnalysisMeta:
    """Parsed contents of one meta.toml file."""
    analysis: str
    description: str
    owners: list
    experiments: list
    tomls: list                  # list[TomlEntry]
    defaults: dict = field(default_factory=dict)


@dataclass
class ProjectUnit:
    """One (analysis, role, experiment) combination ready for batch creation."""
    analysis: str
    role: str
    experiment: str
    toml_path: str
    enable_keys: list = field(default_factory=list)
    batch_size: int = 50


# ---------------------------------------------------------------------------
# Database helpers
# ---------------------------------------------------------------------------

SCHEMA_CAMPAIGN_META = """
CREATE TABLE IF NOT EXISTS campaign_meta (
    name TEXT NOT NULL,
    tag  TEXT NOT NULL,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
"""

SCHEMA_PROJECTS = """
CREATE TABLE IF NOT EXISTS projects (
    project_id INTEGER PRIMARY KEY AUTOINCREMENT,
    analysis TEXT NOT NULL,
    role TEXT NOT NULL,
    experiment TEXT NOT NULL,
    toml_file TEXT NOT NULL,
    project_dir TEXT NOT NULL,
    batch_size INTEGER NOT NULL,
    n_jobs INTEGER DEFAULT 0,
    n_completed INTEGER DEFAULT 0,
    status TEXT DEFAULT 'created',
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    submitted_at TIMESTAMP,
    completed_at TIMESTAMP,
    UNIQUE(analysis, role, experiment)
);
"""


@contextlib.contextmanager
def _open_db(campaign_dir):
    """Context manager: copy campaign.db to a local temp file, yield (conn,
    cursor), then copy back.  Avoids /pnfs SQLite locking failures."""
    db_path = Path(campaign_dir) / 'campaign.db'
    if not db_path.exists():
        raise FileNotFoundError(f"Campaign database not found: {db_path}")
    with tempfile.NamedTemporaryFile(
        suffix='.db', prefix='medulla_campaign_', delete=False
    ) as f:
        tmp = Path(f.name)
    shutil.copy2(db_path, tmp)
    conn = sqlite3.connect(tmp)
    conn.row_factory = sqlite3.Row
    curs = conn.cursor()
    # Migration: add n_completed column if this DB was created before Stage 5.
    curs.execute("PRAGMA table_info(projects)")
    if 'n_completed' not in {row[1] for row in curs.fetchall()}:
        curs.execute(
            "ALTER TABLE projects ADD COLUMN n_completed INTEGER DEFAULT 0"
        )
        conn.commit()
    try:
        yield conn, curs
    finally:
        conn.close()
        # /pnfs does not allow overwriting files in place; delete first.
        db_path.unlink(missing_ok=True)
        shutil.copy2(tmp, db_path)
        tmp.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Discovery helpers
# ---------------------------------------------------------------------------

def discover_analyses(toml_root):
    """
    Walk toml_root/*/ looking for meta.toml files.
    Return one AnalysisMeta per directory that has one.

    Parameters
    ----------
    toml_root : str | Path
        Root directory to search (e.g. selection/toml/).

    Returns
    -------
    list[AnalysisMeta]
    """
    analyses = []
    for meta_path in sorted(Path(toml_root).glob('*/meta.toml')):
        meta = toml.load(meta_path)
        top_experiments = meta['meta'].get('experiments', [])
        tomls = []
        for t in meta.get('toml', []):
            enable_raw = t.get('enable', {})
            enable = {exp: block.get('keys', []) for exp, block in enable_raw.items()}
            tomls.append(TomlEntry(
                role=t['role'],
                file=str(meta_path.parent / t['file']),
                experiments=t.get('experiments', top_experiments),
                enable=enable,
            ))
        analyses.append(AnalysisMeta(
            analysis=meta['meta']['analysis'],
            description=meta['meta'].get('description', ''),
            owners=meta['meta'].get('owners', []),
            experiments=top_experiments,
            tomls=tomls,
            defaults=meta.get('defaults', {}),
        ))
    return analyses


def expand_campaign(analyses, catalog_path,
                    experiment=None, roles=None, analysis_filter=None):
    """
    Expand a list of AnalysisMeta into one ProjectUnit per
    (analysis, role, experiment) combination, filtered by the caller's
    constraints.

    Parameters
    ----------
    analyses : list[AnalysisMeta]
    catalog_path : str | Path
        Path to the sample catalog (passed through; used by create_campaign).
    experiment : str | None
        Keep only units for this experiment.
    roles : list[str] | None
        Keep only units whose role is in this list.
    analysis_filter : list[str] | None
        Keep only units whose analysis name is in this list.

    Returns
    -------
    list[ProjectUnit]
    """
    units = []
    for a in analyses:
        if analysis_filter and a.analysis not in analysis_filter:
            continue
        for t in a.tomls:
            if roles and t.role not in roles:
                continue
            for exp in t.experiments:
                if experiment and exp != experiment:
                    continue
                units.append(ProjectUnit(
                    analysis=a.analysis,
                    role=t.role,
                    experiment=exp,
                    toml_path=t.file,
                    enable_keys=t.enable.get(exp, []),
                    batch_size=a.defaults.get('batch_size', 50),
                ))
    return units


def create_campaign(campaign_dir, project_units, catalog_path,
                    name, tag,
                    batch_size_override=None, campaign_cfg=None):
    """
    Create a campaign directory with campaign.db and one project sub-directory
    per ProjectUnit.

    Parameters
    ----------
    campaign_dir : str | Path
        Directory to create.  Raises FileExistsError if campaign.db already
        exists there.
    project_units : list[ProjectUnit]
    catalog_path : str | Path
        Sample catalog passed to create_new_project for sample resolution.
    name : str
        Human-readable campaign name stored in campaign.db.
    tag : str
        Git ref / version tag stored in campaign.db.
    batch_size_override : int | None
        When set, overrides every project's batch size.
    campaign_cfg : dict | None
        Optional config dict.  ``campaign_cfg["overrides"]`` is a list of
        dicts with keys (analysis, role, experiment, batch_size) that provide
        per-project batch size overrides, taking precedence over
        batch_size_override.

    Raises
    ------
    FileExistsError
        If campaign.db already exists in campaign_dir.
    """
    campaign_dir = Path(campaign_dir)
    db_path = campaign_dir / 'campaign.db'
    if db_path.exists():
        raise FileExistsError(f"Campaign database already exists: {db_path}")

    # Build the entire campaign tree in a local staging directory to avoid
    # SQLite locking failures on /pnfs, then copy everything to the final
    # destination in one shot.
    staging = Path(tempfile.mkdtemp(prefix='medulla_campaign_'))
    try:
        conn = sqlite3.connect(staging / 'campaign.db')
        conn.row_factory = sqlite3.Row
        curs = conn.cursor()
        curs.executescript(SCHEMA_CAMPAIGN_META + SCHEMA_PROJECTS)
        curs.execute("INSERT INTO campaign_meta (name, tag) VALUES (?, ?)", (name, tag))
        conn.commit()

        # Build per-project batch_size override lookup from campaign_cfg.
        cfg_overrides = {}
        if campaign_cfg:
            for ov in campaign_cfg.get('overrides', []):
                key = (ov['analysis'], ov['role'], ov['experiment'])
                cfg_overrides[key] = ov

        for u in project_units:
            # Record the final destination path in the DB so submit.sh finds it.
            proj_dir_final = campaign_dir / f"{u.analysis}_{u.role}_{u.experiment}"
            proj_dir_stage = staging   / f"{u.analysis}_{u.role}_{u.experiment}"

            # Precedence: per-project cfg override > global override > unit default.
            batch_size = u.batch_size
            if batch_size_override is not None:
                batch_size = batch_size_override
            ov = cfg_overrides.get((u.analysis, u.role, u.experiment), {})
            if 'batch_size' in ov:
                batch_size = ov['batch_size']

            create_new_project(
                project_dir=proj_dir_stage,
                tml=u.toml_path,
                batch_size=batch_size,
                catalog_path=catalog_path,
                enable_keys=u.enable_keys,
            )

            # Read back the job count from the locally-built project.db.
            proj_conn = sqlite3.connect(proj_dir_stage / 'project.db')
            (n_jobs,) = proj_conn.execute("SELECT COUNT(*) FROM jobs").fetchone()
            proj_conn.close()

            curs.execute(
                "INSERT INTO projects "
                "(analysis, role, experiment, toml_file, project_dir, batch_size, n_jobs) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                (u.analysis, u.role, u.experiment,
                 u.toml_path, str(proj_dir_final), batch_size, n_jobs),
            )
            conn.commit()

        conn.close()

        # Copy the complete staging tree to the final destination.
        campaign_dir.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(staging, campaign_dir)
    finally:
        shutil.rmtree(staging, ignore_errors=True)


# ---------------------------------------------------------------------------
# Sync helpers
# ---------------------------------------------------------------------------

def _sync_project_status(project_dir):
    """
    Inspect output files for a single project, update its project.db, and
    return completion counts.

    A job is considered complete when its output file exists and is at least
    1 KB — the same criterion used by utilities.check_project_status.

    Parameters
    ----------
    project_dir : str | Path

    Returns
    -------
    dict with keys ``n_jobs`` and ``n_completed``, or None if project.db
    does not exist.
    """
    project_dir = Path(project_dir)
    db_path = project_dir / 'project.db'
    if not db_path.exists():
        return None

    output_files = glob(str(project_dir / 'output' / 'output_jobid*.root'))
    completed_ids = [
        int(Path(f).stem.split('jobid')[-1])
        for f in output_files
        if Path(f).stat().st_size >= 1024
    ]

    # Copy to a local temp file to avoid /pnfs locking failures.
    with tempfile.NamedTemporaryFile(
        suffix='.db', prefix='medulla_proj_', delete=False
    ) as f:
        tmp = Path(f.name)
    try:
        shutil.copy2(db_path, tmp)
        conn = sqlite3.connect(tmp)
        curs = conn.cursor()
        if completed_ids:
            curs.executemany(
                "UPDATE jobs SET status = 'completed' WHERE jobid = ?",
                [(jid,) for jid in completed_ids],
            )
            conn.commit()
        curs.execute("SELECT status, COUNT(*) FROM jobs GROUP BY status")
        counts = dict(curs.fetchall())
        conn.close()
        # /pnfs does not allow overwriting files in place; delete first.
        db_path.unlink(missing_ok=True)
        shutil.copy2(tmp, db_path)
    finally:
        tmp.unlink(missing_ok=True)

    return {
        'n_jobs':      sum(counts.values()),
        'n_completed': counts.get('completed', 0),
    }


# ---------------------------------------------------------------------------
# `sync` subcommand
# ---------------------------------------------------------------------------

def cmd_sync(args):
    """Sync per-project job status into the campaign database."""
    campaign_dir = Path(args.campaign)
    W_AN, W_RO, W_EX, W_ST, W_JB = 28, 18, 10, 12, 10
    synced = 0

    with _open_db(campaign_dir) as (conn, curs):
        curs.execute(
            "SELECT project_id, analysis, role, experiment, project_dir, status "
            "FROM projects ORDER BY analysis, role, experiment"
        )
        rows = list(curs.fetchall())

        if not rows:
            print("[SYNC] No projects in this campaign.")
            return

        _print_table(
            ['Analysis', 'Role', 'Experiment', 'Status', 'Jobs'],
            [W_AN, W_RO, W_EX, W_ST, W_JB],
            [],
        )

        for row in rows:
            result = _sync_project_status(Path(row['project_dir']))
            if result is None:
                continue

            n_jobs      = result['n_jobs']
            n_completed = result['n_completed']

            if n_jobs == 0:
                new_status = row['status']
            elif n_completed == n_jobs:
                new_status = 'completed'
            elif n_completed > 0:
                new_status = 'partial'
            else:
                new_status = row['status']

            curs.execute(
                "UPDATE projects "
                "SET n_jobs = ?, n_completed = ?, status = ? "
                "WHERE project_id = ?",
                (n_jobs, n_completed, new_status, row['project_id']),
            )
            conn.commit()
            synced += 1

            job_str = f"{n_completed}/{n_jobs}" if n_jobs else '—'
            job_color = (
                _A.GREEN  if n_jobs and n_completed == n_jobs else
                _A.YELLOW if n_completed > 0                  else
                _A.RED    if n_jobs                           else None
            )
            print('  '.join([
                _cell(row['analysis'],  W_AN, color=_A.MAGENTA),
                _cell(row['role'],      W_RO),
                _cell(row['experiment'],W_EX),
                _cell(new_status,       W_ST, color=_STATUS_COLOR.get(new_status)),
                _cell(job_str,          W_JB, color=job_color),
            ]))

    print(f"\n[SYNC] Synced {synced} project(s).")


def _discover_meta_files():
    """Return a list of all meta.toml paths under TOML_DIR."""
    return sorted(TOML_DIR.glob('*/meta.toml'))


def _expand_projects(meta_path, filters_exp=None, filters_roles=None):
    """
    Expand a meta.toml into a list of (analysis, role, experiment, toml_file,
    enable_keys, batch_size) tuples.
    """
    meta = toml.load(meta_path)
    analysis = meta['meta']['analysis']
    experiments = meta['meta'].get('experiments', [])
    batch_size = meta.get('defaults', {}).get('batch_size', 50)

    units = []
    for toml_entry in meta.get('toml', []):
        role = toml_entry['role']
        toml_file = meta_path.parent / toml_entry['file']
        # Roles supported by this entry (may override top-level experiments).
        entry_exps = toml_entry.get('experiments', experiments)

        for exp in entry_exps:
            if filters_exp and exp not in filters_exp:
                continue
            if filters_roles and role not in filters_roles:
                continue

            enable_block = toml_entry.get('enable', {}).get(exp, {})
            enable_keys = enable_block.get('keys', [])
            units.append({
                'analysis': analysis,
                'role': role,
                'experiment': exp,
                'toml_file': str(toml_file),
                'enable_keys': enable_keys,
                'batch_size': batch_size,
            })
    return units


# ---------------------------------------------------------------------------
# `create` subcommand
# ---------------------------------------------------------------------------

def cmd_create(args):
    """Discover analyses, expand project units, and create the campaign."""
    output_base = Path(args.output) if args.output else None

    filters_exp      = set(args.experiment) if args.experiment else None
    filters_roles    = list(args.roles) if args.roles else None
    filters_analyses = list(args.analyses) if args.analyses else None

    # Discover and expand all project units.
    catalog_path = TOML_DIR / 'common' / 'samples.toml'
    analyses     = discover_analyses(TOML_DIR)
    # expand_campaign takes a single experiment string; handle multi-experiment
    # filters by expanding without the filter then pruning the result.
    single_exp = next(iter(filters_exp)) if filters_exp and len(filters_exp) == 1 else None
    all_units  = expand_campaign(
        analyses, catalog_path,
        experiment=single_exp,
        roles=filters_roles,
        analysis_filter=filters_analyses,
    )
    if filters_exp and len(filters_exp) > 1:
        all_units = [u for u in all_units if u.experiment in filters_exp]

    if not all_units:
        print("[CAMPAIGN] No project units matched the given filters.")
        return

    # Print summary table.
    W_AN, W_RO, W_EX = 28, 18, 10
    _print_table(
        ['Analysis', 'Role', 'Experiment'],
        [W_AN, W_RO, W_EX],
        [
            [
                _cell(u.analysis, W_AN, color=_A.MAGENTA),
                _cell(u.role, W_RO),
                _cell(u.experiment, W_EX),
            ]
            for u in all_units
        ],
    )
    print(f"\n[CAMPAIGN] Total: {_A.CYAN}{len(all_units)}{_A.RESET} project unit(s).")

    if args.dry_run:
        print("[CAMPAIGN] Dry-run mode: no directories or databases created.")
        return

    if output_base is None:
        print("[CAMPAIGN] Error: --output is required when not using --dry-run.")
        sys.exit(1)

    # Confirm with the user before creating anything.
    resp = input("\n[CAMPAIGN] Proceed with campaign creation? [Y/N] ")
    if resp.strip().lower() != 'y':
        print("[CAMPAIGN] Aborted.")
        return

    tag = args.tag
    ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')
    campaign_name = f"campaign_{tag}_{ts}"
    campaign_dir = output_base / campaign_name

    manifest_cfg = toml.load(args.manifest) if args.manifest else None

    create_campaign(
        campaign_dir=campaign_dir,
        project_units=all_units,
        catalog_path=catalog_path,
        name=campaign_name,
        tag=tag,
        campaign_cfg=manifest_cfg,
    )

    # Write a CLI-level manifest snapshot alongside campaign.db.
    ts_str = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')
    manifest_snapshot = {
        'campaign': {'name': campaign_name, 'tag': tag, 'created_at': ts_str},
        'projects': [
            {'analysis': u.analysis, 'role': u.role, 'experiment': u.experiment,
             'toml_path': u.toml_path, 'batch_size': u.batch_size}
            for u in all_units
        ],
    }
    with open(campaign_dir / 'campaign_manifest.toml', 'w') as f:
        toml.dump(manifest_snapshot, f)

    print(f"\n[CAMPAIGN] Campaign '{campaign_name}' created at {campaign_dir}")
    print(f"[CAMPAIGN] {len(all_units)} project(s) created.")


# ---------------------------------------------------------------------------
# `status` subcommand
# ---------------------------------------------------------------------------

def cmd_status(args):
    """Print aggregated status for all projects in the campaign."""
    with _open_db(args.campaign) as (conn, curs):
        curs.execute("SELECT * FROM projects ORDER BY analysis, role, experiment")
        rows = list(curs.fetchall())

    if not rows:
        print("[CAMPAIGN] No projects registered in this campaign.")
        return

    W_AN, W_RO, W_EX, W_ST, W_JB = 28, 18, 10, 12, 10
    _print_table(
        ['Analysis', 'Role', 'Experiment', 'Status', 'Jobs'],
        [W_AN, W_RO, W_EX, W_ST, W_JB],
        [
            [
                _cell(row['analysis'],   W_AN, color=_A.MAGENTA),
                _cell(row['role'],       W_RO),
                _cell(row['experiment'], W_EX),
                _cell(row['status'],     W_ST, color=_STATUS_COLOR.get(row['status'])),
                _cell(
                    f"{row['n_completed']}/{row['n_jobs']}" if row['n_jobs'] else '—',
                    W_JB,
                    color=(
                        _A.GREEN  if row['n_jobs'] and row['n_completed'] == row['n_jobs'] else
                        _A.YELLOW if row['n_completed']                                     else
                        None
                    ),
                ),
            ]
            for row in rows
        ],
    )
    print(f"\n[CAMPAIGN] {len(rows)} project(s) total.")


# ---------------------------------------------------------------------------
# `list` subcommand
# ---------------------------------------------------------------------------

def cmd_list(args):
    """Discover and display all analyses and their expanded project units."""
    toml_root    = Path(args.toml_root) if args.toml_root else TOML_DIR
    catalog_path = TOML_DIR / 'common' / 'samples.toml'
    analyses     = discover_analyses(toml_root)

    if not analyses:
        print(f"[LIST] No analyses found under {toml_root}")
        return

    # Column widths for the per-analysis role table.
    W_ROLE = 18
    W_EXPS = 16
    W_KEYS = 50

    # Column widths for the full expansion table.
    W_AN = 28
    W_RO = 18
    W_EX = 10

    rule_width = W_ROLE + W_EXPS + W_KEYS + 4  # +4 for two sep pairs
    rule = _A.GREY + '─' * rule_width + _A.RESET

    print(f"\n{_A.BOLD}Analyses under:{_A.RESET} {toml_root}")
    print(rule)

    for a in analyses:
        owners = ', '.join(a.owners) if a.owners else '—'
        batch  = a.defaults.get('batch_size', '(default)')
        print(
            f"\n  {_A.MAGENTA}{_A.BOLD}{a.analysis}{_A.RESET}"
            f"  {_A.DIM}owners: {owners}  batch_size: {batch}{_A.RESET}"
        )
        _print_table(
            ['Role', 'Experiments', 'Enable keys'],
            [W_ROLE, W_EXPS, W_KEYS],
            [
                [
                    _cell(t.role, W_ROLE),
                    _cell(', '.join(t.experiments), W_EXPS),
                    _cell(
                        '  '.join(
                            f"{exp}: [{', '.join(keys)}]"
                            for exp, keys in t.enable.items()
                        ) or '(none)',
                        W_KEYS,
                    ),
                ]
                for t in a.tomls
            ],
            sep='  ',
        )

    # ── Full expansion table ──────────────────────────────────────────────
    units = expand_campaign(analyses, catalog_path)
    print(f"\n{rule}")
    print(
        f"  {_A.BOLD}Full expansion:{_A.RESET} "
        f"{_A.CYAN}{len(units)}{_A.RESET} project unit(s)\n"
    )
    _print_table(
        ['Analysis', 'Role', 'Experiment'],
        [W_AN, W_RO, W_EX],
        [
            [
                _cell(u.analysis, W_AN, color=_A.MAGENTA),
                _cell(u.role, W_RO),
                _cell(u.experiment, W_EX),
            ]
            for u in units
        ],
    )
    print()


# ---------------------------------------------------------------------------
# `launch` subcommand
# ---------------------------------------------------------------------------

def cmd_launch(args):
    """Authenticate per experiment and launch pending projects."""
    W_AN, W_RO, W_EX = 28, 18, 10

    # Resolve how many jobs to submit per project.
    if args.test:
        njobs = 1
        mode_label = " (test: 1 job per project)"
    elif args.njobs is not None:
        njobs = args.njobs
        mode_label = f" (--njobs {njobs} per project)"
    else:
        njobs = -1  # launch_jobsub default: all pending
        mode_label = ""

    with _open_db(args.campaign) as (conn, curs):
        # Build query (optionally filtered by experiment).
        # --relaunch also picks up 'submitted' projects (e.g. after an auth failure
        # that caused jobsub to silently succeed on the campaign side but not the grid).
        statuses = "('created', 'partial', 'submitted')" if args.relaunch else "('created', 'partial')"
        if args.experiment:
            curs.execute(
                f"SELECT * FROM projects WHERE experiment = ? AND status IN {statuses}",
                (args.experiment,),
            )
        else:
            curs.execute(f"SELECT * FROM projects WHERE status IN {statuses}")

        rows = list(curs.fetchall())

        if not rows:
            print("[CAMPAIGN] No projects ready for launch.")
            return

        # Show what will be launched and confirm once.
        _print_table(
            ['Analysis', 'Role', 'Experiment'],
            [W_AN, W_RO, W_EX],
            [
                [
                    _cell(row['analysis'],   W_AN, color=_A.MAGENTA),
                    _cell(row['role'],       W_RO),
                    _cell(row['experiment'], W_EX),
                ]
                for row in rows
            ],
        )
        print(f"\n[CAMPAIGN] {len(rows)} project(s) to launch{mode_label}.")
        resp = input("\n[CAMPAIGN] Proceed with job submission? [Y/N] ")
        if resp.strip().lower() != 'y':
            print("[CAMPAIGN] Aborted.")
            return

        # Group by experiment so we authenticate once per experiment.
        by_exp = {}
        for row in rows:
            by_exp.setdefault(row['experiment'], []).append(row)

        for exp, projects in by_exp.items():
            print(f"\n[CAMPAIGN] Authenticating for experiment: {exp}")
            ok = authenticate(exp)
            if not ok:
                print(f"[CAMPAIGN] Skipping {exp} (authentication declined).")
                continue

            for row in projects:
                proj_dir = Path(row['project_dir'])
                print(f"[CAMPAIGN] Launching: {row['analysis']}/{row['role']}_{row['experiment']}")
                try:
                    ok = launch_jobsub(proj_dir, exp=exp, njobs=njobs, confirm=False)
                except Exception as e:
                    print(f"[CAMPAIGN] Launch failed for {proj_dir}: {e}")
                    ok = False
                if ok:
                    curs.execute(
                        "UPDATE projects SET status = 'submitted', submitted_at = ? WHERE project_id = ?",
                        (datetime.now(timezone.utc).isoformat(), row['project_id']),
                    )
                    conn.commit()

    print("\n[CAMPAIGN] Launch complete.")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="medulla campaign management tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest='command', required=True)

    # -- create ---------------------------------------------------------------
    p_create = sub.add_parser('create', help='Create a new campaign')
    p_create.add_argument('--tag', default='develop',
                          help='Git ref for grid nodes (default: develop)')
    p_create.add_argument('--output', metavar='PATH',
                          help='Base output directory for the campaign')
    p_create.add_argument('--manifest', metavar='PATH',
                          help='Path to a campaign.toml override file')
    p_create.add_argument('--dry-run', action='store_true',
                          help='Print expansion table without creating anything')
    p_create.add_argument('--experiment', metavar='EXP', action='append',
                          help='Filter by experiment (repeatable)')
    p_create.add_argument('--roles', metavar='ROLE', nargs='+',
                          help='Filter by role')
    p_create.add_argument('--analyses', metavar='NAME', nargs='+',
                          help='Filter by analysis name')

    # -- list -----------------------------------------------------------------
    p_list = sub.add_parser('list', help='List discovered analyses and their project units')
    p_list.add_argument('--toml-root', metavar='PATH', dest='toml_root', default=None,
                        help=f'Root directory to search (default: {TOML_DIR})')

    # -- status ---------------------------------------------------------------
    p_status = sub.add_parser('status', help='Show campaign status')
    p_status.add_argument('--campaign', metavar='PATH', required=True,
                          help='Path to the campaign directory')

    # -- sync -----------------------------------------------------------------
    p_sync = sub.add_parser('sync', help='Sync job completion status into campaign.db')
    p_sync.add_argument('--campaign', metavar='PATH', required=True,
                        help='Path to the campaign directory')

    # -- launch ---------------------------------------------------------------
    p_launch = sub.add_parser('launch', help='Launch pending campaign jobs')
    p_launch.add_argument('--campaign', metavar='PATH', required=True,
                          help='Path to the campaign directory')
    p_launch.add_argument('--experiment', metavar='EXP',
                          help='Restrict launch to one experiment')
    p_launch.add_argument('--relaunch', action='store_true',
                          help='Also relaunch projects already marked submitted')
    launch_grp = p_launch.add_mutually_exclusive_group()
    launch_grp.add_argument('--test', action='store_true',
                            help='Submit one job per project to verify setup')
    launch_grp.add_argument('--njobs', type=int, metavar='N',
                            help='Submit at most N jobs per project')

    args = parser.parse_args()

    if args.command == 'list':
        cmd_list(args)
    elif args.command == 'create':
        cmd_create(args)
    elif args.command == 'status':
        cmd_status(args)
    elif args.command == 'sync':
        cmd_sync(args)
    elif args.command == 'launch':
        cmd_launch(args)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("\n[CAMPAIGN] Interrupted.")
        sys.exit(130)
