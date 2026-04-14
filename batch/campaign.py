"""Campaign management CLI for medulla batch processing."""
import argparse
import sqlite3
import sys
import toml
from datetime import datetime, timezone
from pathlib import Path

from auth import authenticate
from catalog import resolve_samples
from utilities import create_new_project, check_project_status, launch_jobsub

# Repo root is two levels above this script (batch/ -> repo root).
REPO_ROOT = Path(__file__).resolve().parent.parent
TOML_DIR = REPO_ROOT / 'selection' / 'toml'

# ---------------------------------------------------------------------------
# Database helpers
# ---------------------------------------------------------------------------

SCHEMA_CAMPAIGN_META = """
CREATE TABLE IF NOT EXISTS campaign_meta (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
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
    status TEXT DEFAULT 'created',
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    submitted_at TIMESTAMP,
    completed_at TIMESTAMP,
    UNIQUE(analysis, role, experiment)
);
"""


def _open_db(campaign_dir):
    """Open the campaign SQLite database and return (conn, cursor)."""
    db_path = Path(campaign_dir) / 'campaign.db'
    if not db_path.exists():
        raise FileNotFoundError(f"Campaign database not found: {db_path}")
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    return conn, conn.cursor()


# ---------------------------------------------------------------------------
# Discovery helpers
# ---------------------------------------------------------------------------

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
    manifest_override = {}
    if args.manifest:
        manifest_override = toml.load(args.manifest)

    filters_exp = set(args.experiment) if args.experiment else None
    filters_roles = set(args.roles) if args.roles else None
    filters_analyses = set(args.analyses) if args.analyses else None

    # Discover and expand all project units.
    all_units = []
    for meta_path in _discover_meta_files():
        if filters_analyses:
            meta = toml.load(meta_path)
            if meta['meta']['analysis'] not in filters_analyses:
                continue
        units = _expand_projects(meta_path, filters_exp, filters_roles)
        all_units.extend(units)

    if not all_units:
        print("[CAMPAIGN] No project units matched the given filters.")
        return

    # Print summary table.
    col_w = [20, 18, 10]
    header = f"{'Analysis':<{col_w[0]}} {'Role':<{col_w[1]}} {'Experiment':<{col_w[2]}}"
    print(header)
    print('-' * sum(col_w))
    for u in all_units:
        print(f"{u['analysis']:<{col_w[0]}} {u['role']:<{col_w[1]}} {u['experiment']:<{col_w[2]}}")
    print(f"\n[CAMPAIGN] Total: {len(all_units)} project unit(s).")

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

    # Build a campaign name from the tag and timestamp.
    tag = args.tag
    ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')
    campaign_name = f"campaign_{tag}_{ts}"
    campaign_dir = output_base / campaign_name
    campaign_dir.mkdir(parents=True, exist_ok=True)

    # Initialise the campaign database.
    conn = sqlite3.connect(campaign_dir / 'campaign.db')
    conn.row_factory = sqlite3.Row
    curs = conn.cursor()
    curs.executescript(SCHEMA_CAMPAIGN_META + SCHEMA_PROJECTS)
    curs.execute("INSERT OR REPLACE INTO campaign_meta VALUES (?, ?)", ('tag', tag))
    curs.execute("INSERT OR REPLACE INTO campaign_meta VALUES (?, ?)", ('name', campaign_name))
    curs.execute("INSERT OR REPLACE INTO campaign_meta VALUES (?, ?)", ('created_at', ts))
    conn.commit()

    created = []
    for u in all_units:
        proj_dir = campaign_dir / u['analysis'] / f"{u['role']}_{u['experiment']}"
        try:
            create_new_project(
                project_dir=proj_dir,
                tml=u['toml_file'],
                batch_size=u['batch_size'],
            )
            status = 'created'
        except FileExistsError:
            print(f"[CAMPAIGN] Skipping existing project: {proj_dir}")
            status = 'exists'
        except Exception as e:
            print(f"[CAMPAIGN] Failed to create {proj_dir}: {e}")
            status = 'failed'

        curs.execute(
            "INSERT OR IGNORE INTO projects "
            "(analysis, role, experiment, toml_file, project_dir, batch_size, status) "
            "VALUES (?, ?, ?, ?, ?, ?, ?)",
            (u['analysis'], u['role'], u['experiment'],
             u['toml_file'], str(proj_dir), u['batch_size'], status),
        )
        conn.commit()
        created.append((u, str(proj_dir), status))

    conn.close()

    # Write a manifest snapshot.
    manifest_snapshot = {
        'campaign': {'name': campaign_name, 'tag': tag, 'created_at': ts},
        'projects': [
            {k: v for k, v in u.items() if k != 'enable_keys'}
            for u in all_units
        ],
    }
    with open(campaign_dir / 'campaign_manifest.toml', 'w') as f:
        toml.dump(manifest_snapshot, f)

    print(f"\n[CAMPAIGN] Campaign '{campaign_name}' created at {campaign_dir}")
    print(f"[CAMPAIGN] {sum(1 for _, _, s in created if s == 'created')} project(s) created successfully.")


# ---------------------------------------------------------------------------
# `status` subcommand
# ---------------------------------------------------------------------------

def cmd_status(args):
    """Print aggregated status for all projects in the campaign."""
    conn, curs = _open_db(args.campaign)
    curs.execute("SELECT * FROM projects ORDER BY analysis, role, experiment")
    rows = list(curs.fetchall())
    conn.close()

    if not rows:
        print("[CAMPAIGN] No projects registered in this campaign.")
        return

    col_w = [20, 18, 10, 12]
    header = (f"{'Analysis':<{col_w[0]}} {'Role':<{col_w[1]}} "
              f"{'Experiment':<{col_w[2]}} {'Status':<{col_w[3]}}")
    print(header)
    print('-' * sum(col_w))
    for row in rows:
        print(
            f"{row['analysis']:<{col_w[0]}} {row['role']:<{col_w[1]}} "
            f"{row['experiment']:<{col_w[2]}} {row['status']:<{col_w[3]}}"
        )
    print(f"\n[CAMPAIGN] {len(rows)} project(s) total.")


# ---------------------------------------------------------------------------
# `launch` subcommand
# ---------------------------------------------------------------------------

def cmd_launch(args):
    """Authenticate per experiment and launch pending projects."""
    conn, curs = _open_db(args.campaign)

    # Build query (optionally filtered by experiment).
    if args.experiment:
        curs.execute(
            "SELECT * FROM projects WHERE experiment = ? AND status = 'created'",
            (args.experiment,),
        )
    else:
        curs.execute("SELECT * FROM projects WHERE status = 'created'")

    rows = list(curs.fetchall())
    conn.close()

    if not rows:
        print("[CAMPAIGN] No projects ready for launch.")
        return

    # Group by experiment.
    by_exp = {}
    for row in rows:
        by_exp.setdefault(row['experiment'], []).append(row)

    for exp, projects in by_exp.items():
        print(f"\n[CAMPAIGN] Authenticating for experiment: {exp}")
        ok = authenticate(exp)
        if not ok:
            print(f"[CAMPAIGN] Skipping {exp} (authentication declined).")
            continue

        conn = sqlite3.connect(Path(args.campaign) / 'campaign.db')
        conn.row_factory = sqlite3.Row
        curs = conn.cursor()

        for row in projects:
            proj_dir = Path(row['project_dir'])
            print(f"[CAMPAIGN] Launching: {row['analysis']}/{row['role']}_{row['experiment']}")
            try:
                launch_jobsub(proj_dir, exp=exp)
                curs.execute(
                    "UPDATE projects SET status = 'submitted', submitted_at = ? WHERE project_id = ?",
                    (datetime.now(timezone.utc).isoformat(), row['project_id']),
                )
            except Exception as e:
                print(f"[CAMPAIGN] Launch failed for {proj_dir}: {e}")
            conn.commit()

        conn.close()

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

    # -- status ---------------------------------------------------------------
    p_status = sub.add_parser('status', help='Show campaign status')
    p_status.add_argument('--campaign', metavar='PATH', required=True,
                          help='Path to the campaign directory')

    # -- launch ---------------------------------------------------------------
    p_launch = sub.add_parser('launch', help='Launch pending campaign jobs')
    p_launch.add_argument('--campaign', metavar='PATH', required=True,
                          help='Path to the campaign directory')
    p_launch.add_argument('--experiment', metavar='EXP',
                          help='Restrict launch to one experiment')

    args = parser.parse_args()

    if args.command == 'create':
        cmd_create(args)
    elif args.command == 'status':
        cmd_status(args)
    elif args.command == 'launch':
        cmd_launch(args)


if __name__ == '__main__':
    main()
