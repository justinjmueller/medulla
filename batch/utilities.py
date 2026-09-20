# Utilities for batch processing in medulla projects using jobsub
import os
import re
import shutil
import sqlite3
import tempfile
import time
import toml
from catalog import resolve_samples
from glob import glob
import subprocess
from pathlib import Path
from typing import Optional

# Prefix used in a sample's 'path' config to indicate that files should
# be located via a SAMWeb dataset definition rather than a filesystem
# glob pattern, e.g. path = "defname:my_definition_name". The SAMWeb
# station defaults to the job's --experiment, but can be overridden per
# sample for definitions registered under a different station (e.g. the
# combined 'sbn' station rather than 'sbnd' or 'icarus') by prefixing the
# definition name with "<station>:", e.g. path = "defname:sbn:my_definition_name".
DEFNAME_PREFIX = 'defname:'

# ANSI helpers (no third-party dependency)
_INFO     = '\033[1m\033[94m[INFO]\033[0m'      # bold blue
_ERROR    = '\033[1m\033[91m[ERROR]\033[0m'     # bold red
_CAMPAIGN = '\033[1m\033[96m[CAMPAIGN]\033[0m'  # bold cyan

# SQL schema for the configuration table for storing job configurations
SCHEMA_CONFIGURATION = """
CREATE TABLE IF NOT EXISTS configuration (
    jobid INTEGER PRIMARY KEY,
    cfg TEXT NOT NULL
);
"""

# SQL schema for the jobs table for tracking job statuses
SCHEMA_JOBS = """
CREATE TABLE IF NOT EXISTS jobs (
    jobid INTEGER PRIMARY KEY,
    status TEXT,
    sample TEXT,
    catalog_experiment TEXT,
    FOREIGN KEY (jobid) REFERENCES configuration(jobid)
);
"""

def command(
    curs : sqlite3.Cursor,
    comm : str,
    vals : tuple = None
):
    """
    Execute a command defined in a string using the provided SQLite 
    cursor. Multiple values can be executed if provided as a list.

    Parameters
    ----------
    curs : sqlite3.Cursor
        The SQLite cursor handle.
    comm : str
        The base command.
    vals : tuple
        Values to use as arguments for the sql command (tuple).

    Returns
    -------
    None.
    """
    try:
        if isinstance(vals, list):
            curs.executemany(comm, vals)
        elif vals:
            curs.execute(comm, vals)
        else:
            curs.execute(comm)
    except Exception as e:
        print(e)

def safe_copy(src, dst):
    """
    Copy a file, working around dCache access quirks that break both
    shutil.copy2()/GNU cp and even a plain buffered read.

    shutil.copy2() (and GNU cp) try an os.sendfile()/copy_file_range()-based
    zero-copy fast path by default on Linux, which dCache's NFS layer has
    been observed to reject with EPERM for some files (e.g. a large
    project.db) even though ordinary buffered reads of the same file
    succeed -- and Python's shutil does not treat EPERM as recoverable for
    its automatic fallback, it just crashes. This function avoids that fast
    path by doing a plain buffered read/write first.

    That alone is not always sufficient, though: a plain buffered read can
    also fail partway through (observed: EPERM after the first ~16KB) for
    some large files on dCache, independent of sendfile. When the buffered
    copy fails, this function falls back to `ifdh cp`, which negotiates
    dCache's access protocol directly instead of relying on the NFS mount,
    and is already a hard dependency of this codebase's grid-submission
    path (submit.sh, auth.py). Any partial output from the failed attempt
    is removed first so ifdh cp starts clean.

    Every project.db/campaign.db copy in this codebase goes through this
    function instead of shutil.copy2/subprocess cp for that reason.

    Parameters
    ----------
    src : str | Path
        Source file path.
    dst : str | Path
        Destination file path.

    Returns
    -------
    None.
    """
    dst = Path(dst)
    try:
        with open(src, 'rb') as fsrc, open(dst, 'wb') as fdst:
            shutil.copyfileobj(fsrc, fdst)
    except OSError:
        dst.unlink(missing_ok=True)
        subprocess.run(['ifdh', 'cp', str(src), str(dst)], check=True)

def safe_write_text(dst, content):
    """
    Write text to a file, working around dCache quirks: /pnfs does not
    allow overwriting a file in place (a repeat write to an existing
    destination fails, confirmed via gfal-copy's "DESTINATION OVERWRITE"
    error), and even a fresh direct open()/write() can fail with EPERM the
    way safe_copy() already works around for existing-file copies.

    Deletes the destination first (matching the same "delete before
    write" precaution already used for project.db/campaign.db elsewhere
    in this codebase), then tries a plain direct write; on OSError, falls
    back to writing a local temp file and moving it into place with
    `ifdh cp`.

    Parameters
    ----------
    dst : str | Path
        Destination file path.
    content : str
        Text content to write.

    Returns
    -------
    None.
    """
    dst = Path(dst)
    dst.unlink(missing_ok=True)
    try:
        dst.write_text(content)
    except OSError:
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tmp') as f:
            f.write(content)
            tmp = Path(f.name)
        try:
            subprocess.run(['ifdh', 'cp', str(tmp), str(dst)], check=True)
        finally:
            tmp.unlink(missing_ok=True)

def locate_samweb_files(defname : str, experiment : str) -> list:
    """
    Query SAMWeb for the files belonging to a dataset definition and
    return their full paths on disk.

    SAMWeb's listFiles() only returns bare file names, so the location
    of each file is looked up separately (in batches of 1000, SAMWeb's
    practical limit for a single metadata query) via
    getMultipleMetadata(locations=True).

    Parameters
    ----------
    defname : str
        Name of the SAMWeb dataset definition.
    experiment : str
        Experiment name used to configure the SAMWeb client (e.g.
        'sbnd' or 'icarus').

    Returns
    -------
    list[str]
        Sorted, full paths (directory + file name) of the files
        belonging to the definition.
    """
    import samweb_client as sam
    samweb = sam.SAMWebClient(experiment=experiment)
    files = samweb.listFiles(defname=defname)
    if not files:
        return []

    paths = []
    for i in range(0, len(files), 1000):
        chunk = files[i:i+1000]
        res = samweb.getMultipleMetadata(chunk, locations=True)
        for x in res:
            loc = x['locations'][0]['location'].split(':', 1)[-1]
            paths.append(loc + '/' + x['file_name'])

    return sorted(paths)

def parse_defname_path(sample_path : str, default_station : str) -> tuple:
    """
    Parse a 'defname:...' sample path into a (station, defname) pair.

    Accepts 'defname:<name>', which uses default_station as the SAMWeb
    station, or 'defname:<station>:<name>', which looks up the
    definition under the named station instead (e.g. 'sbn' for
    definitions shared across experiments rather than registered under
    'sbnd' or 'icarus').

    Parameters
    ----------
    sample_path : str
        The sample's 'path' value; must start with DEFNAME_PREFIX.
    default_station : str
        SAMWeb station to use when the path does not specify one.

    Returns
    -------
    tuple[str, str]
        (station, defname).
    """
    rest = sample_path[len(DEFNAME_PREFIX):]
    head, sep, tail = rest.partition(':')
    if sep:
        return head, tail
    return default_station, head

def get_samples(
    tml : str,
    batch_size : int,
    catalog_path = None,
    enable_keys = None,
    experiment : str = 'sbnd',
):
    """
    Get the list of samples from the TOML file after filtering the list
    for samples that have been disabled. The batch size is used to
    split samples into multiple separate samples if requested (i.e. for
    processing large samples in smaller chunks).

    Parameters
    ----------
    tml : str
        Path to the TOML file.
    batch_size : int
        Number of files to include in each batch. If <= 0, no batching
        is performed.
    catalog_path : str | Path | None
        Path to the sample catalog.  Passed to resolve_samples.
    enable_keys : list[str] | None
        Sample keys to enable.  Passed to resolve_samples.
    experiment : str
        Default SAMWeb station used to resolve a sample's path when it
        is a SAMWeb definition (see DEFNAME_PREFIX), unless the sample
        overrides the station itself.

    Returns
    -------
    samples : list[dict]
        List of samples that are enabled.
    """
    # Get the initial list of samples from the TOML file that have not
    # been disabled.
    cfg = toml.load(tml)
    cfg = resolve_samples(cfg, catalog_path=catalog_path, enable_keys=enable_keys)
    samples = cfg.get('sample', [])
    enabled_samples = [s for s in samples if not s.get('disable', False)]

    # Process the samples and batch them if requested.
    batches = []
    for sample in enabled_samples:
        sample_path = sample['path']
        if isinstance(sample_path, str) and sample_path.startswith(DEFNAME_PREFIX):
            station, defname = parse_defname_path(sample_path, experiment)
            paths = locate_samweb_files(defname, station)
        else:
            paths = glob(sample_path)
        if len(paths) == 0:
            raise FileNotFoundError(f"No files found for sample {sample.get('name', '<unknown>')} with path {sample['path']}")
        if batch_size is None or batch_size <= 0:
            batches.append(sample)
        else:
            for i in range(0, len(paths), batch_size):
                batch_paths = paths[i:i+batch_size]
                if len(batch_paths) == 0:
                    continue
                new_sample = sample.copy()
                new_sample['path'] = batch_paths
                batches.append(new_sample)

    # Return the list of enabled samples.
    return batches

def create_systematics_cfg(
    base_cfg : dict,
    trees : list[dict],
    samples : list[dict],
):
    """
    Create a TOML configuration for running systematics on the given
    samples. The configuration is based on the provided base
    configuration file, which must implement all systematics. Each pair
    of selection sample and selection tree configuration blocks
    represents a unique output in the final systematics output file.
    Systematics are only valid for MC samples, and must specifically be
    requested in the tree configuration block.

    Parameters
    ----------
    base_cfg : dict
        Base configuration dictionary.
    trees : list[dict]
        List of tree configurations in the selection configuration.
    samples : list[dict]
        List of sample configurations in the selection configuration.
    
    Returns
    -------
    syst_cfg : list[dict]
        List of configuration dictionaries for each sample.
    """
    # Loop over each tree and sample combination. If the sample is data
    # or the tree does not have systematics enabled, skip it. There are
    # some sanity checks as well to ensure that the proper branches are
    # present in the tree configuration.
    syst_trees = {}
    for tree in trees:
        for sample in samples:
            # Check if this combination is already configured. If so,
            # skip it (this can happen due to the expansion of samples
            # into batches).
            key = f"events/{sample['name']}/{tree['name']}"
            if key in syst_trees:
                continue

            # Data samples and samples not requesting systematics are
            # configured with a "copy" action that just copies the
            # selected events to the output without applying any
            # systematics.
            if not sample['ismc'] or not tree.get('add_systematics', False):
                syst_trees[key] = {
                    'origin' : key,
                    'destination' : f'events/{sample["name"]}/',
                    'name' : tree['name'],
                    'action' : 'copy',
                }
            # If the sample is MC and the tree requests systematics, do
            # some additional checking and then configure it with a
            # "add_weights" action.
            else:
                # We need to check that the tree configuration includes
                # both a "neutrino_id" branch and a "neutrino_energy"
                # branch (if the systematics template has
                # "use_additional_hash" set to true). These are used by
                # the systematics code and must be present. Better to
                # catch it here than have the job fail later.
                branch_variables = [(b['name'], b['type']) for b in tree['branch']]
                if ('neutrino_id', 'true') not in branch_variables:
                    raise ValueError(f"Tree {tree['name']} for sample {sample['name']} requests systematics but does not define a 'neutrino_id' branch.")
                if base_cfg.get('input.use_additional_hash', False) and ('neutrino_energy', 'mctruth') not in branch_variables:
                    raise ValueError(f"Tree {tree['name']} for sample {sample['name']} requests systematics but does not define a 'neutrino_energy' branch.")
                syst_trees[key] = {
                    'origin' : key,
                    'destination' : f'events/{sample["name"]}/',
                    'name' : tree['name'],
                    'action' : 'add_weights',
                    'table_types': ['multisim', 'multisigma']
                }

    # Create a new configuration dictionary based on the base
    # configuration. For grid submission purposes, we always set the
    # following:
    # - input.path = "output.root"
    # - input.weights = "data/*flat*.root"
    # - output.path = "output_sys.root"
    # - tree = list of syst_trees values
    syst_cfg = base_cfg.copy()
    syst_cfg['input']['path'] = 'output.root'
    syst_cfg['input']['weights'] = 'data/*flat*.root'
    syst_cfg['output']['path'] = 'output_sys.root'
    syst_cfg['tree'] = list(syst_trees.values())
    return syst_cfg

# Field separator for validation_manifest.txt. Chosen because it cannot
# appear in a ROOT directory or tree name, so the worker-side parser in
# validate_pair.C can split on it without quoting rules.
MANIFEST_SEP = '|'

def write_validation_manifest(
    syst_cfg : dict,
    dst : Path,
):
    """
    Write the manifest that validate_pair.C checks a job's output pair
    against, describing what the systematics step is *expected* to have
    produced for every tree.

    This is derived from the same syst_cfg object that is dumped to
    systematics.toml, rather than recomputed, so the manifest and the
    configuration the job actually runs cannot drift apart.

    The manifest covers every sample in the project, not just one job's;
    each line is keyed by an 'origin' of the form
    "events/<sample>/<tree>", and the worker selects the lines matching
    its own sample. One line per tree:

        origin|name|action|table_types

    where table_types is a comma-separated list (empty for a 'copy'
    action, which produces no systematics trees).

    Parameters
    ----------
    syst_cfg : dict
        The systematics configuration as returned by
        create_systematics_cfg.
    dst : Path
        Path to write the manifest to.

    Returns
    -------
    None.
    """
    lines = []
    for tree in syst_cfg.get('tree', []):
        table_types = ','.join(tree.get('table_types', []))
        lines.append(MANIFEST_SEP.join([
            tree['origin'],
            tree['name'],
            tree['action'],
            table_types,
        ]))
    Path(dst).write_text('\n'.join(lines) + ('\n' if lines else ''))

def create_new_project(
    project_dir : Path,
    tml : str,
    batch_size : int,
    sys : str = None,
    catalog_path = None,
    enable_keys = None,
    experiment : str = 'sbnd',
):
    """
    Create a new project directory with the necessary subdirectories
    and a SQLite database to manage the project. Each sample in the
    TOML file is added as a separate job in the database, with the
    configuration modified to include only that sample.

    Parameters
    ----------
    project_dir : Path
        Path to the base directory for the job directory.
    tml : str
        Path to the TOML file containing the configuration.
    batch_size : int
        Number of files to process in each batch.
    sys : str
        Path to the TOML file containing the systematics configuration
        template. If not provided, the default template in the
        medulla/batch directory is used.
    experiment : str
        Experiment name used to configure the SAMWeb client for samples
        whose path is a SAMWeb definition (see DEFNAME_PREFIX in
        get_samples).

    Returns
    -------
    None.
    """
    # Create the project directory and a subdirectory for job output,
    # if they do not already exist.
    os.makedirs(project_dir, exist_ok=True)
    os.makedirs(project_dir / 'output', exist_ok=True)

    # Connect to the project database. If the database does not exist,
    # it will be created. If the project database does already exist,
    # throw an error because we do not want to overwrite an existing
    # project.
    if (project_dir / 'project.db').exists():
        raise FileExistsError(f"Project database {project_dir / 'project.db'} already exists.")
    conn = sqlite3.connect(project_dir / 'project.db')
    curs = conn.cursor()
    command(curs, SCHEMA_CONFIGURATION)
    command(curs, SCHEMA_JOBS)
    conn.commit()

    # Load the TOML file and get the samples.
    cfg = toml.load(tml)
    cfg = resolve_samples(cfg, catalog_path=catalog_path, enable_keys=enable_keys)
    samples = get_samples(tml, batch_size, catalog_path=catalog_path, enable_keys=enable_keys, experiment='sbnd')

    # Create a systematics configuration based on the selection
    # configuration. This will be used by each job to run systematics
    # after the selection step.
    if sys is None:
        sys = Path(__file__).resolve().parent / 'sys_template.toml'
    sys = create_systematics_cfg(toml.load(sys), cfg.get('tree', []), samples)
    with open(project_dir / 'systematics.toml', 'w') as f:
        toml.dump(sys, f)

    # Write the companion manifest describing what each job's systematics
    # output must contain, so the grid node can validate its output pair
    # before transferring anything back. See write_validation_manifest.
    write_validation_manifest(sys, project_dir / 'validation_manifest.txt')

    # Form a "batch" config for each sample: i.e., each sample gets a
    # copy of the TOML configuration with the [[tree]] list preserved,
    # the [general] section modified to set the 'output' key to its
    # base name plus a batch suffix, and the singular [[sample]]
    # section corresponding to the sample.
    base = cfg['general']['output']
    ins_configurations = []
    ins_jobs = []
    for si, sample in enumerate(samples):
        job_tml = cfg.copy()
        job_tml['general']['output'] = 'output'
        job_tml['sample'] = [sample,]

        ins_configurations.append((si, toml.dumps(job_tml),))
        ins_jobs.append((si, 'pending', sample['name'], sample.get('experiment')))

    # Insert the job configuration into the database.
    command(curs, "INSERT INTO configuration (jobid, cfg) VALUES (?, ?)", ins_configurations)
    command(curs, "INSERT INTO jobs (jobid, status, sample, catalog_experiment) VALUES (?, ?, ?, ?)", ins_jobs)
    conn.commit()
    conn.close()

# Minimum size for an output file to count as real rather than a stub
# left behind by a job that failed after creating its output.
MIN_OUTPUT_BYTES = 1024

def survey_project_output(
    project_dir : Path,
):
    """
    Classify a project's job output by *pair*: a job is complete only when
    both its selection and its systematics output are present and
    non-stub.

    Keying completion on the selection file alone -- which is what this
    replaces -- means a job whose systematics step failed or silently
    produced nothing is still marked 'completed'. Such a job never returns
    to 'pending', so it is never resubmitted, and the gap only surfaces at
    merge time when the two file sets are cross-checked. submit.sh now
    refuses to transfer a half-good pair at all, but projects predating
    that still contain orphans, and this is what finds them.

    Parameters
    ----------
    project_dir : Path
        Path to the project directory (containing output/).

    Returns
    -------
    dict with keys:
        completed : list[int]
            Job IDs with both outputs present and >= MIN_OUTPUT_BYTES.
        orphaned : list[int]
            Job IDs with a good selection output but no usable
            systematics partner. These must be reverted to 'pending'.
        stub_files : list[Path]
            Output files (either kind) below MIN_OUTPUT_BYTES.
    """
    project_dir = Path(project_dir)
    out_dir = project_dir / 'output'

    def _by_id(pattern):
        found = {}
        for f in glob(str(out_dir / pattern)):
            path = Path(f)
            try:
                jid = int(path.stem.split('jobid')[-1])
            except ValueError:
                continue
            found[jid] = path
        return found

    # 'output_jobid*' cannot match 'output_systematics_jobid*', so the two
    # globs stay disjoint without further filtering.
    nosyst = _by_id('output_jobid*.root')
    wsyst = _by_id('output_systematics_jobid*.root')

    def _ok(path):
        return path is not None and path.stat().st_size >= MIN_OUTPUT_BYTES

    completed, orphaned = [], []
    for jid, path in nosyst.items():
        if not _ok(path):
            continue
        if _ok(wsyst.get(jid)):
            completed.append(jid)
        else:
            orphaned.append(jid)

    stub_files = [p for p in list(nosyst.values()) + list(wsyst.values())
                  if p.stat().st_size < MIN_OUTPUT_BYTES]

    return {
        'completed':  sorted(completed),
        'orphaned':   sorted(orphaned),
        'stub_files': sorted(stub_files),
    }

def check_project_status(
    project_dir : str,
):
    """
    Check the status of the project by inspecting the job output in the
    project directory.

    Parameters
    ----------
    project_dir : str
        Path to the base directory for the job directory.

    Returns
    -------
    None.
    """
    # Check if the project database exists.
    if not (project_dir / 'project.db').exists():
        raise FileNotFoundError(f"Project database {project_dir / 'project.db'} does not exist.")
    
    # Copy the project database locally to dodge dcache issues.
    safe_copy(project_dir / 'project.db', './project.db')
    conn = sqlite3.connect('./project.db')
    curs = conn.cursor()

    # Classify the output by pair: a job counts as complete only when both
    # its selection and systematics outputs are present and non-stub. The
    # size floor avoids marking a job complete on an empty file left
    # behind by a failure.
    survey = survey_project_output(project_dir)
    completed_jobs = survey['completed']
    orphaned_jobs = survey['orphaned']

    ins = [('completed', jid) for jid in completed_jobs]
    command(curs, "UPDATE jobs SET status = ? WHERE jobid = ?", ins)

    # A job previously marked complete on its selection output alone, but
    # with no usable systematics partner, has to move backward so it is
    # resubmitted. This is the only place status regresses.
    if orphaned_jobs:
        command(curs, "UPDATE jobs SET status = 'pending' WHERE jobid = ?",
                [(jid,) for jid in orphaned_jobs])
    conn.commit()
    conn.close()

    stub_files = survey['stub_files']
    if stub_files:
        resp = input(
            f"[INFO] -- Found {len(stub_files)} stub output file(s) <"
            f" 1024 bytes.\nDelete these stub outputs? [Y/N] "
        )
        if resp.strip().lower() != 'y':
            print(
                "[INFO] -- Keeping stub output files. Please check"
                " these files manually to determine if they are valid"
                " outputs or if the jobs need to be resubmitted."
            )
        else:
            for stub_file in stub_files:
                if stub_file.exists():
                    stub_file.unlink()
            print(f"[INFO] -- Deleted {len(stub_files)} stub output file(s).")

    # Replace the project database copy with the updated version.
    subprocess.run(['mv', './project.db', project_dir / 'project.db'], check=True)

    print(f"[INFO] -- Found {len(completed_jobs)} completed jobs.")
    if orphaned_jobs:
        print(f"[INFO] -- Reverted {len(orphaned_jobs)} job(s) to pending: "
              f"selection output present but no usable systematics output.")

def _submit_jobsub_once(
    cmd : list,
    njobs : int,
    confirm : bool,
    verbose : bool,
    label : str = '',
):
    """
    Run a single jobsub_submit invocation and report the result. Handles
    the transient vault-credential race and expired-token failure modes,
    retrying once for the former. This is a helper for launch_jobsub,
    factored out so it can be called once per project (the default) or
    once per sample (when njobs_per_sample is used).

    Parameters
    ----------
    cmd : list
        The full jobsub_submit command/argument list to run.
    njobs : int
        The number of jobs requested by this particular invocation, used
        only for reporting.
    confirm : bool
        Whether this is a single-project, user-facing launch (affects how
        much output is shown).
    verbose : bool
        Whether to show full jobsub_submit stdout/stderr on success.
    label : str
        Optional suffix describing what this invocation is for (e.g.
        " for sample 'sbnd_mc'"), appended to the printed messages.

    Returns
    -------
    bool
        True if the submission succeeded, False otherwise.
    """
    # Launch the jobs. If the command raises an "ExpiredSignatureError"
    # exception, it likely means that the user's token has expired and
    # they need to run `htgettoken` to refresh it. The exception is
    # printed to stdout by jobsub, so we just need to catch it and
    # print a more user-friendly message.
    #
    # A separate, transient failure mode has been observed when multiple
    # jobsub_submit calls run in quick succession: HTCondor's vault
    # credential manager (condor_vault_storer) can race against a
    # still-in-progress credential write from a previous submission and
    # refuse to proceed ("Credentials exist that do not match the
    # request"). The requested scopes/handle are unchanged in this case
    # (no real credential problem), so it is safe to retry once after a
    # short delay rather than failing outright.
    max_attempts = 2
    retry_delay = 5  # seconds
    for attempt in range(1, max_attempts + 1):
        try:
            out = subprocess.run(cmd, check=True, capture_output=True, text=True)
            break
        except subprocess.CalledProcessError as e:
            if 'condor_vault_storer' in e.stderr and attempt < max_attempts:
                print(f"{_ERROR} -- Transient vault credential conflict detected, "
                      f"retrying in {retry_delay}s...")
                time.sleep(retry_delay)
                continue
            if 'ExpiredSignatureError' in (output := e.stderr.strip()):
                print(f"{_ERROR} -- Job submission failed{label} due to expired token. Please run `htgettoken` to refresh your token and try again.")
            else:
                print(f"{_ERROR} -- Job submission failed{label} with error: {output}")
            if verbose:
                print(f"{_ERROR} -- Full stdout:\n{e.stdout}")
                print(f"{_ERROR} -- Full stderr:\n{e.stderr}")
            return False

    if confirm:
        # Single-project workflow: show full output so the user can verify.
        stdout = out.stdout.strip()
        print('\n'.join(stdout.split('\n')[-4:]))
        print(f"{_INFO} -- Launched {njobs} jobs{label}.")
    elif verbose:
        # Campaign workflow with verbose requested: jobsub_submit can exit 0
        # while still failing to submit some individual jobs, so show the
        # full output rather than just the one-line summary.
        print(f"{_INFO} -- Full jobsub_submit stdout:\n{out.stdout.strip()}")
        if out.stderr.strip():
            print(f"{_INFO} -- Full jobsub_submit stderr:\n{out.stderr.strip()}")
        match = re.search(r'job id\s+(\S+)', out.stdout)
        job_id = match.group(1) if match else 'unknown'
        print(f"{_CAMPAIGN} Submitted {njobs} job(s){label}. Job ID: {job_id}")
    else:
        # Campaign workflow: one clean line per project.
        match = re.search(r'job id\s+(\S+)', out.stdout)
        job_id = match.group(1) if match else 'unknown'
        print(f"{_CAMPAIGN} Submitted {njobs} job(s){label}. Job ID: {job_id}")
    return True

# The worker-side output validator. It is transferred with every job
# (jobsub_submit -f dropbox://) rather than read out of the checkout the
# job builds, because that checkout is whatever --tag names: a release tag
# that predates the validator would otherwise fail validation on every job
# and copy nothing back. Shipping it from here keeps the validator in step
# with the submit.sh it is shipped alongside.
VALIDATE_MACRO_PATH = Path(__file__).resolve().parent / 'validate_pair.C'

def launch_jobsub(
    project_dir : str,
    exp : str = 'sbnd',
    njobs : int = -1,
    njobs_per_sample : Optional[int] = None,
    confirm : bool = True,
    tag : str = 'develop',
    memory : int = 1800,
    disk : Optional[int] = None,
    lifetime : str = '1h',
    verbose : bool = False,
    force : bool = False,
    jobs_per_process : int = 1,
):
    """
    Launch jobs using jobsub for the given project directory. If njobs
    is provided, only that many jobs will be launched in total. If
    njobs_per_sample is provided instead, up to that many jobs will be
    launched for *each* sample in the project (as a separate jobsub_submit
    call per sample), which ensures every sample gets some jobs launched
    even when the project has an uneven mix of sample sizes.

    Parameters
    ----------
    project_dir : str
        Path to the base directory for the job directory.
    exp : str
        Experiment name (default: sbnd).
    njobs : int
        Number of jobs to launch in total. If -1 (default), launch all
        pending jobs. Mutually exclusive with njobs_per_sample.
    njobs_per_sample : int
        Number of jobs to launch per sample. If provided, njobs must be
        left at its default (-1). Requires a project database created
        with per-sample job tracking (i.e. created after this feature was
        added); older projects should be recreated to use this option.
    confirm : bool
        If True (default), prompt the user before submitting.  Pass
        False when the caller has already obtained confirmation (e.g.
        campaign launch confirms once for all projects).
    tag : str
        Git ref passed to submit.sh as --tag (default: develop).
    memory : int
        Amount of memory to request for each job in MB. If None, use default.
    disk : int
        Amount of disk to request for each job in GB. If None, use default.
    lifetime : str
        Expected lifetime of each job (e.g., '1h', '30m'). If None, use default.
    verbose : bool
        If True, print the full jobsub_submit command and its complete
        stdout/stderr, even on a successful submission. jobsub_submit can
        exit 0 while still failing to submit some individual jobs, and
        those failures are otherwise only visible in the full output,
        which is normally discarded down to a one-line summary.
    force : bool
        If True, pass --force through to submit.sh so that a job's output
        copy-back overwrites any pre-existing file at the destination
        (e.g. left behind by a prior failed or resubmitted attempt at the
        same job) instead of failing. Off by default, since silently
        overwriting could mask two jobs unexpectedly racing to write the
        same output.
    jobs_per_process : int
        Number of job IDs each grid process takes on (default 1). njobs
        and njobs_per_sample still count job IDs; the number of processes
        submitted is that count divided by jobs_per_process, rounded up.
        submit.sh is also told the job-ID total, so the last process stops
        at the requested count instead of claiming a full jobs_per_process.
        The build and environment setup are paid once per process, so this
        trades submission overhead against the lifetime each process needs:
        raise lifetime to cover jobs_per_process job IDs.

    Returns
    -------
    bool
        True if at least one jobsub_submit invocation succeeded, False
        otherwise.
    """
    if njobs_per_sample is not None and njobs != -1:
        raise ValueError("njobs and njobs_per_sample are mutually exclusive.")
    if jobs_per_process < 1:
        raise ValueError("jobs_per_process must be at least 1.")

    # Refuse before touching anything if the validator cannot be shipped.
    # Discovering it on the grid would cost every job its full event loop.
    if not VALIDATE_MACRO_PATH.is_file():
        raise FileNotFoundError(
            f"Output validator {VALIDATE_MACRO_PATH} not found; it is transferred "
            f"with every job and must exist to launch."
        )

    # Check if the project database exists.
    if not (project_dir / 'project.db').exists():
        raise FileNotFoundError(f"Project database {project_dir / 'project.db'} does not exist.")

    # Copy the project database locally to dodge dcache issues.
    safe_copy(project_dir / 'project.db', './project.db')
    conn = sqlite3.connect('./project.db')
    curs = conn.cursor()

    # Build the list of submission targets. Each target is a (sample, count)
    # pair, where sample is None for a plain project-wide launch (today's
    # default behavior) or a sample name when njobs_per_sample is used.
    if njobs_per_sample is not None:
        command(curs, "PRAGMA table_info(jobs)")
        columns = {row[1] for row in curs.fetchall()}
        if 'sample' not in columns:
            conn.close()
            raise RuntimeError(
                "This project database does not have per-sample job "
                "tracking (missing 'jobs.sample' column). Create a new "
                "project to use njobs_per_sample."
            )
        command(
            curs,
            "SELECT sample, COUNT(*) FROM jobs WHERE status = 'pending' "
            "GROUP BY sample ORDER BY MIN(jobid)",
        )
        sample_counts = curs.fetchall()
        conn.close()

        if len(sample_counts) == 0:
            if confirm:
                print(f"{_INFO} -- No pending jobs to launch.")
            return False

        targets = [
            (sample, pending if njobs_per_sample <= 0 else min(njobs_per_sample, pending))
            for sample, pending in sample_counts
        ]

        if confirm:
            print(f"{_INFO} -- Found pending jobs for {len(targets)} sample(s):")
            for sample, count in targets:
                print(f"{_INFO} --   {sample}: {count} job(s)")
    else:
        # Get the list of pending jobs.
        command(curs, "SELECT jobid FROM jobs WHERE status = 'pending'")
        pending_jobs = [row[0] for row in curs.fetchall()]
        conn.close()

        # Do some checking that the request is sane. Naturally, if there
        # are no pending jobs, there is nothing to launch. Similarly, if
        # the user requested more jobs than are pending, just launch all
        # of the pending jobs.
        if len(pending_jobs) == 0:
            if confirm:
                print(f"{_INFO} -- No pending jobs to launch.")
            return False
        if njobs > len(pending_jobs):
            njobs = len(pending_jobs)
            if confirm:
                print(f"{_INFO} -- Requested number of jobs exceeds pending jobs. Preparing {njobs} jobs instead.")
        if njobs == -1:
            njobs = len(pending_jobs)

        if confirm:
            print(f"{_INFO} -- Found {len(pending_jobs)} pending jobs.")

        targets = [(None, njobs)]

    # Determine the disk request.
    if disk is not None:
        disk_flag = f'--disk={disk}GB'
    elif exp == 'sbnd':
        disk_flag = '--disk=10GB'
    else:
        disk_flag = '--disk=25GB'

    # Form the jobsub command(s) to launch the jobs, one per target. count
    # is a number of job IDs, while -N is a number of grid processes, each
    # of which takes up to jobs_per_process of them.
    submissions = []
    for sample, count in targets:
        n_processes = -(-count // jobs_per_process)  # ceil(count / jobs_per_process)
        cmd = [
            'jobsub_submit',
            '-G', exp,
            '-N', str(n_processes),
            f'--memory={memory}MB',
            disk_flag,
            f'--expected-lifetime={lifetime}',
            '--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC',
            '--site=FermiGrid',
            "--append_condor_requirements='(TARGET.HAS_Singularity==true)'",
            '--singularity-image=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest',
            # Ship the validator with the job; see VALIDATE_MACRO_PATH. The
            # batch directory is not grid-accessible, hence dropbox://.
            '-f', f'dropbox://{VALIDATE_MACRO_PATH}',
            f'file://{Path(__file__).resolve().parent / "submit.sh"}',
            '--',
            f'--project={project_dir.resolve()}',
            f'--tag={tag}',
            f'--jobs-per-process={jobs_per_process}',
            # Lets the last process stop at the requested count.
            f'--total-jobs={count}',
        ]
        if sample is not None:
            cmd.append(f'--sample={sample}')
        if force:
            cmd.append('--force')
        label = f" for sample '{sample}'" if sample is not None else ""
        if jobs_per_process > 1:
            label += f" ({n_processes} process(es) of up to {jobs_per_process} job IDs)"
        submissions.append((cmd, count, label))

    # Show what will be launched and confirm once for the whole set.
    if confirm or verbose:
        for cmd, count, label in submissions:
            print(f"{_INFO} -- Launching {count} jobs{label} with command: {' '.join(cmd)}")
    if confirm:
        resp = input("Confirm job launch? [Y/N] ")
        if resp.lower() != 'y':
            print(f"{_INFO} -- User aborted job launch.")
            return False

    any_success = False
    for cmd, count, label in submissions:
        ok = _submit_jobsub_once(cmd, count, confirm, verbose, label=label)
        any_success = any_success or ok

    return any_success

def check_git_branch(
    branch : str,
    repo_url : str = 'https://github.com/justinjmueller/medulla',
):
    """
    Check if the specified branch or tag exists in the given Git 
    repository. First checks for branches, then tags if not found.

    Parameters
    ----------
    branch : str
        Branch or tag name to check for existence.
    repo_url : str
        URL to the Git repository.

    Returns
    -------
    bool
        True if the branch or tag exists, False otherwise.
    """
    # Check if it exists as a branch
    result = subprocess.run(
        ["git", "ls-remote", "--exit-code", "--heads", repo_url, branch],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if result.returncode == 0:
        return True
    
    # If not a branch, check if it exists as a tag
    result = subprocess.run(
        ["git", "ls-remote", "--exit-code", "--tags", repo_url, branch],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return result.returncode == 0
