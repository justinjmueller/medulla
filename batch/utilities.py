# Utilities for batch processing in medulla projects using jobsub
import os
import re
import sqlite3
import toml
from catalog import resolve_samples
from glob import glob
import subprocess
from pathlib import Path
from typing import Optional

# ANSI helpers (no third-party dependency)
_INFO     = '\033[1m\033[94m[INFO]\033[0m'      # bold blue
_ERROR    = '\033[1m\033[91m[ERROR]\033[0m'     # bold red
_CAMPAIGN = '\033[1m\033[96m[CAMPAIGN]\033[0m'  # bold cyan

def _ifdh_cp(src: str, dest: str):
    """Copy src to dest via ifdh, removing dest first if it already exists."""
    subprocess.run(['ifdh', 'rm', dest], capture_output=True)
    subprocess.run(['ifdh', 'cp', src, dest], check=True)

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

def get_samples(
    tml : str,
    batch_size : int,
    catalog_path = None,
    enable_keys = None,
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
        paths = glob(sample['path'])
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
                table_types = ['multisim', 'multisigma']
                if 'variations' in base_cfg:
                    table_types.append('variation')
                syst_trees[key] = {
                    'origin' : key,
                    'destination' : f'events/{sample["name"]}/',
                    'name' : tree['name'],
                    'action' : 'add_weights',
                    'table_types': table_types
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

def create_new_project(
    project_dir : Path,
    tml : str,
    batch_size : int,
    sys : str = None,
    catalog_path = None,
    enable_keys = None,
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
    samples = get_samples(tml, batch_size, catalog_path=catalog_path, enable_keys=enable_keys)

    # Create a systematics configuration based on the selection
    # configuration. This will be used by each job to run systematics
    # after the selection step.
    if sys is None:
        sys = Path(__file__).resolve().parent / 'sys_template.toml'
    sys = create_systematics_cfg(toml.load(sys), cfg.get('tree', []), samples)
    with open(project_dir / 'systematics.toml', 'w') as f:
        toml.dump(sys, f)

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
        ins_jobs.append((si, 'pending'))

    # Insert the job configuration into the database.
    command(curs, "INSERT INTO configuration (jobid, cfg) VALUES (?, ?)", ins_configurations)
    command(curs, "INSERT INTO jobs (jobid, status) VALUES (?, ?)", ins_jobs)
    conn.commit()
    conn.close()

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
    subprocess.run(['cp', project_dir / 'project.db', './project.db'], check=True)
    conn = sqlite3.connect('./project.db')
    curs = conn.cursor()

    # Get the list of job outputs in the output directory. We require
    # that the output file be at least 1 KB in size to be considered
    # complete. This helps avoid marking jobs as complete if they
    # failed and produced an empty output file.
    output_files = glob(str(project_dir / 'output' / 'output_jobid*.root'))
    output_jobids = {
        int(Path(f).stem.split('jobid')[-1])
        for f in output_files
    }

    # Find jobs tracked in project.db that do not have a corresponding
    # output file in output/.
    command(curs, "SELECT jobid, status FROM jobs")
    db_rows = curs.fetchall()
    db_jobids = {row[0] for row in db_rows}
    status_by_jobid = {row[0]: row[1] for row in db_rows}

    missing_output_jobs = sorted(db_jobids - output_jobids)
    missing_output_nonpending = [
        jid for jid in missing_output_jobs
        if status_by_jobid.get(jid) != 'pending'
    ]

    completed_jobs = [
        int(Path(f).stem.split('jobid')[-1])
        for f in output_files if Path(f).stat().st_size >= 1024
    ]
    ins = [('completed', jid) for jid in completed_jobs]
    command(curs, "UPDATE jobs SET status = ? WHERE jobid = ?", ins)
    conn.commit()

    stub_jobs = [
        int(Path(f).stem.split("jobid")[-1])
        for f in output_files
        if Path(f).stat().st_size < 1024
    ]
    if stub_jobs:
        resp = input(
            f"[INFO] -- Found {len(stub_jobs)} stub output file(s) <"
            f" 1024 bytes.\nDelete these stub outputs? [Y/N] "
        )
        if resp.strip().lower() != 'y':
            print(
                "[INFO] -- Keeping stub output files. Please check"
                " these files manually to determine if they are valid"
                " outputs or if the jobs need to be resubmitted."
            )
        else:
            for jid in stub_jobs:
                stub_file = project_dir / 'output' / f'output_jobid{jid:04d}.root'
                if stub_file.exists():
                    stub_file.unlink()
            print(f"[INFO] -- Deleted {len(stub_jobs)} stub output file(s).")

    if missing_output_jobs:
        print(
            f"[INFO] -- Found {len(missing_output_jobs)} job(s) present in project.db"
            " but missing output files in output/."
        )
        print(f"[INFO] -- Missing output job IDs: {missing_output_jobs}")
    else:
        print("[INFO] -- No jobs are missing output files (project.db and output/ are in sync).")

    if missing_output_nonpending:
        print(
            f"[WARN] -- {len(missing_output_nonpending)} missing-output job(s) are"
            " not pending; these may need investigation/resubmission."
        )
        print(f"[WARN] -- Non-pending missing output job IDs: {missing_output_nonpending}")

    conn.close()

    # Replace the project database copy with the updated version.
    subprocess.run(['mv', './project.db', project_dir / 'project.db'], check=True)

    print(f"[INFO] -- Found {len(completed_jobs)} completed jobs.")
def launch_jobsub(
    project_dir : str,
    exp : str = 'sbnd',
    njobs : int = -1,
    tag : str = 'develop',
    memory : int = 1800,
    disk : Optional[int] = None,
    lifetime : str = '1h',
    relaunch_missing : bool = False,
    confirm : bool = True,
):
    """
    Launch jobs using jobsub for the given project directory. If njobs
    is provided, only that many jobs will be launched.

    Parameters
    ----------
    project_dir : str
        Path to the base directory for the job directory.
    exp : str
        Experiment name (default: sbnd).
    njobs : int
        Number of jobs to launch. If None, launch all pending jobs.
    tag : str
        Git ref passed to submit.sh as --tag (default: develop).
    memory : int | None
        Amount of memory to request for each job in MB. If None, use default.
    disk : int | None
        Amount of disk to request for each job in GB. If None, use default.
    lifetime : str | None
        Expected lifetime of each job (e.g., '1h', '30m'). If None, use default.
    relaunch_missing : bool
        If True, jobs found in project.db but missing output files in output/
        are reset to 'pending' so they can be relaunched.
    confirm : bool
        If True (default), prompt the user before submitting.  Pass
        False when the caller has already obtained confirmation (e.g.
        campaign launch confirms once for all projects).

    Returns
    -------
    None.
    """
    # Check if the project database exists.
    if not (project_dir / 'project.db').exists():
        raise FileNotFoundError(f"Project database {project_dir / 'project.db'} does not exist.")

    # Copy the project database locally to dodge dcache issues.
    subprocess.run(['cp', project_dir / 'project.db', './project.db'], check=True)
    conn = sqlite3.connect('./project.db')
    curs = conn.cursor()

    # Get the list of jobs from the DB and derive pending jobs.
    command(curs, "SELECT jobid, status FROM jobs")
    db_rows = curs.fetchall()
    db_jobids = {row[0] for row in db_rows}
    status_by_jobid = {row[0]: row[1] for row in db_rows}
    pending_jobs = [jid for jid, status in db_rows if status == 'pending']

    # Optionally relaunch jobs tracked in project.db but missing output files.
    if relaunch_missing:
        output_files = glob(str(project_dir / 'output' / 'output_jobid*.root'))
        output_jobids = {
            int(Path(f).stem.split('jobid')[-1])
            for f in output_files
        }
        missing_output_jobs = sorted(db_jobids - output_jobids)
        relaunch_jobs = [
            jid for jid in missing_output_jobs
            if status_by_jobid.get(jid) != 'pending'
        ]
        if relaunch_jobs:
            ins = [('pending', jid) for jid in relaunch_jobs]
            command(curs, "UPDATE jobs SET status = ? WHERE jobid = ?", ins)
            conn.commit()
            pending_jobs = sorted(set(pending_jobs).union(relaunch_jobs))
            print(
                f"[INFO] -- Marked {len(relaunch_jobs)} missing-output job(s) as pending for relaunch: {relaunch_jobs}"
            )
        else:
            print("[INFO] -- No additional missing-output jobs needed relaunch.")

    conn.close()
    subprocess.run(['mv', './project.db', project_dir / 'project.db'], check=True)

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
    if njobs > 10000:
        raise ValueError(f"Requested number of jobs ({njobs}) exceeds reasonable limits. Please check the project status and reduce the number of jobs to launch.")

    if confirm:
        print(f"{_INFO} -- Found {len(pending_jobs)} pending jobs.")

    disk_size = '10GB' if exp == 'sbnd' else '25GB'
    if disk is not None:
        disk_size = f'{disk}GB'

    # Form the jobsub command to launch the jobs.
    cmd = [
        'jobsub_submit',
        '-G', exp,
        '-N', str(njobs),
        f'--memory={memory}MB',
        f'--expected-lifetime={lifetime}',
        '--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE',
        "--append_condor_requirements='(TARGET.HAS_Singularity==true)'",
        '--singularity-image=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest',
        f'file://{Path(__file__).resolve().parent / "submit.sh"}',
        '--',
        f'--project={project_dir.resolve()}',
        f'--branch={branch}',
    ]

    if disk is not None:
        cmd.append(f'--disk={disk}GB')
    elif exp == 'sbnd':
        cmd.append(f'--disk=10GB')
    else:
        cmd.append(f'--disk=25GB')

    # Query the user to confirm that they want to launch the jobs.
    if confirm:
        print(f"{_INFO} -- Launching {njobs} jobs with command: {' '.join(cmd)}")
        resp = input("Confirm job launch? [Y/N] ")
        if resp.lower() != 'y':
            print(f"{_INFO} -- User aborted job launch.")
            return False

    # Launch the jobs. If the command raises an "ExpiredSignatureError"
    # exception, it likely means that the user's token has expired and
    # they need to run `htgettoken` to refresh it. The exception is
    # printed to stdout by jobsub, so we just need to catch it and
    # print a more user-friendly message.
    try:
        out = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        if 'ExpiredSignatureError' in (output := e.stderr.strip()):
            print(f"{_ERROR} -- Job submission failed due to expired token. Please run `htgettoken` to refresh your token and try again.")
        else:
            print(f"[ERROR] -- Job submission failed with error: {output}")
        return False

    if confirm:
        # Single-project workflow: show full output so the user can verify.
        stdout = out.stdout.strip()
        print('\n'.join(stdout.split('\n')[-4:]))
        print(f"{_INFO} -- Launched {njobs} jobs.")
    else:
        # Campaign workflow: one clean line per project.
        match = re.search(r'job id\s+(\S+)', out.stdout)
        job_id = match.group(1) if match else 'unknown'
        print(f"{_CAMPAIGN} Submitted {njobs} job(s). Job ID: {job_id}")
    return True

def run_variation_phase1_interactive(
    project_dir: Path,
    variation_input: str,
    variation_toml: str,
):
    """
    Run Phase 1 (spline building) interactively on the current node.

    Loads the variation TOML, strips [[tree]] blocks, rewrites the
    [input] and [output] paths, runs run_systematics locally, then
    stages the resulting splines file to project_dir.

    Parameters
    ----------
    project_dir : Path
        Path to the project directory (PNFS or local).
    variation_input : str
        Path or xrootd URL to the merged variation+CV ROOT file.
    variation_toml : str
        Local path to the variation systematics TOML.

    Returns
    -------
    bool
        True if successful, False otherwise.
    """
    import tempfile

    binary = Path(__file__).resolve().parent.parent / 'build' / 'systematics' / 'run_systematics'
    if not binary.exists():
        print(f"{_ERROR} -- run_systematics binary not found at {binary}. Build medulla first.")
        return False

    variation_toml_path = Path(variation_toml).resolve()
    if not variation_toml_path.exists():
        print(f"{_ERROR} -- Variation systematics TOML not found: {variation_toml_path}")
        return False

    cfg = toml.load(str(variation_toml_path))
    cfg['input']['path'] = variation_input
    cfg['output']['path'] = 'variation_splines.root'
    cfg.pop('tree', None)

    print(f"{_INFO} -- Running Phase 1 interactively:")
    print(f"{_INFO} --   TOML:   {variation_toml_path}")
    print(f"{_INFO} --   Input:  {variation_input}")
    print(f"{_INFO} --   Output: {project_dir}/variation_splines.root")

    resp = input("Confirm Phase 1 interactive run? [Y/N] ")
    if resp.strip().lower() != 'y':
        print(f"{_INFO} -- User aborted.")
        return False

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        toml_path = tmpdir / 'variation_systematics_phase1.toml'
        with open(toml_path, 'w') as f:
            toml.dump(cfg, f)

        try:
            subprocess.run([str(binary), str(toml_path)], check=True, cwd=str(tmpdir))
        except subprocess.CalledProcessError as e:
            print(f"{_ERROR} -- run_systematics failed: {e}")
            return False

        splines_local = tmpdir / 'variation_splines.root'
        if not splines_local.exists():
            print(f"{_ERROR} -- Expected output not found: {splines_local}")
            return False

        splines_dest = str(project_dir / 'variation_splines.root')
        print(f"{_INFO} -- Staging splines to {splines_dest}")
        _ifdh_cp(str(splines_local), splines_dest)

    print(f"{_INFO} -- Phase 1 complete. Splines at {project_dir}/variation_splines.root")
    return True


def launch_variation_phase1_jobsub(
    project_dir: Path,
    variation_input: str,
    variation_toml: str,
    exp: str = 'sbnd',
    branch: str = 'develop',
    memory: int = 8000,
    disk: Optional[int] = None,
    lifetime: str = '8h',
):
    """
    Submit a single batch job to run Phase 1 (spline building) on the
    grid. Stages the variation TOML to the project directory before
    submission so the grid worker can retrieve it.

    Parameters
    ----------
    project_dir : Path
        Path to the project directory (PNFS).
    variation_input : str
        PNFS or xrootd path to the merged variation+CV ROOT file.
    variation_toml : str
        Local path to the variation systematics TOML.
    exp : str
        Experiment name (default: sbnd).
    branch : str
        Medulla git branch (default: develop).
    memory : int
        Memory to request in MB (default: 8000).
    disk : int | None
        Disk to request in GB. If None, defaults to 80 GB.
    lifetime : str
        Expected job lifetime (default: '8h').

    Returns
    -------
    bool
        True if submission succeeded, False otherwise.
    """
    runner = Path(__file__).resolve().parent / 'submit_variation_phase1.sh'
    if not runner.exists():
        print(f"{_ERROR} -- Runner script not found: {runner}")
        return False

    variation_toml_path = Path(variation_toml).resolve()
    if not variation_toml_path.exists():
        print(f"{_ERROR} -- Variation systematics TOML not found: {variation_toml_path}")
        return False

    dest_toml = str(project_dir / 'variation_systematics.toml')
    print(f"{_INFO} -- Copying {variation_toml_path} -> {dest_toml}")
    _ifdh_cp(str(variation_toml_path), dest_toml)

    disk_size = f'{disk}GB' if disk is not None else '80GB'

    cmd = [
        'jobsub_submit',
        '-G', exp,
        '-N', '1',
        f'--memory={memory}MB',
        f'--disk={disk_size}',
        f'--expected-lifetime={lifetime}',
        '--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE',
        "--append_condor_requirements='(TARGET.HAS_Singularity==true)'",
        '--singularity-image=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest',
        f'file://{runner}',
        '--',
        f'--project={project_dir.resolve()}',
        f'--branch={branch}',
        f'--variation-input={variation_input}',
    ]

    print(f"{_INFO} -- Submitting Phase 1 variation systematics job:")
    print(f"{_INFO} --   TOML:   {variation_toml_path}")
    print(f"{_INFO} --   Input:  {variation_input}")
    print(f"{_INFO} --   Output: {project_dir}/variation_splines.root")
    print(f"{_INFO} --   Command: {' '.join(cmd)}")

    resp = input("Confirm Phase 1 batch submission? [Y/N] ")
    if resp.strip().lower() != 'y':
        print(f"{_INFO} -- User aborted.")
        return False

    try:
        out = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        if 'ExpiredSignatureError' in (output := e.stderr.strip()):
            print(f"{_ERROR} -- Submission failed: expired token. Run `htgettoken` and retry.")
        else:
            print(f"{_ERROR} -- Submission failed: {output}")
        return False

    stdout = out.stdout.strip()
    print('\n'.join(stdout.split('\n')[-4:]))
    print(f"{_INFO} -- Phase 1 job submitted.")
    return True


def launch_variation_phase2_jobsub(
    project_dir: Path,
    variation_toml: str,
    exp: str = 'sbnd',
    branch: str = 'develop',
    memory: int = 4000,
    disk: Optional[int] = None,
    lifetime: str = '4h',
    njobs: int = -1,
):
    """
    Submit one batch job per completed selection output file to apply
    pre-built detector-variation splines (Phase 2). Files that already
    have a corresponding output_varsys_jobid<NNNN>.root are skipped.

    A manifest (variation_phase2_manifest.txt) listing the job IDs to
    process is staged to the project directory so each grid worker can
    look up its own input file via $PROCESS.

    Parameters
    ----------
    project_dir : Path
        Path to the project directory (PNFS).
    variation_toml : str
        Local path to the variation systematics TOML.
    exp : str
        Experiment name (default: sbnd).
    branch : str
        Medulla git branch (default: develop).
    memory : int
        Memory to request in MB (default: 4000).
    disk : int | None
        Disk to request in GB. If None, defaults to 25 GB.
    lifetime : str
        Expected job lifetime (default: '4h').
    njobs : int
        Maximum number of jobs to submit. -1 submits all pending files.

    Returns
    -------
    bool
        True if submission succeeded, False otherwise.
    """
    runner = Path(__file__).resolve().parent / 'submit_variation_phase2.sh'
    if not runner.exists():
        print(f"{_ERROR} -- Runner script not found: {runner}")
        return False

    variation_toml_path = Path(variation_toml).resolve()
    if not variation_toml_path.exists():
        print(f"{_ERROR} -- Variation systematics TOML not found: {variation_toml_path}")
        return False

    splines_path = project_dir / 'variation_splines.root'
    if not splines_path.exists():
        print(f"{_ERROR} -- Phase 1 splines file not found: {splines_path}")
        print(f"{_ERROR} -- Run Phase 1 first (--variation-phase1).")
        return False

    all_output_files = sorted(glob(str(project_dir / 'output' / 'output_jobid*.root')))
    if not all_output_files:
        print(f"{_ERROR} -- No selection output files found in {project_dir}/output/")
        return False

    already_done = {
        int(Path(f).stem.split('jobid')[-1])
        for f in glob(str(project_dir / 'output' / 'output_varsys_jobid*.root'))
    }

    pending_files = [
        f for f in all_output_files
        if int(Path(f).stem.split('jobid')[-1]) not in already_done
    ]

    if not pending_files:
        print(f"{_INFO} -- All selection outputs already have Phase 2 results. Nothing to do.")
        return False

    jobids = [int(Path(f).stem.split('jobid')[-1]) for f in pending_files]
    if njobs > 0:
        jobids = jobids[:njobs]
    njobs = len(jobids)

    # Build the Phase 2 TOML in Python so the grid job doesn't need awk
    # rewriting. The only substitution left for the grid is replacing the
    # input path placeholder with the per-job staged file name.
    cfg = toml.load(str(variation_toml_path))
    cfg['input']['path'] = '__INPUT_FILE__'
    cfg['output']['path'] = 'output_varsys.root'
    if 'variations' in cfg:
        cfg['variations']['splines_file'] = 'variation_splines.root'
    for tree in cfg.get('tree', []):
        if tree.get('action') == 'add_weights':
            tree['action'] = 'add_detsys_weights'
    phase2_toml_local = Path('variation_systematics_phase2.toml')
    with open(phase2_toml_local, 'w') as f:
        toml.dump(cfg, f)
    _ifdh_cp(str(phase2_toml_local), str(project_dir / 'variation_systematics_phase2.toml'))
    phase2_toml_local.unlink()

    manifest_local = Path('variation_phase2_manifest.txt')
    manifest_local.write_text('\n'.join(str(jid) for jid in jobids) + '\n')
    manifest_dest = str(project_dir / 'variation_phase2_manifest.txt')
    _ifdh_cp(str(manifest_local), manifest_dest)
    manifest_local.unlink()

    disk_size = f'{disk}GB' if disk is not None else '25GB'

    cmd = [
        'jobsub_submit',
        '-G', exp,
        '-N', str(njobs),
        f'--memory={memory}MB',
        f'--disk={disk_size}',
        f'--expected-lifetime={lifetime}',
        '--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE',
        "--append_condor_requirements='(TARGET.HAS_Singularity==true)'",
        '--singularity-image=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest',
        f'file://{runner}',
        '--',
        f'--project={project_dir.resolve()}',
        f'--branch={branch}',
    ]

    print(f"{_INFO} -- Submitting {njobs} Phase 2 variation systematics jobs:")
    print(f"{_INFO} --   TOML:    {variation_toml_path}")
    print(f"{_INFO} --   Splines: {splines_path}")
    print(f"{_INFO} --   Jobs:    {njobs}")
    print(f"{_INFO} --   Output:  {project_dir}/output/output_varsys_jobid<NNNN>.root")
    print(f"{_INFO} --   Command: {' '.join(cmd)}")

    resp = input("Confirm Phase 2 batch submission? [Y/N] ")
    if resp.strip().lower() != 'y':
        print(f"{_INFO} -- User aborted.")
        return False

    try:
        out = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        if 'ExpiredSignatureError' in (output := e.stderr.strip()):
            print(f"{_ERROR} -- Submission failed: expired token. Run `htgettoken` and retry.")
        else:
            print(f"{_ERROR} -- Submission failed: {output}")
        return False

    stdout = out.stdout.strip()
    print('\n'.join(stdout.split('\n')[-4:]))
    print(f"{_INFO} -- Submitted {njobs} Phase 2 jobs.")
    return True


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
