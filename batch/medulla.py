#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from utilities import (
    create_new_project, check_project_status, launch_jobsub, check_git_branch,
    run_variation_phase1_interactive, launch_variation_phase1_jobsub,
    launch_variation_phase2_jobsub,
)
from typing import Optional

def main(
    project_dir : str,
    experiment : str,
    create_project : bool,
    launch_jobs : int = None,
    test_job : bool = False,
    tml : str = None,
    batch_size : int = None,
    systematic : str = None,
    tag : str = 'develop',
    memory : int = 1800,
    disk : Optional[int] = None,
    lifetime : str = '1h',
    relaunch_missing : bool = False,
    variation_systematics : str = None,
    variation_input : str = None,
    variation_phase1 : bool = False,
    variation_phase2 : bool = False,
    variation_interactive : bool = False,
    variation_test : bool = False,
):
    """
    Main function to run the medulla script.

    Parameters
    ----------
    project_dir : str
        Path to the base directory for the job directory.
    experiment : str
        Experiment name (default: sbnd).
    create_project : bool
        Whether to create a new project. If this is True, then tml and
        batch_size must be provided.
    launch_jobs : int
        Number of jobs to launch. If None, launch all pending jobs.
    test_job : bool
        Whether to run a single test job to verify the configuration.
    tml : str
        Path to the TOML file containing the configuration.
    batch_size : int
        Number of files to process in each batch.
    systematic : str
        Path to the systematic template file to use. If None, use the
        default template file in the batch directory.
    tag : str
        Tag to use for the medulla repository (defaults to develop).
    memory : int
        Amount of memory to request for each job in MB. 
    disk : int | None
        Amount of disk to request for each job in GB. If None, use default.
    lifetime : str
        Expected lifetime of each job (e.g., '1h', '30m').
    relaunch_missing : bool
        If True, relaunch jobs that are in project.db but missing output files.
    variation_systematics : str
        Path to the variation systematics TOML. Used with --variation-phase1
        and --variation-phase2.
    variation_input : str
        PNFS or xrootd path to the merged variation+CV ROOT file. Required
        for --variation-phase1.
    variation_phase1 : bool
        Build detector-variation splines from the merged variation+CV input.
        Runs interactively if variation_interactive is True, otherwise submits
        a single batch job.
    variation_phase2 : bool
        Apply pre-built splines to each individual selection output file via
        parallel batch jobs. Requires variation_splines.root to already exist
        in the project directory (i.e. Phase 1 must have completed first).
    variation_interactive : bool
        When True and variation_phase1 is set, runs Phase 1 on the current
        node instead of submitting a batch job.

    Returns
    -------
    None.
    """
    # Check if the project file already exists. This will be used to
    # inform logic later on in the function.
    project_dir = Path(project_dir)
    project_exists = (project_dir / 'project.db').exists()

    # Create a new project if requested.
    if create_project:
        if tml is None:
            raise ValueError("TOML file must be provided when creating a new project.")
        if batch_size is None:
            raise ValueError("Batch size must be provided when creating a new project.")
        if project_exists:
            raise FileExistsError(f"Project database {project_dir / 'project.db'} already exists.")
        create_new_project(project_dir, tml, batch_size, systematic)
        print(f"[INFO] -- Created new project in {project_dir}")

    # If the project exists, do a check of the project status.
    if project_exists:
        check_project_status(project_dir)

    # If the user requested a test job, run a single job to test the
    # configuration.
    if test_job:
        if not project_exists:
            raise FileNotFoundError(f"Project database {project_dir / 'project.db'} does not exist. Please create a new project first.")
        launch_jobsub(project_dir, experiment, njobs=1, tag=tag, memory=memory, disk=disk, lifetime=lifetime, relaunch_missing=relaunch_missing)

    # If the user requested to launch jobs, do so.
    if launch_jobs is not None:
        if not project_exists:
            raise FileNotFoundError(f"Project database {project_dir / 'project.db'} does not exist. Please create a new project first.")
        launch_jobsub(project_dir, experiment, njobs=launch_jobs, tag=tag, memory=memory, disk=disk, lifetime=lifetime, relaunch_missing=relaunch_missing)

    # Phase 1: build detector-variation splines.
    if variation_phase1:
        if not project_exists:
            raise FileNotFoundError(f"Project database {project_dir / 'project.db'} does not exist. Please create a new project first.")
        if variation_interactive:
            run_variation_phase1_interactive(
                project_dir=project_dir,
                variation_input=variation_input,
                variation_toml=variation_systematics,
            )
        else:
            launch_variation_phase1_jobsub(
                project_dir=project_dir,
                variation_input=variation_input,
                variation_toml=variation_systematics,
                exp=experiment,
                branch=tag,
                memory=memory,
                disk=disk,
                lifetime=lifetime,
            )

    # Phase 2: apply pre-built splines to individual selection outputs.
    if variation_phase2:
        if not project_exists:
            raise FileNotFoundError(f"Project database {project_dir / 'project.db'} does not exist. Please create a new project first.")
        launch_variation_phase2_jobsub(
            project_dir=project_dir,
            variation_toml=variation_systematics,
            exp=experiment,
            branch=tag,
            memory=memory,
            disk=disk,
            lifetime=lifetime,
            njobs=1 if variation_test else -1,
        )

if __name__ == '__main__':
    p = ArgumentParser(description='Run medulla.')

    # The project directory is always required.
    p.add_argument(
        '--project-dir', '-p', type=str, required=True,
        help='Path to the base directory for the job directory.'
    )

    # The experiment is always required.
    p.add_argument(
        '--experiment', '-e', type=str, required=False, default='sbnd',
        help='Experiment name (default: sbnd).'
    )

    # The --create-project flag indicates that a new project should be
    # created. If this flag is set, then --toml and --batch-size are
    # required.
    p.add_argument(
        '--create-project', '-c', action='store_true',
        help='Create a new project for medulla job submission.'
    )
    p.add_argument(
        '--toml', '-t', type=str,
        help='Path to the TOML file containing the configuration.'
    )
    p.add_argument(
        '--batch-size', '-b', type=int,
        help='Number of files to process in each batch (only used when creating a new project).'
    )

    # The --systematic flag indicates the systematic template file to
    # use. The default is the template file in the batch directory.
    # This is used when creating a new project.
    p.add_argument(
        '--systematic', '-s', type=str, default=None,
        help='Path to the systematic template file to use.'
    )

    # The --test-job flag indicates that a test job should be run. That
    # is, launch a single job to test the configuration.
    p.add_argument(
        '--test-job', '-T', action='store_true',
        help='Run a single test job to verify the configuration.'
    )

    # The --launch-jobs flag allows an (optional) number of jobs to be
    # launched. If this flag is provided without a number, then all
    # pending jobs will be launched.
    p.add_argument(
        '--launch-jobs', '-l', type=int, nargs='?', const=-1,
        help='Launch the specified number of jobs. If no number is provided, '
             'then all pending jobs will be launched.'
    )

    p.add_argument(
        '--relaunch-missing', '-r', action='store_true',
        help='Also relaunch jobs tracked in project.db but missing output files.'
    )

    p.add_argument(
        '--tag', type=str, default='develop',
        help='Tag to use for the medulla repository (defaults to develop).'
    )

    p.add_argument(
        '--memory', '-m', type=int, default=1800,
        help='Amount of memory to request for each job in MB (default: 1800).'
    )

    p.add_argument(
        '--disk', '-d', type=int, default=None,
        help='Amount of disk to request for each job in GB (default: None).'
    )

    p.add_argument(
        '--lifetime', '-f', type=str, default='1h',
        help="Expected lifetime of each job (e.g., '1h', '30m') (default: '1h')."
    )

    p.add_argument(
        '--variation-systematics', '-V', type=str, default=None,
        help='Path to the variation systematics TOML. Required with '
             '--variation-phase1 and --variation-phase2.'
    )

    p.add_argument(
        '--variation-input', '-I', type=str, default=None,
        help='PNFS or xrootd path to the merged variation+CV ROOT file '
             '(required with --variation-phase1).'
    )

    p.add_argument(
        '--variation-phase1', action='store_true',
        help='Build detector-variation splines from the merged variation+CV '
             'input (requires --variation-systematics and --variation-input). '
             'Submits a single batch job unless --variation-interactive is set.'
    )

    p.add_argument(
        '--variation-phase2', action='store_true',
        help='Apply pre-built splines to each individual selection output '
             'file via parallel batch jobs (requires --variation-systematics '
             'and a completed Phase 1 run).'
    )

    p.add_argument(
        '--variation-interactive', action='store_true',
        help='Run Phase 1 on the current node instead of submitting a batch '
             'job (only valid with --variation-phase1).'
    )

    p.add_argument(
        '--variation-test', action='store_true',
        help='Submit only one Phase 2 job as a test (only valid with '
             '--variation-phase2).'
    )

    args = p.parse_args()

    # Requirement: the experiment must be sbnd or icarus.
    if args.experiment not in ['sbnd', 'icarus']:
        p.error("Experiment must be either 'sbnd' or 'icarus'.")

    # Conditional requirement: if --create-project is set, then --toml
    # and --batch-size are required.
    if args.create_project and args.toml is None:
        p.error('--toml is required when --create-project is set.')
    if args.create_project and args.batch_size is None:
        p.error('--batch-size is required when --create-project is set.')

    # Conditional requirement: flags --test-job and --launch-jobs are
    # mutually exclusive.
    if args.test_job and args.launch_jobs is not None:
        p.error('--test-job and --launch-jobs are mutually exclusive.')

    # Conditional requirements for variation phases.
    if args.variation_phase1 and args.variation_phase2:
        p.error('--variation-phase1 and --variation-phase2 are mutually exclusive.')
    if args.variation_phase1 and args.variation_systematics is None:
        p.error('--variation-systematics is required with --variation-phase1.')
    if args.variation_phase1 and args.variation_input is None:
        p.error('--variation-input is required with --variation-phase1.')
    if args.variation_phase2 and args.variation_systematics is None:
        p.error('--variation-systematics is required with --variation-phase2.')
    if args.variation_interactive and not args.variation_phase1:
        p.error('--variation-interactive is only valid with --variation-phase1.')
    if args.variation_test and not args.variation_phase2:
        p.error('--variation-test is only valid with --variation-phase2.')

    if args.tag != 'develop':
        print(f"[INFO] -- Using tag '{args.tag}' for medulla repository.")
        if not check_git_branch(args.tag):
            p.error(f"Tag '{args.tag}' does not exist in the medulla repository.")

    # Run the main function.
    main(
        project_dir=args.project_dir,
        experiment=args.experiment,
        create_project=args.create_project,
        launch_jobs=args.launch_jobs,
        test_job=args.test_job,
        tml=args.toml,
        batch_size=args.batch_size,
        systematic=args.systematic,
        tag=args.tag,
        memory=args.memory,
        disk=args.disk,
        lifetime=args.lifetime,
        relaunch_missing=args.relaunch_missing,
        variation_systematics=args.variation_systematics,
        variation_input=args.variation_input,
        variation_phase1=args.variation_phase1,
        variation_phase2=args.variation_phase2,
        variation_interactive=args.variation_interactive,
        variation_test=args.variation_test,
    )