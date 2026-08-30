#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from utilities import create_new_project, check_project_status, launch_jobsub, check_git_branch
from typing import Optional

def main(
    project_dir : str,
    experiment : str,
    create_project : bool,
    launch_jobs : int = None,
    njobs_per_sample : int = None,
    test_job : bool = False,
    tml : str = None,
    batch_size : int = None,
    systematic : str = None,
    tag : str = 'develop',
    memory : int = 1800,
    disk : Optional[int] = None,
    lifetime : str = '1h',
    force : bool = False,
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
    njobs_per_sample : int
        Number of jobs to launch per sample. Mutually exclusive with
        launch_jobs and test_job. Requires a project database created
        with per-sample job tracking -- recreate the project if it
        predates this option.
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
    force : bool
        If True, overwrite pre-existing output files at the destination
        (e.g. left behind by a prior failed or resubmitted attempt at the
        same job) instead of failing the copy-back. Off by default, since
        silently overwriting could mask two jobs unexpectedly racing to
        write the same output.

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
        create_new_project(project_dir, tml, batch_size, systematic, experiment=experiment)
        print(f"[INFO] -- Created new project in {project_dir}")

    # If the project exists, do a check of the project status.
    if project_exists:
        check_project_status(project_dir)

    # If the user requested a test job, run a single job to test the
    # configuration.
    if test_job:
        if not project_exists:
            raise FileNotFoundError(f"Project database {project_dir / 'project.db'} does not exist. Please create a new project first.")
        launch_jobsub(project_dir, experiment, njobs=1, tag=tag, memory=memory, disk=disk, lifetime=lifetime, force=force)

    # If the user requested to launch jobs, do so.
    elif launch_jobs is not None:
        if not project_exists:
            raise FileNotFoundError(f"Project database {project_dir / 'project.db'} does not exist. Please create a new project first.")
        launch_jobsub(project_dir, experiment, njobs=launch_jobs, tag=tag, memory=memory, disk=disk, lifetime=lifetime, force=force)

    # If the user requested to launch jobs per sample, do so.
    elif njobs_per_sample is not None:
        if not project_exists:
            raise FileNotFoundError(f"Project database {project_dir / 'project.db'} does not exist. Please create a new project first.")
        launch_jobsub(project_dir, experiment, njobs_per_sample=njobs_per_sample, tag=tag, memory=memory, disk=disk, lifetime=lifetime, force=force)

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

    # The --njobs-per-sample flag launches up to N jobs per sample rather
    # than N jobs total, ensuring every sample gets jobs launched even
    # when the project has an uneven sample mix.
    p.add_argument(
        '--njobs-per-sample', type=int, default=None, metavar='N',
        help='Launch at most N jobs per sample (ensures every sample gets jobs launched '
             'even when the project has an uneven sample mix). Mutually exclusive with '
             '--launch-jobs/--test-job. Requires a project database created with '
             'per-sample job tracking -- recreate the project if it predates this option.'
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
        '--force', action='store_true',
        help='Overwrite pre-existing output files at the destination (e.g. left behind by a '
             'prior failed or resubmitted attempt at the same job) instead of failing the '
             'copy-back. Off by default, since silently overwriting could mask two jobs '
             'unexpectedly racing to write the same output.'
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

    # Conditional requirement: flags --test-job, --launch-jobs, and
    # --njobs-per-sample are mutually exclusive.
    exclusive_flags = [args.test_job, args.launch_jobs is not None, args.njobs_per_sample is not None]
    if sum(bool(x) for x in exclusive_flags) > 1:
        p.error('--test-job, --launch-jobs, and --njobs-per-sample are mutually exclusive.')

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
        njobs_per_sample=args.njobs_per_sample,
        test_job=args.test_job,
        tml=args.toml,
        batch_size=args.batch_size,
        systematic=args.systematic,
        tag=args.tag,
        memory=args.memory,
        disk=args.disk,
        lifetime=args.lifetime,
        force=args.force,
    )
