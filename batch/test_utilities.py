"""
Tests for medulla/batch/utilities.py job creation and per-sample launch
logic that isn't already covered by test_campaign.py's campaign-level
tests.
"""
import subprocess
import sqlite3
import textwrap
from pathlib import Path
from unittest import mock

import pytest

from utilities import create_new_project, launch_jobsub


def _write_selection_toml(path: Path):
    path.write_text(textwrap.dedent("""\
        [general]
        output = "test_output"

        [[sample]]
        name = "placeholder"
        path = "/fake/placeholder.root"
        ismc = true

        [[tree]]
        name = "selected"
        sim_only = false
        mode = "reco"
        cut = []
        branch = []
    """))
    return path


@pytest.fixture
def uneven_project(tmp_path):
    """A project with two samples of very different pending job counts:
    5 chunks of 'sbnd_mc' and 1 chunk of 'sbnd_offbeam'."""
    tml = _write_selection_toml(tmp_path / "selection.toml")
    fake_samples = (
        [{"name": "sbnd_mc", "path": [f"/fake/mc_{i}.root"], "ismc": True, "disable": False}
         for i in range(5)]
        + [{"name": "sbnd_offbeam", "path": ["/fake/offbeam_0.root"], "ismc": False, "disable": False}]
    )
    project_dir = tmp_path / "project"
    with mock.patch("utilities.get_samples", return_value=fake_samples):
        create_new_project(project_dir, str(tml), batch_size=1)
    return project_dir


class TestLaunchJobsubPerSample:
    """launch_jobsub(njobs_per_sample=...) submits one jobsub_submit call
    per sample, each capped/clamped independently."""

    def test_one_call_per_sample_with_clamped_counts(self, uneven_project, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        real_run = subprocess.run
        calls = []

        def fake_run(cmd, **kwargs):
            if cmd[0] == "jobsub_submit":
                calls.append(cmd)
                return subprocess.CompletedProcess(cmd, 0, stdout="job id 12345.0@fnal.gov\n", stderr="")
            return real_run(cmd, **kwargs)

        with mock.patch("subprocess.run", side_effect=fake_run):
            ok = launch_jobsub(uneven_project, njobs_per_sample=2, confirm=False)

        assert ok is True
        assert len(calls) == 2
        by_sample = {}
        for cmd in calls:
            sample_arg = [a for a in cmd if a.startswith("--sample=")][0].split("=", 1)[1]
            n_idx = cmd.index("-N")
            by_sample[sample_arg] = int(cmd[n_idx + 1])
        # sbnd_mc has 5 pending chunks, clamped to njobs_per_sample=2;
        # sbnd_offbeam only has 1 pending chunk, so it's clamped to 1.
        assert by_sample == {"sbnd_mc": 2, "sbnd_offbeam": 1}

    def test_sample_with_no_pending_jobs_is_skipped(self, uneven_project, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        conn = sqlite3.connect(str(uneven_project / "project.db"))
        conn.execute("UPDATE jobs SET status = 'completed' WHERE sample = 'sbnd_offbeam'")
        conn.commit()
        conn.close()

        real_run = subprocess.run
        calls = []

        def fake_run(cmd, **kwargs):
            if cmd[0] == "jobsub_submit":
                calls.append(cmd)
                return subprocess.CompletedProcess(cmd, 0, stdout="job id 12345.0@fnal.gov\n", stderr="")
            return real_run(cmd, **kwargs)

        with mock.patch("subprocess.run", side_effect=fake_run):
            ok = launch_jobsub(uneven_project, njobs_per_sample=2, confirm=False)

        assert ok is True
        assert len(calls) == 1
        assert any(a == "--sample=sbnd_mc" for a in calls[0])

    def test_partial_failure_still_returns_true(self, uneven_project, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        real_run = subprocess.run

        def fake_run(cmd, **kwargs):
            if cmd[0] == "jobsub_submit":
                if "--sample=sbnd_offbeam" in cmd:
                    raise subprocess.CalledProcessError(1, cmd, output="", stderr="boom")
                return subprocess.CompletedProcess(cmd, 0, stdout="job id 12345.0@fnal.gov\n", stderr="")
            return real_run(cmd, **kwargs)

        with mock.patch("subprocess.run", side_effect=fake_run):
            ok = launch_jobsub(uneven_project, njobs_per_sample=2, confirm=False)

        # sbnd_mc succeeded even though sbnd_offbeam failed.
        assert ok is True

    def test_missing_sample_column_raises(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        project_dir = tmp_path / "legacy_project"
        project_dir.mkdir()
        conn = sqlite3.connect(str(project_dir / "project.db"))
        conn.execute("CREATE TABLE configuration (jobid INTEGER PRIMARY KEY, cfg TEXT NOT NULL)")
        conn.execute("CREATE TABLE jobs (jobid INTEGER PRIMARY KEY, status TEXT)")
        conn.execute("INSERT INTO configuration (jobid, cfg) VALUES (0, 'sample = [{name=\"x\"}]')")
        conn.execute("INSERT INTO jobs (jobid, status) VALUES (0, 'pending')")
        conn.commit()
        conn.close()

        with pytest.raises(RuntimeError):
            launch_jobsub(project_dir, njobs_per_sample=1, confirm=False)

    def test_njobs_and_njobs_per_sample_mutually_exclusive(self, uneven_project):
        with pytest.raises(ValueError):
            launch_jobsub(uneven_project, njobs=5, njobs_per_sample=2, confirm=False)


def _capture_jobsub(project_dir, **kwargs):
    """Run launch_jobsub with jobsub_submit intercepted, returning the
    result and every jobsub_submit command it would have run."""
    real_run = subprocess.run
    calls = []

    def fake_run(cmd, **kw):
        if cmd[0] == "jobsub_submit":
            calls.append(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout="job id 12345.0@fnal.gov\n", stderr="")
        return real_run(cmd, **kw)

    with mock.patch("subprocess.run", side_effect=fake_run):
        ok = launch_jobsub(project_dir, confirm=False, **kwargs)
    return ok, calls


def _n_processes(cmd):
    return int(cmd[cmd.index("-N") + 1])


def _script_arg(cmd, name):
    """Value of submit.sh's own --name=value argument. Only the arguments
    after the '--' separator belong to submit.sh; everything before it is
    jobsub_submit's."""
    tail = cmd[cmd.index("--") + 1:]
    found = [a.split("=", 1)[1] for a in tail if a.startswith(f"--{name}=")]
    return found[0] if found else None


class TestLaunchJobsubJobsPerProcess:
    """launch_jobsub(jobs_per_process=M) packs M job IDs into each grid
    process. -N then counts processes, while submit.sh is told both M and
    the total number of job IDs, so the last process stops at the requested
    count instead of claiming job IDs nobody asked for."""

    def test_default_is_one_jobid_per_process(self, uneven_project, tmp_path, monkeypatch):
        """Unchanged behaviour unless asked: one process per job ID."""
        monkeypatch.chdir(tmp_path)
        ok, calls = _capture_jobsub(uneven_project, njobs=4)

        assert ok is True
        assert len(calls) == 1
        assert _n_processes(calls[0]) == 4
        assert _script_arg(calls[0], "jobs-per-process") == "1"
        assert _script_arg(calls[0], "total-jobs") == "4"

    def test_process_count_rounds_up(self, uneven_project, tmp_path, monkeypatch):
        """5 job IDs at 2 per process need 3 processes. The third holds
        only one job ID, which submit.sh enforces via --total-jobs."""
        monkeypatch.chdir(tmp_path)
        ok, calls = _capture_jobsub(uneven_project, njobs=5, jobs_per_process=2)

        assert ok is True
        assert _n_processes(calls[0]) == 3
        assert _script_arg(calls[0], "jobs-per-process") == "2"
        assert _script_arg(calls[0], "total-jobs") == "5"

    def test_all_pending_with_jobs_per_process(self, uneven_project, tmp_path, monkeypatch):
        """The project has 6 pending job IDs; at 4 per process that is 2."""
        monkeypatch.chdir(tmp_path)
        ok, calls = _capture_jobsub(uneven_project, jobs_per_process=4)

        assert _n_processes(calls[0]) == 2
        assert _script_arg(calls[0], "total-jobs") == "6"

    def test_per_sample_counts_processes_per_sample(self, uneven_project, tmp_path, monkeypatch):
        """Each sample's own job-ID count is divided independently."""
        monkeypatch.chdir(tmp_path)
        ok, calls = _capture_jobsub(uneven_project, njobs_per_sample=5, jobs_per_process=2)

        by_sample = {
            _script_arg(c, "sample"): (_n_processes(c), _script_arg(c, "total-jobs"))
            for c in calls
        }
        # sbnd_mc: 5 job IDs -> 3 processes; sbnd_offbeam: 1 job ID -> 1.
        assert by_sample == {"sbnd_mc": (3, "5"), "sbnd_offbeam": (1, "1")}

    @pytest.mark.parametrize("bad", [0, -1])
    def test_nonpositive_jobs_per_process_raises(self, uneven_project, bad):
        with pytest.raises(ValueError):
            launch_jobsub(uneven_project, jobs_per_process=bad, confirm=False)


class TestLaunchJobsubShipsMacros:
    """The job's ROOT macros travel with the job instead of being read out
    of the tagged checkout the job builds."""

    def test_macros_are_transferred_with_the_job(self, uneven_project, tmp_path, monkeypatch):
        """A release tag that predates these macros would otherwise fail every
        job: validation would copy nothing back, and the input check would have
        nothing to run."""
        monkeypatch.chdir(tmp_path)
        ok, calls = _capture_jobsub(uneven_project, njobs=1)
        cmd = calls[0]

        f_values = [cmd[i + 1] for i, a in enumerate(cmd) if a == "-f"]
        shipped = {Path(v[len("dropbox://"):]).name
                   for v in f_values if v.startswith("dropbox://")}
        assert {"validate_pair.C", "check_inputs.C"} <= shipped
        for v in f_values:
            assert Path(v[len("dropbox://"):]).is_file()

        # jobsub_submit only honours options that precede the executable.
        exe = next(i for i, a in enumerate(cmd) if a.startswith("file://") and a.endswith("submit.sh"))
        assert all(i < exe for i, a in enumerate(cmd) if a == "-f")

    def test_missing_macro_fails_before_submitting(self, uneven_project, tmp_path, monkeypatch):
        """Better to refuse at launch than to discover on the grid, after
        every job has spent its event loop, that a macro is missing."""
        monkeypatch.chdir(tmp_path)
        with mock.patch("utilities.JOB_MACRO_PATHS", (tmp_path / "missing.C",)):
            with pytest.raises(FileNotFoundError):
                _capture_jobsub(uneven_project, njobs=1)
