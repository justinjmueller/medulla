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
