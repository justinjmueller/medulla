"""
Tests for medulla/batch/utilities.py job creation and per-sample launch
logic that isn't already covered by test_campaign.py's campaign-level
tests.
"""
import gzip
import shutil
import subprocess
import sqlite3
import textwrap
from pathlib import Path
from unittest import mock

import pytest

from utilities import (
    SUBMISSION_DB_NAME, create_new_project, create_systematics_cfg,
    launch_jobsub, write_submission_db,
)


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


def _dropbox_paths(cmd):
    return [Path(cmd[i + 1][len("dropbox://"):])
            for i, a in enumerate(cmd)
            if a == "-f" and cmd[i + 1].startswith("dropbox://")]


def _capture_jobsub(project_dir, shipped=None, **kwargs):
    """Run launch_jobsub with jobsub_submit intercepted, returning the
    result and every jobsub_submit command it would have run.

    Dropbox files only exist while jobsub_submit runs -- launch_jobsub
    removes the job database snapshots afterwards -- so if `shipped` is a
    directory, each call's shipped database is copied there as
    <call index>.db.gz for inspection, and each call's dropbox paths are
    checked to exist at the moment of submission."""
    real_run = subprocess.run
    calls = []

    def fake_run(cmd, **kw):
        if cmd[0] == "jobsub_submit":
            for path in _dropbox_paths(cmd):
                assert path.is_file(), f"{path} missing at submission time"
                if shipped is not None and path.name == SUBMISSION_DB_NAME:
                    shutil.copy(path, Path(shipped) / f"{len(calls)}.db.gz")
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

        # _capture_jobsub checks each one existed at submission time.
        shipped = {p.name for p in _dropbox_paths(cmd)}
        assert {"validate_pair.C", "check_inputs.C"} <= shipped

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


def _unpack(gz_path, tmp_path, name):
    """Decompress a shipped job database the way submit.sh does."""
    out = tmp_path / name
    with gzip.open(gz_path, "rb") as fsrc, open(out, "wb") as fdst:
        shutil.copyfileobj(fsrc, fdst)
    return out


def _claims(db, n_processes, m, total, sample=None):
    """Every process's job IDs under submit.sh's claim query: process P takes
    ranks P*M .. P*M+M-1 of the non-completed job IDs, capped at total."""
    conn = sqlite3.connect(db)
    out = []
    for p in range(n_processes):
        offset = p * m
        limit = min(m, total - offset)
        where = "status != 'completed'"
        params = []
        if sample is not None:
            where += " AND sample = ?"
            params.append(sample)
        rows = conn.execute(
            f"SELECT jobid FROM jobs WHERE {where} ORDER BY jobid LIMIT ? OFFSET ?",
            (*params, limit, offset),
        ).fetchall()
        out.append([r[0] for r in rows])
    conn.close()
    return out


class TestShippedJobDatabase:
    """Each submission ships only the rows of project.db its processes can
    claim, instead of every process copying the whole database (2.7 GB for
    a large project) from dCache. The snapshot must give every process
    exactly the job IDs the full database would have."""

    def _complete(self, project_dir, jobids):
        conn = sqlite3.connect(project_dir / "project.db")
        conn.executemany("UPDATE jobs SET status = 'completed' WHERE jobid = ?",
                         [(j,) for j in jobids])
        conn.commit()
        conn.close()

    def test_submission_ships_the_database_and_says_so(self, uneven_project, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        ok, calls = _capture_jobsub(uneven_project, njobs=3)

        assert ok is True
        cmd = calls[0]
        assert SUBMISSION_DB_NAME in {p.name for p in _dropbox_paths(cmd)}
        assert "--shipped-db" in cmd[cmd.index("--") + 1:]

    def test_snapshot_is_removed_after_launch(self, uneven_project, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        ok, calls = _capture_jobsub(uneven_project, njobs=3)

        shipped = next(p for p in _dropbox_paths(calls[0]) if p.name == SUBMISSION_DB_NAME)
        assert not shipped.exists()
        assert not list(tmp_path.glob("medulla_submission_*"))

    def test_every_process_claims_what_the_full_database_gives_it(
            self, uneven_project, tmp_path, monkeypatch):
        """The property that matters: identical claims, process by process,
        with completed job IDs interleaved so the ranks are not just the job
        IDs themselves."""
        monkeypatch.chdir(tmp_path)
        self._complete(uneven_project, [0, 2])
        shipped = tmp_path / "shipped"
        shipped.mkdir()
        ok, calls = _capture_jobsub(uneven_project, shipped=shipped, njobs=3, jobs_per_process=2)

        db = _unpack(shipped / "0.db.gz", tmp_path, "snap.db")
        n = _n_processes(calls[0])
        full = _claims(uneven_project / "project.db", n, 2, 3)
        assert _claims(db, n, 2, 3) == full
        assert full == [[1, 3], [4]]

    def test_snapshot_holds_only_claimable_rows_with_their_configs(
            self, uneven_project, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        self._complete(uneven_project, [0])
        shipped = tmp_path / "shipped"
        shipped.mkdir()
        _capture_jobsub(uneven_project, shipped=shipped, njobs=2)

        db = _unpack(shipped / "0.db.gz", tmp_path, "snap.db")
        snap = sqlite3.connect(db)
        full = sqlite3.connect(uneven_project / "project.db")
        jobs = snap.execute("SELECT jobid, status, sample FROM jobs ORDER BY jobid").fetchall()
        assert jobs == full.execute(
            "SELECT jobid, status, sample FROM jobs WHERE jobid IN (1, 2) ORDER BY jobid").fetchall()
        # Each job's configuration is carried over byte for byte: it is what
        # submit.sh writes out as the job's config.
        for (jobid, _, _) in jobs:
            q = "SELECT cfg FROM configuration WHERE jobid = ?"
            assert snap.execute(q, (jobid,)).fetchone() == full.execute(q, (jobid,)).fetchone()
        assert snap.execute("SELECT COUNT(*) FROM configuration").fetchone()[0] == 2

    def test_per_sample_submissions_each_ship_their_own_sample(
            self, uneven_project, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        shipped = tmp_path / "shipped"
        shipped.mkdir()
        ok, calls = _capture_jobsub(uneven_project, shipped=shipped, njobs_per_sample=5, jobs_per_process=2)

        paths = [next(p for p in _dropbox_paths(c) if p.name == SUBMISSION_DB_NAME) for c in calls]
        assert len(set(paths)) == len(calls)
        for i, cmd in enumerate(calls):
            sample = _script_arg(cmd, "sample")
            total = int(_script_arg(cmd, "total-jobs"))
            db = _unpack(shipped / f"{i}.db.gz", tmp_path, f"snap{i}.db")
            samples = {r[0] for r in sqlite3.connect(db).execute("SELECT sample FROM jobs")}
            assert samples == {sample}
            n = _n_processes(cmd)
            assert (_claims(db, n, 2, total, sample)
                    == _claims(uneven_project / "project.db", n, 2, total, sample))

    def test_legacy_jobs_table_without_sample_column(self, tmp_path):
        """A project predating the sample column still snapshots; the
        missing columns are left NULL."""
        src = tmp_path / "legacy.db"
        conn = sqlite3.connect(src)
        conn.execute("CREATE TABLE configuration (jobid INTEGER PRIMARY KEY, cfg TEXT NOT NULL)")
        conn.execute("CREATE TABLE jobs (jobid INTEGER PRIMARY KEY, status TEXT)")
        conn.executemany("INSERT INTO configuration VALUES (?, ?)", [(i, f"cfg{i}") for i in range(4)])
        conn.executemany("INSERT INTO jobs VALUES (?, ?)",
                         [(0, "completed"), (1, "pending"), (2, "pending"), (3, "pending")])
        conn.commit()
        conn.close()

        dst = tmp_path / SUBMISSION_DB_NAME
        assert write_submission_db(src, dst, 2) == 2
        db = _unpack(dst, tmp_path, "snap.db")
        rows = sqlite3.connect(db).execute(
            "SELECT j.jobid, j.status, j.sample, c.cfg FROM jobs j "
            "JOIN configuration c USING (jobid) ORDER BY jobid").fetchall()
        assert rows == [(1, "pending", None, "cfg1"), (2, "pending", None, "cfg2")]
        assert not (tmp_path / (SUBMISSION_DB_NAME + ".tmp")).exists()


class TestSystematicsSwitch:
    """create_new_project(systematics=False) copies every tree through the
    systematics step, and the validation manifest -- generated from the same
    configuration -- checks them as copies."""

    TREES = [
        {"name": "selected", "add_systematics": True,
         "branch": [{"name": "neutrino_id", "type": "true"},
                    {"name": "neutrino_energy", "type": "mctruth"}]},
        {"name": "signal", "branch": []},
    ]
    SAMPLES = [{"name": "sbnd_mc", "ismc": True}, {"name": "sbnd_offbeam", "ismc": False}]

    def _actions(self, cfg):
        return {t["origin"]: t["action"] for t in cfg["tree"]}

    def _base(self):
        return {"input": {}, "output": {}}

    def test_default_reweights_mc_trees_that_ask_for_it(self):
        cfg = create_systematics_cfg(self._base(), self.TREES, self.SAMPLES)
        assert self._actions(cfg) == {
            "events/sbnd_mc/selected": "add_weights",
            "events/sbnd_mc/signal": "copy",
            "events/sbnd_offbeam/selected": "copy",
            "events/sbnd_offbeam/signal": "copy",
        }

    def test_off_copies_every_tree(self):
        cfg = create_systematics_cfg(self._base(), self.TREES, self.SAMPLES, add_weights=False)
        assert set(self._actions(cfg).values()) == {"copy"}
        assert all("table_types" not in t for t in cfg["tree"])

    def test_project_manifest_follows(self, tmp_path):
        tml = tmp_path / "selection.toml"
        tml.write_text(textwrap.dedent("""\
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
            add_systematics = true
            cut = []

            [[tree.branch]]
            name = "neutrino_id"
            type = "true"

            [[tree.branch]]
            name = "neutrino_energy"
            type = "mctruth"
        """))
        fake = [{"name": "sbnd_detvar_x", "path": ["/fake/x.root"], "ismc": True, "disable": False}]

        manifests = {}
        for on in (True, False):
            project_dir = tmp_path / f"project_{on}"
            with mock.patch("utilities.get_samples", return_value=fake):
                create_new_project(project_dir, str(tml), batch_size=1, systematics=on)
            manifests[on] = (project_dir / "validation_manifest.txt").read_text().split()

        assert manifests[True] == ["events/sbnd_detvar_x/selected|selected|add_weights|multisim,multisigma"]
        assert manifests[False] == ["events/sbnd_detvar_x/selected|selected|copy|"]
