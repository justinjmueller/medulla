"""
Tests for the bucketed output layout (utilities.output_bucket and friends)
and for every reader handling both it and the legacy flat layout.

The layout matters because a flat output/ directory does not survive a large
campaign, and because the change has to be invisible to projects written
before it: a reader that looked only in the new place would report finished
jobs as missing and resubmit them.
"""
import sqlite3
import threading
from pathlib import Path
from unittest import mock

import pytest

from utilities import (output_bucket, output_path_for, find_output_file,
                       iter_output_files, output_globs, jobid_of,
                       survey_project_output, OUTPUT_BUCKET_SIZE)
from campaign import (_resolve_hadd_files, _run_one_scan, cmd_scan,
                      SCHEMA_CAMPAIGN_META, SCHEMA_PROJECTS)


def _touch(path: Path, size: int = 2048) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"x" * size)
    return path


def _flat(proj: Path, jid: int, kind: str) -> Path:
    name = {"nosyst": f"output_jobid{jid:04d}.root",
            "wsyst": f"output_systematics_jobid{jid:04d}.root",
            "record": f"validation_jobid{jid:04d}.txt"}[kind]
    return proj / "output" / ("val" if kind == "record" else "") / name


# ---------------------------------------------------------------------------
# The layout itself
# ---------------------------------------------------------------------------

class TestBucket:

    @pytest.mark.parametrize("jid, bucket", [
        (0, "000"), (999, "000"), (1000, "001"), (1999, "001"), (123456, "123"),
    ])
    def test_bucket_is_jobid_over_bucket_size(self, jid, bucket):
        assert OUTPUT_BUCKET_SIZE == 1000
        assert output_bucket(jid) == bucket

    def test_submit_sh_uses_the_same_bucket_arithmetic(self):
        """submit.sh computes the bucket independently in bash; if the two
        ever disagree, jobs write where no reader looks."""
        text = (Path(__file__).parent / "submit.sh").read_text()
        assert 'printf -v BUCKET "%03d" $(( JOBID / 1000 ))' in text

    def test_paths(self, tmp_path):
        assert output_path_for(tmp_path, 42, "nosyst") == tmp_path / "output/000/output_jobid0042.root"
        assert output_path_for(tmp_path, 1042, "wsyst") == tmp_path / "output/001/output_systematics_jobid1042.root"
        assert output_path_for(tmp_path, 42, "record") == tmp_path / "output/val/000/validation_jobid0042.txt"

    def test_five_digit_jobids_keep_their_full_number(self, tmp_path):
        """%04d is a minimum width; job 12345 must not be truncated."""
        assert output_path_for(tmp_path, 12345, "nosyst").name == "output_jobid12345.root"
        assert jobid_of("output_jobid12345.root") == 12345

    def test_unknown_kind_raises(self, tmp_path):
        with pytest.raises(ValueError):
            output_path_for(tmp_path, 1, "nonsense")


class TestReadersSeeBothLayouts:

    def test_iter_finds_flat_and_bucketed(self, tmp_path):
        _touch(_flat(tmp_path, 7, "nosyst"))
        _touch(output_path_for(tmp_path, 1500, "nosyst"))
        found = iter_output_files(tmp_path, "nosyst")
        assert set(found) == {7, 1500}

    def test_bucketed_wins_when_both_exist(self, tmp_path):
        flat = _touch(_flat(tmp_path, 7, "nosyst"))
        nested = _touch(output_path_for(tmp_path, 7, "nosyst"))
        assert iter_output_files(tmp_path, "nosyst")[7] == nested
        assert find_output_file(tmp_path, 7, "nosyst") == nested
        assert flat.exists()  # nothing was deleted

    def test_find_falls_back_to_flat_then_none(self, tmp_path):
        flat = _touch(_flat(tmp_path, 7, "wsyst"))
        assert find_output_file(tmp_path, 7, "wsyst") == flat
        assert find_output_file(tmp_path, 8, "wsyst") is None

    def test_kinds_do_not_leak_into_each_other(self, tmp_path):
        """'output_jobid*' must not pick up systematics files, and neither
        must pick up records sitting under output/val/."""
        _touch(output_path_for(tmp_path, 3, "nosyst"))
        _touch(output_path_for(tmp_path, 3, "wsyst"))
        _touch(output_path_for(tmp_path, 3, "record"))
        assert [p.name for p in iter_output_files(tmp_path, "nosyst").values()] == ["output_jobid0003.root"]
        assert [p.name for p in iter_output_files(tmp_path, "wsyst").values()] == ["output_systematics_jobid0003.root"]
        assert [p.name for p in iter_output_files(tmp_path, "record").values()] == ["validation_jobid0003.txt"]

    def test_records_in_both_layouts(self, tmp_path):
        _touch(_flat(tmp_path, 1, "record"), 10)
        _touch(output_path_for(tmp_path, 2001, "record"), 10)
        assert set(iter_output_files(tmp_path, "record")) == {1, 2001}


class TestSurveyAcrossLayouts:

    def test_mixed_project_classifies_by_pair(self, tmp_path):
        """A project mid-transition: an old flat pair, a new bucketed pair,
        a bucketed orphan, and a pair split across the two layouts."""
        _touch(_flat(tmp_path, 0, "nosyst")); _touch(_flat(tmp_path, 0, "wsyst"))
        _touch(output_path_for(tmp_path, 1500, "nosyst")); _touch(output_path_for(tmp_path, 1500, "wsyst"))
        _touch(output_path_for(tmp_path, 2000, "nosyst"))
        _touch(_flat(tmp_path, 3000, "nosyst")); _touch(output_path_for(tmp_path, 3000, "wsyst"))

        s = survey_project_output(tmp_path)

        assert s["completed"] == [0, 1500, 3000]
        assert s["orphaned"] == [2000]


# ---------------------------------------------------------------------------
# finalize: hadd input resolution
# ---------------------------------------------------------------------------

class TestHaddAcrossLayouts:

    def test_glob_list_spans_both_layouts(self, tmp_path):
        a = _touch(_flat(tmp_path, 0, "nosyst"))
        b = _touch(output_path_for(tmp_path, 1500, "nosyst"))
        row = {"project_dir": str(tmp_path)}
        files = _resolve_hadd_files(row, output_globs(tmp_path, "nosyst"), "nosyst")
        assert files == sorted([str(a), str(b)])

    def test_single_glob_string_still_works(self, tmp_path):
        a = _touch(_flat(tmp_path, 0, "nosyst"))
        row = {"project_dir": str(tmp_path)}
        assert _resolve_hadd_files(row, str(tmp_path / "output" / "output_jobid*.root"), "nosyst") == [str(a)]

    def test_job_ids_resolve_in_either_layout(self, tmp_path):
        a = _touch(_flat(tmp_path, 0, "wsyst"))
        b = _touch(output_path_for(tmp_path, 1500, "wsyst"))
        row = {"project_dir": str(tmp_path)}
        files = _resolve_hadd_files(row, None, "wsyst", job_ids=[0, 1500, 9999])
        assert files == sorted([str(a), str(b)])


# ---------------------------------------------------------------------------
# scan: both halves checked, partner removed with a corrupt file
# ---------------------------------------------------------------------------

def _fake_root_marking(bad_suffix):
    """A stand-in for the batched ROOT check that reports as corrupt every
    listed file whose name ends with bad_suffix."""
    def run(cmd, **kwargs):
        call = cmd[-1]
        args = call[call.index("(") + 1:call.rindex(")")]
        filelist, report = (a.strip().strip('"') for a in args.split(","))
        lines = Path(filelist).read_text().splitlines()
        Path(report).write_text("".join(l + "\n" for l in lines if l.endswith(bad_suffix)))
        return mock.Mock(returncode=0, stderr="")
    return run


class TestScanAcrossLayouts:

    def test_systematics_files_are_checked_too(self, tmp_path):
        proj = tmp_path / "proj"
        for jid in (0, 1500):
            _touch(output_path_for(proj, jid, "nosyst"))
            _touch(output_path_for(proj, jid, "wsyst"))
        with mock.patch("subprocess.run", side_effect=_fake_root_marking("output_systematics_jobid1500.root")):
            r = _run_one_scan(tmp_path / "campaign", "proj", proj, threading.Lock())
        assert r["n_checked"] == 4
        assert [p.name for p in r["bad_files"]] == ["output_systematics_jobid1500.root"]

    def _campaign(self, tmp_path):
        campaign_dir = tmp_path / "campaign"
        campaign_dir.mkdir()
        proj = campaign_dir / "eps_primary_sbnd"
        conn = sqlite3.connect(campaign_dir / "campaign.db")
        conn.executescript(SCHEMA_CAMPAIGN_META + SCHEMA_PROJECTS)
        conn.execute("INSERT INTO campaign_meta (name, tag) VALUES ('c1', 'v1')")
        conn.execute(
            "INSERT INTO projects (analysis, role, experiment, toml_file, project_dir, "
            "batch_size, n_jobs, n_completed, status) "
            "VALUES ('eps', 'primary', 'sbnd', 'x.toml', ?, 25, 2, 2, 'completed')", (str(proj),))
        conn.commit(); conn.close()
        (proj / "output").mkdir(parents=True)
        p = sqlite3.connect(proj / "project.db")
        p.execute("CREATE TABLE configuration (jobid INTEGER PRIMARY KEY, cfg TEXT NOT NULL)")
        p.execute("CREATE TABLE jobs (jobid INTEGER PRIMARY KEY, status TEXT, sample TEXT)")
        for jid in (0, 1500):
            p.execute("INSERT INTO configuration VALUES (?, '')", (jid,))
            p.execute("INSERT INTO jobs VALUES (?, 'completed', 's')", (jid,))
        p.commit(); p.close()
        _touch(_flat(proj, 0, "nosyst")); _touch(_flat(proj, 0, "wsyst"))
        _touch(output_path_for(proj, 1500, "nosyst")); _touch(output_path_for(proj, 1500, "wsyst"))
        return campaign_dir, proj

    class _Args:
        name = None; experiment = None; dry_run = False; workers = 1

    def test_corrupt_systematics_file_removes_its_pair_and_reverts(self, tmp_path):
        campaign_dir, proj = self._campaign(tmp_path)
        args = self._Args(); args.campaign = str(campaign_dir)
        with mock.patch("subprocess.run", side_effect=_fake_root_marking("output_systematics_jobid1500.root")), \
             mock.patch("shutil.which", return_value="/usr/bin/root"), \
             mock.patch("builtins.input", return_value="y"):
            cmd_scan(args)

        # The corrupt file and its partner are both gone; the other pair is intact.
        assert find_output_file(proj, 1500, "wsyst") is None
        assert find_output_file(proj, 1500, "nosyst") is None
        assert find_output_file(proj, 0, "nosyst") is not None
        assert find_output_file(proj, 0, "wsyst") is not None

        conn = sqlite3.connect(proj / "project.db")
        statuses = dict(conn.execute("SELECT jobid, status FROM jobs").fetchall())
        conn.close()
        assert statuses == {0: "completed", 1500: "pending"}
