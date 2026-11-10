"""
Tests for utilities.reconcile_project and the reconcile.py entry point.

Reconciliation is the tool that is supposed to make it impossible for part of
a dataset to go missing unnoticed, so the tests build one project engineered
to hit every category at once and check that each job ID lands in exactly one
place, that every input file is accounted for, and that the anomalies nothing
else reports are the ones that surface.
"""
import sqlite3
from pathlib import Path

import pytest
import toml

from utilities import (reconcile_project, parse_validation_record, format_reconcile_report,
                       output_path_for, RECONCILE_CATEGORIES)
import reconcile as reconcile_cli


def _cfg(jid, n_files=3):
    """A stored job configuration shaped like the real ones: the whole
    selection plus one [[sample]] with a path list."""
    return toml.dumps({
        "general": {"output": "output"},
        "tree": [{"name": "selected", "branch": [{"name": "neutrino_id", "type": "true"}]}],
        "sample": [{"name": "sbnd", "ismc": True,
                    "path": [f"/pnfs/fake/job{jid}_file{i}.flat.caf.root" for i in range(n_files)]}],
    })


def _record(jid, status, cluster, process, start, claimed,
            n_inputs=3, n_dropped=0, input_events=30, events_written=None, pot=1e16):
    events_written = input_events if events_written is None else events_written
    return "\n".join([
        f"JOB ({jid},{process}) VALIDATION",
        "TAG=test", "COMMIT=abc", f"CLUSTER={cluster}", f"PROCESS={process}",
        f"PROCESS_START={start}", f"CLAIMED_JOBIDS={claimed}", "SAMPLE=sbnd",
        f"N_INPUTS={n_inputs}", f"N_DROPPED_INPUTS={n_dropped}", f"INPUT_EVENTS={input_events}",
        "SAMPLE=sbnd",                       # the validator repeats it; first must win
        f"POT={pot}", "LIVETIME=0",
        f"TREE_selected=10,10,0", f"TREE_events={events_written},{events_written},-1",
        f"STATUS={status}", "EXIT_CODE=0", "",
    ])


def _write(path: Path, content, binary=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    (path.write_bytes if binary else path.write_text)(content)


@pytest.fixture
def project(tmp_path):
    """
    Job IDs 0-7, three input files each, arranged to cover every category:

      0  complete               bucketed pair, ok record
      1  complete               legacy flat pair, match_partial, 1 dropped input
      2  complete_unrecorded    bucketed pair, no record
      3  complete_stale_record  pair, but the record is from a failed attempt
      4  failed                 no pair, match_unknown record
      5  output_lost            no pair, but an ok record
      6  no_trace               nothing
      7  complete               ok record, but events written != events read

    plus a record for job 99, which is not in project.db.

    Claims: job 1 is claimed by two processes of cluster 100 (duplicate work);
    job 5 by cluster 200 and later cluster 100 (a retry).
    """
    proj = tmp_path / "proj"
    (proj / "output").mkdir(parents=True)
    conn = sqlite3.connect(proj / "project.db")
    conn.execute("CREATE TABLE configuration (jobid INTEGER PRIMARY KEY, cfg TEXT NOT NULL)")
    conn.execute("CREATE TABLE jobs (jobid INTEGER PRIMARY KEY, status TEXT, sample TEXT)")
    for jid in range(8):
        conn.execute("INSERT INTO configuration VALUES (?, ?)", (jid, _cfg(jid)))
        conn.execute("INSERT INTO jobs VALUES (?, 'pending', 'sbnd')", (jid,))
    conn.commit(); conn.close()

    big = b"x" * 2048
    for jid in (0, 2, 3, 7):
        _write(output_path_for(proj, jid, "nosyst"), big, binary=True)
        _write(output_path_for(proj, jid, "wsyst"), big, binary=True)
    _write(proj / "output" / "output_jobid0001.root", big, binary=True)              # legacy flat
    _write(proj / "output" / "output_systematics_jobid0001.root", big, binary=True)
    _write(output_path_for(proj, 1, "badfiles"), "/pnfs/fake/job1_file2.flat.caf.root\n")

    recs = {
        0: _record(0, "ok",            100, 0, 1000, "0,1"),
        1: _record(1, "match_partial", 100, 0, 1000, "0,1", n_inputs=2, n_dropped=1, input_events=20),
        3: _record(3, "sys_crash",     100, 1, 1001, "1,3"),
        4: _record(4, "match_unknown", 200, 0, 2000, "4,5"),
        5: _record(5, "ok",            100, 2, 1002, "5,7"),
        7: _record(7, "ok",            100, 2, 1002, "5,7", input_events=50, events_written=48),
        99: _record(99, "ok",          100, 3, 1003, "99"),
    }
    for jid, text in recs.items():
        path = (proj / "output" / "val" / f"validation_jobid{jid:04d}.txt") if jid == 1 \
            else output_path_for(proj, jid, "record")
        _write(path, text)
    return proj


class TestCategories:

    def test_every_job_id_in_exactly_one_category(self, project):
        r = reconcile_project(project)
        placed = [j for c in RECONCILE_CATEGORIES for j in r["categories"][c]]
        assert sorted(placed) == list(range(8))
        assert len(placed) == len(set(placed))
        assert r["n_jobs"] == 8

    def test_each_category(self, project):
        c = reconcile_project(project)["categories"]
        assert c["complete"] == [0, 1, 7]
        assert c["complete_unrecorded"] == [2]
        assert c["complete_stale_record"] == [3]
        assert c["failed"] == [4]
        assert c["output_lost"] == [5]
        assert c["no_trace"] == [6]

    def test_failed_split_by_status(self, project):
        assert reconcile_project(project)["failed_by_status"] == {"match_unknown": [4]}

    def test_record_for_unknown_job_id_is_reported(self, project):
        assert reconcile_project(project)["unknown_records"] == [99]


class TestInputFiles:

    def test_every_file_accounted_for(self, project):
        f = reconcile_project(project)["files"]
        assert f["expected"] == 24                      # 8 jobs x 3 files
        # complete: 0 (3) + 1 (2, one dropped) + 2 (3, unverified) + 3 (3) + 7 (3)
        assert f["processed"] == 14
        assert f["dropped"] == 1
        assert f["in_incomplete"] == 9                  # jobs 4, 5, 6
        assert f["unverified"] == 3                     # job 2 has no record
        assert f["balanced"] is True

    def test_dropped_inputs_are_named(self, project):
        assert reconcile_project(project)["dropped_files"] == ["/pnfs/fake/job1_file2.flat.caf.root"]

    def test_imbalance_is_detected(self, project):
        """A record claiming more inputs than the job owns must not balance."""
        _write(output_path_for(project, 0, "record"),
               _record(0, "ok", 100, 0, 1000, "0,1", n_inputs=5))
        assert reconcile_project(project)["files"]["balanced"] is False


class TestEventsAndClaims:

    def test_events_read_vs_written(self, project):
        e = reconcile_project(project)["events"]
        assert e["checked"] == 3                        # complete, recorded: 0, 1, 7
        assert e["mismatched"] == [7]
        assert e["unknown"] == []

    def test_unknown_event_count_is_separate_from_a_mismatch(self, project):
        _write(output_path_for(project, 0, "record"),
               _record(0, "ok", 100, 0, 1000, "0,1", input_events=-1, events_written=30))
        e = reconcile_project(project)["events"]
        assert e["unknown"] == [0]
        assert 0 not in e["mismatched"]

    def test_duplicate_within_a_launch_vs_retry_across_launches(self, project):
        """Job 1: two processes of cluster 100 -- duplicate work, a bug.
        Job 5: clusters 200 then 100 -- a retry, expected."""
        r = reconcile_project(project)
        assert list(r["duplicate_claims"]) == [1]
        assert list(r["retried"]) == [5]

    def test_pot_sums_complete_recorded_jobs(self, project):
        assert reconcile_project(project)["pot"] == pytest.approx(3e16)


class TestParsing:

    def test_first_occurrence_wins_and_lists_collect(self):
        rec = parse_validation_record(
            "JOB (4,2) VALIDATION\nSAMPLE=a\nSAMPLE=b\nERROR=x\nERROR=y\nWARN=w\n"
            "TREE_selected=33,32,1\nSTATUS=ok\n")
        assert (rec["jobid"], rec["process"]) == (4, 2)
        assert rec["SAMPLE"] == "a"
        assert rec["errors"] == ["x", "y"] and rec["warns"] == ["w"]
        assert rec["trees"]["selected"] == (33, 32, 1)

    def test_report_mentions_every_category(self, project):
        text = format_reconcile_report("proj", reconcile_project(project))
        for c in RECONCILE_CATEGORIES:
            assert c in text


class TestCli:

    def test_exit_code_flags_gaps(self, project, capsys):
        """The fixture has output lost, a stale record, duplicate work and an
        event mismatch, so the CLI must fail -- it is meant to gate a merge."""
        assert reconcile_cli.main([str(project)]) == 1
        assert "GAPS in proj" in capsys.readouterr().out

    def test_clean_project_exits_zero(self, tmp_path):
        proj = tmp_path / "clean"
        (proj / "output").mkdir(parents=True)
        conn = sqlite3.connect(proj / "project.db")
        conn.execute("CREATE TABLE configuration (jobid INTEGER PRIMARY KEY, cfg TEXT NOT NULL)")
        conn.execute("CREATE TABLE jobs (jobid INTEGER PRIMARY KEY, status TEXT, sample TEXT)")
        conn.execute("INSERT INTO configuration VALUES (0, ?)", (_cfg(0),))
        conn.execute("INSERT INTO jobs VALUES (0, 'completed', 'sbnd')")
        conn.commit(); conn.close()
        _write(output_path_for(proj, 0, "nosyst"), b"x" * 2048, binary=True)
        _write(output_path_for(proj, 0, "wsyst"), b"x" * 2048, binary=True)
        _write(output_path_for(proj, 0, "record"), _record(0, "ok", 1, 0, 1, "0"))
        assert reconcile_cli.main([str(proj)]) == 0
