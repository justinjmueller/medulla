"""
Tests for medulla/batch/validate_pair.C, the worker-side check that gates
a job's copy-back on its (selection, systematics) output pair being good.

These build real ROOT files (via make_fixture.C) and run the macro against
them as a subprocess, exactly as submit.sh does, so what is tested is the
thing that actually runs on the grid -- including its exit code, which is
what submit.sh branches on.

Fixtures are written with ROOT rather than uproot deliberately: uproot
pulls in pyarrow, which in the SL7 analysis container is linked against a
newer glibc than the container provides and fails to import. Depending
only on ROOT keeps the tests runnable wherever the jobs themselves are.

The fixtures reproduce the failure modes seen in production rather than
hypothetical ones: a systematics file that opens cleanly but contains
nothing (systematics/src/main.cc skips missing trees and still returns 0),
one missing a single tree, one whose weight tables are short, and one
whose POT disagrees with the selection file.

Skipped when ROOT is unavailable, so the suite still runs outside the
analysis environment.
"""
import subprocess
import shutil
from pathlib import Path

import pytest

pytestmark = pytest.mark.skipif(
    shutil.which("root") is None,
    reason="ROOT is not on PATH; set up the analysis environment to run these.",
)

HERE = Path(__file__).resolve().parent
MACRO = HERE / "validate_pair.C"
FIXTURE_MACRO = HERE / "make_fixture.C"

SAMPLE = "sbnd"
TREE = "selectedNu"

# Exit codes from validate_pair.C, mirrored here so a change to the
# contract has to be made deliberately in both places.
RC_OK = 0
RC_MANIFEST = 1
RC_NOSYST = 2
RC_SYST = 3
RC_EXPOSURE = 4
RC_TREES = 5


def _make_fixture(path, sample=SAMPLE, tree=TREE, n_main=10, n_non=0,
                  n_table=None, pot=1.0e19, livetime=0.0, tables=(),
                  include_tree=True, include_dir=True):
    """Build one output file by invoking make_fixture.C."""
    n_table = n_main if n_table is None else n_table
    args = ', '.join([
        f'"{path}"', f'"{sample}"', f'"{tree}"',
        str(n_main), str(n_non), str(n_table),
        repr(pot), repr(livetime),
        f'"{",".join(tables)}"',
        'true' if include_tree else 'false',
        'true' if include_dir else 'false',
    ])
    proc = subprocess.run(
        ["root", "-l", "-b", "-q", f"{FIXTURE_MACRO}({args})"],
        capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, f"fixture build failed:\n{proc.stdout}\n{proc.stderr}"
    assert Path(path).exists(), f"fixture not created:\n{proc.stdout}\n{proc.stderr}"


def _write_nosyst(path, n_events=10, pot=1.0e19, livetime=0.0, tree=TREE):
    _make_fixture(path, tree=tree, n_main=n_events, pot=pot, livetime=livetime)


def _write_syst(path, n_sel=10, n_non=0, pot=1.0e19, livetime=0.0,
                tables=("multisim", "multisigma"), n_table=None,
                include_tree=True, include_dir=True):
    """
    Write a systematics output. The knobs correspond one-to-one with the
    ways this file has been observed to go wrong in production.
    """
    _make_fixture(path, n_main=n_sel, n_non=n_non, n_table=n_table,
                  pot=pot, livetime=livetime, tables=tables,
                  include_tree=include_tree, include_dir=include_dir)


def _write_manifest(path, action="add_weights",
                    table_types="multisim,multisigma", sample=SAMPLE):
    path.write_text(
        f"events/{sample}/{TREE}|{TREE}|{action}|{table_types}\n"
        # A second sample's entry, to confirm it is filtered out rather
        # than demanded of this job's output.
        f"events/{sample}_offbeam/{TREE}|{TREE}|copy|\n"
    )


def _run(tmp_path, nosyst="nosyst.root", syst="syst.root",
         manifest="manifest.txt", sample=SAMPLE):
    report = tmp_path / "report.txt"
    call = (f'{MACRO}("{tmp_path / nosyst}","{tmp_path / syst}",'
            f'"{sample}","{tmp_path / manifest}","{report}")')
    proc = subprocess.run(["root", "-l", "-b", "-q", call],
                          capture_output=True, text=True, timeout=300)
    text = report.read_text() if report.exists() else ""
    return proc.returncode, text


def _status(report):
    for line in report.splitlines():
        if line.startswith("STATUS="):
            return line.split("=", 1)[1]
    return None


# ---------------------------------------------------------------------------
# The good case
# ---------------------------------------------------------------------------

def test_matched_pair_passes(tmp_path):
    _write_nosyst(tmp_path / "nosyst.root", n_events=10)
    _write_syst(tmp_path / "syst.root", n_sel=10)
    _write_manifest(tmp_path / "manifest.txt")

    rc, report = _run(tmp_path)
    assert rc == RC_OK, report
    assert _status(report) == "ok"
    assert f"TREE_{TREE}=10,10,-1" in report


def test_nonmatched_events_are_accounted_not_rejected(tmp_path):
    """Matched + non-matched must sum to the input; a small unmatched
    fraction is normal and must not block the transfer."""
    _write_nosyst(tmp_path / "nosyst.root", n_events=1000)
    _write_syst(tmp_path / "syst.root", n_sel=999, n_non=1, n_table=999)
    _write_manifest(tmp_path / "manifest.txt")

    rc, report = _run(tmp_path)
    assert rc == RC_OK, report
    assert _status(report) == "ok"


def test_large_nonmatched_fraction_warns_but_passes(tmp_path):
    _write_nosyst(tmp_path / "nosyst.root", n_events=100)
    _write_syst(tmp_path / "syst.root", n_sel=50, n_non=50, n_table=50)
    _write_manifest(tmp_path / "manifest.txt")

    rc, report = _run(tmp_path)
    assert rc == RC_OK, report
    assert "WARN=match_partial" in report
    assert _status(report) == "match_partial"


def test_offbeam_sample_with_zero_pot_passes(tmp_path):
    """Offbeam/intime samples carry no POT and are bookkept by livetime;
    requiring POT > 0 would reject every one of them."""
    _write_nosyst(tmp_path / "nosyst.root", n_events=5, pot=0.0, livetime=123.0)
    _write_syst(tmp_path / "syst.root", n_sel=5, pot=0.0, livetime=123.0,
                tables=())
    _write_manifest(tmp_path / "manifest.txt", action="copy", table_types="")

    rc, report = _run(tmp_path)
    assert rc == RC_OK, report


# ---------------------------------------------------------------------------
# The failure this whole mechanism exists for
# ---------------------------------------------------------------------------

def test_empty_systematics_output_is_rejected(tmp_path):
    """The production failure: output.root fine, output_sys.root opens
    cleanly but has no sample directory because every tree was skipped."""
    _write_nosyst(tmp_path / "nosyst.root", n_events=10)
    _write_syst(tmp_path / "syst.root", include_dir=False)
    _write_manifest(tmp_path / "manifest.txt")

    rc, report = _run(tmp_path)
    assert rc == RC_SYST, report
    assert _status(report) == "match_empty"


def test_missing_systematics_file_is_rejected(tmp_path):
    _write_nosyst(tmp_path / "nosyst.root", n_events=10)
    _write_manifest(tmp_path / "manifest.txt")

    rc, report = _run(tmp_path)
    assert rc == RC_SYST, report
    # The report must survive the early exit -- it is the only per-job
    # signal, and gSystem->Exit() does not unwind local destructors.
    assert "ERROR=syst_unopenable" in report


def test_missing_tree_in_systematics_is_rejected(tmp_path):
    _write_nosyst(tmp_path / "nosyst.root", n_events=10)
    _write_syst(tmp_path / "syst.root", n_sel=10, include_tree=False)
    _write_manifest(tmp_path / "manifest.txt")

    rc, report = _run(tmp_path)
    assert rc == RC_TREES, report
    assert f"ERROR=syst_missing_tree:{TREE}" in report


def test_missing_weight_table_is_rejected(tmp_path):
    """A tree copied across without its universe weights -- the file looks
    populated, but the systematics are simply absent."""
    _write_nosyst(tmp_path / "nosyst.root", n_events=10)
    _write_syst(tmp_path / "syst.root", n_sel=10, tables=("multisim",))
    _write_manifest(tmp_path / "manifest.txt")

    rc, report = _run(tmp_path)
    assert rc == RC_TREES, report
    assert f"ERROR=missing_weight_table:{TREE}_multisigmaTree" in report


def test_short_weight_table_is_rejected(tmp_path):
    _write_nosyst(tmp_path / "nosyst.root", n_events=10)
    _write_syst(tmp_path / "syst.root", n_sel=10, n_table=7)
    _write_manifest(tmp_path / "manifest.txt")

    rc, report = _run(tmp_path)
    assert rc == RC_TREES, report
    assert "ERROR=weight_table_entries" in report


def test_event_accounting_mismatch_is_rejected(tmp_path):
    """Matched + non-matched must equal the input; silent event loss here
    is invisible to every exit-code-based check."""
    _write_nosyst(tmp_path / "nosyst.root", n_events=10)
    _write_syst(tmp_path / "syst.root", n_sel=6, n_non=1, n_table=6)
    _write_manifest(tmp_path / "manifest.txt")

    rc, report = _run(tmp_path)
    assert rc == RC_TREES, report
    assert f"ERROR=match_accounting:{TREE}" in report


def test_all_events_unmatched_is_rejected(tmp_path):
    _write_nosyst(tmp_path / "nosyst.root", n_events=10)
    _write_syst(tmp_path / "syst.root", n_sel=0, n_non=10, n_table=0)
    _write_manifest(tmp_path / "manifest.txt")

    rc, report = _run(tmp_path)
    assert rc == RC_TREES, report
    assert f"ERROR=match_empty:{TREE}" in report


def test_copy_action_entry_mismatch_is_rejected(tmp_path):
    _write_nosyst(tmp_path / "nosyst.root", n_events=10)
    _write_syst(tmp_path / "syst.root", n_sel=9, tables=())
    _write_manifest(tmp_path / "manifest.txt", action="copy", table_types="")

    rc, report = _run(tmp_path)
    assert rc == RC_TREES, report
    assert "ERROR=copy_entries" in report


# ---------------------------------------------------------------------------
# Exposure
# ---------------------------------------------------------------------------

def test_pot_mismatch_is_rejected(tmp_path):
    """The check collect_good_syst_files.py does offline, done here before
    anything is transferred."""
    _write_nosyst(tmp_path / "nosyst.root", n_events=10, pot=1.0e19)
    _write_syst(tmp_path / "syst.root", n_sel=10, pot=2.0e19)
    _write_manifest(tmp_path / "manifest.txt")

    rc, report = _run(tmp_path)
    assert rc == RC_EXPOSURE, report
    assert "ERROR=exposure_mismatch" in report


def test_zero_exposure_is_rejected(tmp_path):
    _write_nosyst(tmp_path / "nosyst.root", n_events=10, pot=0.0, livetime=0.0)
    _write_syst(tmp_path / "syst.root", n_sel=10, pot=0.0, livetime=0.0)
    _write_manifest(tmp_path / "manifest.txt")

    rc, report = _run(tmp_path)
    assert rc == RC_NOSYST, report
    assert "ERROR=nosyst_zero_exposure" in report


# ---------------------------------------------------------------------------
# Selection output and manifest handling
# ---------------------------------------------------------------------------

def test_missing_selection_file_is_rejected(tmp_path):
    _write_syst(tmp_path / "syst.root", n_sel=10)
    _write_manifest(tmp_path / "manifest.txt")

    rc, report = _run(tmp_path)
    assert rc == RC_NOSYST, report
    assert "ERROR=nosyst_unopenable" in report


def test_wrong_sample_directory_is_rejected(tmp_path):
    """A job whose output is keyed by a different sample than the one it
    was assigned is a real, silent mix-up -- the offline checker cannot
    see it because it infers the sample from the file itself."""
    _write_nosyst(tmp_path / "nosyst.root", n_events=10)
    _write_syst(tmp_path / "syst.root", n_sel=10)
    _write_manifest(tmp_path / "manifest.txt", sample="sbnd_dirt")

    rc, report = _run(tmp_path, sample="sbnd_dirt")
    assert rc == RC_NOSYST, report
    assert "ERROR=nosyst_missing_sample_dir" in report


def test_sample_absent_from_manifest_is_rejected(tmp_path):
    _write_nosyst(tmp_path / "nosyst.root", n_events=10)
    _write_syst(tmp_path / "syst.root", n_sel=10)
    _write_manifest(tmp_path / "manifest.txt", sample="icarus_run2")

    rc, report = _run(tmp_path)
    assert rc == RC_MANIFEST, report
    assert "ERROR=no_manifest_entries_for_sample" in report


def test_manifest_prefix_match_is_exact(tmp_path):
    """'sbnd' must not pick up 'sbnd_offbeam' entries: the prefix includes
    the trailing slash precisely so that one sample name cannot swallow
    another that extends it."""
    _write_nosyst(tmp_path / "nosyst.root", n_events=10)
    _write_syst(tmp_path / "syst.root", n_sel=10)
    _write_manifest(tmp_path / "manifest.txt")

    rc, report = _run(tmp_path)
    assert rc == RC_OK, report
    assert "N_EXPECTED_TREES=1" in report


def test_unreadable_manifest_is_rejected(tmp_path):
    _write_nosyst(tmp_path / "nosyst.root", n_events=10)
    _write_syst(tmp_path / "syst.root", n_sel=10)

    rc, report = _run(tmp_path)
    assert rc == RC_MANIFEST, report
    assert "ERROR=manifest_unreadable" in report
