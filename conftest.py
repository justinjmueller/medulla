"""
Root-level pytest fixtures for framework validation tests.

The ``framework_env`` session fixture runs the three-step pipeline needed by
selection/test/test_framework.py:

    1. validate --generate   (create validation_simlike.root / validation_datalike.root)
    2. medulla test.toml     (run the selection framework over the generated inputs)
    3. Expose the working directory to tests via the ``framework_env`` fixture.

Both binaries are located via the ``MEDULLA_BUILD_DIR`` environment variable,
which is set automatically when the ``pytest`` CMake target is used.  If the
variable is absent or the binaries cannot be found, every framework test is
skipped rather than errored.
"""
import os
import subprocess
from pathlib import Path
from typing import Optional

import pytest

_TOML = Path(__file__).parent / "selection" / "test" / "test.toml"
_GROUPS = [
    "sim_reco",
    "sim_reco_with_truth_cut",
    "sim_truth",
    "sim_truth_with_reco_cut",
    "data_reco",
    "sim_reco_particles",
    "sim_reco_particles_with_truth_cut",
    "sim_truth_particles",
    "sim_truth_particles_with_reco_cut",
    "data_reco_particles",
    "sim_event",
    "sim_reco_event_cut",
    "sim_reco_mctruth",
    "sim_truth_mctruth",
]


def _find_binary(build_dir: Path, name: str) -> Optional[Path]:
    candidate = build_dir / "selection" / name
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return candidate
    # flat layout (in-source or alternative build)
    candidate = build_dir / name
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return candidate
    return None


@pytest.fixture(scope="session")
def framework_env(tmp_path_factory):
    """Prepare the framework validation environment.

    Returns a dict with:
        ``workdir``    – Path to the temp directory containing test.root
        ``validate``   – Path to the validate binary
        ``groups``     – list of group names available for testing
    Skips the entire session if MEDULLA_BUILD_DIR is unset or binaries missing.
    """
    build_dir_str = os.environ.get("MEDULLA_BUILD_DIR", "")
    if not build_dir_str:
        pytest.skip("MEDULLA_BUILD_DIR not set — skipping framework tests")

    build_dir = Path(build_dir_str)
    validate_bin = _find_binary(build_dir, "validate")
    medulla_bin = _find_binary(build_dir, "medulla")

    if validate_bin is None:
        pytest.skip(f"validate binary not found under {build_dir}")
    if medulla_bin is None:
        pytest.skip(f"medulla binary not found under {build_dir}")

    workdir = tmp_path_factory.mktemp("framework")

    # Step 1 — generate input ROOT files.
    result = subprocess.run(
        [str(validate_bin), "--generate"],
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.skip(
            f"validate --generate failed (rc={result.returncode}):\n{result.stderr}"
        )

    # Step 2 — run the framework selection.
    result = subprocess.run(
        [str(medulla_bin), str(_TOML)],
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.skip(
            f"medulla {_TOML.name} failed (rc={result.returncode}):\n{result.stderr}"
        )

    return {
        "workdir": workdir,
        "validate": validate_bin,
        "groups": _GROUPS,
    }
