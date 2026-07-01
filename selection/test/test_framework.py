"""
Framework validation tests.

Each test calls ``validate --validate --group <name>`` in the pre-built
working directory prepared by the ``framework_env`` session fixture
(defined in the root conftest.py).  A non-zero exit code from the binary
means at least one condition in that group failed.
"""
import subprocess

import pytest


@pytest.mark.framework
@pytest.mark.parametrize("group", [
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
    "data_reco_track_selectors",
])
def test_framework_group(framework_env, group):
    """Run validate --validate --group <group> and assert exit code 0."""
    result = subprocess.run(
        [str(framework_env["validate"]), "--validate", "--group", group],
        cwd=framework_env["workdir"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"validate --validate --group {group} reported failures:\n"
        f"{result.stdout}\n{result.stderr}"
    )
