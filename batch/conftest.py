"""
Shared pytest fixtures for medulla/batch tests.
"""
import textwrap
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Low-level write helpers
# ---------------------------------------------------------------------------

def _write(path: Path, content: str) -> Path:
    """Write dedented *content* to *path*, creating parent dirs as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(content))
    return path


def _write_toml(directory: Path, filename: str, content: str) -> Path:
    """Write dedented *content* into *directory*/*filename* and return the path."""
    return _write(directory / filename, content)


# ---------------------------------------------------------------------------
# catalog fixture — used by test_catalog.py
#
# Four entries: two experiments × {MC, data}.  Kept minimal so that
# test_catalog.py can assert exact counts without surprises.
# ---------------------------------------------------------------------------

@pytest.fixture
def catalog(tmp_path):
    """Four-entry sample catalog covering both experiments and both data types."""
    return _write_toml(tmp_path, "samples.toml", """\
        [[sample]]
        key = "sbnd_mc_nominal"
        name = "sbnd"
        path = "/pnfs/sbnd/mc/nominal/*.flat.root"
        ismc = true
        experiment = "sbnd"

        [[sample]]
        key = "sbnd_offbeam"
        name = "sbnd_offbeam"
        path = "/pnfs/sbnd/data/offbeam/*.flat.root"
        ismc = false
        experiment = "sbnd"

        [[sample]]
        key = "icarus_mc_nominal"
        name = "icarus"
        path = "/pnfs/icarus/mc/nominal/*.flat.root"
        ismc = true
        experiment = "icarus"

        [[sample]]
        key = "icarus_onbeam"
        name = "icarus_onbeam"
        path = "/pnfs/icarus/data/onbeam/*.flat.root"
        ismc = false
        experiment = "icarus"
    """)


# ---------------------------------------------------------------------------
# workspace fixture — used by test_campaign.py
#
# Builds a complete mock workspace under tmp_path that mirrors the real
# selection/toml/ directory layout.  Contains four analyses:
#
#   alpha_2026  — two experiments, two roles (primary + data_blind_safe)
#   beta_2026   — SBND only, one role (primary)
#   gamma_legacy — no meta.toml; should be silently skipped by discovery
#   delta_2026  — SBND only, inline [[sample]] blocks (no include_samples)
#
# The catalog has six entries to cover all enable-key lists used in the
# analysis TOMLs (two more than the four-entry catalog above).
# ---------------------------------------------------------------------------

def _make_campaign_catalog(root: Path) -> Path:
    """Create common/samples.toml under *root* with six entries."""
    return _write(root / "common" / "samples.toml", """\
        [[sample]]
        key = "sbnd_mc_nominal"
        name = "sbnd"
        path = "/pnfs/sbnd/mc/nominal/*.flat.root"
        ismc = true
        experiment = "sbnd"

        [[sample]]
        key = "sbnd_offbeam"
        name = "sbnd_offbeam"
        path = "/pnfs/sbnd/data/offbeam/*.flat.root"
        ismc = false
        experiment = "sbnd"

        [[sample]]
        key = "sbnd_bnblight"
        name = "sbnd_bnblight"
        path = "/pnfs/sbnd/data/devsample/*.flat.root"
        ismc = false
        experiment = "sbnd"

        [[sample]]
        key = "icarus_mc_nominal"
        name = "icarus"
        path = "/pnfs/icarus/mc/nominal/*.flat.root"
        ismc = true
        experiment = "icarus"

        [[sample]]
        key = "icarus_offbeam"
        name = "icarus_offbeam"
        path = "/pnfs/icarus/data/offbeam/*.flat.root"
        ismc = false
        experiment = "icarus"

        [[sample]]
        key = "icarus_onbeam"
        name = "icarus_onbeam"
        path = "/pnfs/icarus/data/onbeam/*.flat.root"
        ismc = false
        experiment = "icarus"
    """)


def _make_analysis_alpha(toml_root: Path):
    """Two-experiment analysis with primary and data_blind_safe roles."""
    d = toml_root / "alpha_2026"
    _write(d / "meta.toml", """\
        [meta]
        analysis = "alpha"
        description = "Alpha analysis"
        owners = ["alice", "bob"]
        experiments = ["sbnd", "icarus"]

        [defaults]
        batch_size = 50
        systematics_template = "batch/sys_template.toml"

        [[toml]]
        role = "primary"
        file = "selection.toml"
        experiments = ["sbnd", "icarus"]
        [toml.enable.sbnd]
        keys = ["sbnd_mc_nominal", "sbnd_offbeam"]
        [toml.enable.icarus]
        keys = ["icarus_mc_nominal", "icarus_offbeam"]

        [[toml]]
        role = "data_blind_safe"
        file = "data_blind_safe.toml"
        experiments = ["sbnd", "icarus"]
        [toml.enable.sbnd]
        keys = ["sbnd_bnblight"]
        [toml.enable.icarus]
        keys = ["icarus_onbeam"]
    """)
    _write(d / "selection.toml", """\
        [general]
        output = "alpha_2026"

        [[include_samples]]
        keys = ["sbnd_mc_nominal", "sbnd_offbeam", "icarus_mc_nominal", "icarus_offbeam"]

        [[tree]]
        name = "selected"
        sim_only = false
        mode = "reco"
        cut = []
        branch = []
    """)
    _write(d / "data_blind_safe.toml", """\
        [general]
        output = "alpha_2026"

        [[include_samples]]
        keys = ["sbnd_bnblight", "icarus_onbeam"]

        [[tree]]
        name = "selected"
        sim_only = false
        mode = "reco"
        cut = []
        branch = []
    """)


def _make_analysis_beta(toml_root: Path):
    """Single-experiment (SBND) analysis with one primary role."""
    d = toml_root / "beta_2026"
    _write(d / "meta.toml", """\
        [meta]
        analysis = "beta"
        description = "Beta analysis (SBND only)"
        owners = ["charlie"]
        experiments = ["sbnd"]

        [defaults]
        batch_size = 25

        [[toml]]
        role = "primary"
        file = "selection.toml"
        experiments = ["sbnd"]
        [toml.enable.sbnd]
        keys = ["sbnd_mc_nominal"]
    """)
    _write(d / "selection.toml", """\
        [general]
        output = "beta_2026"

        [[include_samples]]
        keys = ["sbnd_mc_nominal"]

        [[tree]]
        name = "selected"
        sim_only = false
        mode = "reco"
        cut = []
        branch = []
    """)


def _make_analysis_gamma_no_meta(toml_root: Path):
    """Legacy analysis directory without meta.toml; must be skipped by discovery."""
    d = toml_root / "gamma_legacy"
    _write(d / "selection.toml", """\
        [general]
        output = "gamma_legacy"

        [[sample]]
        name = "inline_only"
        path = "/pnfs/sbnd/mc/*.root"
        ismc = true

        [[tree]]
        name = "selected"
        sim_only = false
        mode = "reco"
        cut = []
        branch = []
    """)


def _make_analysis_delta_inline(toml_root: Path):
    """SBND analysis whose TOML uses inline [[sample]] blocks (no include_samples)."""
    d = toml_root / "delta_2026"
    _write(d / "meta.toml", """\
        [meta]
        analysis = "delta"
        description = "Delta analysis with inline samples"
        owners = ["eve"]
        experiments = ["sbnd"]

        [defaults]
        batch_size = 10

        [[toml]]
        role = "data_quality"
        file = "data_quality.toml"
        experiments = ["sbnd"]
    """)
    _write(d / "data_quality.toml", """\
        [general]
        output = "data_quality"

        [[sample]]
        name = "offbeamlight"
        path = "/pnfs/sbnd/data/offbeamlight/input*.flat.root"
        ismc = false
        disable = false

        [[tree]]
        name = "noCut_nreco"
        sim_only = false
        mode = "event"
        cut = []
        branch = []
    """)


@pytest.fixture
def workspace(tmp_path):
    """Complete mock workspace with four analyses and a six-entry catalog.

    Layout under tmp_path::

        selection/toml/
            common/samples.toml          ← six-entry catalog
            alpha_2026/{meta,selection,data_blind_safe}.toml
            beta_2026/{meta,selection}.toml
            gamma_legacy/selection.toml  ← no meta.toml; skipped by discovery
            delta_2026/{meta,data_quality}.toml

    Returns a dict with keys ``root``, ``toml_root``, and ``catalog``.
    """
    toml_root = tmp_path / "selection" / "toml"
    toml_root.mkdir(parents=True)

    catalog = _make_campaign_catalog(toml_root)
    _make_analysis_alpha(toml_root)
    _make_analysis_beta(toml_root)
    _make_analysis_gamma_no_meta(toml_root)
    _make_analysis_delta_inline(toml_root)

    return {
        "root": tmp_path,
        "toml_root": toml_root,
        "catalog": catalog,
    }
