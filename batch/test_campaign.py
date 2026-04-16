"""
Integration tests for the campaign management system.
======================================================

Tests the full campaign orchestration layer in campaign.py and its
integration with the sample catalog resolution in catalog.py.

These tests describe the *intended* interface for three functions that
are not yet in campaign.py.  All test classes are decorated with
``@skip_campaign`` and will be automatically skipped until the
implementation is added.  The tests serve as a specification: once
``discover_analyses``, ``expand_campaign``, and ``create_campaign``
exist, the suite should pass without modification.

Target API
----------
campaign.discover_analyses(toml_root) -> list[AnalysisMeta]
    Walk selection/toml/*/ looking for meta.toml files.
    Return one AnalysisMeta per directory that has one.

campaign.expand_campaign(
    analyses, catalog_path,
    experiment=None, roles=None, analysis_filter=None,
) -> list[ProjectUnit]
    Expand analyses into one ProjectUnit per
    (analysis, role, experiment) combination, filtered by the
    caller's constraints.  Each unit carries: analysis name, role,
    experiment, TOML path, and resolved enable_keys.

campaign.create_campaign(
    campaign_dir, project_units, catalog_path,
    name, tag,
    batch_size_override=None, campaign_cfg=None,
) -> None
    Create campaign_dir with campaign.db (metadata + projects table).
    For each ProjectUnit, call create_new_project() in a subdirectory.

The meta.toml contract
----------------------
::

    [meta]
    analysis = "<name>"
    description = "..."
    owners = ["person_a", "person_b"]
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

Run with:
    pytest               # from medulla/ (uses pytest.ini testpaths)
    pytest batch/        # from medulla/ explicitly
    make pytest          # from the build directory
"""

import sqlite3
import textwrap
from pathlib import Path
from unittest import mock

import toml
import pytest


# ---------------------------------------------------------------------------
# Import the modules under test.  Each function is imported individually so
# that test classes can be skipped at the granularity of the staged
# implementation:
#   Stage 1 — discover_analyses  → TestDiscovery
#   Stage 2 — expand_campaign    → TestExpansionFilters
#   Stage 3 — create_campaign    → everything else
# ---------------------------------------------------------------------------
try:
    from campaign import discover_analyses
    _DISCOVER_AVAILABLE = True
except ImportError:
    _DISCOVER_AVAILABLE = False

try:
    from campaign import expand_campaign
    _EXPAND_AVAILABLE = True
except ImportError:
    _EXPAND_AVAILABLE = False

try:
    from campaign import create_campaign
    _CREATE_AVAILABLE = True
except ImportError:
    _CREATE_AVAILABLE = False

try:
    from auth import authenticate
    _AUTH_AVAILABLE = True
except ImportError:
    _AUTH_AVAILABLE = False

skip_discover = pytest.mark.skipif(
    not _DISCOVER_AVAILABLE,
    reason="discover_analyses not yet implemented in campaign.py",
)
skip_expand = pytest.mark.skipif(
    not (_DISCOVER_AVAILABLE and _EXPAND_AVAILABLE),
    reason="discover_analyses/expand_campaign not yet implemented in campaign.py",
)
skip_campaign = pytest.mark.skipif(
    not (_DISCOVER_AVAILABLE and _EXPAND_AVAILABLE and _CREATE_AVAILABLE),
    reason="discover_analyses/expand_campaign/create_campaign not yet implemented in campaign.py",
)
skip_auth = pytest.mark.skipif(
    not _AUTH_AVAILABLE,
    reason="auth module not yet implemented",
)


# ===================================================================
# T3.1 — meta.toml discovery
# ===================================================================

@skip_discover
class TestDiscovery:
    """discover_analyses() finds all directories with meta.toml."""

    def test_discovers_all_analyses(self, workspace):
        analyses = discover_analyses(workspace["toml_root"])
        names = {a.analysis for a in analyses}
        # alpha, beta, and delta have meta.toml; gamma does not
        assert names == {"alpha", "beta", "delta"}

    def test_skips_directories_without_meta(self, workspace):
        analyses = discover_analyses(workspace["toml_root"])
        names = {a.analysis for a in analyses}
        assert "gamma_legacy" not in names

    def test_analysis_has_expected_fields(self, workspace):
        analyses = discover_analyses(workspace["toml_root"])
        alpha = [a for a in analyses if a.analysis == "alpha"][0]

        assert alpha.description == "Alpha analysis"
        assert set(alpha.owners) == {"alice", "bob"}
        assert set(alpha.experiments) == {"sbnd", "icarus"}

    def test_analysis_toml_roles(self, workspace):
        analyses = discover_analyses(workspace["toml_root"])
        alpha = [a for a in analyses if a.analysis == "alpha"][0]

        roles = {t.role for t in alpha.tomls}
        assert roles == {"primary", "data_blind_safe"}

    def test_analysis_toml_experiments(self, workspace):
        analyses = discover_analyses(workspace["toml_root"])
        alpha = [a for a in analyses if a.analysis == "alpha"][0]
        primary = [t for t in alpha.tomls if t.role == "primary"][0]

        assert set(primary.experiments) == {"sbnd", "icarus"}

    def test_analysis_toml_enable_keys(self, workspace):
        analyses = discover_analyses(workspace["toml_root"])
        alpha = [a for a in analyses if a.analysis == "alpha"][0]
        primary = [t for t in alpha.tomls if t.role == "primary"][0]

        assert set(primary.enable["sbnd"]) == {"sbnd_mc_nominal", "sbnd_offbeam"}
        assert set(primary.enable["icarus"]) == {"icarus_mc_nominal", "icarus_offbeam"}

    def test_defaults_batch_size(self, workspace):
        analyses = discover_analyses(workspace["toml_root"])
        alpha = [a for a in analyses if a.analysis == "alpha"][0]
        assert alpha.defaults["batch_size"] == 50

        beta = [a for a in analyses if a.analysis == "beta"][0]
        assert beta.defaults["batch_size"] == 25

    def test_empty_directory_returns_empty(self, tmp_path):
        empty = tmp_path / "empty_toml_root"
        empty.mkdir()
        analyses = discover_analyses(empty)
        assert analyses == []


# ===================================================================
# T3.2 — CLI filters in expand_campaign()
# ===================================================================

@skip_expand
class TestExpansionFilters:
    """expand_campaign() correctly filters by experiment, role, and analysis."""

    def _analyses(self, workspace):
        return discover_analyses(workspace["toml_root"])

    # -- Experiment filter ---

    def test_filter_experiment_sbnd(self, workspace):
        units = expand_campaign(
            self._analyses(workspace), workspace["catalog"], experiment="sbnd"
        )
        assert {u.experiment for u in units} == {"sbnd"}

    def test_filter_experiment_icarus(self, workspace):
        units = expand_campaign(
            self._analyses(workspace), workspace["catalog"], experiment="icarus"
        )
        assert {u.experiment for u in units} == {"icarus"}
        # beta is sbnd-only, delta is sbnd-only: neither should appear
        analysis_names = {u.analysis for u in units}
        assert "beta" not in analysis_names
        assert "delta" not in analysis_names

    # -- Role filter ---

    def test_filter_role_primary(self, workspace):
        units = expand_campaign(
            self._analyses(workspace), workspace["catalog"], roles=["primary"]
        )
        assert {u.role for u in units} == {"primary"}

    def test_filter_role_data_blind_safe(self, workspace):
        units = expand_campaign(
            self._analyses(workspace), workspace["catalog"], roles=["data_blind_safe"]
        )
        assert {u.role for u in units} == {"data_blind_safe"}
        assert {u.analysis for u in units} == {"alpha"}

    def test_filter_role_data_quality(self, workspace):
        units = expand_campaign(
            self._analyses(workspace), workspace["catalog"], roles=["data_quality"]
        )
        assert {u.role for u in units} == {"data_quality"}
        assert {u.analysis for u in units} == {"delta"}

    # -- Analysis filter ---

    def test_filter_analysis_alpha(self, workspace):
        units = expand_campaign(
            self._analyses(workspace), workspace["catalog"], analysis_filter=["alpha"]
        )
        assert {u.analysis for u in units} == {"alpha"}

    def test_filter_analysis_beta(self, workspace):
        units = expand_campaign(
            self._analyses(workspace), workspace["catalog"], analysis_filter=["beta"]
        )
        assert {u.analysis for u in units} == {"beta"}

    # -- Combination filters ---

    def test_filter_experiment_and_role(self, workspace):
        units = expand_campaign(
            self._analyses(workspace), workspace["catalog"],
            experiment="icarus", roles=["primary"],
        )
        assert len(units) >= 1
        for u in units:
            assert u.experiment == "icarus"
            assert u.role == "primary"

    def test_filter_analysis_and_experiment(self, workspace):
        units = expand_campaign(
            self._analyses(workspace), workspace["catalog"],
            experiment="sbnd", analysis_filter=["alpha"],
        )
        for u in units:
            assert u.analysis == "alpha"
            assert u.experiment == "sbnd"

    def test_triple_filter(self, workspace):
        units = expand_campaign(
            self._analyses(workspace), workspace["catalog"],
            experiment="sbnd", roles=["primary"], analysis_filter=["beta"],
        )
        assert len(units) == 1
        u = units[0]
        assert u.analysis == "beta"
        assert u.role == "primary"
        assert u.experiment == "sbnd"

    # -- No filter: full expansion ---

    def test_no_filter_full_expansion(self, workspace):
        """Full combinatorial expansion:
        alpha: primary×{sbnd,icarus} + data_blind_safe×{sbnd,icarus} = 4
        beta:  primary×{sbnd} = 1
        delta: data_quality×{sbnd} = 1
        Total = 6
        """
        units = expand_campaign(self._analyses(workspace), workspace["catalog"])
        assert len(units) == 6

    # -- Enable keys are correct per experiment ---

    def test_enable_keys_per_experiment(self, workspace):
        units = expand_campaign(
            self._analyses(workspace), workspace["catalog"],
            analysis_filter=["alpha"], roles=["primary"],
        )
        sbnd_units   = [u for u in units if u.experiment == "sbnd"]
        icarus_units = [u for u in units if u.experiment == "icarus"]

        assert len(sbnd_units) == 1
        assert set(sbnd_units[0].enable_keys) == {"sbnd_mc_nominal", "sbnd_offbeam"}

        assert len(icarus_units) == 1
        assert set(icarus_units[0].enable_keys) == {"icarus_mc_nominal", "icarus_offbeam"}

    # -- Edge: filter matches nothing ---

    def test_filter_nonexistent_analysis(self, workspace):
        units = expand_campaign(
            self._analyses(workspace), workspace["catalog"],
            analysis_filter=["nonexistent"],
        )
        assert units == []

    def test_filter_nonexistent_role(self, workspace):
        units = expand_campaign(
            self._analyses(workspace), workspace["catalog"],
            roles=["nonexistent_role"],
        )
        assert units == []


# ===================================================================
# T3.3 — Campaign manifest overrides
# ===================================================================

@skip_campaign
class TestCampaignOverrides:
    """A campaign config dict can override per-project defaults."""

    def test_batch_size_override_global(self, workspace):
        """A global batch_size_override replaces the per-analysis default."""
        analyses = discover_analyses(workspace["toml_root"])
        units = expand_campaign(analyses, workspace["catalog"])
        campaign_dir = workspace["root"] / "campaign_override"

        fake_sample = {"name": "fake", "path": ["/fake/file.root"], "ismc": True, "disable": False}
        with mock.patch("utilities.get_samples", return_value=[fake_sample]):
            create_campaign(
                campaign_dir=campaign_dir,
                project_units=units,
                catalog_path=workspace["catalog"],
                name="test_campaign",
                tag="v0.1.0",
                batch_size_override=999,
            )

        conn = sqlite3.connect(str(campaign_dir / "campaign.db"))
        curs = conn.cursor()
        curs.execute("SELECT batch_size FROM projects")
        rows = curs.fetchall()
        conn.close()

        for (bs,) in rows:
            assert bs == 999

    def test_per_project_override_in_campaign_cfg(self, workspace):
        """campaign_cfg can override settings for a specific project."""
        analyses = discover_analyses(workspace["toml_root"])
        units = expand_campaign(
            analyses, workspace["catalog"],
            analysis_filter=["alpha"], roles=["primary"],
        )
        campaign_dir = workspace["root"] / "campaign_per_proj"
        campaign_cfg = {
            "overrides": [{
                "analysis": "alpha", "role": "primary",
                "experiment": "sbnd", "batch_size": 100,
            }]
        }

        fake_sample = {"name": "fake", "path": ["/fake/file.root"], "ismc": True, "disable": False}
        with mock.patch("utilities.get_samples", return_value=[fake_sample]):
            create_campaign(
                campaign_dir=campaign_dir,
                project_units=units,
                catalog_path=workspace["catalog"],
                name="test_campaign",
                tag="v0.1.0",
                campaign_cfg=campaign_cfg,
            )

        conn = sqlite3.connect(str(campaign_dir / "campaign.db"))
        curs = conn.cursor()
        curs.execute(
            "SELECT batch_size FROM projects "
            "WHERE analysis = 'alpha' AND role = 'primary' AND experiment = 'sbnd'"
        )
        rows = curs.fetchall()
        conn.close()

        assert len(rows) == 1
        assert rows[0][0] == 100


# ===================================================================
# T3.4 — create_campaign() with resolved samples
# ===================================================================

@skip_campaign
class TestCreateCampaignResolved:
    """create_campaign() creates the expected directory structure and
    project databases for include_samples-based TOMLs."""

    def test_campaign_directory_structure(self, workspace):
        analyses = discover_analyses(workspace["toml_root"])
        units = expand_campaign(
            analyses, workspace["catalog"], analysis_filter=["beta"]
        )
        campaign_dir = workspace["root"] / "campaign_struct"

        fake_sample = {"name": "fake", "path": ["/fake/file.root"], "ismc": True, "disable": False}
        with mock.patch("utilities.get_samples", return_value=[fake_sample]):
            create_campaign(
                campaign_dir=campaign_dir,
                project_units=units,
                catalog_path=workspace["catalog"],
                name="beta_test",
                tag="v1.0.0",
            )

        assert (campaign_dir / "campaign.db").exists()
        project_dirs = [
            p for p in campaign_dir.iterdir()
            if p.is_dir() and (p / "project.db").exists()
        ]
        assert len(project_dirs) == 1  # beta has 1 unit

    def test_project_db_has_resolved_samples(self, workspace):
        """The configuration stored in project.db must have resolved
        [[sample]] blocks, not [[include_samples]]."""
        analyses = discover_analyses(workspace["toml_root"])
        units = expand_campaign(
            analyses, workspace["catalog"], analysis_filter=["beta"]
        )
        campaign_dir = workspace["root"] / "campaign_resolved"

        fake_sample = {
            "name": "sbnd",
            "path": ["/pnfs/sbnd/mc/nominal/file.flat.root"],
            "ismc": True,
            "disable": False,
        }
        with mock.patch("utilities.get_samples", return_value=[fake_sample]):
            create_campaign(
                campaign_dir=campaign_dir,
                project_units=units,
                catalog_path=workspace["catalog"],
                name="beta_resolved",
                tag="v1.0.0",
            )

        project_dirs = [
            p for p in campaign_dir.iterdir()
            if p.is_dir() and (p / "project.db").exists()
        ]
        assert len(project_dirs) == 1

        conn = sqlite3.connect(str(project_dirs[0] / "project.db"))
        curs = conn.cursor()
        curs.execute("SELECT cfg FROM configuration")
        rows = curs.fetchall()
        conn.close()

        for (cfg_text,) in rows:
            cfg = toml.loads(cfg_text)
            assert "include_samples" not in cfg, (
                "include_samples should be resolved before storage"
            )
            assert "sample" in cfg
            assert len(cfg["sample"]) >= 1

    def test_project_job_count(self, workspace):
        """Number of jobs should match the number of batched samples."""
        analyses = discover_analyses(workspace["toml_root"])
        units = expand_campaign(
            analyses, workspace["catalog"], analysis_filter=["beta"]
        )
        campaign_dir = workspace["root"] / "campaign_jobcount"

        fake_samples = [
            {"name": "sbnd", "path": [f"/fake/{i}.root"], "ismc": True, "disable": False}
            for i in range(3)
        ]
        with mock.patch("utilities.get_samples", return_value=fake_samples):
            create_campaign(
                campaign_dir=campaign_dir,
                project_units=units,
                catalog_path=workspace["catalog"],
                name="beta_count",
                tag="v1.0.0",
            )

        project_dirs = [
            p for p in campaign_dir.iterdir()
            if p.is_dir() and (p / "project.db").exists()
        ]
        conn = sqlite3.connect(str(project_dirs[0] / "project.db"))
        curs = conn.cursor()
        curs.execute("SELECT COUNT(*) FROM jobs")
        (count,) = curs.fetchone()
        conn.close()

        assert count == 3


# ===================================================================
# T3.5 — create_campaign() backward compatibility (inline samples)
# ===================================================================

@skip_campaign
class TestCreateCampaignInline:
    """Analyses whose TOMLs use inline [[sample]] blocks (no include_samples)
    work seamlessly through the campaign tool."""

    def test_inline_samples_project_created(self, workspace):
        analyses = discover_analyses(workspace["toml_root"])
        units = expand_campaign(
            analyses, workspace["catalog"], analysis_filter=["delta"]
        )
        campaign_dir = workspace["root"] / "campaign_inline"

        fake_sample = {
            "name": "offbeamlight",
            "path": ["/pnfs/sbnd/data/offbeamlight/input0.flat.root"],
            "ismc": False,
            "disable": False,
        }
        with mock.patch("utilities.get_samples", return_value=[fake_sample]):
            create_campaign(
                campaign_dir=campaign_dir,
                project_units=units,
                catalog_path=workspace["catalog"],
                name="delta_inline",
                tag="v1.0.0",
            )

        assert (campaign_dir / "campaign.db").exists()
        project_dirs = [
            p for p in campaign_dir.iterdir()
            if p.is_dir() and (p / "project.db").exists()
        ]
        assert len(project_dirs) == 1

    def test_inline_cfg_preserved(self, workspace):
        analyses = discover_analyses(workspace["toml_root"])
        units = expand_campaign(
            analyses, workspace["catalog"], analysis_filter=["delta"]
        )
        campaign_dir = workspace["root"] / "campaign_inline_cfg"

        fake_sample = {
            "name": "offbeamlight",
            "path": ["/pnfs/sbnd/data/offbeamlight/input0.flat.root"],
            "ismc": False,
            "disable": False,
        }
        with mock.patch("utilities.get_samples", return_value=[fake_sample]):
            create_campaign(
                campaign_dir=campaign_dir,
                project_units=units,
                catalog_path=workspace["catalog"],
                name="delta_inline_cfg",
                tag="v1.0.0",
            )

        project_dirs = [
            p for p in campaign_dir.iterdir()
            if p.is_dir() and (p / "project.db").exists()
        ]
        conn = sqlite3.connect(str(project_dirs[0] / "project.db"))
        curs = conn.cursor()
        curs.execute("SELECT cfg FROM configuration LIMIT 1")
        (cfg_text,) = curs.fetchone()
        conn.close()

        cfg = toml.loads(cfg_text)
        assert cfg["sample"][0]["name"] == "offbeamlight"


# ===================================================================
# T3.6 — Campaign database integrity
# ===================================================================

@skip_campaign
class TestCampaignDatabase:
    """The campaign.db contains correct metadata and project rows."""

    def _create_full_campaign(self, workspace, name="full"):
        analyses = discover_analyses(workspace["toml_root"])
        units = expand_campaign(analyses, workspace["catalog"])
        campaign_dir = workspace["root"] / f"campaign_{name}"

        fake_sample = {"name": "fake", "path": ["/fake/file.root"], "ismc": True, "disable": False}
        with mock.patch("utilities.get_samples", return_value=[fake_sample]):
            create_campaign(
                campaign_dir=campaign_dir,
                project_units=units,
                catalog_path=workspace["catalog"],
                name="full_campaign",
                tag="v2.0.0",
            )
        return campaign_dir

    def test_campaign_meta_table(self, workspace):
        campaign_dir = self._create_full_campaign(workspace)
        conn = sqlite3.connect(str(campaign_dir / "campaign.db"))
        curs = conn.cursor()
        curs.execute("SELECT name, tag FROM campaign_meta")
        row = curs.fetchone()
        conn.close()

        assert row is not None
        assert row[0] == "full_campaign"
        assert row[1] == "v2.0.0"

    def test_projects_table_row_count(self, workspace):
        """6 project units expected (see TestExpansionFilters.test_no_filter_full_expansion)."""
        campaign_dir = self._create_full_campaign(workspace, name="rowcount")
        conn = sqlite3.connect(str(campaign_dir / "campaign.db"))
        curs = conn.cursor()
        curs.execute("SELECT COUNT(*) FROM projects")
        (count,) = curs.fetchone()
        conn.close()

        assert count == 6

    def test_projects_table_columns(self, workspace):
        campaign_dir = self._create_full_campaign(workspace, name="cols")
        conn = sqlite3.connect(str(campaign_dir / "campaign.db"))
        curs = conn.cursor()
        curs.execute("SELECT analysis, role, experiment, status FROM projects")
        rows = curs.fetchall()
        conn.close()

        for analysis, role, experiment, status in rows:
            assert analysis in {"alpha", "beta", "delta"}
            assert role in {"primary", "data_blind_safe", "data_quality"}
            assert experiment in {"sbnd", "icarus"}
            assert status == "created"

    def test_projects_unique_combinations(self, workspace):
        """Each (analysis, role, experiment) tuple must be unique."""
        campaign_dir = self._create_full_campaign(workspace, name="unique")
        conn = sqlite3.connect(str(campaign_dir / "campaign.db"))
        curs = conn.cursor()
        curs.execute("SELECT analysis, role, experiment FROM projects")
        rows = curs.fetchall()
        conn.close()

        assert len(rows) == len(set(rows)), "Duplicate project entries found"

    def test_campaign_already_exists_error(self, workspace):
        """Creating a campaign in a directory that already has campaign.db
        must raise."""
        campaign_dir = self._create_full_campaign(workspace, name="exists")
        analyses = discover_analyses(workspace["toml_root"])
        units = expand_campaign(analyses, workspace["catalog"])

        fake_sample = {"name": "fake", "path": ["/fake/file.root"], "ismc": True, "disable": False}
        with mock.patch("utilities.get_samples", return_value=[fake_sample]):
            with pytest.raises((FileExistsError, RuntimeError)):
                create_campaign(
                    campaign_dir=campaign_dir,
                    project_units=units,
                    catalog_path=workspace["catalog"],
                    name="duplicate",
                    tag="v3.0.0",
                )


# ===================================================================
# T3.7 — auth.py mock authentication
# ===================================================================

@skip_auth
class TestAuthentication:
    """Mock-based tests for the authenticate() helper in auth.py."""

    def test_htgettoken_not_found(self):
        """If htgettoken is not on PATH, authenticate() returns False."""
        with mock.patch("shutil.which", return_value=None):
            with mock.patch("builtins.input", return_value="n"):
                result = authenticate("sbnd")
        assert result is False

    def test_htgettoken_success(self):
        """If htgettoken exits 0, authenticate() returns True."""
        with mock.patch("shutil.which", return_value="/usr/bin/htgettoken"):
            with mock.patch("subprocess.run", return_value=mock.Mock(returncode=0)):
                result = authenticate("sbnd")
        assert result is True

    def test_htgettoken_failure(self):
        """If htgettoken exits non-zero, authenticate() returns False."""
        with mock.patch("shutil.which", return_value="/usr/bin/htgettoken"):
            with mock.patch("subprocess.run", return_value=mock.Mock(returncode=1)):
                result = authenticate("sbnd")
        assert result is False

    def test_different_experiments_are_independent(self):
        """authenticate() is called per experiment; each call includes the
        experiment name in the htgettoken arguments."""
        calls = []

        def fake_run(cmd, **kwargs):
            calls.append(cmd)
            return mock.Mock(returncode=0)

        with mock.patch("shutil.which", return_value="/usr/bin/htgettoken"):
            with mock.patch("subprocess.run", side_effect=fake_run):
                authenticate("sbnd")
                authenticate("icarus")

        assert len(calls) == 2
        assert "sbnd"   in " ".join(str(c) for c in calls[0])
        assert "icarus" in " ".join(str(c) for c in calls[1])


# ===================================================================
# T3.8 — Status synchronisation
# ===================================================================

try:
    from campaign import _sync_project_status, cmd_sync
    _SYNC_AVAILABLE = _CREATE_AVAILABLE
except ImportError:
    _SYNC_AVAILABLE = False

skip_sync = pytest.mark.skipif(
    not _SYNC_AVAILABLE,
    reason="sync helpers not yet implemented in campaign.py",
)


@skip_sync
class TestStatusSync:
    """n_jobs is populated at creation; sync updates completion counts and
    campaign status from job output files."""

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _make_campaign(self, workspace, name="sync", n_fake_samples=1):
        """Create a campaign for the beta analysis with *n_fake_samples* jobs."""
        analyses = discover_analyses(workspace["toml_root"])
        units = expand_campaign(
            analyses, workspace["catalog"], analysis_filter=["beta"]
        )
        campaign_dir = workspace["root"] / f"campaign_{name}"
        fake_samples = [
            {"name": f"sbnd_{i}", "path": [f"/fake/{i}.root"],
             "ismc": True, "disable": False}
            for i in range(n_fake_samples)
        ]
        with mock.patch("utilities.get_samples", return_value=fake_samples):
            create_campaign(
                campaign_dir=campaign_dir,
                project_units=units,
                catalog_path=workspace["catalog"],
                name=f"campaign_{name}",
                tag="v1.0",
            )
        return campaign_dir

    def _project_dir(self, campaign_dir):
        """Return the single project directory inside *campaign_dir*."""
        dirs = [p for p in campaign_dir.iterdir()
                if p.is_dir() and (p / "project.db").exists()]
        assert len(dirs) == 1
        return dirs[0]

    def _write_output(self, proj_dir, jobid, size=2048):
        """Create a fake output file for *jobid* with *size* bytes."""
        out = proj_dir / "output" / f"output_jobid{jobid:04d}.root"
        out.write_bytes(b"x" * size)

    # ------------------------------------------------------------------
    # n_jobs at creation
    # ------------------------------------------------------------------

    def test_n_jobs_populated_at_creation(self, workspace):
        """n_jobs in campaign.db must equal the number of rows in project.db
        jobs table immediately after creation."""
        campaign_dir = self._make_campaign(workspace, name="njobs", n_fake_samples=3)

        conn = sqlite3.connect(str(campaign_dir / "campaign.db"))
        curs = conn.cursor()
        curs.execute("SELECT n_jobs FROM projects")
        rows = curs.fetchall()
        conn.close()

        assert len(rows) == 1
        assert rows[0][0] == 3

    def test_n_completed_zero_at_creation(self, workspace):
        campaign_dir = self._make_campaign(workspace, name="ncomp0")

        conn = sqlite3.connect(str(campaign_dir / "campaign.db"))
        curs = conn.cursor()
        curs.execute("SELECT n_completed FROM projects")
        (n_completed,) = curs.fetchone()
        conn.close()

        assert n_completed == 0

    # ------------------------------------------------------------------
    # _sync_project_status
    # ------------------------------------------------------------------

    def test_sync_returns_none_for_missing_db(self, tmp_path):
        result = _sync_project_status(tmp_path / "nonexistent")
        assert result is None

    def test_sync_no_output_files(self, workspace):
        campaign_dir = self._make_campaign(workspace, name="noout")
        proj_dir = self._project_dir(campaign_dir)

        result = _sync_project_status(proj_dir)

        assert result["n_jobs"] == 1
        assert result["n_completed"] == 0

    def test_sync_one_complete(self, workspace):
        campaign_dir = self._make_campaign(workspace, name="onecomp")
        proj_dir = self._project_dir(campaign_dir)
        self._write_output(proj_dir, 0)

        result = _sync_project_status(proj_dir)

        assert result["n_completed"] == 1
        assert result["n_jobs"] == 1

    def test_sync_stub_file_ignored(self, workspace):
        """Output files smaller than 1 KB must not count as completed."""
        campaign_dir = self._make_campaign(workspace, name="stub", n_fake_samples=2)
        proj_dir = self._project_dir(campaign_dir)
        self._write_output(proj_dir, 0, size=512)   # stub — too small
        self._write_output(proj_dir, 1, size=2048)  # real

        result = _sync_project_status(proj_dir)

        assert result["n_completed"] == 1
        assert result["n_jobs"] == 2

    def test_sync_updates_project_db(self, workspace):
        """project.db jobs table must be updated to 'completed' by _sync."""
        campaign_dir = self._make_campaign(workspace, name="projdb")
        proj_dir = self._project_dir(campaign_dir)
        self._write_output(proj_dir, 0)

        _sync_project_status(proj_dir)

        conn = sqlite3.connect(str(proj_dir / "project.db"))
        curs = conn.cursor()
        curs.execute("SELECT status FROM jobs WHERE jobid = 0")
        (status,) = curs.fetchone()
        conn.close()
        assert status == "completed"

    # ------------------------------------------------------------------
    # cmd_sync status transitions
    # ------------------------------------------------------------------

    def _db_row(self, campaign_dir, field):
        conn = sqlite3.connect(str(campaign_dir / "campaign.db"))
        curs = conn.cursor()
        curs.execute(f"SELECT {field} FROM projects LIMIT 1")
        val = curs.fetchone()[0]
        conn.close()
        return val

    def test_cmd_sync_all_complete(self, workspace):
        campaign_dir = self._make_campaign(workspace, name="allcomp")
        proj_dir = self._project_dir(campaign_dir)
        self._write_output(proj_dir, 0)

        class _Args:
            campaign = str(campaign_dir)

        cmd_sync(_Args())

        assert self._db_row(campaign_dir, "status") == "completed"
        assert self._db_row(campaign_dir, "n_completed") == 1

    def test_cmd_sync_partial(self, workspace):
        campaign_dir = self._make_campaign(workspace, name="partial", n_fake_samples=2)
        proj_dir = self._project_dir(campaign_dir)
        self._write_output(proj_dir, 0)  # only job 0 done

        class _Args:
            campaign = str(campaign_dir)

        cmd_sync(_Args())

        assert self._db_row(campaign_dir, "status") == "partial"
        assert self._db_row(campaign_dir, "n_completed") == 1

    def test_cmd_sync_no_output_preserves_status(self, workspace):
        campaign_dir = self._make_campaign(workspace, name="nostatus")

        class _Args:
            campaign = str(campaign_dir)

        cmd_sync(_Args())

        # Status must remain 'created' — no output files means no change.
        assert self._db_row(campaign_dir, "status") == "created"
