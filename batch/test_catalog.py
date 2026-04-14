"""
Unit tests for catalog.resolve_samples()
=========================================

Tests the [[include_samples]] directive resolution logic in catalog.py.
resolve_samples() reads a shared sample catalog TOML and expands
[[include_samples]] directives in a selection config into concrete
[[sample]] blocks, applying enable/disable filtering along the way.

The function under test — resolve_samples(cfg, catalog_path, enable_keys=None) —
has the following contract:

    Parameters
    ----------
    cfg : dict
        A parsed TOML configuration dictionary (e.g. from toml.load()).
        May contain an 'include_samples' key (list of directives).
    catalog_path : str | Path
        Path to the shared sample catalog TOML file.
    enable_keys : list[str] | None
        If provided, overrides each directive's own 'enable' list.
        Only samples whose catalog key is in this list will have
        disable = False. All others will have disable = True.

    Returns
    -------
    cfg : dict
        The configuration dictionary with:
        - 'include_samples' removed.
        - Resolved [[sample]] entries appended to cfg['sample'].
        - Each resolved sample has: name, path, ismc, disable.
        - Existing inline [[sample]] entries are preserved.

    Raises
    ------
    KeyError
        If a requested key is not found in the catalog.
    FileNotFoundError
        If the catalog file does not exist.

The catalog TOML format:
    [[sample]]
    key = "unique_id"
    name = "human_name"
    path = "/some/glob/path/*.root"
    ismc = true
    experiment = "sbnd"

Run with:
    pytest               # from medulla/ (uses pytest.ini testpaths)
    pytest batch/        # from medulla/ explicitly
    make pytest          # from the build directory
"""

import textwrap
from pathlib import Path

import toml
import pytest

# ---------------------------------------------------------------------------
# Import the function under test.  When the module doesn't exist yet the
# entire test file is skipped with a clear message.
# ---------------------------------------------------------------------------
try:
    from catalog import resolve_samples
except ImportError:
    pytestmark = pytest.mark.skip(
        reason="catalog.py not found or resolve_samples not implemented"
    )
    # Define a stub so the rest of the file parses without NameError.
    def resolve_samples(cfg, catalog_path, enable_keys=None):
        raise NotImplementedError


# ---------------------------------------------------------------------------
# Helpers — local factories (shared fixtures live in conftest.py)
# ---------------------------------------------------------------------------

def _write_toml(directory: Path, filename: str, content: str) -> Path:
    """Write a TOML string into *directory*/*filename* and return its path."""
    p = directory / filename
    p.write_text(textwrap.dedent(content))
    return p


def _make_config_with_includes(
    tmp: Path,
    keys: list,
    enable: list = None,
    inline_samples: str = "",
) -> dict:
    """
    Build and return a parsed config dict that contains an
    ``[[include_samples]]`` directive.  Optionally add inline ``[[sample]]``
    blocks and an ``enable`` list inside the directive.
    """
    enable_line = ""
    if enable is not None:
        enable_line = f'enable = {enable}'

    keys_line = f'keys = {keys}'

    content = f"""\
        [general]
        output = "test_output"

        {inline_samples}

        [[include_samples]]
        {keys_line}
        {enable_line}
    """
    path = _write_toml(tmp, "config.toml", content)
    return toml.load(str(path))


# ===================================================================
# T1.1 — Basic resolution: keys are resolved to [[sample]] entries
# ===================================================================

class TestBasicResolution:
    """resolve_samples() correctly resolves catalog keys into sample dicts."""

    def test_single_key_resolved(self, tmp_path, catalog):
        cfg = _make_config_with_includes(tmp_path, keys=["sbnd_mc_nominal"])
        result = resolve_samples(cfg, catalog)

        samples = result.get("sample", [])
        assert len(samples) == 1
        s = samples[0]
        assert s["name"] == "sbnd"
        assert s["path"] == "/pnfs/sbnd/mc/nominal/*.flat.root"
        assert s["ismc"] is True

    def test_multiple_keys_resolved(self, tmp_path, catalog):
        cfg = _make_config_with_includes(
            tmp_path, keys=["sbnd_mc_nominal", "icarus_mc_nominal"]
        )
        result = resolve_samples(cfg, catalog)

        samples = result.get("sample", [])
        assert len(samples) == 2
        names = {s["name"] for s in samples}
        assert names == {"sbnd", "icarus"}

    def test_include_samples_key_removed(self, tmp_path, catalog):
        cfg = _make_config_with_includes(tmp_path, keys=["sbnd_mc_nominal"])
        result = resolve_samples(cfg, catalog)
        assert "include_samples" not in result

    def test_catalog_only_fields_excluded(self, tmp_path, catalog):
        """Catalog-internal fields (key, experiment, metadata) must not
        leak into the resolved sample entries."""
        cfg = _make_config_with_includes(tmp_path, keys=["sbnd_mc_nominal"])
        result = resolve_samples(cfg, catalog)

        for s in result["sample"]:
            assert "key" not in s
            assert "experiment" not in s


# ===================================================================
# T1.2 — Enable / disable logic
# ===================================================================

class TestEnableDisable:
    """The enable mechanism correctly sets the disable flag."""

    @pytest.mark.parametrize("directive_enable,fn_enable_keys,expect_disabled", [
        # No enable anywhere → everything disabled
        (None,                              None,                       {"sbnd", "icarus"}),
        # Directive enable list → only named key enabled
        (["sbnd_mc_nominal"],               None,                       {"icarus"}),
        # Directive enable list covers all keys → nothing disabled
        (["sbnd_mc_nominal", "icarus_mc_nominal"], None,                set()),
        # enable_keys argument overrides directive's own enable list
        (["sbnd_mc_nominal"],               ["icarus_mc_nominal"],      {"sbnd"}),
    ], ids=["no_enable", "directive_enables_sbnd", "directive_enables_all", "fn_arg_overrides"])
    def test_enable_logic(
        self, tmp_path, catalog,
        directive_enable, fn_enable_keys, expect_disabled,
    ):
        cfg = _make_config_with_includes(
            tmp_path,
            keys=["sbnd_mc_nominal", "icarus_mc_nominal"],
            enable=directive_enable,
        )
        result = resolve_samples(cfg, catalog, enable_keys=fn_enable_keys)

        disabled_names = {s["name"] for s in result["sample"] if s.get("disable", False)}
        enabled_names  = {s["name"] for s in result["sample"] if not s.get("disable", False)}
        assert disabled_names == expect_disabled
        assert enabled_names  == ({"sbnd", "icarus"} - expect_disabled)


# ===================================================================
# T1.3 — Coexistence with inline [[sample]] blocks
# ===================================================================

class TestInlineCoexistence:
    """Inline samples survive resolution unchanged."""

    INLINE = """\
        [[sample]]
        name = "local_test"
        path = "/tmp/test/*.root"
        ismc = false
        disable = false
    """

    def test_inline_preserved(self, tmp_path, catalog):
        cfg = _make_config_with_includes(
            tmp_path,
            keys=["sbnd_mc_nominal"],
            enable=["sbnd_mc_nominal"],
            inline_samples=self.INLINE,
        )
        result = resolve_samples(cfg, catalog)

        samples = result.get("sample", [])
        names = [s["name"] for s in samples]
        assert "local_test" in names, "Inline sample must be preserved"
        assert "sbnd" in names, "Resolved sample must also be present"

    def test_inline_not_modified(self, tmp_path, catalog):
        cfg = _make_config_with_includes(
            tmp_path,
            keys=["sbnd_mc_nominal"],
            enable=["sbnd_mc_nominal"],
            inline_samples=self.INLINE,
        )
        result = resolve_samples(cfg, catalog)

        inline = [s for s in result["sample"] if s["name"] == "local_test"]
        assert len(inline) == 1
        s = inline[0]
        assert s["path"] == "/tmp/test/*.root"
        assert s["ismc"] is False
        assert s.get("disable", False) is False

    def test_total_count(self, tmp_path, catalog):
        cfg = _make_config_with_includes(
            tmp_path,
            keys=["sbnd_mc_nominal", "icarus_mc_nominal"],
            inline_samples=self.INLINE,
        )
        result = resolve_samples(cfg, catalog)
        assert len(result["sample"]) == 3  # 1 inline + 2 resolved


# ===================================================================
# T1.4 — Missing key raises KeyError
# ===================================================================

class TestMissingKey:
    """Requesting a catalog key that doesn't exist must fail loudly."""

    def test_missing_key_raises(self, tmp_path, catalog):
        cfg = _make_config_with_includes(tmp_path, keys=["nonexistent_key"])
        with pytest.raises(KeyError):
            resolve_samples(cfg, catalog)

    def test_partial_missing_key_raises(self, tmp_path, catalog):
        """Even if some keys are valid, one bad key should raise."""
        cfg = _make_config_with_includes(
            tmp_path, keys=["sbnd_mc_nominal", "nonexistent_key"]
        )
        with pytest.raises(KeyError):
            resolve_samples(cfg, catalog)


# ===================================================================
# T1.5 — No [[include_samples]] is a no-op
# ===================================================================

class TestNoIncludeDirective:
    """Configs without include_samples are returned unchanged."""

    def test_no_include_noop(self, tmp_path, catalog):
        content = """\
            [general]
            output = "test_output"

            [[sample]]
            name = "inline_only"
            path = "/some/path/*.root"
            ismc = true
        """
        path = _write_toml(tmp_path, "config.toml", content)
        cfg = toml.load(str(path))

        result = resolve_samples(cfg, catalog)

        assert result["sample"] == cfg["sample"]
        assert "include_samples" not in result

    def test_no_samples_at_all(self, tmp_path, catalog):
        """A config with neither [[sample]] nor [[include_samples]]."""
        content = """\
            [general]
            output = "test_output"
        """
        path = _write_toml(tmp_path, "config.toml", content)
        cfg = toml.load(str(path))

        result = resolve_samples(cfg, catalog)

        assert result.get("sample", []) == []


# ===================================================================
# T1.6 — Empty keys list
# ===================================================================

class TestEmptyKeys:
    """Edge case: include_samples with an empty keys list."""

    def test_empty_keys_resolves_nothing(self, tmp_path, catalog):
        cfg = _make_config_with_includes(tmp_path, keys=[])
        result = resolve_samples(cfg, catalog)

        assert result.get("sample", []) == []
        assert "include_samples" not in result


# ===================================================================
# T1.7 — Missing catalog file raises FileNotFoundError
# ===================================================================

class TestMissingCatalog:
    """A missing catalog file must fail with FileNotFoundError."""

    def test_missing_catalog_file(self, tmp_path):
        cfg = _make_config_with_includes(tmp_path, keys=["sbnd_mc_nominal"])
        bogus_path = tmp_path / "does_not_exist.toml"
        with pytest.raises(FileNotFoundError):
            resolve_samples(cfg, bogus_path)


# ===================================================================
# T1.8 — Multiple [[include_samples]] directives
# ===================================================================

class TestMultipleDirectives:
    """A config can have more than one [[include_samples]] block,
    e.g. one for SBND samples and one for ICARUS samples with
    different enable lists."""

    def test_two_directives(self, tmp_path, catalog):
        content = """\
            [general]
            output = "test_output"

            [[include_samples]]
            keys = ["sbnd_mc_nominal", "sbnd_offbeam"]
            enable = ["sbnd_mc_nominal"]

            [[include_samples]]
            keys = ["icarus_mc_nominal", "icarus_onbeam"]
            enable = ["icarus_onbeam"]
        """
        path = _write_toml(tmp_path, "config.toml", content)
        cfg = toml.load(str(path))
        result = resolve_samples(cfg, catalog)

        samples = result.get("sample", [])
        assert len(samples) == 4

        sbnd_mc  = [s for s in samples if s["path"] == "/pnfs/sbnd/mc/nominal/*.flat.root"]
        sbnd_off = [s for s in samples if s["path"] == "/pnfs/sbnd/data/offbeam/*.flat.root"]
        ic_mc    = [s for s in samples if s["path"] == "/pnfs/icarus/mc/nominal/*.flat.root"]
        ic_on    = [s for s in samples if s["path"] == "/pnfs/icarus/data/onbeam/*.flat.root"]

        assert len(sbnd_mc) == 1  and sbnd_mc[0].get("disable", False)  is False
        assert len(sbnd_off) == 1 and sbnd_off[0].get("disable", False) is True
        assert len(ic_mc) == 1    and ic_mc[0].get("disable", False)    is True
        assert len(ic_on) == 1    and ic_on[0].get("disable", False)    is False


# ===================================================================
# T1.9 — Duplicate keys across directives
# ===================================================================

class TestDuplicateKeys:
    """The same key appearing in two separate directives produces two
    entries (no cross-directive deduplication)."""

    def test_duplicate_key_across_directives(self, tmp_path, catalog):
        content = """\
            [general]
            output = "test_output"

            [[include_samples]]
            keys = ["sbnd_mc_nominal"]
            enable = ["sbnd_mc_nominal"]

            [[include_samples]]
            keys = ["sbnd_mc_nominal"]
            enable = ["sbnd_mc_nominal"]
        """
        path = _write_toml(tmp_path, "config.toml", content)
        cfg = toml.load(str(path))
        result = resolve_samples(cfg, catalog)

        sbnd_entries = [
            s for s in result["sample"]
            if s["path"] == "/pnfs/sbnd/mc/nominal/*.flat.root"
        ]
        assert len(sbnd_entries) == 2, (
            "Each directive resolves independently; duplicate keys are not deduplicated"
        )


# ===================================================================
# T1.10 — Mutation contract
# ===================================================================

class TestMutationContract:
    """resolve_samples() mutates cfg in place (pops include_samples,
    sets sample).  Other fields must not be corrupted."""

    def test_non_include_fields_unchanged(self, tmp_path, catalog):
        cfg = _make_config_with_includes(tmp_path, keys=["sbnd_mc_nominal"])
        general_before = dict(cfg["general"])

        resolve_samples(cfg, catalog)

        assert cfg["general"] == general_before, "Unrelated fields must survive resolution"
        assert "include_samples" not in cfg, "include_samples must be consumed, not left in cfg"


# ===================================================================
# T1.11 — Integration with get_samples() contract
# ===================================================================

class TestGetSamplesIntegration:
    """After resolve_samples(), the resulting config should be
    consumable by the existing get_samples() pipeline: i.e., every
    resolved sample has the fields that get_samples() expects
    (name, path, ismc, and optionally disable)."""

    REQUIRED_FIELDS = {"name", "path", "ismc"}

    def test_resolved_samples_have_required_fields(self, tmp_path, catalog):
        all_keys = [
            "sbnd_mc_nominal", "sbnd_offbeam",
            "icarus_mc_nominal", "icarus_onbeam",
        ]
        cfg = _make_config_with_includes(tmp_path, keys=all_keys, enable=all_keys)
        result = resolve_samples(cfg, catalog)

        for s in result["sample"]:
            missing = self.REQUIRED_FIELDS - set(s.keys())
            assert not missing, (
                f"Sample '{s.get('name', '?')}' is missing fields: {missing}"
            )

    def test_disable_is_bool(self, tmp_path, catalog):
        cfg = _make_config_with_includes(
            tmp_path,
            keys=["sbnd_mc_nominal", "icarus_mc_nominal"],
            enable=["sbnd_mc_nominal"],
        )
        result = resolve_samples(cfg, catalog)

        for s in result["sample"]:
            if "disable" in s:
                assert isinstance(s["disable"], bool), (
                    f"Sample '{s['name']}' disable field must be bool, "
                    f"got {type(s['disable'])}"
                )
