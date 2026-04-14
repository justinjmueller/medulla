"""
Shared pytest fixtures for medulla/batch tests.
"""
import textwrap
from pathlib import Path

import pytest


def _write_toml(directory: Path, filename: str, content: str) -> Path:
    """Write *content* (dedented) into *directory*/*filename* and return the path."""
    p = directory / filename
    p.write_text(textwrap.dedent(content))
    return p


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
