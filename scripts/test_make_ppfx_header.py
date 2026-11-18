"""
Tests for scripts/make_ppfx_header.py and the generated
selection/include/ppfx_weights.h (the selection's true_ppfx_cv_weight).

Run: python3 -m pytest scripts/test_make_ppfx_header.py -v   (host python; g++ for the C++ test)
"""
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))
import make_ppfx_header as mph  # noqa: E402

HEADER = mph.DEFAULT_OUT
have_flux = Path(mph.DEFAULT_FLUX).exists()

FAKE = {
    'numu': (np.array([0.0, 1.0, 2.0]), np.array([0.5, 0.6])),
    'numubar': (np.array([0.0, 1.0, 2.0]), np.array([0.7, 0.8])),
    'nue': (np.array([0.0, 0.5, 1.0, 4.0]), np.array([1.1, 1.2, 1.3])),
    'nuebar': (np.array([0.0, 4.0]), np.array([0.9])),
}


def compile_and_eval(header_dir, queries, tmp_path):
    """Evaluate ppfx::cv_weight(pdg, E) for each (pdg, E) with a tiny C++ program."""
    src = tmp_path / 'eval.cc'
    src.write_text('#include <cstdio>\n#include <cstdlib>\n#include "ppfx_weights.h"\n'
                   'int main(int argc, char ** argv) {\n'
                   '  for(int i = 1; i + 1 < argc; i += 2)\n'
                   '    std::printf("%.17g\\n", ppfx::cv_weight(std::atoi(argv[i]), std::strtod(argv[i + 1], nullptr)));\n'
                   '}\n')
    exe = tmp_path / 'eval'
    subprocess.run(['g++', '-std=c++17', '-Wall', '-Werror', '-I', str(header_dir), str(src), '-o', str(exe)], check=True)
    args = [str(x) for q in queries for x in (int(q[0]), repr(float(q[1])))]
    out = subprocess.run([str(exe)] + args, check=True, capture_output=True, text=True).stdout.split()
    return [float(x) for x in out]


def stage5_lookup(tables, pdg, e):
    """add_numi_flux_sys.py's lookup (values[searchsorted(edges, e) - 1]), NaN for other flavors."""
    name = {14: 'numu', -14: 'numubar', 12: 'nue', -12: 'nuebar'}.get(pdg)
    if name is None:
        return float('nan')
    edges, values = tables[name]
    return values[np.searchsorted(edges, e) - 1]


@pytest.mark.skipif(shutil.which('g++') is None, reason='needs g++')
def test_lookup_bins_and_edges(tmp_path):
    (tmp_path / 'ppfx_weights.h').write_text(mph.render(FAKE, 'fake.root', 'fhc', '0'))
    q = [(14, 0.5), (14, 1.0), (14, 1.5), (14, 0.0),        # [low, high) bins
         (14, 2.0), (14, 25.0), (14, -0.1), (14, float('nan')),  # outside / NaN -> 1
         (-14, 1.9), (12, 0.75), (12, 3.9), (12, 4.0), (-12, 2.0),
         (16, 1.0), (0, 1.0), (2212, 1.0)]                   # other PDG -> 1
    assert compile_and_eval(tmp_path, q, tmp_path) == \
        [0.5, 0.6, 0.6, 0.5, 1.0, 1.0, 1.0, 1.0, 0.8, 1.2, 1.3, 1.0, 0.9, 1.0, 1.0, 1.0]


@pytest.mark.skipif(not have_flux, reason='needs flux.root')
def test_header_up_to_date():
    assert mph.main(['--check']) == 0, 'rerun scripts/make_ppfx_header.py'


@pytest.mark.skipif(not have_flux or shutil.which('g++') is None, reason='needs flux.root and g++')
def test_header_matches_stage5_lookup(tmp_path):
    """The committed header agrees with add_numi_flux_sys.py's lookup inside each bin."""
    tables = mph.read_tables(mph.DEFAULT_FLUX, 'fhc')
    rng = np.random.default_rng(1)
    q = []
    for pdg, name in ((14, 'numu'), (-14, 'numubar'), (12, 'nue'), (-12, 'nuebar')):
        edges = tables[name][0]
        q += [(pdg, e) for e in rng.uniform(edges[0], edges[-1], 200)]
        q += [(pdg, 0.5 * (lo + hi)) for lo, hi in zip(edges[:-1], edges[1:])]
    got = compile_and_eval(HEADER.parent, q, tmp_path)
    want = [stage5_lookup(tables, p, e) for p, e in q]
    assert got == want


def test_check_binning():
    mph.check_binning('numu', [0.0, 1.0, 2.0], [1.0, 1.0])
    with pytest.raises(ValueError):
        mph.check_binning('numu', [0.0, 1.0, 1.0], [1.0, 1.0])   # empty bin
    with pytest.raises(ValueError):
        mph.check_binning('numu', [0.0, 1.0], [1.0, 1.0])        # count mismatch
