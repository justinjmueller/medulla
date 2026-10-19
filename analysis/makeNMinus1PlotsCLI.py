#!/usr/bin/env python3
"""
Command-line version of makeNMinus1Plots.py: make the N-1 cut plots for one TOML
and one input file. The original (hardcoded) script is left untouched.

Usage:
    python3 analysis/makeNMinus1PlotsCLI.py --toml analysis/nminus1_plots.toml \
        --input /exp/icarus/data/users/micarrig/merged_v30.root \
        --output-dir plots_v30/nminus1
    python3 analysis/makeNMinus1PlotsCLI.py ... --figures nminus1_flash_cut --dry-run

Each figure is saved as plot_<name>.png (300 dpi) and plot_<name>.pkl.
"""
from plot_cli_common import main

# The default plotVars from makeNMinus1Plots.py.
DEFAULT_FIGURES = [
    'nminus1_fiducial_x',
    'nminus1_fiducial_y',
    'nminus1_fiducial_z',
    'nminus1_flash_cut',
    'nminus1_containment_x',
    'nminus1_containment_y',
    'nminus1_containment_z',
    'nminus1_particle_multiplicity',
    'nminus1_electron_softmax',
    'nminus1_primary_softmax',
    'nminus1_no_muons',
    'nminus1_vertex_distance',
    'nminus1_electron_energy',
    'nminus1_leading_muon_ke',
]


if __name__ == '__main__':
    main(__doc__, DEFAULT_FIGURES)
