#!/usr/bin/env python3
"""
Command-line version of makeSelectionPlots.py: make the selection plots for one
TOML and one input file. The original (hardcoded) script is left untouched.

Usage:
    python3 analysis/makeSelectionPlotsCLI.py --toml analysis/nue_plots.toml \
        --input /exp/icarus/data/users/micarrig/merged_v30.root \
        --output-dir plots_v30/selection
    python3 analysis/makeSelectionPlotsCLI.py ... --figures reco_nu true_nu --dry-run

Each figure is saved as plot_<name>.png (300 dpi) and plot_<name>.pkl.
"""
from plot_cli_common import main

# The default plotVars from makeSelectionPlots.py (the uncommented list).
DEFAULT_FIGURES = [
    'reco_leading_electron_energy',
    'reco_ele_beam_open_angle',
    'reco_leading_electron_start_dedx',
    'reco_leading_electron_axial_spread',
    'reco_leading_electron_directional_spread',
    'reco_leading_electron_start_x',
    'reco_leading_electron_start_y',
    'reco_leading_electron_start_z',
    'reco_leading_electron_start_dir_x',
    'reco_leading_electron_start_dir_y',
    'reco_leading_electron_start_dir_z',
    'reco_leading_electron_mip_softmax',
    'reco_leading_electron_photon_softmax',
    'reco_leading_electron_muon_softmax',
    'reco_leading_electron_proton_softmax',
    'reco_leading_electron_pion_softmax',
    'reco_leading_electron_electron_softmax',
    'reco_leading_electron_primary_softmax',
    'reco_leading_electron_secondary_softmax',
    'reco_leading_electron_vertex_distance',
    'reco_muon_multiplicity',
    'reco_photon_multiplicity',
    'reco_pion_multiplicity',
    'reco_proton_multiplicity',
    'reco_electron_multiplicity',
    'reco_nu',
    'true_neutrino_energy',
    'true_leading_electron_energy',
    'true_ele_beam_open_angle',
    'true_leading_electron_start_x',
    'true_leading_electron_start_y',
    'true_leading_electron_start_z',
    'true_leading_electron_start_dir_x',
    'true_leading_electron_start_dir_y',
    'true_leading_electron_start_dir_z',
    'true_muon_multiplicity',
    'true_photon_multiplicity',
    'true_pion_multiplicity',
    'true_proton_multiplicity',
    'true_electron_multiplicity',
    'true_momentum_transfer',
    'true_hadronic_invariant_mass',
    'true_nu',
    'nu_diff',
    'true_vs_reco_leading_electron_energy',
    'true_vs_reco_leading_electron_angle',
    'true_vs_reco_nu',
    'ele_diff',
    'angle_diff',
]


if __name__ == '__main__':
    main(__doc__, DEFAULT_FIGURES)
