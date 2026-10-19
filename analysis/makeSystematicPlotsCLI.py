#!/usr/bin/env python3
"""
Command-line version of makeSystematicPlots.py: make the systematic-band plots
for one TOML and one input file. The original (hardcoded) script is left
untouched.

The signal-only / background-only variants are separate named figures in the
TOML (the `_signal` / `_background` suffixes, set via each artist's category
`groups`), so selecting them is just a matter of the --figures subset.

--systematics-set names a list in analysis/systematic_sets.toml (xsec, flux,
detector, total, ...) and overrides the `systematics` of every SpineSystematics
artist, so the one TOML produces each breakdown in its own run - no hand-editing
the commented `systematics = ...` lines.

Usage:
    python3 analysis/makeSystematicPlotsCLI.py --toml analysis/systematics.toml \
        --input /exp/icarus/data/users/micarrig/merged_v30.root \
        --output-dir plots_v30/systematics_xsec --systematics-set xsec
    python3 analysis/makeSystematicPlotsCLI.py ... --figures sys_nuE sys_nuE_signal --dry-run

Each figure is saved as plot_<name>.png (300 dpi) and plot_<name>.pkl.
"""
from plot_cli_common import main

# The default plotVars from makeSystematicPlots.py.
DEFAULT_FIGURES = [
    'sys_eleE',
    'sys_eleE_signal',
    'sys_eleE_background',
    'sys_eleCosThetaNUMI',
    'sys_eleCosThetaNUMI_signal',
    'sys_eleCosThetaNUMI_background',
    'sys_nuE',
    'sys_nuE_signal',
    'sys_nuE_background',
    'sys_trueEleE',
    'sys_trueEleE_signal',
    'sys_trueEleE_background',
    'sys_trueEleCosThetaNUMI',
    'sys_trueEleCosThetaNUMI_signal',
    'sys_trueEleCosThetaNUMI_background',
]


if __name__ == '__main__':
    main(__doc__, DEFAULT_FIGURES, with_sys_set=True)
