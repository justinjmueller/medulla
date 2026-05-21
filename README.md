# Medulla

**Medulla** is a modular C++ framework for physics selection development, systematic uncertainty evaluation, and analysis plotting in LArTPC-based neutrino experiments. It is the companion analysis layer to the [SPINE](https://github.com/DeepLearnPhysics/spine) deep-learning reconstruction framework, providing a declarative TOML-based configuration interface designed for rapid iteration across SBND and ICARUS.

---

## Features

- **Declarative selection configuration** — entire event selections, including cuts, branch variables, and truth-level categories, are expressed in human-readable TOML files without writing C++ per analysis.
- **Flexible tree modes** — selections operate at the event level (`event`), reconstructed-interaction level (`reco`), or truth-interaction level (`true`), with optional particle-level loops and two-particle (bivariate) observables.
- **Cross-type logic** — reco-mode trees can apply truth-level cuts (and vice versa) by following the reco↔truth match, with NaN fill-values for unmatched objects.
- **MCTruth integration** — generator-level cuts and variables (`type = "mctruth"`) operate directly on `SRTrueInteraction` objects, enabling GENIE-level signal definitions without touching SPINE truth objects.
- **Batch processing** — a Python batch layer handles HTCondor job creation, submission, monitoring, and output sync for large-scale grid campaigns at Fermilab.
- **Campaign coordination** — the campaign layer automates multi-analysis, multi-experiment, multi-role batch runs with a local registry and SQLite-backed state tracking.
- **Systematic uncertainties** — a dedicated systematics layer produces reweighted output trees alongside the nominal selection.
- **Plotting** — the `spineplot` Python library provides publication-quality 1D/2D spectra, efficiency curves, confusion matrices, and ROC curves built on top of the sbruce tree format.

---

## Requirements

Medulla requires the SBN software stack available via CVMFS:

| Dependency | Version |
|---|---|
| `sbnana` | `v10_01_02_01 -q e26:prof` |
| `cmake` | `≥ 3.14` |
| Python | `≥ 3.9` (interpreter only; no install needed) |

Grid submission additionally requires `jobsub_client` and valid OIDC tokens via `htgettoken`.

---

## Installation

```bash
# Source the experiment environment (choose one):
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh    # SBND
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh  # ICARUS

# Set up dependencies:
setup sbnana v10_01_02_01 -q e26:prof
setup cmake v3_27_4

# Clone and build:
git clone https://github.com/justinjmueller/medulla.git medulla
cd medulla && git checkout v1.0.3
mkdir build && cd build
cmake .. && make -j4
```

The build produces two executables under `build/selection/`:

| Binary | Purpose |
|---|---|
| `medulla` | Run a TOML-configured selection over CAF input files |
| `validate` | Generate synthetic test inputs and validate framework output |

---

## Quick Start

```bash
# Run the example CCQE-like selection (requires XRootD token):
htgettoken -a htvaultprod.fnal.gov --vaulttokenttl=1d -i sbnd
./build/selection/medulla tutorial/examples/example01_ccqe.toml
```

The output is a ROOT file containing one `TDirectory` per sample, each holding the configured `TTree`s in `events/<sample>/<tree>` layout.

For a complete walkthrough of the configuration format, all cut and variable types, batch submission, and the campaign layer, see the **[tutorial](tutorial/tutorial.md)**.

---

## Repository Layout

```
medulla/
├── selection/
│   ├── include/          # C++ headers: cuts, variables, particles, MCTruth, framework
│   ├── src/              # Framework, main executable, and validation binary
│   ├── test/             # Validation TOML (test.toml) and pytest module
│   └── toml/             # Analysis selection TOMLs and shared sample catalog
├── batch/
│   ├── medulla.py        # Single-project batch CLI
│   ├── campaign.py       # Multi-analysis campaign CLI
│   ├── utilities.py      # Shared submission and sync utilities
│   ├── submit.sh         # Grid worker script
│   └── test_*.py         # Python unit tests for the batch layer
├── systematics/          # Systematic uncertainty evaluation
├── spineplot/            # Python plotting library
├── shared/               # Shared C++ library (configuration, utilities)
├── cmake/                # CMake modules (environment checks, path setup)
├── tutorial/             # Tutorial document and example TOML files
├── conftest.py           # Root pytest fixture (framework validation environment)
└── pytest.ini            # Test discovery configuration
```

---

## Testing

The full test suite is run via CMake after a successful build:

```bash
make pytest
```

This rebuilds `medulla` and `validate` if needed, injects `MEDULLA_BUILD_DIR` into the environment, and runs both test suites:

- **Python unit tests** (`batch/`) — cover the catalog parser, campaign creation, and job state transitions using synthetic fixtures; no grid credentials or ROOT files required.
- **Framework validation tests** (`selection/test/`) — run the `validate` binary to generate synthetic CAF inputs, process them with `medulla`, and assert 124 named conditions across 14 validation groups covering every combination of tree mode, sample type, and cut variant.

To run only the batch tests (no build required):

```bash
python3 -m pytest batch/ -v
```

To run the framework tests against an existing build:

```bash
MEDULLA_BUILD_DIR=/path/to/build/selection python3 -m pytest selection/test/ -v
```

Framework tests are automatically skipped if `MEDULLA_BUILD_DIR` is unset. Full details of the validation design — including the group table and manual invocation of `validate` — are in the [Testing and Validation section of the tutorial](tutorial/tutorial.md#testing-and-validation).

---

## Documentation

The [tutorial](tutorial/tutorial.md) covers:

- Full TOML configuration reference (parameters, samples, trees, cuts, variables, categories)
- Tree modes, cross-type cut logic, particle-level and bivariate branches, MCTruth integration
- Running a selection locally and authenticating for XRootD access
- Single-project batch workflow (`medulla.py`)
- Multi-analysis campaign workflow (`campaign.py`)
- Test suite design and manual validation procedures

---

## Contact

Questions and discussion: **#medulla** in the SBN Slack workspace.

Bug reports and feature requests: open an issue on [GitHub](https://github.com/justinjmueller/medulla).