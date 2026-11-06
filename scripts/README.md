# scripts/

Helper scripts for the ICARUS NuMI nue CC inclusive analysis. They fall into five groups:

1. Tooling for Claude Code and the SL7 container.
2. Workflow helpers.
3. Studies of the detector-variation (det-sys) samples.
4. Background and sideband studies.
5. Event lookup.

Unless a script says otherwise:

- Python scripts run in the pixi env: `spineplot/.pixi/envs/default/bin/python scripts/<name>.py ...`
- The analysis scripts are read-only on ROOT files. They write only CSV, PNG, PDF and NPZ files to `--outdir` / `--output`.
- Every script has a full description at the top of the file and in `--help`.
- Input ROOT files are hard-coded near the top of each script, in `FILE`, `FILES` or `SETS`. They currently point at v31 files in `/exp/icarus/data/users/micarrig/`:
  - `merged_final_v31e_flux.root`: the main project, after stage 5 and the sideband graft.
  - `merged_nue_v31.root`: the nue-enhanced project.

  Edit those constants to point at a new production.
- Many of these scripts load whole trees from files of 30 GB or more. The GPVM has about 11 GB of RAM, so run the heavy ones on the bigger-memory server. The scans that read `events/NuMIFull/all` read it in chunks or cache it to an `.npz` file.

## 1. Tooling (Claude Code guard and container)

| File | What it does |
|------|--------------|
| `guard_bash.py` | Claude Code PreToolUse hook for the Bash tool. It is wired in `.claude/settings.local.json`. It splits each command into its simple parts and checks each part against the CLAUDE.md rules.<br>**deny:** `pkill`/`killall`, raw `jobsub_*`, `setup_icarus_env.sh`, deleting or moving files under `/pnfs` or a project dir, and writes to `guard_allow.txt`.<br>**ask:** grid submissions, `medulla.py --yes`, plotting on a low-memory node, `hadd` onto an existing file, `rm`/`mv`, `kill` and `htgettoken`.<br>In auto mode every "ask" becomes "deny". It never prints "allow". Tests: `batch/test_guard_bash.py`. |
| `guard_allow.txt` | The owner's allow list for `guard_bash.py`: one Python regex per line, full-matched against a single simple command. A match skips the ask rules but never the deny rules. **Only the owner edits this file.** |
| `in_sl7.sh` | Runs one command non-interactively inside the SL7 apptainer container (`fnal-dev-sl7:jsl`). It first sources `~/scripts/startSpine.sh` and `build/setup_medulla.sh`, then runs the command from the directory you called it from. It never runs `htgettoken` and has no jobsub inside. Usage: `scripts/in_sl7.sh '<command>'`. Use it for `hadd`, `run_systematics`, `to_gundam` and PyROOT. Do not use it for `jobsub_submit`/`jobsub_q`. |

## 2. Workflow helpers

| File | What it does |
|------|--------------|
| `graft_tree.py` | Copies a big merged file and adds one TTree from a small file to the copy. Example: adding `events/NuMIFull/nc_sideband_variationTree` to make `merged_final_v31e_flux.root`, as in the "nc_sideband detector systematics" section of CLAUDE.md.<br>Before writing, it checks that:<br>• the tree is absent from `--base`;<br>• it has the same entry count as `--ref`;<br>• Run/Subrun/Evt match entry by entry.<br>It never writes to `--base` or `--extra`, refuses an existing `--output`, and writes `<output>.partial` first. Supports `--dry-run`. Needs PyROOT, so run it through `in_sl7.sh`. |

## 3. Detector-variation (det-sys) studies

These scripts compare each det-sys variation with its CV sample, for two sets:

- **main:** CV = `cv`.
- **nue-enhanced:** CV = `NuMI_nue`.

The scripts depend on each other in this order:

```
run_corrected_pot.py (runs analysis/compute_corrected_pot_list.C per chunk)
    -> detsys_corrected_pot.csv            POT only (pot_th1d, failed_pot, corrected_pot, ...)
detsys_event_counts.py
    -> detsys_pot_and_counts.csv           + n_all, n_selected, n_sel_signal, ...
                                             (used by detsys_cutflow, detsys_resolution, plot_pot_table)
no script in the repo (made by hand)
    -> detsys_pot_and_counts_norm2p5e20.csv        + scale_to_2p5e20, n_*_2p5e20
    -> detsys_pot_and_counts_norm2p5e20_pctcv.csv  + cv_reference, n_*_pct_vs_cv
                                             (used by plot_detsys_counts, plot_detsys_multiplicity, plot_detsys_softmax)
```

`analysis/compute_corrected_pot.C` is the older, single-directory version of the macro. It only prints its results and is not part of this chain.

If `detsys_pot_and_counts.csv` changes, redo the two hand-made files.

Outputs default to `plots_v31/pot/`. Unless a script says otherwise, every count is scaled to 2.5e20 POT.

| File | What it does |
|------|--------------|
| `run_corrected_pot.py` | Runs `analysis/compute_corrected_pot_list.C` on every `detector_variation` sample in `selection.toml` and `selection_nue.toml`, plus `NuMI_nue`. It lists the files over xrootd, runs the macro in parallel chunks, and sums the results into `<out>/detsys_corrected_pot.csv`. Unreadable files go to `failed_files.txt`. Needs `$BEARER_TOKEN_FILE` and ROOT with xrootd. Supports `--dry-run` and `--only`. |
| `detsys_event_counts.py` | Adds event counts to the POT CSV for each sample: `all`-tree entries, selected, selected signal, background, non-nue, and uncategorized. Writes a new CSV and refuses an existing one. |
| `plot_pot_table.py` | Table images (`pot_table_<set>.png`) and `pot_table.csv` comparing uncorrected, failed and corrected POT, with file, subrun and event counts per sample. |
| `plot_detsys_counts.py` | Percent change from the CV in the selected, signal, background and non-nue counts, with statistical error bars (`detsys_pct_vs_cv_<set>.png`). Also writes a summary table image (`detsys_table_<set>.png`). |
| `detsys_cutflow.py` | Cutflow tables for one variation, or `all`, against the CV, using the per-cut flags in the `all` tree. Each table shows step and cumulative efficiency, the change from the CV in percentage points, the POT-scaled N and purity. Writes PNGs, CSVs and a PDF per variation. |
| `detsys_resolution.py` | Electron energy and cos(theta_beam) resolution for the CV and each nue-enhanced variation. Each resolution is fitted with a double-sided crystal ball. It has two modes:<br>• `splines`: applies the +1 sigma phase-1 spline weights to v31 NuMIFull signal events.<br>• `direct`: uses the variation sample itself.<br>Writes per-variation plots, fit tables, CSVs and a PDF. |
| `plot_detsys_multiplicity.py` | Reco particle multiplicities (e, p, pi, gamma, mu) and truth category of selected events, CV vs each variation, with ratio panels. Also plots the % change per reco interaction type (1e1p0pi, ...). Writes one PDF per set. |
| `plot_detsys_softmax.py` | Particle-ID softmax and primary softmax of each leading particle in selected events, CV vs every variation, with ratio panels. |
| `plot_detsys_matched_energy.py` | Matches the same true neutrinos between the CV and a variation on truth: E_nu, Q2, true electron energy and mode. Plots the migration of reco electron energy and cos(theta_beam), as 2D histograms normalized per column. Also writes `shared_events_*.png/.csv`: how many signal neutrinos the two samples share. |
| `plot_detsys_matched_position.py` | Position dependence of a variation, using the matched neutrinos above. Plots the reco electron start position in x, y and z (1D and 2D), CV vs variation, with ratio. Imports its constants from `plot_detsys_matched_energy.py`. |

## 4. Background and sideband studies (NuMIFull nominal MC)

All of these use raw, unweighted counts with no POT scaling and no flux weights.

| File | What it does |
|------|--------------|
| `background_topologies.py` | Particle content of each `true_category` in `selected` (or in `nc_sideband` with `--tree nc_sideband`). It reports, reco and true side by side:<br>• the fraction of events with at least one particle of each type;<br>• mean multiplicities;<br>• the most common (Ne, Np, Npi, Ngamma) combinations;<br>• reco vs true electron counts.<br>Writes CSVs and table PNGs. |
| `compare_region_kinematics.py` | Compares each background category in the signal region with the same category in the NC pi0 sideband. Variables: true E_nu, W, Q2, reco and true visible energy, and interaction mode. Writes one figure per category and a summary CSV. |
| `plot_nc_electron_energy.py` | Selected NC events vs signal: reco electron energy, and reco E_e / true photon energy, to tell single photons from two merged photons. Also plots against the photon opening angle and the vertex distance. |
| `plot_nc_electron_softmax.py` | Electron, photon and primary softmax of the leading reco electron, selected NC vs signal, area normalized. |
| `inverted_primary_sideband.py` | Emulates a candidate sideband from the `all` tree: the signal cuts with the leading electron's primary softmax cut inverted. It also makes variants with the energy cut and/or the muon veto put back. Reports each region's composition, the true particle behind the reco "electron" and its kinematics. Checks the rebuilt signal region against `selected`. |
| `photon_sideband_scan.py` | Scans candidate sidebands that copy the signal-region background: photon, exactly 1 gamma, at least 2 gamma, and low electron score. For each it reports composition, the true identity of the leading shower, kinematics, and KS distances to the SR background. Caches the preselection in `scan_events.npz`; rerun from the cache with `--from-npz`. |
| `plot_photon_sideband.py` | Plots the photon sideband vs the SR background, with the NC pi0 sideband as a reference: visible and shower energy, W and Q2, stacks by interaction mode, per-mode shapes, and the sideband reweighted in E_vis. Needs `scan_events.npz` from `photon_sideband_scan.py`. Options: `--threshold`, `--cats`. |
| `photon_sideband_eveto_scan.py` | Scans the photon sideband's electron-veto threshold (0 to 300 MeV): does allowing soft primary reco electrons help? Reports signal %, background makeup, DIS fraction per E_vis fit bin, and overlap with the SR. Caches `eveto_events.npz`; rerun from it with `--from-npz`. |

## 5. Event lookup (for event displays)

| File | What it does |
|------|--------------|
| `find_caf_events.py` | Finds which flat CAF files hold a list of (Run, Evt) events, given a CSV. It reads only the `rec.hdr` run/subrun/evt branches and stops once every event is found. `--files` takes globs, paths, `root://` URLs or `@list.txt`. Supports `--dry-run`. The output CSV is rewritten after each match, so an interrupted scan keeps what it found. |
| `find_raw_events.py` | For each event, uses SAM queries and xrootd to find:<br>• the raw compressed art file holding it, found by reading `EventAuxiliary`;<br>• its stage0 and stage1 children.<br>Default SAM definition: `run2_compression_production_..._numimajority_compressed_data`. Needs `$BEARER_TOKEN_FILE`. Writes a new CSV. |

`__pycache__/` exists because `plot_photon_sideband.py` and `plot_detsys_matched_position.py` import other scripts in this directory.
