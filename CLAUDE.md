# Medulla: ICARUS nue CC inclusive working copy

Medulla is a C++/TOML selection and systematics framework for LArTPC analyses
(SBND and ICARUS), with a Python batch layer for Fermilab grid jobs and a
`spineplot` library for plotting. This working copy is used for the ICARUS NuMI
nue CC inclusive cross-section measurement (analysis name
`icarus_nueCC_inclusive_2026`). Owner: micarrig. Upstream: justinjmueller/medulla.

See `README.md` and `tutorial/tutorial.md` for the full framework docs. This file
is the short operational version for agents: what runs where, what is safe, and
the order of the pipeline.

## Rules for agents (read first)

1. Never submit grid jobs without explicit approval in the current conversation.
   This covers `medulla.py --launch-jobs`, `--test-job`, `--variation-phase1`,
   `--variation-phase2`, and any raw `jobsub_submit`. Approval is per stage: an
   OK to launch stage 2 is not an OK to start stage 3 or to relaunch failures.
   Before asking, run the same command with `--dry-run` and show the owner its
   output: project, number of jobs, `--tag` and the git SHA it resolves to,
   memory/disk, and the exact `jobsub_submit` command. The owner does
   not use `campaign.py`; do not use it unless asked.
2. Never delete outputs on your own initiative. `check_project_status` (called by
   `medulla.py` whenever a project exists) and `--check-zombies` delete stub
   (< 1 KB) and zombie/recovered ROOT files and reset their jobs to `pending`.
   `medulla.py --yes` / `-y` does this without prompting, after printing a table
   of every file (name, size, creation time). Use `--yes` only when the user has
   approved that specific run in the current conversation. Deletion is
   irreversible, so run with stdout captured (e.g. `tee`) so the table is kept.
   Corrupt-file detection beyond the size check needs `-Z`.
3. `medulla.py` without `--yes`, and `campaign.py` (which has no `--yes` yet),
   block on `input("... [Y/N]")`. Do not pipe `yes` into them. If a command stops
   at a prompt you were not asked to answer, stop and report.
4. Never handle credentials. `htgettoken` needs the user in a browser the first
   time. If a command fails on auth, report it and let the user re-authenticate.
5. Do not push, force-push, or change the remote without being asked. Grid jobs
   clone `--tag` from GitHub, so unpushed commits do not reach the grid, and a
   moving branch such as `develop` can change between jobs in one campaign. For
   anything beyond a smoke test, ask which tag or SHA to use.
6. Never read or print `.root` or `.df` files, and never `cat` a whole `*.toml`
   under `selection/toml/` (they are about 1400 lines). Use grep or targeted
   `sed -n` ranges.
7. Do not edit anything under `.history/` (VS Code Local History, not source).

8. Never run a script, binary, or command you do not know. "Know" means it is
   named in this file or you have read it in this session and understand what
   it submits, deletes, and writes. Otherwise read it, tell the owner what it
   does, and ask before running it. This applies most to anything near jobs
   (`*.sh`, `submit*`, `sandbox.sh`, `gensys.py`, notebooks, one-off scripts in
   the repo root or `analysis/`). Scripts that only look are fine to read, never
   to run blindly. Do not source env scripts without reading them (see Token).
9. While running or preparing jobs, never delete, move, overwrite, or edit any
   script, TOML, or source file without asking, even to "fix" a problem. Grid
   jobs ship `batch/*.sh` and the project's frozen configs, so a quiet edit
   changes what runs. If a change is needed, show the diff you propose and wait.
   For anything new (plot configs, helpers), create a new file with a new name
   and leave the original untouched. Edits the owner asked for in this
   conversation are fine.

10. Ask before running anything heavy, including scripts you wrote yourself in
   this session. Plotting scripts load whole trees into memory (the plotting
   inputs are 30 GB+ files) and have nearly crashed the shared node. Before any
   such run: use `--dry-run` if it has one, say what it will read, how much
   memory it may use, and where it writes, and wait for a yes. Never start these
   with `&`/`nohup` and walk away. Do not `pkill -f` by script name from Bash
   (the pattern matches the shell running the command); use
   `pgrep -fa "[n]ame"` to look and kill by PID.

Enforcement and shortcuts:
- `scripts/guard_bash.py` is a PreToolUse hook for Bash (wired in
  `.claude/settings.local.json`). It denies `pkill`/`killall`, `setup_icarus_env.sh`,
  raw `jobsub_*` and deleting under `/pnfs` or a project dir, and asks for
  `medulla.py` without `--dry-run`, plotting scripts, merges, stages 5 and 6, and
  `htgettoken`, and `kill` (even by PID). In auto mode (and `bypassPermissions`/
  `dontAsk`) every "ask" becomes "deny", because those modes would approve it
  without the owner seeing a prompt; the owner switches out of auto mode
  (shift+tab) to approve such a command. Exception: a plotting run (alone, not
  chained with another guarded command) stays "ask" in auto mode, and
  `autoMode.soft_deny` in `.claude/settings.local.json` lets it through only when
  the owner's latest message explicitly approves that run, so they can approve
  from a phone. Grid/merge/flux/gundam stages and `kill` are in
  `autoMode.hard_deny` as a second layer. If a command
  is denied or asked, do not rephrase it to get around the guard; tell the owner.
  Tests: `batch/test_guard_bash.py`.
- The skill `/medulla-stage <stage> <project>` (`.claude/skills/medulla-stage/`) is
  the order of operations for one stage: preflight, dry-run, ask, run, wait,
  verify, record in `<proj>/STATE.md`, notify. Use it for every stage.
- The skill `/medulla-dev` (`.claude/skills/medulla-dev/`) is the local dev loop:
  compile medulla, run a selection or `run_systematics`, and run the test suites,
  via `scripts/in_sl7.sh`. Use it for build/run/test; it is not for grid stages.

## Where things run (environments)

Everything runs on a Fermilab ICARUS GPVM (`/pnfs/icarus/...`, `/exp/icarus/...`).
There are three environments; use the right one per step.

1. Host shell, NOT in the container: grid submission. `jobsub_submit` and
   `jobsub_q` must be run outside the SL7 container and need nothing but a valid
   token. This covers `medulla.py` creating a project, `--launch-jobs`,
   `--test-job`, non-interactive `--variation-phase1`, and `--variation-phase2`.
2. SL7 container, for anything that needs ROOT or the compiled binaries
   (merges with `hadd`, interactive phase 1 = `run_systematics`, stage 5
   `add_numi_flux_sys.py`, stage 6 `to_gundam`, rebuilding):
   - enter with `sh /exp/sbnd/data/users/vito/podman/start_SL7dev_jsl.sh`
     (apptainer, `fnal-dev-sl7:jsl`; it binds `/exp`, `/pnfs`, `/cvmfs`, `/nashome`;
     on exit it deletes `/run/user/$UID/bt_u$UID`)
   - then `source /nashome/m/micarrig/scripts/startSpine.sh` (sources the ICARUS
     CVMFS setup on icarus hosts, then `setup sbnana v10_01_02_01 -q e26:prof` and
     `setup cmake v3_27_4`)
   - then `source /exp/icarus/app/users/micarrig/nuESpine/medulla/build/setup_medulla.sh`
     (sets `MEDULLA_*_DIR`, puts `build/selection` and `build/systematics` on
     `PATH`, adds `build/{selection,systematics,shared}` to `LD_LIBRARY_PATH`)
   - Do not use `~/scripts/setup_icarus_env.sh` for this: it also runs
     `htgettoken` and sets `set -euo pipefail`.
   - From a tool call use the wrapper `scripts/in_sl7.sh '<one command>'` (run it
     from the repo root). It does the three steps above non-interactively, in the
     directory you call it from, with the same image and binds. It does not run
     `htgettoken`, does not source `~/.bash_profile`, has no jobsub, and skips the
     start script's token cleanup on exit. Quote the whole command; chain with
     `&&` inside the quotes. Verified: SL7.9, Python 3.9.15 with `toml`, `uproot`,
     `pandas`, `tqdm`, `ROOT`, `numpy`, `awkward`; `hadd`, `root`, `run_systematics`
     and `to_gundam` on PATH; no missing libraries for `run_systematics`.
   - Never use it for `jobsub_submit`/`jobsub_q` (host only).
3. Plotting: the pixi environment from `spineplot/pyproject.toml` (python >=3.10,
   numpy, scipy, pandas, matplotlib, uproot, toml, scikit-learn, root,
   memory_profiler). The env is at `spineplot/.pixi/envs/default`; the `pixi`
   binary is not on PATH. `/usr/bin/python3` also has the needed packages, but
   prefer the pixi env. Plots are produced by the owner on a bigger-memory
   server: do NOT run the real plotting here (this node has about 11 GB of RAM);
   `--dry-run` only.

The Mac copy at `~/Desktop/medulla` is for editing TOMLs and code and running the
Python unit tests; it cannot run `jobsub_submit`, `ifdh`, or the C++ binaries.

## The analysis config

- Analysis dir: `selection/toml/icarus_nueCC_inclusive_2026/`
  - `meta.toml`: experiment `icarus`, one role `primary` -> `selection.toml`,
    `batch_size = 20`. It lists `enable` keys (`icarus_nominal`, `icarus_cvext`,
    `icarus_offbeam`) but the selection file uses an inline `[[sample]]` instead.
  - `selection.toml` and `selection_nue.toml` are near-duplicates (1426/1137
    lines; line counts drift as the owner edits live, re-`wc -l` before trusting
    either number). As of this writing the only difference besides the samples
    below is the `NuMI_nue` sample path: `selection.toml` lists 10 explicit
    `root://fndcadoor...` files (a subset), while `selection_nue.toml` uses the
    glob `/pnfs/icarus/persistent/users/dcarber/spine/NuMI_Nue_CV_flat_cafs/*.flat.root`.
    `meta.toml` points at `selection.toml`, so a campaign built from it runs the
    10-file subset. Production selection is `selection.toml` (owner confirmed); `selection_nue.toml` is only for the nue-enhanced variation samples.
  - `selection_g4.toml`: same `[parameters]`/`[general]`/`[[tree]]`/`[[category]]`
    as `selection.toml`, but with the `[[sample]]` list trimmed to just the
    `g4` sample (tag `g4`, path under
    `/pnfs/icarus/persistent/users/kyjung/G4Reweight/NuMI/...`). It exists
    because the g4 sample's `rec.mc.nu.wgt` layout only matches the nominal
    production's up to index ~114; indices 115-120 hold the six G4
    reinteraction weight groups there instead of the Flux groups nominal has at
    those slots, so it cannot share a project or a systematics TOML with
    `selection.toml`. Used only with `NuMI_g4_gundam.toml` (see below).
- Selection TOML layout: `[parameters]` (thresholds, referenced as `"@name"`),
  `[general]`, `[[sample]]`, `[[tree]]` (6: event/reco/true modes with `[[cut]]`
  and `[[branch]]`, including `nc_sideband`, a reco-mode NC control region), `[[category]]`
  (6 signal/background categories). Samples carry a `tag` (`nominal`,
  `nue_enhanced`, `data`, `detector_variation`, `g4`) that later steps use to
  pick output filenames. Reference: tutorial.md, "The medulla Selection
  Configuration File".
  - `nc_sideband`'s `no_electrons` cut uses parameter `@electron_threshold_sb`.
    It used to be misspelled `@electron_sb_threshold` in `selection.toml`
    (the parameter `electron_threshold_sb` is defined but `electron_sb_threshold`
    is not), which would have thrown a `ConfigurationError` the first time
    anything ran the `nc_sideband` tree. Fixed in `selection.toml`,
    `selection_nue.toml`, and `selection_g4.toml` as of 2026-10-19.
- Systematics TOMLs: `systematics/toml/NuMI_nueEnhanced_gundam.toml`,
  `NuMI_nue_gundam.toml`, and `NuMI_g4_gundam.toml` (GUNDAM-oriented, `[[tree]]`
  blocks with `action = 'add_weights'`, `[variations]` keys such as
  `induction_high`, `wire_gap`, `recombination`). They contain hard-coded input
  and output paths under `/exp/icarus/...`; check them before reuse.
  `NuMI_g4_gundam.toml` only has `[[tree]]` blocks for `events/g4/...` and only
  six `[[sys]]` entries (`reinteractions_{kminus,kplus,neutron,piminus,piplus,
  proton}_Geant4`, all `type = 'multisim'`, indices 115-120); every tree in it
  has `gundam_store_syst = false` because stage 6's TGraph conversion only
  reads `_multisigmaTree`/`_variationTree`/`_NuMIfluxsimTree`, none of which
  this all-multisim file produces. Do not add GENIE/flux `[[sys]]` entries to
  it, and do not add the G4 entries to `NuMI_nue_gundam.toml` -- same indices,
  different physics in the two productions' `rec.mc.nu.wgt`.
- `batch/sys_template.toml` is the default per-job systematics config used when
  a project is created (`add_systematics = true` on a tree turns it on).

## The workflow (grid stages via `batch/medulla.py`)

A full analysis needs THREE projects. The main and nue-enhanced ones exist
because the detector-systematic splines come from two different samples; the
g4 one is a dedicated single-sample project because its `rec.mc.nu.wgt` index
layout collides with the nominal production's (see "The analysis config"
above) and so cannot share a project, a per-job systematics template, or a
`to_gundam` run with the other two:

| Project | Selection TOML | Systematics TOML | Provides |
|---------|----------------|------------------|----------|
| main `<proj>` | `selection.toml` | `NuMI_nue_gundam.toml` | nominal splines; the final merged output |
| nue-enhanced `<proj_nue>` | `selection_nue.toml` | `NuMI_nueEnhanced_gundam.toml` | nue splines |
| g4 `<proj_g4>` | `selection_g4.toml` | `NuMI_g4_gundam.toml` | G4 reinteraction systematics for the `g4` sample only; combined into the main project's gundam output by hand at the end (stage G) |

In phase 2 (main project only), every true electron-neutrino event uses the
nue-enhanced splines and everything else uses the nominal splines. Phase 2
therefore needs both `<proj>/variation_splines.root` and the nue splines, so
`<proj_nue>` must have finished phase 1 before the main project's phase 2.
The g4 project is independent of all of this: it carries no detector-variation
samples, so it never touches phase 1/2 or stage 5, and nothing in the main or
nue-enhanced projects depends on it finishing first (or at all -- it can run
whenever, and stage G can wait for it).
Selection TOML paths are under `selection/toml/icarus_nueCC_inclusive_2026/`,
systematics TOMLs under `systematics/toml/`.

Run from the repo root on the GPVM. Each Bash call is a fresh shell, so
environment setup does not carry over; see "Where things run" for which stages
run on the host (grid submission, outside the container) and which in the SL7
container (merges A and B, interactive phase 1, stages 5, 6 and G). Project dirs are on
PNFS scratch. Merged files and logs live in the project directory (logs
in `<proj>/logs/`). Capture output:
`... 2>&1 < /dev/null | tee -a <proj>/logs/<stage>.log` (`< /dev/null` makes an
unexpected prompt fail with EOFError instead of hanging).

`medulla.py --dry-run` (works for every stage flag) prints what it would create,
mark `completed`, delete, copy into the project dir, and submit, including the
full `jobsub_submit` command, and does none of it; `project.db` and the project
dir are untouched. It beats `--yes` if both are given. With `-Z` it still opens
the output ROOT files to find corrupt ones. Create (`-c`) stops after printing.

Standard inputs (use these unless the owner says otherwise):
- Git ref (`--tag`): `feature/carriganm95_working`. It is a moving branch, so
  resolve it to a SHA with `git ls-remote origin feature/carriganm95_working`,
  tell the owner the SHA, and confirm it is pushed and equals local HEAD.
- Batch size `-b 30` (files per job), `--memory 10000`, `--disk 10` (GB) for the
  job-submitting stages. Do not pass `--lifetime`: the owner is fine with the
  default. Note that default is `medulla.py`'s own `1h`, which it always sends as
  `--expected-lifetime`; it is not jobsub's default. If jobs die at the lifetime
  limit, tell the owner rather than guessing a larger value.
- Project dirs: `/pnfs/icarus/scratch/users/micarrig/nueCCInclusiveSpine<X>` (main),
  `.../nueCCInclusiveSpine_nue_<X>` (nue-enhanced), and `.../nueCCInclusiveSpine_g4_<X>`
  (g4), `<X>` = version number. The next version is 31, so `nueCCInclusiveSpine31`,
  `nueCCInclusiveSpine_nue_31`, and `nueCCInclusiveSpine_g4_31`.
  Version 30 is the previous one; never reuse or modify it. Always create a new
  project dir (`-c` refuses an existing `project.db`).
- The right systematics TOML must be passed at EVERY stage that takes one, or the
  wrong systematics are produced. The flag differs by stage: `-s/--systematic`
  at project creation (stage 1; the per-job systematics template, frozen into the
  project) and `-V/--variation-systematics` for phases 1 and 2 (main/nue-enhanced
  only; the g4 project never passes `-V`, it has no phase 1/2). `--systematic` is
  ignored after creation, so stage 2 does not take it. If a project was created
  without `-s` it used `batch/sys_template.toml`; stop and tell the owner.
- Never mix the three sets of TOMLs (main/nue/g4) within one project, and never
  point the g4 project's `-s`/`to_gundam` config at a non-`g4`-tagged sample or
  vice versa -- see the index-collision note in "The analysis config".

Token: for now, ASK the owner before running `htgettoken`, at the start and
after every stage ("should I run htgettoken?"); never run it unprompted. The
command, for when the owner says yes:
`export BEARER_TOKEN_FILE=/tmp/bt_u$(id -u) && htgettoken -a htvaultprod.fnal.gov -i icarus`
(it may need a browser, rule 4). Sourcing `setup_icarus_env.sh` also runs it, so
that counts too. Set `BEARER_TOKEN_FILE` in every command that submits or runs
`jobsub_q`, since the shell does not persist.

```
main:  1 create -> 2 selection -> [A merge] -> 3 phase 1 -----------\
nue:   1 create -> 2 selection -> [A merge] -> 3 phase 1 -> nue splines -> 4 phase 2 (main)
                                                                            -> [B merge] -> 5 flux -> 6 gundam -\
g4:    1 create -> 2 selection -------------------------------------------[B merge] -> 6 gundam (g4 TOML) -----+-> G merge -> merged_gundam_final.root
```

Stages 1 to 3 run once per project (they can run side by side if the owner
approves both). Stages 4 to 6 run on the main project. The nue-enhanced project
needs nothing after phase 1 (no phase 2, merge B, or stages 5 and 6). The g4
project needs only 1, 2, its own merge B (now that `merge_list.py` knows the
`g4` tag), and its own stage 6 with `NuMI_g4_gundam.toml` -- no phase 1/2 (it
has no detector-variation samples) and no stage 5 (its tree blocks are all
`gundam_store_syst = false`, so there is no `_NuMIfluxsimTree` for flux
systematics to populate). Stage G (new, below) combines the main project's
and the g4 project's gundam outputs; it can run whenever both are ready and
does not block or get blocked by phase 2/stage 4-5 on the main project.

`<S>` = `selection/toml/icarus_nueCC_inclusive_2026`, `<T>` = `systematics/toml`
(`<T>/NuMI_nue_gundam.toml` for main, `<T>/NuMI_nueEnhanced_gundam.toml` for nue,
`<T>/NuMI_g4_gundam.toml` for g4).
`<G>` = `--tag <SHA or branch> --memory 10000 --disk 10`.

| # | Stage | Command | Done when |
|---|-------|---------|-----------|
| 1 | Create project | `python3 batch/medulla.py -p <proj> -e icarus -c -t <S>/selection.toml -s <T>/NuMI_nue_gundam.toml -b 30` (nue: `selection_nue.toml`, `NuMI_nueEnhanced_gundam.toml`; g4: `selection_g4.toml`, `NuMI_g4_gundam.toml`) | `<proj>/project.db` exists; nothing submitted |
| 2 | Selection jobs | `python3 batch/medulla.py -p <proj> -e icarus -Z -y --launch-jobs [N] <G>` | `jobsub_q --user micarrig` shows nothing for it, and every project.db job is `completed` (see Failures) |
| A | Merge for phase 1 (main/nue only) | `python3 batch/merge_list.py <proj> --tags detector_variation --no-varsys --output <proj>/merge_list_stage0.txt`, then `~/scripts/mergeFiles.sh -d -l <proj>/merge_list_stage0.txt <proj>/merged_stage0.root` | `<proj>/merged_stage0.root` exists; "N of N files are readable" |
| 3 | Phase 1 splines (interactive, no grid; main/nue only) | `python3 batch/medulla.py -p <proj> -e icarus -y --variation-phase1 --variation-interactive -V <T>/NuMI_nue_gundam.toml -I <proj>/merged_stage0.root` (nue project: its own TOML) | `<proj>/variation_splines.root` exists |
| 4 | Phase 2 apply (main only) | `python3 batch/medulla.py -p <proj> -e icarus -Z -y --variation-phase2 -V <T>/NuMI_nue_gundam.toml [--dataset-tag nominal] <G>` (`--test-job` first for 1 job) | one `output/output_varsys_jobid<NNNN>.root` per `output_systematics_jobid<NNNN>.root` of the selected tag |
| B | Final merge | `python3 batch/merge_list.py <proj> --output <proj>/merge_list_final.txt` (all tags plus `output_varsys*.root`), then `~/scripts/mergeFiles.sh -d -l <proj>/merge_list_final.txt <proj>/merged_final.root` | `<proj>/merged_final.root` exists; "N of N files are readable" (for g4, `<proj>` = `<proj_g4>`, and there is no `output_varsys*.root` to include) |
| 5 | Flux systematics + PPFX (main only) | `python3 systematics/add_numi_flux_sys.py <proj>/merged_final.root events/NuMIFull all_signal all_signal_nonmatched selected selected_nonmatched` | log shows each tree processed; `events/NuMIFull/<tree>_NuMIfluxsimTree` and the `ppfx_cv_weight` branch exist |
| 6 | GUNDAM format | `mkdir -p <proj>/gundam_tmp && cd <proj>/gundam_tmp && <repo>/build/systematics/to_gundam <repo>/<T>/NuMI_nue_gundam.toml <proj>/merged_final.root`, then `mv output.root <proj>/merged_gundam.root` (g4: `<T>/NuMI_g4_gundam.toml` and `<proj_g4>/merged_final.root`) | `<proj>/merged_gundam.root` exists; log has a `Processing events/...` line per `[[tree]]` |
| G | Final gundam merge (main + g4) | `hadd <proj>/merged_gundam_final.root <proj>/merged_gundam.root <proj_g4>/merged_gundam.root` | `<proj>/merged_gundam_final.root` exists; `events/NuMIFull/...` and `events/g4/...` are both present with no duplicate `(Run, Subrun, Evt)` within either (check with `checkGundamDuplicates.py`, per directory -- the two inputs populate disjoint `events/` subdirectories so plain `hadd` is safe, there is nothing to deduplicate across them) |

Where each stage runs: stages 1, 2 and 4 (and `jobsub_q`) on the host. Merge A,
phase 1 (interactive), merge B, stage 5, stage 6, and stage G inside the
container: wrap the whole command, e.g. `scripts/in_sl7.sh 'python3 batch/medulla.py -p <proj> -e icarus -y --variation-phase1 --variation-interactive -V ... -I ...'`
or `scripts/in_sl7.sh 'python3 batch/merge_list.py ... && ~/scripts/mergeFiles.sh -d -l ... ...'`.
The g4 project only ever uses stages 1, 2, B and 6, so it only ever needs the
host (1, 2) and the container (B, 6); it never touches phase 1/2 or stage 5.

Phase 1 runs on the node itself (`--variation-interactive`), not on the grid:
it takes no `--tag`/`--memory`/`--disk`, runs `build/systematics/run_systematics`
(needs the SL7 container + medulla setup), and copies its TOML to the hard-coded
`/exp/icarus/app/users/micarrig/nuESpine/variation_systematics_phase1.toml`.
It is quick, so run it in the background and notify when it exits.

Plots (`analysis/makeSelectionPlots.py`, `makeSystematicPlots.py`,
`plot_gundam_graphs.py`) are separate from this sequence.

Notes:
- The script is `batch/merge_list.py` (underscore). `--tags T [T ...]` keeps only
  jobs with those sample tags and `--no-varsys` drops the varsys files. Merge A
  uses `--tags detector_variation`, which covers the variation samples and the
  `cv` sample (all tagged `detector_variation`). Without flags it picks
  `output_jobid*.root` for tags `data`/`detector_variation`/`rock` and
  `output_systematics_jobid*.root` for `nominal`/`nue_enhanced`/`g4` (g4 is MC
  with `add_weights` trees, same as nominal, so it needs the systematics-merged
  file too), plus any `output_varsys*.root`. Missing files are only a `[WARN]`
  line; read them.
- `~/scripts/mergeFiles.sh` (= `/nashome/m/micarrig/scripts/mergeFiles.sh`):
  `-l/--file-list` merges exactly the listed files, `-d` opens each with ROOT to
  drop corrupt ones. It never fails on bad files, it skips them and prints
  `[WARN] Skipping N ...` and `[INFO] X of Y files are readable`. If X < Y, tell
  the owner before treating the merge as good. It writes `hadd` temp dirs next to
  the output, so the output directory must already exist.
- Each selection job runs `medulla` then `run_systematics` (`batch/submit.sh`),
  producing `output_jobid<NNNN>.root` and `output_systematics_jobid<NNNN>.root`.
- `-Z` adds the PyROOT zombie/recovered check on top of the 1 KB size check.
  Stage 2 job completion in `project.db` is only updated when `medulla.py` runs
  and sees the output files, it is not updated by the grid.

- Phase 2 reads the nue splines from `[variations] splines_file_nue` in the
  `-V` TOML and copies that file into the jobs. `NuMI_nue_gundam.toml` currently
  points it at an old project (`/pnfs/icarus/scratch/users/micarrig/nueCCInclusiveSpine30/`).
  Before stage 4 it must be `<proj_nue>/variation_splines.root` (for version 31,
  `/pnfs/icarus/scratch/users/micarrig/nueCCInclusiveSpine_nue_31/variation_splines.root`); show the owner
  the line and the proposed change and wait before editing a shared TOML. The
  main splines come from `<proj>/variation_splines.root`, not from `splines_file`
  (that key is overwritten). `[input]`/`[output]` paths in these TOMLs are
  hard-coded to `/exp/icarus/...` and are also rewritten per job by phases 1 and 2.
- Stage 5 (`add_numi_flux_sys.py`): needs `uproot`, `pandas`, `tqdm` and PyROOT
  (setup script), and reads `/nashome/m/micarrig/icarus/nuESpine/flux.root`
  (hard-coded, horn current `fhc`). It opens the file in `UPDATE` mode, so it
  edits `merged_final.root` in place; if it dies midway, redo merge B. Always
  pass the tree names: with none it processes every TTree in the directory, which
  includes the `*_multisigmaTree`/`*_variationTree`/`*_NuMIfluxsimTree` trees and
  will fail or write junk. It loops over events in Python, so it is slow: run it
  in the background with a log. It reads `Run`, `Subrun`, `Evt`, `true_pdg`,
  `true_parent_pdg`, `true_neutrino_energy` from each tree.
- Stage 6 (`to_gundam`): the argument order is `<config.toml> [input.root]`; if
  the input is omitted it uses `output.path` from the TOML (a debug file), so
  always pass it. It writes `output.root` in the CURRENT directory with
  `RECREATE`, hence the empty `gundam_tmp` directory. For a `[[tree]]` block
  with `gundam_store_syst = true` it needs that tree's `*_multisigmaTree`,
  `*_variationTree` and `*_NuMIfluxsimTree` (the last one is what stage 5
  adds), so stage 6 before stage 5 gives a broken file for those trees -- this
  is why `NuMI_nue_gundam.toml`'s trees need stage 5 first. `NuMI_g4_gundam.toml`
  sets `gundam_store_syst = false` on every tree (it only has `multisim`-type
  `[[sys]]` entries, which this conversion pass does not touch), so the g4
  project's stage 6 does not need stage 5 and can run right after its merge B.
  The binary is `build/systematics/to_gundam`; rebuild with `make` in `build/`
  after any C++ change (ask first).
- Stage G (final gundam merge): plain `hadd` of `<proj>/merged_gundam.root` and
  `<proj_g4>/merged_gundam.root`. Safe as a flat `hadd` (no need for
  `merge_list.py`/`mergeFiles.sh`'s per-file corruption skipping, since there
  are only two whole, already-verified inputs) because each `to_gundam` run
  only ever writes the destinations in its own `[[tree]]` list --
  `NuMI_nue_gundam.toml` writes `events/{NuMIFull,data,dataBeamOff,rock}/`,
  `NuMI_g4_gundam.toml` writes only `events/g4/` -- so the two files have no
  overlapping tree paths to deduplicate. Run `checkGundamDuplicates.py` on the
  result anyway (per the owner's `events/g4/selected` and `events/NuMIFull/selected`
  trees) as a sanity check before calling this stage done.

## Failures and relaunching

Stages 2 and 4 both detect failed jobs the same way: on every `medulla.py` run
the stub/zombie/missing outputs are listed, and (with `-y`) the bad files are
deleted and the jobs are resubmitted. Do not run `-y` to just look.
- Stage 2 (`--launch-jobs`): stubs/zombies are reset to `pending` and relaunched
  together with any other pending jobs. Jobs that are tracked but have no output
  at all (not `pending`) are only picked up with `-r/--relaunch-missing`, and
  that can duplicate jobs still running; check the queue first.
- Stage 4 (`--variation-phase2`): it skips jobs that already have a good
  `output_varsys_*` file, deletes stub/zombie ones (with `-Z -y`), and
  resubmits only the rest. Re-running the same command is the relaunch.
- To inspect without changing anything, run without `-y` and with
  `< /dev/null`: `python3 batch/medulla.py -p <proj> -e icarus -Z < /dev/null`.
  This still updates `project.db` statuses to `completed` for good outputs. If
  it hits a delete prompt it dies on EOF, which is the signal that failed files
  exist; report the job IDs printed above the prompt and do not rerun with `-y`
  until the owner approves.
- Is the queue empty? `jobsub_q --user micarrig` (needs a valid token, see
  above). A stage is only "finished" when the queue shows no jobs
  for it AND the status check has run.

## Stage protocol: notify, then ask (required)

The owner decides a stage at a time whether to get phone notifications.
When asking to start a plotting run or a grid submission (any job-submitting
stage), also ask in the same question: "Push notifications to your phone when
it finishes, so you can answer from there?" If yes, check that Remote Control is
connected (the owner runs `/remote-control`) and send a `PushNotification` at
step 3 below; if no, report in the chat only. The answer covers that run only;
ask again next time. Either way, nothing moves on without the owner's decision.
After submitting a stage (once approved):
1. Start a background watch (Bash `run_in_background` with an `until` loop on
   the queue, or `Monitor`) so you are woken when the jobs leave the queue. For
   the local stages (merges A and B, stages 5 and 6) there is no queue: run the
   command itself with `run_in_background`, logging to `<proj>/logs/`, so you are
   woken when it exits. Do
   not poll in the foreground and do not sleep-loop in a turn.
2. When woken, run the read-only status check above, and ask the owner whether
   to refresh the token before the next submit. For merges, check the
   "X of Y readable" line and that the output file exists.
3. If the owner opted in, send `PushNotification` (under 200 characters, lead
   with the result and the question); otherwise put the same text in the chat.
   For example: `Stage 2 done: 1180/1200 ok, 20 failed. Relaunch the 20, or go on to merge?`
   or `Merge A done, 1200/1200 readable. Start phase 1?`
4. Stop and wait for the owner's reply. Do not start the next stage, relaunch,
   or run `-y` before they answer.

Notification caveats: `PushNotification` reaches the phone only if Remote
Control is connected for this session, and it is skipped when the owner is
actively at the terminal. It cannot receive a reply; the owner answers in the
session. If the tool reports "not sent", say so in the chat and wait.
Notify on both success and failure. A silent failure is worse than a noisy one.

## Things that will bite you

- Project job configs are frozen into `project.db` when the project is created.
  Editing a selection TOML afterwards does not change existing projects; create a
  new project or campaign.
- Grid workers build medulla from GitHub at `--tag` for every job (`git clone`,
  `cmake`, `make -j4`). The C++ they run is whatever is pushed at that ref.
- `project.db` and `campaign.db` live on `/pnfs`, which breaks SQLite locking.
  The code copies them to a local temp file, edits, and copies back. Do not open
  them in place with `sqlite3` or write to them from two shells at once.
- A `--dry-run` of phase 2 on a big project is slow (about 3 minutes for v30: it
  copies and parses the 354 MB `project.db`), and every `medulla.py` run on an
  existing project first runs the status check (and copies `project.db` to the cwd,
  so run from a scratch directory with space).
- Launching does not check for jobs already running, so relaunching can create
  duplicates. Look at the queue first, and be careful with `-r`.
- `merge_list.py` only knows the tags `data`, `detector_variation`, `rock`,
  `nominal`, `nue_enhanced`, `g4`; other tags are skipped with a warning.
- `.gitignore` ignores `*.json`, `*.root`, `*.png`, and `*.log`. New JSON state or
  plan files will not be tracked unless the ignore rule is changed.

## Plotting (spineplot)

- `spineplot/` is a Python library (pixi env in `spineplot/pyproject.toml`,
  Python >= 3.10). The scripts in `analysis/` configure it through TOML files
  (`analysis/systematics.toml`, `nue_plots.toml`, `efficiency_*.toml`,
  `systematic_recipes.toml`).
- The plotting scripts are not CLIs. Input ROOT paths and the `spineplot` path
  are hard-coded near the top of `__main__`, and some TOMLs include absolute
  `/Users/micarrig/Desktop/medulla/...` paths. Edit those, do not guess.
- Systematic bands can be split by `true_category` (signal vs background) via an
  artist's `groups` key; see `analysis/systematics.toml`.

### Efficiency plots (`analysis/plotEfficiency.py`)

- Not a CLI. `main()` hard-codes the TOML, the input ROOT file
  (`/exp/icarus/data/users/micarrig/merged_v30.root`, 32 GB), the output dir and
  the list of figures; switch TOML by editing/commenting lines. Figure names in
  the script must match `[[figure]] name` in the TOML. `efficiency_vs_trueNu` and
  `efficiency_vs_recoFlashTime` are commented out in the script.
- Environment: the pixi env (see "Where things run"); `/usr/bin/python3` also has
  uproot, ROOT, awkward, matplotlib, memory_profiler, scipy, pandas. No SL7
  container or token is needed. The script puts `spineplot/` on `sys.path`.
- `/nashome/m/micarrig/icarus` is a symlink to `/exp/icarus/app/users/micarrig`,
  so the TOMLs' `/nashome/...` include paths are this same tree.
- TOMLs: `efficiency_analysis.toml` (tree `all`, Signal vs Background),
  `efficiency_recoInt.toml` (tree `true_signal`, 6 neutrino-type categories,
  `final_only` curves), `efficiency_withError.toml` (tree `all_signal` plus its
  `_multisigmaTree`, `_NuMIfluxsimTree`, `_variationTree`, and `ppfx_cv_weight`,
  so it needs the file AFTER stage 5; figures pick the error with `draw_error`).
- Efficiency is binomial counts, unweighted (weights are not used), computed per
  sample and per category; each TOML adds a second `NuMIFullNonSel` sample
  (true-signal events with no reco match) to the same ordinate.
- Each figure is saved as a `.png` and a `.pkl`.
- New: `analysis/plotEfficiencyCLI.py --toml T --input F --output-dir D [--dry-run]`
  runs the 8 main figures; plus `efficiency_analysis_allTree.toml` (no
  `NuMIFullNonSel`) and `efficiency_analysis_trueSignalTree.toml`. Output goes to
  `plots_v30/efficiency_<tree>/`. Do not run it for real here: the owner runs the
  plots on a bigger-memory server. `--dry-run` is fine.

## Tests

- After changing anything in `batch/`: `python3 -m pytest batch/ -v` (no build or
  credentials needed).
  `batch/test_dry_run.py` checks that `--dry-run` changes no file (project.db,
  project dir, cwd), never calls `jobsub_submit`, and beats `--yes`; each test has a
  control that does mutate. Rerun it after touching any prompt, delete, or submit path.
- C++ framework validation (needs a build): `MEDULLA_BUILD_DIR=<build>/selection
  python3 -m pytest selection/test/ -v`, or `make pytest` from the build dir.
- Skipped without `MEDULLA_BUILD_DIR`, by design.

## Git

- Working branch: `feature/carriganm95_working`. Remote `origin` is the shared
  upstream `justinjmueller/medulla`; ask before pushing.
- Untracked or modified local files that are not part of any change you were
  asked to make (`.history/`, `plots/`, `test_1file.df`, `output_nueCCNp_mc_5e18.root`,
  notebooks) should be left alone.

## Open items for the owner

- The GUNDAM fit commands themselves (after stage 6).

## Where we left off (2026-10-09)

Read this first in a new session. The previous session was started from
`/exp/icarus/app/users/micarrig`, so it cannot be resumed from this directory.

Built and tested (all `python3 -m pytest batch -q`, 135 tests pass):
- `medulla.py --yes` and `--dry-run` (tests in `batch/test_dry_run.py`); `merge_list.py --tags/--no-varsys`.
- `scripts/in_sl7.sh` (run one command in the SL7 container with the medulla environment).
- `scripts/guard_bash.py` + `.claude/settings.local.json` (hook, allow/ask/deny rules) and the
  `/medulla-stage` skill. The hook was written but had NOT fired yet in the old session (the
  `.claude/` dir did not exist at its start). First thing to do here: type `/hooks`, look for
  PreToolUse > Bash > `guard_bash.py`, then run `echo pkill -f nothing-here`; it must be denied.
- Plotting: `analysis/efficiency_analysis_allTree.toml` (no `NuMIFullNonSel`),
  `analysis/efficiency_analysis_trueSignalTree.toml` (`true_signal` + `true_signal_nonmatched`),
  and `analysis/plotEfficiencyCLI.py --dry-run`. Both dry-runs work. The REAL plotting is run by the
  owner on a bigger-memory server, not here (two earlier attempts here were killed; the node has ~11 GB).
  Leftovers I created: `plots_v30/` (empty dirs and two `.log` files) in `nuESpine/`.
- Dry-runs against the existing v30 project worked (launch, phase 1, phase 2). v30 has one stub
  pair (job 2854: 436 and 423 byte outputs) and one pending job; phase 2 is 2472/2472 complete.
  Creation dry-runs for `nueCCInclusiveSpine31` and `nueCCInclusiveSpine_nue31` worked; neither exists yet.

Not verified / open:
- A container dry-run of `merge_list.py` on v30 was interrupted when the session ended, and the
  node may have frozen around it (cause unknown). The owner was going to reproduce it by hand.
  `merge_list.py` copies the 354 MB `project.db` and stats every output on PNFS; it has no `--dry-run`.
  A leftover `project_merge_tmp.db` may exist in an old scratch directory (do not delete it unasked).
- `NuMI_nue_gundam.toml` `splines_file_nue` still points at the v30 path; it must be the nue31
  project's `variation_splines.root` before stage 4.
- `--lifetime`: `medulla.py` always sends `--expected-lifetime=1h` (its default, not jobsub's).
- `PushNotification` reaching the phone (needs Remote Control) has not been tested.

Next steps, in order:
1. Verify the guard hook is loaded (above).
2. Try `/medulla-stage 1 nueCCInclusiveSpine31` (it should stop at the dry-run and ask).
3. Write the plotting runbook for the bigger server: pixi env at `spineplot/.pixi/envs/default`,
   input `merged_v30.root`, the two TOMLs, output `plots_v30/efficiency_<tree>/`, and a leaner
   loading step. Later: error-bar version (`efficiency_withError.toml`, needs the file after
   stage 5) and the neutrino-type breakdown (`efficiency_recoInt.toml`).
4. Add `--dry-run` to `merge_list.py`; a preflight script (branch pushed, node memory); GUNDAM fit
   commands are still undocumented.
Leave unrelated modified or untracked files alone (see Git).
