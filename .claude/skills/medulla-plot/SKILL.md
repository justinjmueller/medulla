---
name: medulla-plot
description: Make the spineplot plots (selection, systematic, n-1, efficiency, gundam) via the plot_all.py driver and its manifest - validate with --dry-run, then run for real wherever there is enough memory (the built-in gate refuses on the ~11 GB GPVM and runs on a big-memory server). Use when asked to make, prep, re-run, or check plots, e.g. "make all the plots", "make the selection plots", "run plot_all", "prep the nue plots".
---

# medulla-plot

Make the plots with one driver, after validating. Then stop and ask.

The driver is `analysis/plot_all.py`. It reads a manifest (`analysis/plot_jobs.toml`)
of plot jobs and runs each family's CLI as a separate subprocess, one at a time,
so memory is released between families. It has a **built-in memory gate**: a real
run of a job happens only when this node's available memory comfortably exceeds
the input file size. So the SAME command is safe everywhere - on the ~11 GB
ICARUS GPVM it refuses every job and tells you to use the big-memory server; on a
big node it actually makes the plots. `--dry-run` opens no ROOT and writes
nothing. **Never bypass the gate** (do not lower `--mem-factor` to force a run on
a node that cannot hold the input - earlier real runs on the GPVM were killed,
CLAUDE.md rule 10).

Environment and plotting facts live in `CLAUDE.md` ("Plotting (spineplot)" and the
efficiency subsection). Read that first; do not rely on memory of it.

The manifest (the variant matrix lives here, not in code): each `[[job]]` is one
CLI invocation (name, cli, toml, input, output_dir, optional figures, optional
`args` passed verbatim). To plot the same TOML another way (a different figure
subset, a different tree/grouping TOML, or a different input), add another
`[[job]]` - no code change. The systematic breakdowns (xsec / flux / detector /
total / xsec_grouped) are separate jobs that pass `args = ["--systematics-set",
"<name>"]`; the named sets live in `analysis/systematic_sets.toml`, whose
membership was lifted from `systematics.toml` and is physics config - show it to
the owner before a real systematics run (the `total` set uses aggregate names, not
individual knobs; confirm they resolve). Families and CLIs:

| Family | CLI | Default TOML / input |
|--------|-----|----------------------|
| selection | `analysis/makeSelectionPlotsCLI.py` | `nue_plots.toml` / `merged_v30.root` |
| systematic | `analysis/makeSystematicPlotsCLI.py` | `systematics.toml` / `merged_v30.root` |
| n-1 | `analysis/makeNMinus1PlotsCLI.py` | `nminus1_plots.toml` / `merged_v30.root` |
| efficiency | `analysis/plotEfficiencyCLI.py` | `efficiency_*_*.toml` / `merged_v30.root` |
| gundam | `analysis/plot_gundam_graphs.py` | a `to_gundam` output `.root` |

All five share `--toml/--input/--output-dir/--figures/--dry-run`. A plot variant
is selected by (TOML, figure-name subset); there is no `--tree/--groups` flag
because `run_interactively` takes a figure name and the TOML defines each figure's
tree and category `groups`. The hardcoded `make*Plots.py` originals are left in
place (CLAUDE.md rule 9); the `*CLI.py` files are their argument-driven versions.

## 1. Scope (read-only)

What to make: all jobs in the manifest, or a subset (`--only <name> ...`). If a
specific TOML/figures/input is asked for and no job matches, say so - add a
`[[job]]` (show the diff, it is a new data row) or pass a one-off CLI command;
do not guess paths.

## 2. Validate here (cheap, read-only - no data loaded)

`python3 analysis/plot_all.py --dry-run` (add `--only ...` to narrow). It prints,
per job: TOML, input + size, samples/branches, output dir, every figure it would
write, and overwrite warnings; it opens no ROOT file and creates nothing. Read the
output: a missing input is a `[SKIP] ... input missing` line, not an error. Fix
any bad TOML/figure/input before a real run.

## 3. Memory gate (the driver enforces it)

`free -g`. The gate is automatic, but check first so you can tell the owner
whether a real run is even possible here. If free/available memory is well under
the input size (always true here for a 30 GB input), a real run will `[SKIP] ...
skipped-memory` every job - that is the signal to run on the big-memory server
instead, not to force it.

## 4. Run

Only after the owner approves (plotting is heavy; the guard asks for a non
`--dry-run` plot command). In the same question, ask whether they want a push
notification to their phone when it finishes (needs Remote Control,
`/remote-control`); the choice covers this run only. Environment for the real run: the pixi env at
`spineplot/.pixi/envs/default` (the `pixi` binary is not on PATH; `/usr/bin/python3`
also has the packages). No SL7 container, no token. Run the driver in the
background with a log, one process:
`python3 analysis/plot_all.py [--only ...] 2>&1 < /dev/null | tee -a <logdir>/plot_all.log`
(use `run_in_background`; never `nohup`, never a trailing `&`). The driver already
serializes the families, so do not launch several at once.

## 5. Verify (read-only)

Read the `[SUMMARY]` block: how many ran / skipped / failed. For each job that
ran, the output dir should hold one `.png` + one `.pkl` per figure (efficiency and
gundam write their own names/PDF). Report any `failed` job with the `[ERROR]` line
the CLI printed.

## 6. Notify, then stop

If the owner opted in at step 4, `PushNotification`; otherwise the same text in
the chat (under 200 chars, result first), e.g. `Plots done: 5/5 jobs,
312 pngs. Review, or re-run a family?`. If the gate skipped everything here, say
so and give the exact command to run on the big-memory server. Then stop.

## Stop and report instead of continuing when

a TOML does not parse, an input is missing for a job the owner expected, free
memory is low and the owner asks to force a run here, the manifest has a bad/
duplicate job, or a CLI is not one named in `CLAUDE.md`/the manifest.
