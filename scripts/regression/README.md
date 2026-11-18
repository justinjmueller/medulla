# Regression sandbox

Runs the processing chain locally on a few pinned input files, so two code
versions can be compared output by output. Use it before and after any change
to the selection, systematics, batch config logic or the stage scripts: the diff
should contain only what the change was meant to do.

Everything lives under `/exp/icarus/data/users/micarrig/regression/` (`--root`):

| Path | What |
|---|---|
| `inputs/` | pinned CAF files (copied once by `prepare`), the production spline files, `manifest.toml` |
| `<sha>[-dirty][-suffix]/<tier>/` | one run: a dir per sample, `phase1/`, `phase2/`, `mergeB*/`, `stage6*_tmp/`, `stage7/`, `logs/`, `timings.toml` |
| `src_<sha>/`, `build_<sha>/` | git worktree and build of a `--ref` (built once, reused) |

## Use

```bash
# once: copy the inputs listed in scripts/regression/inputs.toml (~1.5 GB)
scripts/in_sl7.sh 'python3 scripts/regression/run_chain.py prepare'

# reference: a commit (built in its own worktree; never touches build/)
scripts/in_sl7.sh 'python3 scripts/regression/run_chain.py run --ref 6cb6044 --tier full'
# your change: this checkout and its build/ (rebuild build/ first)
scripts/in_sl7.sh 'python3 scripts/regression/run_chain.py run --ref WORKTREE --tier full'

# diff (host python is enough); exit 0 = identical
python3 scripts/regression/compare.py \
    /exp/icarus/data/users/micarrig/regression/6cb6044cf1ae/full \
    /exp/icarus/data/users/micarrig/regression/<new>-dirty/full
```

`run` refuses an existing output dir (`--force` replaces it, `--suffix rep2`
makes a second one). By default both refs run with THIS checkout's TOMLs, so
only the code differs; `--configs-from-ref` uses the ref's own TOMLs.
`compare.py --rtol 1e-12` allows rounding differences, `--ignore REGEX` skips
objects or `tree/branch` paths.

## Tiers

- `quick` (~2.5 min): selection + `run_systematics` for NuMIFull (20 files),
  NuMI_nue (10) and data (10).
- `full` (~6 min, 2026-11-17): every input, then
  - phase 1 spline check: `induction_high` + `cv` only (`[variations].keys` and
    `[[sys]]` cut down), seed fixed;
  - phase 2 on the NuMIFull systematics output with the PRODUCTION splines
    (`inputs/splines_{main,nue}.root`, v32), seed fixed; `[[sys]]` blocks the
    production splines lack (today `recombination`) are dropped with a message;
  - merge B (merge_list.py's rule), stage 5, stage 6, the g4 merge/stage 5/stage 6,
    stage 7.

The configs come from the production code: `create_systematics_cfg`,
`make_phase1_cfg`, `make_phase2_cfg` (batch/utilities.py), the same binaries and
`add_numi_flux_sys.py`.

## Coverage, honestly

20 NuMIFull files give ~16k `all` entries but only ~2 `selected`; NuMI_nue gives
~70 `selected`. data/beam-off/rock select ~0-1 events per 10 files, so stage 6's
data trees are nearly empty. Cut logic is well covered (the `all` trees); weights
and the GUNDAM conversion only on a handful of events. Raise `nfiles` in
`inputs.toml` and `prepare` into a new `--root` if that matters.

Phase 1/2 are only reproducible with `[variations] seed` (added 2026-11-17).
Refs older than that draw random universes, so their phase 2 / stage 6 syst
trees differ from run to run: compare those with `--ignore` or use `quick`.

Tests: `python3 -m pytest scripts/regression/ batch/test_phase_cfgs.py`.
