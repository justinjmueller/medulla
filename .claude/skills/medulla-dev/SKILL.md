---
name: medulla-dev
description: Build, run, or test the medulla C++/Python framework locally (the dev loop - compile, run a selection over CAF input, run the test suites). Use when asked to build/rebuild medulla, run a TOML selection, run run_systematics interactively, or run the tests. NOT for grid stages (use medulla-stage) or plotting. Triggers: "build medulla", "rebuild", "run this toml", "run the selection", "make pytest", "run the tests", binary won't link / `sbnana` not found / cmake FATAL_ERROR about SBNANA_VERSION / `setup: command not found` / `perl: ... undefined symbol`.
---

# medulla-dev

The local dev loop: compile medulla, run a selection, run the tests. For the
grid workflow (project create, selection jobs, merges, phases, flux systematics,
to_gundam) use `/medulla-stage` instead.

**Environment facts live in `CLAUDE.md` ("Where things run"). Read that first and
follow it; do not rely on memory.** In short: anything that needs ROOT or the
compiled binaries (`make`, `medulla`, `validate`, `run_systematics`, `to_gundam`,
`hadd`) runs in the **SL7 container**; from a tool call, use the wrapper
`scripts/in_sl7.sh '<one command>'` from the repo root — it sources
`startSpine.sh` + `build/setup_medulla.sh` non-interactively, has no jobsub, and
never runs htgettoken. **Do not** use the container for grid submission (host
EL9 + token) or plotting (`spineplot/` pixi env / host `python3`). Each Bash
call is a fresh shell, so the container+env never persist — one `in_sl7.sh` call
per command.

## Build

`build/` is already configured, so builds are incremental:

```bash
scripts/in_sl7.sh 'make -C build -j4'          # build all targets
scripts/in_sl7.sh 'make -C build medulla'      # one target: medulla, validate,
                                               # run_systematics, to_gundam, shared, ...
```

Clean reconfigure only if the cmake cache is stale or CMakeLists changed:
`scripts/in_sl7.sh 'rm -rf build && mkdir build && cd build && cmake .. && make -j4'`
(deletes/recreates `build/`, so confirm with the owner first per CLAUDE.md).

Binaries: `build/selection/{medulla,validate}`,
`build/systematics/{run_systematics,to_gundam}`.

## Run a selection

```bash
scripts/in_sl7.sh './build/selection/medulla <config.toml>'
```

Output is a ROOT file with one `TDirectory` per sample, trees under
`events/<sample>/<tree>`. Remote (`root://`) inputs need a valid token first —
that is the owner's call (`htgettoken`); local files on disk do not.
`run_systematics` takes one config: `scripts/in_sl7.sh './build/systematics/run_systematics <config.toml>'`.

## Test

```bash
# Batch-layer unit tests — PURE PYTHON, run on the HOST (no container, no build):
python3 -m pytest batch/ -v
# Framework validation — needs the binaries, so run in the container:
scripts/in_sl7.sh 'MEDULLA_BUILD_DIR=$PWD/build/selection python3 -m pytest selection/test/ -v'
# Or the full suite (rebuilds medulla+validate, runs both):
scripts/in_sl7.sh 'make -C build pytest'
```

Framework tests are **silently skipped** if `MEDULLA_BUILD_DIR` is unset — if
`selection/test/` reports everything skipped, that's the cause, not a pass.
Always run the guard tests after touching `scripts/guard_bash.py`:
`python3 -m pytest batch/test_guard_bash.py -q`.

## Gotchas

- Outside the SL7 container → `setup: command not found` / perl symbol errors.
  Use `scripts/in_sl7.sh`.
- Don't rebuild while jobs are preparing without asking: grid jobs clone the C++
  at the pushed `--tag`, but a local rebuild can change what an interactive stage
  (phase 1, merges, to_gundam) runs. See CLAUDE.md.
