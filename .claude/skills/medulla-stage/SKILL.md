---
name: medulla-stage
description: Run one stage of the medulla grid workflow (create project, selection jobs, merges, phase 1, phase 2, flux systematics, to_gundam) with the owner's approval protocol - preflight, dry-run, ask, run, wait, verify, record, notify. Use when asked to run, relaunch, or check a stage for a project such as nueCCInclusiveSpine31.
---

# medulla-stage

Run ONE stage, then stop and ask. Facts (paths, flags, which environment, how
failures work) live in `CLAUDE.md` ("Where things run", "The workflow", "Failures
and relaunching", "Stage protocol"). Read those sections first and follow them;
do not rely on memory of them. This skill is the order of operations.

Input: a stage name and a project. Stages: `1` create, `2` selection jobs, `A`
merge for phase 1, `3` phase 1, `4` phase 2, `B` final merge, `5` flux systematics
(`add_numi_flux_sys.py`), `6` `to_gundam`. Projects: main
`nueCCInclusiveSpine<X>` or nue-enhanced `nueCCInclusiveSpine_nue<X>` under
`/pnfs/icarus/scratch/users/micarrig`. If the stage or project is missing or
ambiguous, ask once; do not guess.

## 1. Preflight (read-only, nothing is changed)

- Which project and which TOMLs apply (main vs nue-enhanced); never mix them.
- The previous stage is done: `<proj>/project.db`, merged files, `variation_splines.root`
  and, for stage 4, `<proj_nue>/variation_splines.root` and the `splines_file_nue` line in
  the `-V` TOML. Show any mismatch and stop.
- Git: `git ls-remote origin feature/carriganm95_working` against `git rev-parse HEAD`, and
  `git status --short batch/ systematics/ selection/`. For anything that submits, the SHA
  must be pushed. Report it; do not push.
- Node health before anything heavy: `free -g`, `uptime`. Stop if free memory is low.
- Token: ask the owner whether to run `htgettoken`; never run it unprompted.

## 2. Show the plan (dry-run)

- Grid and project stages: run the exact command with `--dry-run` and show the
  owner the summary (job count, tag and SHA, memory/disk, the `jobsub_submit`
  command, any stub/zombie files it would delete).
- Local stages with no dry-run (merges, 5, 6): print the command, the input and
  output paths, input sizes, and the memory/disk risk. Do not run them to "test".
- Stage 4 with `-Z`, merges, stage 5: say they are heavy.

## 3. Ask

One question: the exact command(s), the environment (host or `scripts/in_sl7.sh`),
and the log path. For a stage that submits grid jobs, also ask whether the
owner wants a push notification to their phone when it finishes (needs Remote
Control, `/remote-control`). Wait for a yes. A yes covers this stage only, and so
does the notification choice.

## 4. Run

- Create `<proj>/logs/` if missing; log to `<proj>/logs/<stage>.log` with `tee -a`,
  `< /dev/null`.
- Grid stages run on the host (outside the container); ROOT stages inside it via
  `scripts/in_sl7.sh '<one command>'`.
- Long or heavy commands use `run_in_background` (one process). Never `nohup`,
  never a trailing `&`, never `pkill -f`.
- Do not add `--yes` unless the owner approved deleting the listed stub/zombie files.

## 5. Wait

Start one background watch: for grid stages an `until` loop on
`jobsub_q --user micarrig` showing nothing for the project, for local stages the
process exit. Do not poll in the foreground. If the session may end, say so and
point to the log file and `<proj>/STATE.md`.

## 6. Verify (read-only)

- Grid stages: `python3 batch/medulla.py -p <proj> -e icarus -Z < /dev/null` without
  `-y` (a delete prompt dying on EOF means failed files exist: list the job IDs).
- Merges: "X of Y files are readable" and the output file exists. X < Y is reported.
- Stage 5: each of the four trees processed. Stage 6: a `Processing ...` line per tree.

## 7. Record

Append one entry to `<proj>/STATE.md` (create it if missing): date, stage, git SHA,
command, result counts, log path, and what is next. For stage 1 on a project that
does not exist yet, record it after the project is created.

## 8. Notify, then stop

If the owner opted in at step 3, `PushNotification`; otherwise the same text in
the chat (under 200 characters, result first, then the question), for
example `Stage 2 done: 1180/1200 ok, 20 failed. Relaunch the 20, or merge?`. Also
ask whether to refresh the token before the next submit. Then stop. Do not start
the next stage or a relaunch until the owner answers.

## Stop and report instead of continuing when

an unexpected prompt appears, authentication fails, the SHA is not pushed or does
not match, free memory is low, a command is not named in `CLAUDE.md` (read it and ask),
or the dry-run shows anything the owner has not seen.
