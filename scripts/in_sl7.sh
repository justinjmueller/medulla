#!/usr/bin/env bash
#
# Run ONE command inside the SL7 dev container with the medulla environment,
# non-interactively (for agent tool calls). Same container, image and binds as
# /exp/sbnd/data/users/vito/podman/start_SL7dev_jsl.sh, but:
#   - runs `bash -c` instead of an interactive shell (no --rcfile, so ~/.bash_profile
#     and /etc/profile.d/jobsub_lite.sh are NOT sourced: no jobsub inside);
#   - does NOT run that script's exit cleanup (it removes /run/user/$UID/bt_u$UID);
#   - never runs htgettoken.
# Inside, it sources ~/scripts/startSpine.sh and build/setup_medulla.sh, then runs
# the command from the directory you called this from.
#
# Usage:  scripts/in_sl7.sh '<command>'
#   e.g.  scripts/in_sl7.sh 'which hadd && hadd -h | head -3'
#         scripts/in_sl7.sh 'python3 batch/merge_list.py <proj> --tags detector_variation --no-varsys'
#
# Do not use for jobsub_submit / jobsub_q: run those on the host.
set -uo pipefail

if [[ $# -lt 1 ]]; then
    sed -n '2,19p' "${BASH_SOURCE[0]}" >&2
    exit 2
fi

MEDULLA_DIR="/exp/icarus/app/users/micarrig/nuESpine/medulla"
INNER='
source /nashome/m/micarrig/scripts/startSpine.sh > /dev/null 2>&1 || { echo "[in_sl7] startSpine.sh failed" >&2; exit 97; }
source '"${MEDULLA_DIR}"'/build/setup_medulla.sh > /dev/null || { echo "[in_sl7] setup_medulla.sh failed" >&2; exit 98; }
'"$*"'
'

[[ -d '/scratch' ]] && export APPTAINER_BINDPATH+="${APPTAINER_BINDPATH:+,}/scratch"

exec /cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer exec \
    --pid --ipc \
    -B /etc/hosts,/tmp,/opt,/cvmfs,/exp/,/pnfs/,/nashome \
    -B /etc/profile.d/jobsub_lite.sh,/etc/condor/,/etc/grid-security/ \
    --home ~/:"${HOME}" --pwd "${PWD}" \
    /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:jsl \
    /bin/bash -c "${INNER}"
