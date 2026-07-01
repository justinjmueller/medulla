#!/bin/bash

#######################################################################
# Usage: submit_variation_phase1.sh [--project=PROJECT]
#                                   [--tag=BRANCH]
#                                   [--variation-input=FILE]
#
# Phase 1 runner: builds detector-variation splines from the merged
# variation+CV input file and stages them to
# project_dir/variation_splines.root. [[tree]] blocks are stripped so
# only spline construction runs; tree processing happens in Phase 2.
#
# Arguments:
#   --project=PROJECT        : Path to the project directory (PNFS)
#   --tag=BRANCH          : Medulla git branch (default: develop)
#   --variation-input=FILE   : PNFS or xrootd path to the merged
#                              variation+CV ROOT file
#######################################################################

PROJECT=""
BRANCH="develop"
VARIATION_INPUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --project=*)           PROJECT="${1#*=}";        shift ;;
    --project)             PROJECT="$2";             shift 2 ;;
    --tag=*)            BRANCH="${1#*=}";         shift ;;
    --tag)              BRANCH="$2";              shift 2 ;;
    --variation-input=*)   VARIATION_INPUT="${1#*=}"; shift ;;
    --variation-input)     VARIATION_INPUT="$2";     shift 2 ;;
    *) shift ;;
  esac
done

if [[ -z "$PROJECT" ]]; then
  echo "[ERROR] --project is required" >&2
  exit 1
fi
if [[ -z "$VARIATION_INPUT" ]]; then
  echo "[ERROR] --variation-input is required" >&2
  exit 1
fi

#######################################################################
# Initial setup
#######################################################################

export IFDH_CP_MAXRETRIES=0
export IFDH_WEB_TIMEOUT=100
export CAFANA_DISABLE_SNAPSHOTS=1

source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup sbnana v10_01_02_01 -q e26:prof
setup cmake v3_27_4

echo "[INFO] Project:         $PROJECT"
echo "[INFO] Branch:          $BRANCH"
echo "[INFO] Variation input: $VARIATION_INPUT"

#######################################################################
# Build medulla
#######################################################################

git clone https://github.com/justinjmueller/medulla.git
cd medulla
git checkout "$BRANCH"
mkdir build && cd build
export CC=$(which gcc)
export CXX=$(which g++)
cmake .. -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_C_COMPILER=$CC
make -j4

#######################################################################
# Stage input files
#######################################################################

ifdh cp "$PROJECT/variation_systematics.toml" variation_systematics.toml
echo "[INFO] Copied variation_systematics.toml from project directory"

if [[ "$VARIATION_INPUT" == *root://* ]]; then
  local_input="$VARIATION_INPUT"
  echo "[INFO] Using xrootd input: $VARIATION_INPUT"
else
  local_input="$(basename "$VARIATION_INPUT")"
  echo "[INFO] Staging merged input: $VARIATION_INPUT -> $local_input"
  ifdh cp "$VARIATION_INPUT" "$local_input"
fi

#######################################################################
# Rewrite TOML for Phase 1
#
# Pass 1: strip all [[tree]] blocks (spline building only).
# Pass 2: rewrite [input] path and [output] path.
#######################################################################

awk '
  /^\[\[tree\]\]/ { skip=1; next }
  skip && /^\[\[/ { skip=0 }
  skip && /^\[/ && !/^\[\[/ { skip=0 }
  skip { next }
  { print }
' variation_systematics.toml > tmp_phase1.toml

awk -v in_path="$local_input" '
  /^\[input\]/  { in_input=1; in_output=0; print; next }
  /^\[output\]/ { in_input=0; in_output=1; print; next }
  /^\[/         { in_input=0; in_output=0; print; next }
  in_input && /^[[:space:]]*path[[:space:]]*=/ && !input_done {
    print "path = \047" in_path "\047"
    input_done=1; next
  }
  in_output && /^[[:space:]]*path[[:space:]]*=/ && !output_done {
    print "path = \047variation_splines.root\047"
    output_done=1; next
  }
  { print }
' tmp_phase1.toml > variation_systematics_phase1.toml && rm tmp_phase1.toml

echo "[INFO] Rewrote TOML: [input] path  -> $local_input"
echo "[INFO] Rewrote TOML: [output] path -> variation_splines.root"
echo "[INFO] Stripped [[tree]] blocks (Phase 1: spline building only)"

#######################################################################
# Run systematics (Phase 1: build and write splines)
#######################################################################

echo "[INFO] Starting run_systematics (Phase 1)..."
ls -lrth
./systematics/run_systematics variation_systematics_phase1.toml
echo "[INFO] run_systematics completed"
ls -lrth

#######################################################################
# Stage output
#######################################################################

ifdh cp variation_splines.root "$PROJECT/variation_splines.root"
echo "[INFO] Staged splines to: $PROJECT/variation_splines.root"
