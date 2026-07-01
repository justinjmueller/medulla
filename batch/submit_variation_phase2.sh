#!/bin/bash

#######################################################################
# Usage: submit_variation_phase2.sh [--project=PROJECT]
#                                   [--tag=BRANCH]
#
# Phase 2 runner: applies pre-built detector-variation splines to one
# individual selection output file. Uses $PROCESS to look up the job
# ID from the Phase 2 manifest staged in the project directory, then
# stages the corresponding output_jobid<NNNN>.root, applies weights,
# and copies the result to output_varsys_jobid<NNNN>.root.
#
# Arguments:
#   --project=PROJECT   : Path to the project directory (PNFS)
#   --tag=BRANCH     : Medulla git branch (default: develop)
#######################################################################

set -e

PROJECT=""
BRANCH="develop"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --project=*) PROJECT="${1#*=}"; shift ;;
    --project)   PROJECT="$2";      shift 2 ;;
    --tag=*)  BRANCH="${1#*=}";  shift ;;
    --tag)    BRANCH="$2";       shift 2 ;;
    *) shift ;;
  esac
done

if [[ -z "$PROJECT" ]]; then
  echo "[ERROR] --project is required" >&2
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

echo "[INFO] Project: $PROJECT"
echo "[INFO] Branch:  $BRANCH"

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
# Determine this job's input file via the Phase 2 manifest
#######################################################################

ifdh cp "$PROJECT/variation_phase2_manifest.txt" variation_phase2_manifest.txt
JOBID=$(sed -n "$((PROCESS + 1))p" variation_phase2_manifest.txt)
if [[ -z "$JOBID" ]]; then
  echo "[ERROR] Could not determine JOBID for PROCESS=$PROCESS" >&2
  exit 1
fi
printf -v PADDED_JOBID "%04d" "$JOBID"
echo "[INFO] PROCESS=$PROCESS -> JOBID=$JOBID"

#######################################################################
# Stage input files
#######################################################################

ifdh cp "$PROJECT/variation_systematics_phase2.toml" variation_systematics_phase2.toml
echo "[INFO] Copied variation_systematics_phase2.toml"

ifdh cp "$PROJECT/variation_splines.root" variation_splines.root
echo "[INFO] Copied variation_splines.root"

ifdh cp "$PROJECT/output/output_jobid${PADDED_JOBID}.root" input_selection.root
echo "[INFO] Copied output_jobid${PADDED_JOBID}.root -> input_selection.root"

# Replace the per-job input path placeholder written at submission time.
sed -i 's|__INPUT_FILE__|input_selection.root|g' variation_systematics_phase2.toml
echo "[INFO] Set [input] path -> input_selection.root"

#######################################################################
# Run systematics (Phase 2: apply detsys weights from splines)
#######################################################################

echo "[INFO] Starting run_systematics (Phase 2)..."
ls -lrth
./systematics/run_systematics variation_systematics_phase2.toml
echo "[INFO] run_systematics completed"
ls -lrth

#######################################################################
# Stage output
#######################################################################

# Validate that the output file contains the expected variation tree.
# Write the check to a temp macro file because ROOT ignores stdin when a
# .root file is passed on the command line (here-docs are silently dropped).
VALIDATE_MACRO=$(mktemp /tmp/validate_varsys_XXXXXX.C)
cat > "$VALIDATE_MACRO" <<'MACRO'
{
    TFile *f = TFile::Open("output_varsys.root");
    if (!f || f->IsZombie()) {
        std::cerr << "[ERROR] Failed to open output_varsys.root" << std::endl;
        gSystem->Exit(1);
    }
    TTree *t = dynamic_cast<TTree*>(f->Get("events/NuMIFull/selected_variationTree"));
    if (!t) {
        std::cerr << "[ERROR] selected_variationTree not found in output_varsys.root" << std::endl;
        f->Close();
        gSystem->Exit(1);
    }
    std::cout << "[INFO] selected_variationTree has " << t->GetEntries() << " entries" << std::endl;
    f->Close();
}
MACRO
root -l -b -q "$VALIDATE_MACRO"
rm -f "$VALIDATE_MACRO"

printf -v OUTNAME "output_varsys_jobid%04d.root" "$JOBID"
ifdh cp output_varsys.root "$PROJECT/output/$OUTNAME"
echo "[INFO] Staged output to: $PROJECT/output/$OUTNAME"
