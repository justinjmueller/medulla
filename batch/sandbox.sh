#!/bin/bash
#######################################################################
# Usage: sandbox.sh --project=PROJECT --jobid=JOBID [--tag=TAG] [--no-build]
#
# Stages all inputs for a single job into the current directory so
# the job can be run interactively for debugging.
#
# After this script completes, run the job manually with:
#   ./selection/medulla job_config.toml
#   ./systematics/run_systematics systematics.toml
#
# Arguments:
#   --project=PROJECT   Path to the project directory (PNFS)
#   --jobid=JOBID       Job ID to stage
#   --tag=TAG           Medulla git ref to build (default: develop)
#   --no-build          Skip cloning/building medulla (use existing build)
#######################################################################

PROJECT=""
JOBID=""
TAG="develop"
NO_BUILD=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --project=*) PROJECT="${1#*=}"; shift ;;
    --project)   PROJECT="$2";      shift 2 ;;
    --jobid=*)   JOBID="${1#*=}";   shift ;;
    --jobid)     JOBID="$2";        shift 2 ;;
    --tag=*)     TAG="${1#*=}";     shift ;;
    --tag)       TAG="$2";          shift 2 ;;
    --no-build)  NO_BUILD=1;        shift ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$PROJECT" || -z "$JOBID" ]]; then
  echo "Usage: sandbox.sh --project=PROJECT --jobid=JOBID [--tag=TAG] [--no-build]" >&2
  exit 1
fi

export IFDH_CP_MAXRETRIES=2
export IFDH_WEB_TIMEOUT=300

#######################################################################
# Environment setup
# UPS setup commands return non-zero on warnings even when they succeed,
# so run them before enabling set -e.
#######################################################################

source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup sbnana v10_01_02_01 -q e26:prof
setup cmake v3_27_4

set -e

#######################################################################
# Build medulla (unless --no-build)
#######################################################################

if [[ "$NO_BUILD" -eq 0 ]]; then
  echo "[INFO] Cloning medulla (tag: $TAG)..."
  git clone https://github.com/justinjmueller/medulla.git
  cd medulla
  git checkout "$TAG"
  mkdir build && cd build
  export CC=$(which gcc)
  export CXX=$(which g++)
  cmake .. -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_C_COMPILER=$CC
  make -j4
  cd ../..
  echo "[INFO] Build complete."
else
  echo "[INFO] Skipping build (--no-build set)."
  if [[ ! -x "./selection/medulla" || ! -x "./systematics/run_systematics" ]]; then
    echo "[WARN] ./selection/medulla or ./systematics/run_systematics not found — build medulla before running."
  fi
fi

#######################################################################
# Stage project database and extract job config
#######################################################################

echo "[INFO] Copying project.db..."
ifdh cp "$PROJECT/project.db" project.db

echo "[INFO] Extracting config for job $JOBID..."
sqlite3 -noheader -cmd ".mode list" project.db \
  "SELECT cfg FROM configuration WHERE jobid=${JOBID};" > job_config.toml

if [[ ! -s job_config.toml ]]; then
  echo "[ERROR] No configuration found for job ID $JOBID." >&2
  exit 1
fi

#######################################################################
# Stage systematics TOML
#######################################################################

echo "[INFO] Copying systematics.toml..."
ifdh cp "$PROJECT/systematics.toml" systematics.toml

#######################################################################
# Stage input data files
#######################################################################

full_paths=$(grep '"/pnfs' job_config.toml | grep -o '"[^"]*"' | sed 's/"//g')
n_files=$(echo "$full_paths" | grep -c . || true)
echo "[INFO] Staging $n_files input file(s) to data/..."
mkdir -p data

for p in $full_paths; do
  echo "[INFO]   $p"
  ifdh cp "$p" data/
done

# Rewrite paths in job_config.toml to point at local data/
for p in $full_paths; do
  b=$(basename "$p")
  sed -i "s#\"$p\"#\"data/$b\"#g" job_config.toml
done

#######################################################################
# Summary
#######################################################################

echo ""
echo "========================================================"
echo "  Sandbox ready for job $JOBID"
echo "========================================================"
echo ""
echo "  Config   : job_config.toml"
echo "  Systematics TOML : systematics.toml"
echo "  Data files:"
ls data/ | sed 's/^/    /'
echo ""
echo "  To run interactively:"
echo "    ./selection/medulla job_config.toml"
echo "    ./systematics/run_systematics systematics.toml"
echo "========================================================"
