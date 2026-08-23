#!/bin/bash

#######################################################################
# Usage: submit.sh [--project=PROJECT] [--tag=TAG] [--sample=SAMPLE]
#
# Arguments:
#   --project=PROJECT   : Specify the project directory
#   --tag=TAG           : Git ref to checkout on grid nodes (default: develop)
#   --sample=SAMPLE     : Restrict job selection to the given sample
#######################################################################

# Print usage information
usage() {
  echo "Usage: submit.sh [--project=PROJECT] [--tag=TAG] [--sample=SAMPLE]"
  echo ""
  echo "Arguments:"
  echo "  --project=PROJECT   : Specify the project directory"
  echo "  --tag=TAG           : Git ref to checkout on grid nodes (default: develop)"
  echo "  --sample=SAMPLE     : Restrict job selection to the given sample"
}

# Initialize variables
PROJECT=""
TAG="develop"
SAMPLE=""

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project=*)
      PROJECT="${1#*=}"
      shift
      ;;
    --project)
      PROJECT="$2"
      shift 2
      ;;
    --tag=*)
      TAG="${1#*=}"
      shift
      ;;
    --tag)
      TAG="$2"
      shift 2
      ;;
    --sample=*)
      SAMPLE="${1#*=}"
      shift
      ;;
    --sample)
      SAMPLE="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --) # end of options
      shift
      break
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

#######################################################################
# Check for required arguments
#######################################################################
missing_args=()
[[ -z "$PROJECT" ]] && missing_args+=("--project")

if [[ ${#missing_args[@]} -gt 0 ]]; then
    echo "Error: Missing required argument(s): ${missing_args[*]}" >&2
    usage
    exit 1
fi

#######################################################################
# Initial Setup
#######################################################################

# IFDH options
export IFDH_CP_MAXRETRIES=0
export IFDH_WEB_TIMEOUT=100

# Setup CVMFS area
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh

# Setup the required dependencies
setup sbnana v10_01_04 -q e26:prof
setup cmake v3_27_4

ups active

# Build medulla
git clone https://github.com/justinjmueller/medulla.git
cd medulla
git checkout ${TAG}
mkdir build && cd build
export CC=$(which gcc)
export CXX=$(which g++)
cmake .. -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_C_COMPILER=$CC
make -j4

#######################################################################
# Prestage job-specific files
#######################################################################

# Copy the project database
ifdh cp $PROJECT/project.db project.db

# Extract this job's configuration file. First, we get the job ID for this
# process by checking against the list of not-yet-completed jobs in the project
# database. If --sample was given, restrict this to jobs belonging to that
# sample so that a per-sample jobsub_submit call only ever pulls from its
# own sample's pending jobs.
if [[ -n "$SAMPLE" ]]; then
    SAMPLE_ESC="${SAMPLE//\'/\'\'}"
    JOBID=$(sqlite3 -noheader project.db "SELECT jobid FROM jobs WHERE status != 'completed' AND sample = '${SAMPLE_ESC}' ORDER BY jobid LIMIT 1 OFFSET ${PROCESS};")
else
    JOBID=$(sqlite3 -noheader project.db "SELECT jobid FROM jobs WHERE status != 'completed' ORDER BY jobid LIMIT 1 OFFSET ${PROCESS};")
fi
if [[ -z "$JOBID" ]]; then
    echo "Error: could not determine JOBID for PROCESS=${PROCESS}. The project database may be empty or corrupt." >&2
    exit 1
fi
sqlite3 -noheader -cmd ".mode list" project.db "SELECT cfg FROM configuration WHERE jobid=${JOBID};" > job_config.toml

# Copy the systematics TOML file
ifdh cp $PROJECT/systematics.toml systematics.toml

# Copy the input data file(s)
mkdir data

# Extract all paths
full_paths=$(grep '"/pnfs' job_config.toml | grep -o '"[^"]*"' | sed 's/"//g')
echo "Found $(echo "$full_paths" | wc -l) input files to copy."

# Copy input files, validating each staged copy as a well-formed,
# non-corrupt ROOT file. A file that is zero-byte, truncated, or otherwise
# fails to open (TFile::IsZombie()) is dropped from this job's sample list
# rather than aborting the whole job -- this is distinct from a file that
# opens fine but contains zero events, which is a separate failure mode
# handled elsewhere. Dropped files are recorded in bad_files.log and copied
# back alongside this job's output so the exclusion can be reconciled later;
# POT is bookkept at file-read time, so a dropped file is simply excluded
# from the dataset with no further accounting needed here.
mkdir -p data
good_paths=()
bad_files=()
for p in $full_paths; do
    echo "Copying input file: $p"
    ifdh cp "$p" data/
    b=$(basename "$p")
    if root -l -b -q -e "TFile *f = TFile::Open(\"data/$b\"); gSystem->Exit((!f || f->IsZombie()) ? 1 : 0);" >/dev/null 2>&1; then
        good_paths+=("$p")
    else
        echo "Warning: dropping unreadable/corrupt input file: $p" >&2
        bad_files+=("$p")
        rm -f "data/$b"
    fi
done
ls -lrth data/

# Modify the job_config.toml to use local paths for the surviving files, and
# remove the array entries for any files that were dropped above. TOML
# tolerates a trailing comma before the closing bracket, so deleting any one
# entry's line leaves a syntactically valid array.
for p in "${good_paths[@]}"; do
    b=$(basename "$p")
    sed -i "s#\"$p\"#\"data/$b\"#g" job_config.toml
done
for p in "${bad_files[@]}"; do
    sed -i "\#\"$p\"#d" job_config.toml
done

# If every input file for this job was dropped, fail explicitly rather than
# running the selection binary against zero input.
if [[ ${#good_paths[@]} -eq 0 ]]; then
    echo "Error: all input files for this job were unreadable/corrupt." >&2
    exit 1
fi

# Record any dropped files alongside this job's output for bookkeeping.
if [[ ${#bad_files[@]} -gt 0 ]]; then
    printf '%s\n' "${bad_files[@]}" > bad_files.log
    printf -v BADNAME "bad_files_jobid%04d.log" "$JOBID"
    ifdh cp bad_files.log $PROJECT/output/$BADNAME
fi

# Dump some info for debugging
cat job_config.toml
ls -lrth .
ls -lrth data/

#######################################################################
# Run the analysis
#######################################################################

# Run medulla (selection)
./selection/medulla job_config.toml
ls -lrth

# Copy output file to the output directory
printf -v RAWNAME "output_jobid%04d.root" "$JOBID"
ifdh cp output.root $PROJECT/output/$RAWNAME

# Run medulla (systematics)
./systematics/run_systematics systematics.toml
ls -lrth

# Copy output file to the output directory
printf -v SYSTNAME "output_systematics_jobid%04d.root" "$JOBID"
ifdh cp output_sys.root $PROJECT/output/$SYSTNAME