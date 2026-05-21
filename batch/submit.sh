#!/bin/bash

#######################################################################
# Usage: submit.sh [--project=PROJECT] [--branch=BRANCH]
# 
# Arguments:
#   --project=PROJECT   : Specify the project directory
#   --branch=BRANCH     : Git branch to checkout (default: develop)
#######################################################################

# Initialize variables
PROJECT=""
BRANCH="develop"

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
    --branch=*)
      BRANCH="${1#*=}"
      shift
      ;;
    --branch)
      BRANCH="$2"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    --) # end of options
      shift
      break
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
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
setup sbnana v10_01_02_01 -q e26:prof
setup cmake v3_27_4

ups active

# Build medulla
git clone https://github.com/justinjmueller/medulla.git
cd medulla
git checkout "$BRANCH"
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
# database.
JOBID=$(sqlite3 -noheader project.db "SELECT jobid FROM jobs WHERE status != 'completed' ORDER BY jobid LIMIT 1 OFFSET ${PROCESS};")
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

# Copy input files
mkdir -p data
for p in $full_paths; do
    echo "Copying input file: $p"
    ifdh cp "$p" data/
done
ls -lrth data/

# Modify the job_config.toml to use local paths
for p in $full_paths; do
    b=$(basename "$p")
    sed -i "s#\"$p\"#\"data/$b\"#g" job_config.toml
done

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