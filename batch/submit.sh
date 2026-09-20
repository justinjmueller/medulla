#!/bin/bash

#######################################################################
# Usage: submit.sh [--project=PROJECT] [--tag=TAG] [--sample=SAMPLE]
#                  [--force] [--stage-timeout=SECONDS] [--no-validate]
#                  [--jobs-per-process=M] [--total-jobs=N]
#
# Arguments:
#   --project=PROJECT   : Specify the project directory
#   --tag=TAG           : Git ref to checkout on grid nodes (default: develop)
#   --sample=SAMPLE     : Restrict job selection to the given sample
#   --force              : Overwrite pre-existing output files at the
#                          destination (e.g. left behind by a prior failed
#                          or resubmitted attempt at the same job) instead
#                          of failing the copy-back. Off by default so that
#                          two jobs unexpectedly racing to write the same
#                          output do not silently clobber one another.
#   --stage-timeout=SEC  : Wall-clock budget for each of the selection and
#                          systematics stages (default 7200). A stage that
#                          exceeds it is killed and reported as exit 124,
#                          which is recoverable information; being killed
#                          at the job's lifetime wall instead leaves no
#                          record at all.
#   --no-validate        : Skip output validation and copy back
#                          unconditionally (the pre-validation behaviour).
#                          For debugging only -- this re-enables the
#                          asymmetric-output failure mode described below.
#   --jobs-per-process=M : Number of job IDs this grid process takes on
#                          (default 1). See "Multiple job IDs per process".
#   --total-jobs=N       : Total number of job IDs the whole submission
#                          covers. Process P claims pending job IDs
#                          P*M .. P*M+M-1, capped at N, so the last process
#                          stops where the request stops instead of
#                          claiming job IDs nobody asked for. Omit to
#                          always claim a full M.
#
# Output validation
# -----------------
# Both artifacts are validated on the worker node and copied back only if
# the pair is good ("fail closed"). The failure this exists to prevent is
# output.root arriving on dCache without a usable output_sys.root beside
# it: systematics/src/main.cc skips any configured tree missing from its
# input and still returns 0, so a systematics run that produced nothing
# reports success. Because job completion is keyed on the selection file,
# such a job is marked complete and never resubmitted, and the breakage is
# only discovered at merge time.
#
# Every job ID also writes a small validation record to
# $PROJECT/output/val/, whether it passed or failed. That record -- not the
# HTCondor exit code -- is the reliable per-job signal.
#
# Timing and sizing fields in the record
# --------------------------------------
# Absolute times are epoch seconds; durations are whole seconds.
#   PROCESS_START      when this grid process began (shared by every job ID
#                      it runs); join against HTCondor's QDate/JobStartDate
#                      for queue wait
#   JOB_START/JOB_END  when this job ID began and ended; JOB_TIME is the
#                      difference
#   BUILD_TIME         per-process fixed cost, the total of:
#     SETUP_TIME       CVMFS/UPS environment setup
#     CLONE_TIME       git clone + checkout of --tag
#     COMPILE_TIME     cmake + make
#   STAGE_TIME         input staging for this job ID, of which:
#     STAGE_COPY_TIME  ifdh copies
#     STAGE_CHECK_TIME per-file ROOT integrity checks (one ROOT start each)
#   SELECTION_TIME, SYST_TIME, VALIDATE_TIME, COPY_TIME   the later stages
#   INPUT_BYTES        bytes of the inputs actually processed
#   STAGED_BYTES       bytes transferred, dropped inputs included
#   INPUT_EVENTS       recTree entries across the processed inputs; -1 if any
#                      input's count could not be read
# A field appears only once its stage has been reached, so a job ID that
# failed early carries only the fields for the stages it got through.
#
# Multiple job IDs per process
# ----------------------------
# Checking out and building medulla and setting up the environment is a
# fixed cost paid before any physics runs. With --jobs-per-process=M that
# cost is paid once, and the process then runs M job IDs in turn, each with
# its own inputs, outputs, validation, copy-back and record.
#
# The hazard this introduces is state leaking from one job ID into the
# next. The concrete case: a job ID whose selection writes nothing would
# find the previous job ID's output.root still sitting in the working
# directory, run systematics on it, validate it, and copy it back under its
# own name. Two things prevent that:
#   * each job ID runs in a fresh working directory of its own, deleted
#     afterwards, so nothing any stage writes can survive into the next;
#   * each job ID's body runs in a subshell, so its variables, its record,
#     and its exit (finish() exits) end with it.
# A failed job ID does not stop the ones after it. The process exits
# non-zero if any job ID failed, so the HTCondor exit code stays truthful,
# but the records remain the per-job-ID signal.
#######################################################################

# Print usage information
usage() {
  echo "Usage: submit.sh [--project=PROJECT] [--tag=TAG] [--sample=SAMPLE] [--force]"
  echo "                 [--stage-timeout=SECONDS] [--no-validate]"
  echo "                 [--jobs-per-process=M] [--total-jobs=N]"
  echo ""
  echo "Arguments:"
  echo "  --project=PROJECT   : Specify the project directory"
  echo "  --tag=TAG           : Git ref to checkout on grid nodes (default: develop)"
  echo "  --sample=SAMPLE     : Restrict job selection to the given sample"
  echo "  --force              : Overwrite pre-existing output files at the destination"
  echo "  --stage-timeout=SEC  : Per-stage wall-clock budget in seconds (default 7200)"
  echo "  --no-validate        : Skip output validation (debugging only)"
  echo "  --jobs-per-process=M : Job IDs this grid process runs in turn (default 1)"
  echo "  --total-jobs=N       : Job IDs the whole submission covers (caps the last process)"
}

# Initialize variables
PROJECT=""
TAG="develop"
SAMPLE=""
FORCE=""
STAGE_TIMEOUT=7200
VALIDATE=1
JOBS_PER_PROCESS=1
TOTAL_JOBS=""

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
    --force)
      FORCE=1
      shift
      ;;
    --stage-timeout=*)
      STAGE_TIMEOUT="${1#*=}"
      shift
      ;;
    --stage-timeout)
      STAGE_TIMEOUT="$2"
      shift 2
      ;;
    --no-validate)
      VALIDATE=0
      shift
      ;;
    --jobs-per-process=*)
      JOBS_PER_PROCESS="${1#*=}"
      shift
      ;;
    --jobs-per-process)
      JOBS_PER_PROCESS="$2"
      shift 2
      ;;
    --total-jobs=*)
      TOTAL_JOBS="${1#*=}"
      shift
      ;;
    --total-jobs)
      TOTAL_JOBS="$2"
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

# These feed shell arithmetic and SQL LIMIT/OFFSET, so accept digits only.
if [[ ! "$JOBS_PER_PROCESS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: --jobs-per-process must be a positive integer (got '${JOBS_PER_PROCESS}')." >&2
    exit 1
fi
if [[ -n "$TOTAL_JOBS" && ! "$TOTAL_JOBS" =~ ^[0-9]+$ ]]; then
    echo "Error: --total-jobs must be a non-negative integer (got '${TOTAL_JOBS}')." >&2
    exit 1
fi

#######################################################################
# Initial Setup
#######################################################################

# IFDH options
export IFDH_CP_MAXRETRIES=0
export IFDH_WEB_TIMEOUT=100

log_info() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] -- $*"
}

log_error() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] -- $*" >&2
}

# Seconds since the epoch, used to time each stage. The durations end up
# in the validation record, which is what makes it possible to size a
# job's requested lifetime, and how many job IDs one grid process should
# take on, from measurement rather than guesswork.
now_s() { date '+%s'; }

# Copy a local file to the given /pnfs destination. dCache does not allow
# overwriting a file in place, so a stale destination left behind by a
# prior failed or resubmitted attempt at this same job would otherwise
# cause the copy to fail. With --force, remove any existing file at the
# destination first; without it, leave a collision alone and let the copy
# fail loudly (the default, since silently overwriting could mask two
# jobs unexpectedly racing to write the same output).
copy_output() {
    local src="$1"
    local dst="$2"
    if [[ -n "$FORCE" ]]; then
        ifdh rm "$dst" 2>/dev/null
    fi
    ifdh cp "$src" "$dst"
}

# True when something already exists at the given /pnfs path.
dest_exists() {
    ifdh ls "$1" >/dev/null 2>&1
}

# Turn a stage's exit code into a status category. The three cases are
# worth separating because they point in different directions: a timeout
# means the budget or the input size is wrong, a signal means a genuine
# crash (systematics/src/main.cc does not null-check its input TFile, so a
# missing or corrupt output.root lands here), and anything else is an
# application-defined error to be read out of the stage's own log.
classify_stage_rc() {
    local rc="$1"
    local stage="$2"
    if [[ $rc -eq 124 ]]; then
        echo "${stage}_timeout"
    elif [[ $rc -ge 128 ]]; then
        echo "${stage}_crash"
    else
        echo "${stage}_error"
    fi
}

# Accumulates one job ID's validation record as it runs, so that a failure
# at any stage still produces a record describing how far it got. The
# record is copied back unconditionally by finish() below. Each job ID runs
# in its own subshell (see run_jobid), so these start empty for every one.
RECORD_LINES=()
record() { RECORD_LINES+=("$*"); }

# Write the validation record and exit with the given code. Every exit
# after the record is initialized goes through here: a job that dies
# without leaving a record is indistinguishable from one that never
# started, which is the ambiguity this whole mechanism exists to remove.
# Called inside run_jobid's subshell, the exit ends only the current job
# ID, not the process.
STATUS_SET=""
finish() {
    local code="$1"
    local status="$2"

    if [[ -n "$JOBID" ]]; then
        # One clock read for both, so JOB_TIME is exactly JOB_END - JOB_START.
        local job_end
        job_end=$(now_s)
        record "JOB_END=${job_end}"
        [[ -n "$JOB_T0" ]] && record "JOB_TIME=$(( job_end - JOB_T0 ))"
        record "EXIT_CODE=$code"
        [[ -z "$STATUS_SET" ]] && record "STATUS=$status"
        printf -v VALNAME "validation_jobid%04d.txt" "$JOBID"
        printf '%s\n' "${RECORD_LINES[@]}" > "$VALNAME"
        log_info "Validation record:"
        cat "$VALNAME"
        ifdh mkdir_p "$PROJECT/output/val" 2>/dev/null
        # --force semantics apply here too: a resubmitted job must be able
        # to replace its own earlier record.
        if [[ -n "$FORCE" ]]; then
            ifdh rm "$PROJECT/output/val/$VALNAME" 2>/dev/null
        fi
        ifdh cp "$VALNAME" "$PROJECT/output/val/$VALNAME" \
            || log_error "Failed to copy back the validation record."
    fi

    log_info "Job ID ${JOBID:-?} exiting with code $code (status: $status)."
    exit "$code"
}

# An absolute anchor for this process, so every record can be placed on a
# timeline and joined against HTCondor's own timestamps (QDate,
# JobStartDate, CompletionDate) -- durations alone cannot give queue wait,
# concurrency, or throughput over a campaign.
PROCESS_START=$(now_s)

# The fixed per-process cost is timed in three parts, since they call for
# different remedies: slow setup points at CVMFS, slow checkout at the
# network, and slow compilation at shipping a prebuilt tarball instead.

# Setup CVMFS area
t0=$(now_s)
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh

# Setup the required dependencies
setup sbnana v10_01_04 -q e26:prof
setup cmake v3_27_4
setup genie v3_04_02a -q e26:prof

ups active
SETUP_TIME=$(( $(now_s) - t0 ))

# Build medulla
t0=$(now_s)
git clone https://github.com/justinjmueller/medulla.git
cd medulla
git checkout ${TAG}
CLONE_TIME=$(( $(now_s) - t0 ))

# Resolve --tag to the commit this job actually built. A tag name is not
# enough: a branch moves, so two processes of one campaign that cloned either
# side of a push build different code with nothing to tell them apart. That
# happened on 2026-09-20, where jobs launched a minute after a fix was
# committed still ran the code from before it, and the only way to find out
# was to inspect the physics content of their output.
MEDULLA_COMMIT=$(git rev-parse HEAD 2>/dev/null)
[[ -z "$MEDULLA_COMMIT" ]] && MEDULLA_COMMIT=unknown
log_info "Built ${TAG} at commit ${MEDULLA_COMMIT}."

t0=$(now_s)
mkdir build && cd build
export CC=$(which gcc)
export CXX=$(which g++)
cmake .. -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_C_COMPILER=$CC
make -j4
COMPILE_TIME=$(( $(now_s) - t0 ))

# Everything per-process lives here; each job ID works in its own
# subdirectory, so binaries and staged files are addressed absolutely.
BUILD_DIR=$(pwd)
# The total of the three parts above, kept so records stay comparable with
# those written before the split.
BUILD_TIME=$(( $(now_s) - PROCESS_START ))

# Verify the build actually produced both binaries before claiming any job
# IDs. An unchecked build is the single most expensive failure mode on the
# grid: cmake picking up the wrong compiler offsite has previously caused
# every job in a batch to run a binary that was never produced, failing
# with exit 127 after the job had already burned its input staging. Fail
# here instead, while the job IDs are still eligible for resubmission.
SELECTION_BIN="$BUILD_DIR/selection/medulla"
SYSTEMATICS_BIN="$BUILD_DIR/systematics/run_systematics"
if [[ ! -x "$SELECTION_BIN" || ! -x "$SYSTEMATICS_BIN" ]]; then
    log_error "Build did not produce the expected binaries:"
    log_error "  $SELECTION_BIN      exists: $([[ -x $SELECTION_BIN ]] && echo yes || echo no)"
    log_error "  $SYSTEMATICS_BIN  exists: $([[ -x $SYSTEMATICS_BIN ]] && echo yes || echo no)"
    # No job ID claimed yet, so no record can be attributed to one; the
    # nonzero exit is the only signal, and it is an honest one here.
    exit 1
fi
log_info "Build verified in ${BUILD_TIME}s: both binaries present."

# Locate the output validator. launch_jobsub transfers it with the job
# (jobsub_submit -f dropbox://), which lands it in $CONDOR_DIR_INPUT. It is
# deliberately not taken from the checkout above: that is whatever --tag
# names, and a release tag that predates the validator would fail
# validation on every job and copy nothing back. The checkout copy is only
# a fallback for a job submitted by hand without shipping one.
VALIDATE_MACRO=""
if [[ -n "$CONDOR_DIR_INPUT" && -f "$CONDOR_DIR_INPUT/validate_pair.C" ]]; then
    VALIDATE_MACRO="$CONDOR_DIR_INPUT/validate_pair.C"
elif [[ -f "$BUILD_DIR/../batch/validate_pair.C" ]]; then
    VALIDATE_MACRO="$(cd "$BUILD_DIR/../batch" && pwd)/validate_pair.C"
fi
if [[ $VALIDATE -eq 1 ]]; then
    if [[ -z "$VALIDATE_MACRO" ]]; then
        # Refuse before claiming anything: every job ID would otherwise run
        # its full event loop only to fail validation and be discarded.
        log_error "Validation is enabled but validate_pair.C was neither shipped with the job"
        log_error "(\$CONDOR_DIR_INPUT) nor present in the '${TAG}' checkout. Submit through"
        log_error "launch_jobsub, which transfers it, or pass --no-validate."
        exit 1
    fi
    log_info "Using validator: ${VALIDATE_MACRO}"
fi

#######################################################################
# Prestage files shared by every job ID in this process
#######################################################################

DB="$BUILD_DIR/project.db"
SYST_CFG="$BUILD_DIR/systematics.toml"
MANIFEST_PATH="$BUILD_DIR/validation_manifest.txt"

# Copy the project database
if ! ifdh cp $PROJECT/project.db "$DB"; then
    log_error "Failed to stage the project database from $PROJECT."
    exit 1
fi

# Copy the systematics TOML file
if ! ifdh cp $PROJECT/systematics.toml "$SYST_CFG"; then
    log_error "Failed to stage systematics.toml from $PROJECT."
    exit 1
fi

# Copy the validation manifest, which says what each job's systematics
# output must contain. Projects created before the manifest existed do not
# have one; those fall back to skipping validation rather than failing, so
# that this script stays usable against them.
MANIFEST_MISSING=""
if ! ifdh cp $PROJECT/validation_manifest.txt "$MANIFEST_PATH"; then
    if [[ $VALIDATE -eq 1 ]]; then
        log_error "No validation_manifest.txt in $PROJECT -- this project predates output validation."
        log_error "Proceeding without validation. Recreate the project to enable it."
        MANIFEST_MISSING=1
        VALIDATE=0
    fi
fi

#######################################################################
# Claim this process's job IDs
#######################################################################

# Process P takes the M pending job IDs at ranks P*M .. P*M+M-1, so the
# processes of one submission cover disjoint, contiguous slices. With
# --total-jobs, the slice is cut off at the requested total, which keeps
# the last process from claiming job IDs beyond what was asked for.
#
# The ranking is computed against this process's own snapshot of
# project.db. A sync that marks jobs completed while a submission is still
# starting shifts the ranks between snapshots, so two processes can claim
# the same job ID or one can be skipped -- and a shifted rank misaligns a
# whole slice at once. The claimed list goes into every record so that this
# is at least detectable after the fact. Avoid syncing a project while its
# jobs are starting.
OFFSET=$(( ${PROCESS:-0} * JOBS_PER_PROCESS ))
LIMIT=$JOBS_PER_PROCESS
if [[ -n "$TOTAL_JOBS" ]]; then
    REMAINING=$(( TOTAL_JOBS - OFFSET ))
    if (( REMAINING <= 0 )); then
        log_info "Process ${PROCESS} lies beyond the ${TOTAL_JOBS} job ID(s) requested; nothing to do."
        exit 0
    fi
    (( REMAINING < LIMIT )) && LIMIT=$REMAINING
fi

# If --sample was given, restrict this to jobs belonging to that sample so
# that a per-sample jobsub_submit call only ever pulls from its own
# sample's pending jobs.
if [[ -n "$SAMPLE" ]]; then
    SAMPLE_ESC="${SAMPLE//\'/\'\'}"
    CLAIM_QUERY="SELECT jobid FROM jobs WHERE status != 'completed' AND sample = '${SAMPLE_ESC}' ORDER BY jobid LIMIT ${LIMIT} OFFSET ${OFFSET};"
else
    CLAIM_QUERY="SELECT jobid FROM jobs WHERE status != 'completed' ORDER BY jobid LIMIT ${LIMIT} OFFSET ${OFFSET};"
fi
CLAIMED_JOBIDS=()
while IFS= read -r jid; do
    [[ -n "$jid" ]] && CLAIMED_JOBIDS+=("$jid")
done < <(sqlite3 -noheader "$DB" "$CLAIM_QUERY")

if [[ ${#CLAIMED_JOBIDS[@]} -eq 0 ]]; then
    echo "Error: could not determine any JOBID for PROCESS=${PROCESS} (offset ${OFFSET}, limit ${LIMIT}). The project database may be empty or corrupt." >&2
    exit 1
fi
CLAIMED_LIST=$(IFS=,; echo "${CLAIMED_JOBIDS[*]}")
log_info "Process ${PROCESS} claimed ${#CLAIMED_JOBIDS[@]} job ID(s): ${CLAIMED_LIST}"

#######################################################################
# Run one job ID
#######################################################################

# The body of one job ID, run in a subshell (note the parentheses rather
# than braces) inside a fresh working directory. The subshell is what lets
# finish() exit on failure without ending the process, and guarantees that
# nothing one job ID sets -- its record, its status, its variables -- is
# visible to the next.
run_jobid() (
    JOBID="$1"
    SEQ="$2"
    WORKDIR="$3"
    RECORD_LINES=()
    STATUS_SET=""
    JOB_T0=$(now_s)

    if ! cd "$WORKDIR"; then
        log_error "Could not enter working directory $WORKDIR for job ID ${JOBID}."
        exit 1
    fi

    # Start the validation record. The header format matches the one the
    # failure-analysis tooling already parses (failure_ana/sources.py).
    record "JOB (${JOBID},${PROCESS}) VALIDATION"
    record "TAG=${TAG}"
    record "COMMIT=${MEDULLA_COMMIT}"
    record "HOSTNAME=$(hostname -f 2>/dev/null || hostname)"
    record "SITE=${GLIDEIN_Site:-unknown}"
    record "CLUSTER=${CLUSTER:-unknown}"
    record "PROCESS=${PROCESS}"
    record "JOBS_PER_PROCESS=${JOBS_PER_PROCESS}"
    record "SEQ=${SEQ}/${#CLAIMED_JOBIDS[@]}"
    record "CLAIMED_JOBIDS=${CLAIMED_LIST}"
    record "PROCESS_START=${PROCESS_START}"
    record "JOB_START=${JOB_T0}"
    record "BUILD_TIME=${BUILD_TIME}"
    record "SETUP_TIME=${SETUP_TIME}"
    record "CLONE_TIME=${CLONE_TIME}"
    record "COMPILE_TIME=${COMPILE_TIME}"
    [[ -n "$MANIFEST_MISSING" ]] && record "WARN=no_manifest_validation_skipped"

    # The sample this job processes. The validator needs it to address
    # 'events/<sample>' directly rather than inferring it from whatever
    # happens to be in the file, which makes a wrong-sample output a
    # detectable error rather than an accepted one.
    JOB_SAMPLE=$(sqlite3 -noheader "$DB" "SELECT sample FROM jobs WHERE jobid=${JOBID};")
    if [[ -z "$JOB_SAMPLE" ]]; then
        log_error "Could not determine the sample for job ID ${JOBID}."
        finish 1 "db_error"
    fi
    record "SAMPLE=${JOB_SAMPLE}"

    log_info "Job ID ${JOBID} (sample ${JOB_SAMPLE}), process ${PROCESS}, ${SEQ} of ${#CLAIMED_JOBIDS[@]}."

    sqlite3 -noheader -cmd ".mode list" "$DB" "SELECT cfg FROM configuration WHERE jobid=${JOBID};" > job_config.toml

    # Destination names, needed up front for the collision pre-check below.
    printf -v RAWNAME "output_jobid%04d.root" "$JOBID"
    printf -v SYSTNAME "output_systematics_jobid%04d.root" "$JOBID"

    # Check for a pre-existing destination before doing any work. Without
    # --force the copy-back would fail anyway, but only after the job has
    # paid for input staging and both event loops; catching it here costs
    # seconds.
    if [[ -z "$FORCE" ]]; then
        if dest_exists "$PROJECT/output/$RAWNAME" || dest_exists "$PROJECT/output/$SYSTNAME"; then
            log_error "Output already exists at the destination for job ID ${JOBID}."
            log_error "Another job may already have produced it. Re-run with --force to replace it."
            finish 1 "duplicate_record"
        fi
    fi

    # Copy the input data file(s)
    mkdir -p data

    # Extract all paths
    full_paths=$(grep '"/pnfs' job_config.toml | grep -o '"[^"]*"' | sed 's/"//g')
    echo "Found $(echo "$full_paths" | wc -l) input files to copy."

    # Copy input files, validating each staged copy as a well-formed,
    # non-corrupt ROOT file. A file that is zero-byte, truncated, or
    # otherwise fails to open (TFile::IsZombie()) is dropped from this job's
    # sample list rather than aborting the whole job -- this is distinct
    # from a file that opens fine but contains zero events, which is a
    # separate failure mode handled elsewhere. Dropped files are recorded in
    # bad_files.log and copied back alongside this job's output so the
    # exclusion can be reconciled later; POT is bookkept at file-read time,
    # so a dropped file is simply excluded from the dataset with no further
    # accounting needed here.
    #
    # The loop is also where the job's input is sized, so that run times can
    # be normalized per byte and per event: the integrity check already opens
    # every file, so reading recTree's entry count there costs nothing extra
    # and does not depend on the analysis configuration or on validation
    # being reached. Copy and check are timed separately because the check
    # starts a fresh ROOT process per file, a cost that grows with the number
    # of files rather than their size.
    good_paths=()
    bad_files=()
    input_bytes=0       # surviving inputs only: what the job actually processed
    staged_bytes=0      # everything transferred, dropped files included
    input_events=0      # recTree entries across surviving inputs
    events_known=1
    copy_ms=0
    check_ms=0
    stage_t0=$(now_s)
    for p in $full_paths; do
        echo "Copying input file: $p"
        t1=$(date '+%s%3N')
        ifdh cp "$p" data/
        copy_ms=$(( copy_ms + $(date '+%s%3N') - t1 ))
        b=$(basename "$p")
        size=$(stat -c %s "data/$b" 2>/dev/null || echo 0)
        staged_bytes=$(( staged_bytes + size ))

        t1=$(date '+%s%3N')
        check_out=$(root -l -b -q -e "TFile *f = TFile::Open(\"data/$b\"); if(!f || f->IsZombie()) gSystem->Exit(1); TTree *t = dynamic_cast<TTree*>(f->Get(\"recTree\")); std::cout << \"NENTRIES=\" << (t ? t->GetEntries() : -1) << std::endl; gSystem->Exit(0);" 2>/dev/null)
        check_rc=$?
        check_ms=$(( check_ms + $(date '+%s%3N') - t1 ))

        if [[ $check_rc -eq 0 ]]; then
            good_paths+=("$p")
            input_bytes=$(( input_bytes + size ))
            n=$(printf '%s\n' "$check_out" | sed -n 's/^NENTRIES=//p' | tail -1)
            if [[ "$n" =~ ^[0-9]+$ ]]; then
                input_events=$(( input_events + n ))
            else
                # No recTree, or no count printed: a partial sum would be
                # silently wrong, so report the total as unknown instead.
                events_known=0
            fi
        else
            echo "Warning: dropping unreadable/corrupt input file: $p" >&2
            bad_files+=("$p")
            rm -f "data/$b"
        fi
    done
    STAGE_TIME=$(( $(now_s) - stage_t0 ))
    ls -lrth data/

    record "N_INPUTS=${#good_paths[@]}"
    record "N_DROPPED_INPUTS=${#bad_files[@]}"
    record "STAGE_TIME=${STAGE_TIME}"
    record "STAGE_COPY_TIME=$(( (copy_ms + 500) / 1000 ))"
    record "STAGE_CHECK_TIME=$(( (check_ms + 500) / 1000 ))"
    record "INPUT_BYTES=${input_bytes}"
    record "STAGED_BYTES=${staged_bytes}"
    # -1 when any surviving input lacked a readable recTree count.
    record "INPUT_EVENTS=$([[ $events_known -eq 1 ]] && echo "$input_events" || echo -1)"

    # Modify the job_config.toml to use local paths for the surviving files,
    # and remove the array entries for any files that were dropped above.
    # The path array may be serialized on a single line, so a dropped entry
    # must be removed in place (its quoted string plus a trailing comma, if
    # present) -- deleting the whole matching line would risk wiping out
    # every other path sharing that line. TOML tolerates a trailing comma
    # before the closing bracket, so removing just the "path", substring
    # (or, for a trailing entry with no comma of its own, just "path") leaves
    # a syntactically valid array.
    for p in "${good_paths[@]}"; do
        b=$(basename "$p")
        sed -i "s#\"$p\"#\"data/$b\"#g" job_config.toml
    done
    for p in "${bad_files[@]}"; do
        sed -i -e "s#\"$p\",[[:space:]]*##g" -e "s#\"$p\"##g" job_config.toml
    done

    # If every input file for this job was dropped, fail explicitly rather
    # than running the selection binary against zero input.
    if [[ ${#good_paths[@]} -eq 0 ]]; then
        log_error "All input files for this job were unreadable/corrupt."
        finish 1 "transfer_lost"
    fi

    # The systematics configuration reads its universe weights from the
    # fixed glob 'data/*flat*.root'. If the surviving inputs contain no flat
    # CAF, the systematics step would run to completion and silently produce
    # output with no weights applied, which is exactly the kind of
    # quietly-wrong result this validation exists to stop. Check before
    # spending the event loop.
    if [[ $VALIDATE -eq 1 ]] && grep -q 'add_weights' "$MANIFEST_PATH" 2>/dev/null; then
        shopt -s nullglob
        flat_inputs=(data/*flat*.root)
        shopt -u nullglob
        if [[ ${#flat_inputs[@]} -eq 0 ]]; then
            log_error "No 'data/*flat*.root' inputs survived staging, but this project applies weights."
            finish 1 "transfer_lost"
        fi
        log_info "Found ${#flat_inputs[@]} flat CAF input(s) for weights."
    fi

    # Record any dropped files alongside this job's output for bookkeeping.
    if [[ ${#bad_files[@]} -gt 0 ]]; then
        printf '%s\n' "${bad_files[@]}" > bad_files.log
        printf -v BADNAME "bad_files_jobid%04d.log" "$JOBID"
        copy_output bad_files.log $PROJECT/output/$BADNAME
    fi

    # Dump some info for debugging
    cat job_config.toml
    ls -lrth .
    ls -lrth data/

    ###################################################################
    # Run the analysis
    ###################################################################

    # Both stages are wrapped in `timeout` so that a hang produces exit 124
    # and a validation record, rather than the job being killed at its
    # lifetime wall with nothing written back at all.

    # Run medulla (selection)
    log_info "Running selection (timeout ${STAGE_TIMEOUT}s)."
    t0=$(now_s)
    timeout "$STAGE_TIMEOUT" "$SELECTION_BIN" job_config.toml
    SELECTION_RC=$?
    SELECTION_TIME=$(( $(now_s) - t0 ))
    record "SELECTION_EXIT_CODE=${SELECTION_RC}"
    record "SELECTION_TIME=${SELECTION_TIME}"
    log_info "Selection finished with exit code ${SELECTION_RC} in ${SELECTION_TIME}s."
    ls -lrth

    if [[ $SELECTION_RC -ne 0 ]]; then
        log_error "Selection failed; not running systematics and not copying anything back."
        finish "$SELECTION_RC" "$(classify_stage_rc "$SELECTION_RC" selection)"
    fi

    # Run medulla (systematics). Note that this is deliberately run *before*
    # anything is copied back: copying the selection output first is what
    # produced orphaned output_jobid*.root files with no systematics
    # partner, and those get marked complete and never resubmitted.
    log_info "Running systematics (timeout ${STAGE_TIMEOUT}s)."
    t0=$(now_s)
    timeout "$STAGE_TIMEOUT" "$SYSTEMATICS_BIN" "$SYST_CFG"
    SYST_RC=$?
    SYST_TIME=$(( $(now_s) - t0 ))
    record "SYST_EXIT_CODE=${SYST_RC}"
    record "SYST_TIME=${SYST_TIME}"
    log_info "Systematics finished with exit code ${SYST_RC} in ${SYST_TIME}s."
    ls -lrth

    if [[ $SYST_RC -ne 0 ]]; then
        log_error "Systematics failed; not copying anything back."
        finish "$SYST_RC" "$(classify_stage_rc "$SYST_RC" sys)"
    fi

    ###################################################################
    # Validate the output pair
    ###################################################################

    # The systematics binary returns 0 even when it silently skipped every
    # tree, so its exit code alone proves nothing about what it produced.
    # The output has to be inspected.
    if [[ $VALIDATE -eq 1 ]]; then
        log_info "Validating the output pair."
        t0=$(now_s)
        root -l -b -q "${VALIDATE_MACRO}(\"output.root\",\"output_sys.root\",\"${JOB_SAMPLE}\",\"${MANIFEST_PATH}\",\"validation_report.txt\")"
        VALIDATE_RC=$?
        VALIDATE_TIME=$(( $(now_s) - t0 ))
        record "VALIDATE_EXIT_CODE=${VALIDATE_RC}"
        record "VALIDATE_TIME=${VALIDATE_TIME}"

        # Fold the macro's findings into the record. Its last line is a
        # STATUS=, which is more specific than anything this script could
        # infer, so let it stand as the job's status.
        if [[ -f validation_report.txt ]]; then
            log_info "Validation report:"
            cat validation_report.txt
            while IFS= read -r line; do
                [[ -z "$line" ]] && continue
                record "$line"
                [[ "$line" == STATUS=* ]] && STATUS_SET=1
            done < validation_report.txt
        else
            log_error "Validator produced no report; treating as a failure."
            record "ERROR=validator_no_report"
        fi

        if [[ $VALIDATE_RC -ne 0 ]]; then
            log_error "Output validation failed (code ${VALIDATE_RC}); copying nothing back."
            log_error "Job ID ${JOBID} stays pending and is eligible for resubmission."
            finish "$VALIDATE_RC" "sys_error"
        fi
        log_info "Output pair validated in ${VALIDATE_TIME}s."
    else
        record "WARN=validation_disabled"
    fi

    ###################################################################
    # Copy the output back
    ###################################################################

    # Both files or neither: a lone selection output would be counted as a
    # completed job and never revisited. Copied now, as this job ID
    # finishes, rather than at the end of the process, so a process
    # preempted partway through keeps the job IDs it already completed.
    t0=$(now_s)
    if ! copy_output output.root $PROJECT/output/$RAWNAME; then
        log_error "Failed to copy back the selection output."
        finish 1 "transfer_lost"
    fi

    if ! copy_output output_sys.root $PROJECT/output/$SYSTNAME; then
        log_error "Failed to copy back the systematics output."
        # The selection output is already on dCache and would be read as a
        # completed job, so remove it rather than leave the pair
        # half-written.
        log_error "Removing the already-copied selection output to keep the pair atomic."
        ifdh rm $PROJECT/output/$RAWNAME 2>/dev/null \
            || log_error "Could not remove $PROJECT/output/$RAWNAME -- clean it up before the next sync."
        finish 1 "transfer_lost"
    fi
    record "COPY_TIME=$(( $(now_s) - t0 ))"

    log_info "Job ID ${JOBID} completed successfully."
    finish 0 "ok"
)

#######################################################################
# Run every claimed job ID in turn
#######################################################################

N_OK=0
FAILED_JOBIDS=()
SEQ=0
for JOBID in "${CLAIMED_JOBIDS[@]}"; do
    SEQ=$(( SEQ + 1 ))

    # A fresh directory per job ID, removed afterwards: the isolation that
    # stops one job ID's leftovers being picked up by the next, and the
    # reason disk use does not grow with --jobs-per-process.
    printf -v WORKDIR "%s/work_jobid%04d" "$BUILD_DIR" "$JOBID"
    rm -rf "$WORKDIR"
    if ! mkdir -p "$WORKDIR"; then
        log_error "Could not create working directory $WORKDIR; skipping job ID ${JOBID}."
        FAILED_JOBIDS+=("$JOBID")
        continue
    fi

    log_info "===== Job ID ${JOBID} (${SEQ} of ${#CLAIMED_JOBIDS[@]}) ====="
    run_jobid "$JOBID" "$SEQ" "$WORKDIR"
    rc=$?
    rm -rf "$WORKDIR"

    if [[ $rc -eq 0 ]]; then
        N_OK=$(( N_OK + 1 ))
    else
        FAILED_JOBIDS+=("$JOBID")
    fi
done

log_info "Process ${PROCESS} finished: ${N_OK} of ${#CLAIMED_JOBIDS[@]} job ID(s) succeeded."
if [[ ${#FAILED_JOBIDS[@]} -gt 0 ]]; then
    log_error "Failed job ID(s): $(IFS=,; echo "${FAILED_JOBIDS[*]}") -- see their validation records."
    exit 1
fi
exit 0
