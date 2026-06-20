#!/usr/bin/env bash

##############################################################################################################
# MitoHPC Batch Container Script
# 
# This script provides an easy way to run MitoHPC with parallel sample processing using containers.
# It's designed to work similarly to the user's existing apptainer command but with built-in
# parallel processing capabilities.
#
# Usage:
#   mitohpc-batch-container.sh <working_directory> [num_threads] [container_image]
#
# Example:
#   mitohpc-batch-container.sh /path/to/groupCram 4 "docker://ghcr.io/jlanej/mitohpc:main"
#
# This is equivalent to the user's original command but with parallel processing:
#   apptainer exec --bind "$groupCram":"$groupCram" --pwd "$groupCram" \
#     --env HP_ADIR=bams,HP_ODIR=out,HP_IN=in.txt \
#     "docker://ghcr.io/jlanej/mitohpc:main" mitohpc.sh
##############################################################################################################

set -e

# Parse arguments
WORKING_DIR="$1"
NUM_THREADS="${2:-$(nproc 2>/dev/null || echo 4)}"
CONTAINER_IMAGE="${3:-docker://ghcr.io/jlanej/mitohpc:main}"

# Structural-variant (large-deletion) calling is ON by default (HP_SV=callsv); it writes
# additional *.sv.* outputs and never changes the existing deliverables. To disable, set
# HP_SV to empty in the environment:  HP_SV= mitohpc-batch-container.sh <dir> ...
SV_MODE="${HP_SV-callsv}"

# Function to show usage
show_usage() {
    cat << EOF
Usage: $0 <working_directory> [num_threads] [container_image]

Arguments:
    working_directory    Directory containing BAM/CRAM files (equivalent to groupCram)
    num_threads         Number of parallel threads to use (default: auto-detect CPU cores)
    container_image     Container image to use (default: docker://ghcr.io/jlanej/mitohpc:main)

Structural-variant (large-deletion) calling is ON by default (HP_SV=callsv); it adds
*.sv.* outputs without changing existing results, and is SKIPPED with a warning if the
image's python lacks pysam (e.g. the :main image until the SV branch is merged + rebuilt).
Use the :sv-calling image (3rd arg) to run it, or disable with: HP_SV= $0 <dir> ...

Example:
    $0 /data/samples 4 "docker://ghcr.io/jlanej/mitohpc:main"

This script expects the working directory to contain:
    bams/           - Directory with BAM files, OR
    crams/          - Directory with CRAM files
    
Output will be created in:
    out/            - MitoHPC output directory
    in.txt          - Input file list (auto-generated)

The script will:
1. Process all BAM/CRAM files in parallel using the specified number of threads
2. Generate comprehensive results including VCF files, haplogroups, and summaries
3. Provide progress information and error handling

This provides the same functionality as the original MitoHPC container call
but with built-in parallel processing for better performance.
EOF
}

# Validate arguments
if [ $# -lt 1 ] || [ "$WORKING_DIR" = "-h" ] || [ "$WORKING_DIR" = "--help" ]; then
    show_usage
    exit 0
fi

if [ -z "$WORKING_DIR" ]; then
    echo "Error: Working directory is required" >&2
    show_usage
    exit 1
fi

if [ ! -d "$WORKING_DIR" ]; then
    echo "Error: Working directory does not exist: $WORKING_DIR" >&2
    exit 1
fi

# Convert to absolute path
WORKING_DIR=$(cd "$WORKING_DIR" && pwd)

# Validate number of threads
if ! [[ "$NUM_THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: Number of threads must be a positive integer" >&2
    show_usage
    exit 1
fi

# Check if apptainer is available
if ! command -v apptainer > /dev/null 2>&1; then
    echo "Error: apptainer command not found. Please install Apptainer/Singularity" >&2
    exit 1
fi

# Check for input data
if [ ! -d "$WORKING_DIR/bams" ] && [ ! -d "$WORKING_DIR/crams" ]; then
    echo "Error: Neither 'bams' nor 'crams' directory found in $WORKING_DIR" >&2
    exit 1
fi

# Determine data directory: prefer a populated bams/, else a populated crams/. crams/ is a
# documented input, so it must NOT be hardcoded to bams (a crams-only working dir would
# otherwise look at a nonexistent bams/ and find zero inputs).
if [ -d "$WORKING_DIR/bams" ] && find "$WORKING_DIR/bams" -name "*.bam" 2>/dev/null | head -1 | grep -q .; then
    DATA_DIR="bams"
elif [ -d "$WORKING_DIR/crams" ] && find "$WORKING_DIR/crams" -name "*.cram" 2>/dev/null | head -1 | grep -q .; then
    DATA_DIR="crams"
elif [ -d "$WORKING_DIR/bams" ]; then
    DATA_DIR="bams"   # exists but empty — let the input-list step below report the error
else
    DATA_DIR="crams"
fi

echo "MitoHPC Batch Container Processing"
echo "=================================="
echo "Working directory: $WORKING_DIR"
echo "Data directory: $DATA_DIR"
echo "Number of threads: $NUM_THREADS"
echo "Container image: $CONTAINER_IMAGE"
echo "SV calling (HP_SV): ${SV_MODE:-off (disabled)}"
echo

# Count input files (group the -name alternation so any future trailing predicate binds to both)
FILE_COUNT=$(find "$WORKING_DIR/$DATA_DIR" \( -name "*.bam" -o -name "*.cram" \) | wc -l)
echo "Found $FILE_COUNT input files to process"

# Create output directory
mkdir -p "$WORKING_DIR/out"

# Run MitoHPC with parallel processing
# We'll use a modified approach that processes samples in parallel within the container
echo "Starting parallel processing..."

# Create the in-container runner. Use mktemp (atomic O_EXCL create, mode 0600) with a UNIQUE
# name rather than a predictable, world-readable path in the shared working dir: this closes a
# symlink/TOCTOU hijack of the script we then chmod +x and execute, and stops two concurrent
# runs in the same workdir from clobbering each other's runner. The trap removes it on
# exit/signal/failure — not only on the happy path. (It still lives under $WORKING_DIR so the
# existing --bind makes it visible in the container.)
TEMP_SCRIPT=$(mktemp "$WORKING_DIR/.run_parallel_mitohpc.XXXXXX.sh")
trap 'rm -f "$TEMP_SCRIPT"' EXIT INT TERM
cat > "$TEMP_SCRIPT" << 'EOF'
#!/usr/bin/env bash
set -e

NUM_THREADS="$1"

# Source the MitoHPC initialization.
# init.sh unconditionally re-exports HP_SV (to its default) AND HP_ADIR (to $PWD/bams/),
# clobbering the values we passed via --env. Capture both first and restore them afterwards:
# keeps structural-variant calling enabled, and — crucially — keeps HP_ADIR=crams for a
# crams-only run (otherwise init.sh forces bams/ and the input find below returns nothing).
HP_SV_REQ="${HP_SV:-}"
HP_ADIR_REQ="${HP_ADIR:-}"
. $HP_SDIR/init.sh
export HP_SV="$HP_SV_REQ"
[ -n "$HP_ADIR_REQ" ] && export HP_ADIR="$HP_ADIR_REQ"

# Budget per-sample CPU/RAM (HP_P threads, HP_MM sort buffer, HP_JOPT java -Xmx) for the
# NUM_THREADS samples we run concurrently below. resource_budget is cgroup-aware, so inside
# the container it sees the real CPU/RAM limits (not just the host's nproc) and divides them
# across the concurrent samples. Sourced AFTER init.sh so it overrides the per-sample defaults
# without editing init.sh. An explicit --env HP_P (if forwarded) is respected as a cap.
# Guarded: older images (e.g. :main before this branch is merged + rebuilt) lack this file —
# degrade to init.sh's per-sample defaults rather than aborting the run under `set -e`.
if [ -r "$HP_SDIR/resource_budget.sh" ]; then
    . $HP_SDIR/resource_budget.sh
    resource_budget "$NUM_THREADS"
else
    echo "[mitohpc] note: resource_budget.sh not in image; using init.sh per-sample defaults (HP_P=$HP_P)" >&2
fi

# Generate input file
echo "Generating input file list..."
find "$HP_ADIR/" \( -name "*.bam" -o -name "*.cram" \) -readable | ls2in.pl -out "$HP_ODIR" | sort -V > "$HP_IN"

if [ ! -s "$HP_IN" ]; then
    echo "Error: No input files found" >&2
    exit 1
fi

# Security: run.sh emits each per-sample command as SHELL TEXT that parallel/xargs re-parse, so
# a sample name or path containing shell metacharacters (or a space, which also breaks ls2in.pl's
# field split) could inject or corrupt commands. Reject anything outside a conservative safe set
# up front — mtDNA BAM/CRAM names never need these. (tr makes the tab field-separator a newline
# so spaces inside a field are still caught.)
if tr '\t' '\n' < "$HP_IN" | LC_ALL=C grep -qE '[^[:alnum:]._/:+,@%=-]'; then
    echo "Error: $HP_IN has a sample name/path with unsafe characters (shell metacharacters or spaces):" >&2
    tr '\t' '\n' < "$HP_IN" | LC_ALL=C grep -E '[^[:alnum:]._/:+,@%=-]' | sed 's/^/    /' >&2
    echo "Rename the offending file(s) to use only [A-Za-z0-9._/:+,@%=-] and re-run." >&2
    exit 1
fi

TOTAL_SAMPLES=$(wc -l < "$HP_IN")
echo "Processing $TOTAL_SAMPLES samples with $NUM_THREADS threads"

# Generate processing commands
$HP_SDIR/run.sh > run.all.sh

# Extract and run filter commands in parallel
grep "filter.sh" run.all.sh > filter.commands.sh

echo "Running parallel processing..."
if command -v parallel > /dev/null 2>&1; then
    # GNU parallel's --progress draws on /dev/tty; with no controlling terminal (apptainer under
    # a scheduler, nohup, redirected output) it spams "cannot open /dev/tty". Enable it only when
    # a tty is actually writable; always keep --joblog for post-run per-sample status (no tty needed).
    PARALLEL_OPTS="--joblog parallel.log"
    if { : > /dev/tty; } 2>/dev/null; then PARALLEL_OPTS="$PARALLEL_OPTS --progress"; fi
    parallel -j "$NUM_THREADS" $PARALLEL_OPTS < filter.commands.sh
elif command -v xargs > /dev/null 2>&1; then
    cat filter.commands.sh | xargs -I {} -P "$NUM_THREADS" bash -c '{}'
else
    echo "Warning: Running sequentially (parallel/xargs not available)"
    bash filter.commands.sh
fi

# Run summary (getSummary.sh, plus getSVSummary.sh when HP_SV is set; "Summary.sh" matches both)
echo "Generating summary..."
SUMMARY_CMD=$(grep "Summary.sh" run.all.sh || echo "")
if [ -n "$SUMMARY_CMD" ]; then
    eval "$SUMMARY_CMD"
fi

echo "Processing completed successfully"
rm -f run.all.sh filter.commands.sh
EOF

chmod +x "$TEMP_SCRIPT"

# Per-sample resource budgeting (HP_P threads, java -Xmx, samtools sort -m) is computed INSIDE
# the container by scripts/resource_budget.sh, which is cgroup-aware: it reads the container's
# real CPU/RAM limits (not just the host's nproc) and divides them across the NUM_THREADS
# concurrent samples, bounding aggregate CPU and RAM. Here we only forward an explicit user
# HP_P override (respected, but capped to a safe share); otherwise the budget is auto-derived.
HP_P_ENV=""
if [ -n "${HP_P:-}" ]; then
    HP_P_ENV=",HP_P=$HP_P"
    echo "Forwarding user-specified HP_P=$HP_P (resource_budget will cap it to a safe per-sample share)"
else
    echo "Per-sample resources will be auto-budgeted inside the container for $NUM_THREADS concurrent samples"
fi

# Run the container with parallel processing
echo "Executing MitoHPC container..."
# Capture the real exit status: under `set -e` a bare failing apptainer call would abort this
# script before we could read $? (making the error branch below dead). Guard with `|| ...`.
CONTAINER_EXIT_CODE=0
apptainer exec \
    --bind "$WORKING_DIR":"$WORKING_DIR" \
    --pwd "$WORKING_DIR" \
    --env HP_ADIR="$DATA_DIR",HP_ODIR=out,HP_IN=in.txt,HP_SV="$SV_MODE"$HP_P_ENV \
    "$CONTAINER_IMAGE" \
    "./$(basename "$TEMP_SCRIPT")" "$NUM_THREADS" || CONTAINER_EXIT_CODE=$?

# Clean up temporary script
rm -f "$TEMP_SCRIPT"

if [ $CONTAINER_EXIT_CODE -eq 0 ]; then
    echo
    echo "✅ MitoHPC batch processing completed successfully!"
    echo "Results are available in: $WORKING_DIR/out/"
    echo
    echo "Output files include:"
    echo "  - VCF files with variants"
    echo "  - Haplogroup assignments"
    echo "  - Coverage statistics"
    echo "  - Summary reports"
    
    # Show some basic stats if available. Use a nullglob-guarded array: `[ -f "$dir/"*.summary ]`
    # breaks when the glob matches multiple files (too many args to -f) or matches none.
    shopt -s nullglob
    SUMMARY_FILES=("$WORKING_DIR/out/"*.summary)
    shopt -u nullglob
    if [ ${#SUMMARY_FILES[@]} -gt 0 ]; then
        echo
        echo "Summary files created:"
        ls -la "${SUMMARY_FILES[@]}"
    fi
else
    echo "❌ MitoHPC batch processing failed with exit code $CONTAINER_EXIT_CODE" >&2
    exit $CONTAINER_EXIT_CODE
fi
