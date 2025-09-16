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

# Function to show usage
show_usage() {
    cat << EOF
Usage: $0 <working_directory> [num_threads] [container_image]

Arguments:
    working_directory    Directory containing BAM/CRAM files (equivalent to groupCram)
    num_threads         Number of parallel threads to use (default: auto-detect CPU cores)
    container_image     Container image to use (default: docker://ghcr.io/jlanej/mitohpc:main)

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

# Determine data directory
DATA_DIR="bams"
# if [ -d "$WORKING_DIR/bams" ]; then
#     DATA_DIR="bams"
#     if ! find "$WORKING_DIR/bams" -name "*.bam" | head -1 | grep -q .; then
#         echo "Error: No BAM files found in $WORKING_DIR/bams" >&2
#         exit 1
#     fi
# elif [ -d "$WORKING_DIR/crams" ]; then
#     DATA_DIR="crams"
#     if ! find "$WORKING_DIR/crams" -name "*.cram" | head -1 | grep -q .; then
#         echo "Error: No CRAM files found in $WORKING_DIR/crams" >&2
#         exit 1
#     fi
# fi

echo "MitoHPC Batch Container Processing"
echo "=================================="
echo "Working directory: $WORKING_DIR"
echo "Data directory: $DATA_DIR"
echo "Number of threads: $NUM_THREADS"
echo "Container image: $CONTAINER_IMAGE"
echo

# Count input files
FILE_COUNT=$(find "$WORKING_DIR/$DATA_DIR" -name "*.bam" -o -name "*.cram" | wc -l)
echo "Found $FILE_COUNT input files to process"

# Create output directory
mkdir -p "$WORKING_DIR/out"

# Run MitoHPC with parallel processing
# We'll use a modified approach that processes samples in parallel within the container
echo "Starting parallel processing..."

# Create a temporary script that will run inside the container
TEMP_SCRIPT="$WORKING_DIR/run_parallel_mitohpc.sh"
cat > "$TEMP_SCRIPT" << 'EOF'
#!/usr/bin/env bash
set -e

NUM_THREADS="$1"

# Source the MitoHPC initialization
. $HP_SDIR/init.sh

# Generate input file
echo "Generating input file list..."
find $HP_ADIR/ -name "*.bam" -o -name "*.cram" -readable | ls2in.pl -out $HP_ODIR | sort -V > $HP_IN

if [ ! -s "$HP_IN" ]; then
    echo "Error: No input files found" >&2
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
    parallel -j "$NUM_THREADS" --progress < filter.commands.sh
elif command -v xargs > /dev/null 2>&1; then
    cat filter.commands.sh | xargs -I {} -P "$NUM_THREADS" bash -c '{}'
else
    echo "Warning: Running sequentially (parallel/xargs not available)"
    bash filter.commands.sh
fi

# Run summary
echo "Generating summary..."
SUMMARY_CMD=$(grep "getSummary.sh" run.all.sh || echo "")
if [ -n "$SUMMARY_CMD" ]; then
    eval "$SUMMARY_CMD"
fi

echo "Processing completed successfully"
rm -f run.all.sh filter.commands.sh
EOF

chmod +x "$TEMP_SCRIPT"

# Calculate appropriate HP_P value to prevent resource conflicts
# If HP_P is not explicitly set by the user, we should set it to a reasonable value
# to prevent each parallel job from trying to use all available cores
TOTAL_CORES=$(nproc 2>/dev/null || echo 4)
if [ -z "$HP_P" ]; then
    # Calculate threads per sample: max(1, total_cores / num_parallel_samples)
    # This ensures each sample gets a fair share without oversubscription
    HP_P_CALCULATED=$((TOTAL_CORES / NUM_THREADS))
    if [ $HP_P_CALCULATED -lt 1 ]; then
        HP_P_CALCULATED=1
    fi
    HP_P_ENV="HP_P=$HP_P_CALCULATED"
    echo "Setting HP_P=$HP_P_CALCULATED threads per sample (total cores: $TOTAL_CORES, parallel samples: $NUM_THREADS)"
else
    HP_P_ENV="HP_P=$HP_P"
    echo "Using user-specified HP_P=$HP_P threads per sample"
fi

# Run the container with parallel processing
echo "Executing MitoHPC container..."
apptainer exec \
    --bind "$WORKING_DIR":"$WORKING_DIR" \
    --pwd "$WORKING_DIR" \
    --env HP_ADIR="$DATA_DIR",HP_ODIR=out,HP_IN=in.txt,"$HP_P_ENV" \
    "$CONTAINER_IMAGE" \
    ./run_parallel_mitohpc.sh "$NUM_THREADS"

CONTAINER_EXIT_CODE=$?

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
    
    # Show some basic stats if available
    if [ -f "$WORKING_DIR/out/"*.summary ]; then
        echo
        echo "Summary files created:"
        ls -la "$WORKING_DIR/out/"*.summary 2>/dev/null || true
    fi
else
    echo "❌ MitoHPC batch processing failed with exit code $CONTAINER_EXIT_CODE" >&2
    exit $CONTAINER_EXIT_CODE
fi
