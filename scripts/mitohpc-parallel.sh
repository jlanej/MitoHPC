#!/usr/bin/env bash

##############################################################################################################
# MitoHPC Parallel Processing Script
# 
# This script provides an easy way to run MitoHPC with parallel processing of samples.
# It works with the existing MitoHPC container approach but processes samples across multiple threads.
#
# Usage:
#   mitohpc-parallel.sh [num_threads] [container_image]
#
# Example:
#   mitohpc-parallel.sh 4 "docker://ghcr.io/jlanej/mitohpc:main"
#
# This script expects to be run in the same way as the original mitohpc.sh, but will
# process the samples in parallel using the specified number of threads.
##############################################################################################################

set -e

# Get number of threads from command line, default to number of CPU cores
NUM_THREADS=${1:-$(nproc 2>/dev/null || echo 4)}
CONTAINER_IMAGE=${2:-"docker://ghcr.io/jlanej/mitohpc:main"}

# Validate number of threads
if ! [[ "$NUM_THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: Number of threads must be a positive integer" >&2
    echo "Usage: $0 [num_threads] [container_image]" >&2
    echo "Example: $0 4 \"docker://ghcr.io/jlanej/mitohpc:main\"" >&2
    exit 1
fi

echo "Starting MitoHPC parallel processing with $NUM_THREADS threads"
echo "Container image: $CONTAINER_IMAGE"

# Check if apptainer is available
if ! command -v apptainer > /dev/null 2>&1; then
    echo "Error: apptainer command not found. Please install Apptainer/Singularity" >&2
    exit 1
fi

# Source the MitoHPC initialization to get HP_* variables
# This mimics what mitohpc.sh does but we'll handle the parallel execution
if [ -z "$HP_SDIR" ]; then
    echo "Error: HP_SDIR not set. Please ensure MitoHPC environment variables are set." >&2
    echo "This script should be called from within the MitoHPC container or with HP_* variables set." >&2
    exit 1
fi

# Source init.sh to get all HP_ variables
. $HP_SDIR/init.sh

echo "MitoHPC parallel configuration:"
echo "  Per-sample threads (HP_P): $HP_P"
echo "  Parallel sample processing threads: $NUM_THREADS"
echo "  Total potential thread usage: $((HP_P * NUM_THREADS))"

# If total thread usage seems excessive, warn the user
TOTAL_CORES=$(nproc 2>/dev/null || echo 4)
POTENTIAL_THREADS=$((HP_P * NUM_THREADS))
if [ $POTENTIAL_THREADS -gt $((TOTAL_CORES * 2)) ]; then
    echo "Warning: Potential thread oversubscription detected!"
    echo "  Available cores: $TOTAL_CORES"
    echo "  Potential threads: $POTENTIAL_THREADS"
    echo "  Consider reducing NUM_THREADS or setting HP_P to a lower value"
fi

# Generate input file if it doesn't exist
if [ ! -s "$HP_IN" ]; then
    echo "Generating input file: $HP_IN"
    find $HP_ADIR/ -name "*.bam" -o -name "*.cram" -readable | ls2in.pl -out $HP_ODIR | sort -V > $HP_IN
fi

# Check if input file has samples
if [ ! -s "$HP_IN" ]; then
    echo "Error: No BAM/CRAM files found in $HP_ADIR" >&2
    exit 1
fi

# Count total samples
TOTAL_SAMPLES=$(wc -l < "$HP_IN")
echo "Found $TOTAL_SAMPLES samples to process"

# Generate the run script
echo "Generating processing commands..."
$HP_SDIR/run.sh > run.all.sh

# Extract filter.sh commands (these process individual samples)
grep "filter.sh" run.all.sh > filter.commands.sh 2>/dev/null || {
    echo "Error: No filter.sh commands found in run.all.sh" >&2
    exit 1
}

FILTER_COMMANDS=$(wc -l < filter.commands.sh)
echo "Generated $FILTER_COMMANDS processing commands"

# Function to run commands in parallel
run_parallel() {
    if command -v parallel > /dev/null 2>&1; then
        echo "Using GNU parallel to process $FILTER_COMMANDS commands with $NUM_THREADS threads"
        parallel -j "$NUM_THREADS" --joblog parallel.log < filter.commands.sh
    elif command -v xargs > /dev/null 2>&1; then
        echo "Using xargs to process $FILTER_COMMANDS commands with $NUM_THREADS threads"
        cat filter.commands.sh | xargs -I {} -P "$NUM_THREADS" bash -c '{}'
    else
        echo "Warning: Neither GNU parallel nor xargs found. Running sequentially."
        bash filter.commands.sh
    fi
}

# Set up signal handling for cleanup
cleanup() {
    echo "Cleaning up..."
    # Kill any remaining background processes
    jobs -p | xargs -r kill 2>/dev/null || true
    exit 1
}
trap cleanup INT TERM

# Export environment variables that might be needed
export HP_SDIR HP_ODIR HP_IN

# Run the filter commands in parallel
echo "Starting parallel processing..."
START_TIME=$(date +%s)

run_parallel

FILTER_EXIT_CODE=$?

if [ $FILTER_EXIT_CODE -eq 0 ]; then
    echo "Parallel processing completed successfully"
    
    # Run the summary script (this needs to run after all samples are processed)
    echo "Generating summary..."
    SUMMARY_CMD=$(grep "getSummary.sh" run.all.sh)
    if [ -n "$SUMMARY_CMD" ]; then
        eval "$SUMMARY_CMD"
        SUMMARY_EXIT_CODE=$?
        if [ $SUMMARY_EXIT_CODE -eq 0 ]; then
            echo "Summary generation completed successfully"
        else
            echo "Warning: Summary generation failed with exit code $SUMMARY_EXIT_CODE" >&2
        fi
    else
        echo "Warning: No summary command found" >&2
    fi
else
    echo "Error: Parallel processing failed with exit code $FILTER_EXIT_CODE" >&2
    exit $FILTER_EXIT_CODE
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "Total processing time: ${DURATION} seconds"

# Show job log if parallel was used
if [ -f "parallel.log" ]; then
    echo "Parallel job statistics:"
    echo "  Total jobs: $(tail -n +2 parallel.log | wc -l)"
    echo "  Failed jobs: $(tail -n +2 parallel.log | awk '$7 != 0' | wc -l)"
    if [ "$(tail -n +2 parallel.log | awk '$7 != 0' | wc -l)" -gt 0 ]; then
        echo "Failed job details:"
        tail -n +2 parallel.log | awk '$7 != 0 {print "  " $9 " (exit code: " $7 ")"}'
    fi
fi

# Clean up temporary files
rm -f filter.commands.sh run.all.sh

echo "MitoHPC parallel processing completed"