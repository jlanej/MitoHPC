#!/usr/bin/env bash

##############################################################################################################
# MitoHPC Batch Processing Script
# 
# This script allows easy batch processing of multiple sample directories using MitoHPC
# across multiple parallel threads/processes using Apptainer containers.
#
# Usage:
#   mitohpc-batch.sh -d <base_directory> -j <num_jobs> [options]
#
# Example:
#   mitohpc-batch.sh -d /data/samples -j 4 -c "docker://ghcr.io/jlanej/mitohpc:main"
#
##############################################################################################################

set -e

# Default values
NUM_JOBS=1
CONTAINER_IMAGE="docker://ghcr.io/jlanej/mitohpc:main"
BASE_DIR=""
OUTPUT_BASE=""
VERBOSE=false
DRY_RUN=false

# Function to show usage
show_usage() {
    cat << EOF
Usage: $0 -d <base_directory> [-j <num_jobs>] [-c <container_image>] [-o <output_base>] [-v] [-n] [-h]

Options:
    -d <base_directory>    Base directory containing subdirectories with BAM/CRAM files
    -j <num_jobs>         Number of parallel jobs to run (default: 1)
    -c <container_image>  Container image to use (default: docker://ghcr.io/jlanej/mitohpc:main)
    -o <output_base>      Base output directory (default: <base_directory>/batch_output)
    -v                    Verbose output
    -n                    Dry run - show commands without executing
    -h                    Show this help message

Example:
    $0 -d /data/samples -j 4 -c "docker://ghcr.io/jlanej/mitohpc:main"

This will:
1. Find all subdirectories in /data/samples that contain BAM/CRAM files
2. Process them in parallel using 4 jobs
3. Each job runs MitoHPC in a separate container instance
4. Results are collected in /data/samples/batch_output/

Directory structure expected:
    /data/samples/
    ├── sample1/
    │   ├── bams/
    │   │   ├── file1.bam
    │   │   └── file2.bam
    ├── sample2/
    │   ├── bams/
    │   │   └── file3.cram
    └── ...

EOF
}

# Function to log messages
log() {
    if [ "$VERBOSE" = true ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
    fi
}

# Function to find sample directories
find_sample_dirs() {
    local base_dir="$1"
    find "$base_dir" -mindepth 1 -maxdepth 2 -type d -name "bams" -o -name "crams" | \
        xargs -I {} dirname {} | sort -u | \
        while read dir; do
            # Check if directory actually contains BAM/CRAM files
            if find "$dir" -name "*.bam" -o -name "*.cram" 2>/dev/null | head -1 | grep -q .; then
                echo "$dir"
            fi
        done
}

# Function to process a single sample directory
process_sample() {
    local sample_dir="$1"
    local output_dir="$2"
    local container_image="$3"
    
    local sample_name=$(basename "$sample_dir")
    local sample_output="$output_dir/$sample_name"
    
    log "Processing sample: $sample_name"
    log "  Input directory: $sample_dir"
    log "  Output directory: $sample_output"
    
    # Create output directory
    mkdir -p "$sample_output"
    
    # Determine if we have bams or crams directory
    local data_dir=""
    if [ -d "$sample_dir/bams" ]; then
        data_dir="bams"
    elif [ -d "$sample_dir/crams" ]; then
        data_dir="crams"
    else
        echo "Error: No bams or crams directory found in $sample_dir" >&2
        return 1
    fi
    
    # Create input file list
    local input_file="$sample_output/in.txt"
    find "$sample_dir/$data_dir" -name "*.bam" -o -name "*.cram" | sort -V > "$input_file"
    
    if [ ! -s "$input_file" ]; then
        echo "Error: No BAM/CRAM files found in $sample_dir/$data_dir" >&2
        return 1
    fi
    
    # Calculate appropriate HP_P value to prevent resource conflicts
    # Each sample should get a fair share of CPU cores
    local total_cores=$(nproc 2>/dev/null || echo 4)
    local hp_p_env=""
    if [ -z "$HP_P" ]; then
        # Estimate number of parallel jobs from the parent process
        # This is an approximation since we don't have direct access to the parallel job count
        local estimated_parallel_jobs=${NUM_JOBS:-4}
        local hp_p_calculated=$((total_cores / estimated_parallel_jobs))
        if [ $hp_p_calculated -lt 1 ]; then
            hp_p_calculated=1
        fi
        hp_p_env="HP_P=$hp_p_calculated"
        log "Setting HP_P=$hp_p_calculated for sample $sample_name"
    else
        hp_p_env="HP_P=$HP_P"
        log "Using user-specified HP_P=$HP_P for sample $sample_name"
    fi
    
    # Run MitoHPC in container
    local cmd="apptainer exec \
        --bind \"$sample_dir\":\"$sample_dir\" \
        --bind \"$sample_output\":\"$sample_output\" \
        --pwd \"$sample_dir\" \
        --env HP_ADIR=\"$data_dir\",HP_ODIR=\"$sample_output/out\",HP_IN=\"$input_file\",$hp_p_env \
        \"$container_image\" \
        mitohpc.sh"
    
    if [ "$DRY_RUN" = true ]; then
        echo "Would execute: $cmd"
    else
        log "Executing: $cmd"
        eval "$cmd" 2>&1 | while IFS= read -r line; do
            log "[$sample_name] $line"
        done
        
        if [ $? -eq 0 ]; then
            log "Successfully completed processing sample: $sample_name"
        else
            echo "Error processing sample: $sample_name" >&2
            return 1
        fi
    fi
}

# Function to run batch processing
run_batch() {
    local base_dir="$1"
    local num_jobs="$2"
    local container_image="$3"
    local output_base="$4"
    
    log "Starting batch processing"
    log "  Base directory: $base_dir"
    log "  Number of jobs: $num_jobs"
    log "  Container image: $container_image"
    log "  Output base: $output_base"
    
    # Find all sample directories
    local sample_dirs=($(find_sample_dirs "$base_dir"))
    
    if [ ${#sample_dirs[@]} -eq 0 ]; then
        echo "Error: No sample directories with BAM/CRAM files found in $base_dir" >&2
        exit 1
    fi
    
    log "Found ${#sample_dirs[@]} sample directories to process:"
    for dir in "${sample_dirs[@]}"; do
        log "  - $(basename "$dir")"
    done
    
    # Create output base directory
    mkdir -p "$output_base"
    
    # Export function for parallel execution
    export -f process_sample log
    export VERBOSE DRY_RUN NUM_JOBS="$num_jobs"
    
    # Use GNU parallel or xargs for parallel processing
    if command -v parallel > /dev/null 2>&1; then
        log "Using GNU parallel for batch processing"
        printf '%s\n' "${sample_dirs[@]}" | \
            parallel -j "$num_jobs" process_sample {} "$output_base" "$container_image"
    else
        log "GNU parallel not available, using xargs"
        printf '%s\n' "${sample_dirs[@]}" | \
            xargs -I {} -P "$num_jobs" bash -c 'process_sample "$@"' _ {} "$output_base" "$container_image"
    fi
    
    local exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        log "Batch processing completed successfully"
        echo "Results available in: $output_base"
    else
        echo "Batch processing completed with errors" >&2
        exit $exit_code
    fi
}

# Main script
main() {
    # Parse command line arguments
    while getopts "d:j:c:o:vnh" opt; do
        case $opt in
            d)
                BASE_DIR="$OPTARG"
                ;;
            j)
                NUM_JOBS="$OPTARG"
                ;;
            c)
                CONTAINER_IMAGE="$OPTARG"
                ;;
            o)
                OUTPUT_BASE="$OPTARG"
                ;;
            v)
                VERBOSE=true
                ;;
            n)
                DRY_RUN=true
                ;;
            h)
                show_usage
                exit 0
                ;;
            \?)
                echo "Invalid option: -$OPTARG" >&2
                show_usage
                exit 1
                ;;
        esac
    done
    
    # Validate required arguments
    if [ -z "$BASE_DIR" ]; then
        echo "Error: Base directory (-d) is required" >&2
        show_usage
        exit 1
    fi
    
    if [ ! -d "$BASE_DIR" ]; then
        echo "Error: Base directory does not exist: $BASE_DIR" >&2
        exit 1
    fi
    
    # Set default output base if not provided
    if [ -z "$OUTPUT_BASE" ]; then
        OUTPUT_BASE="$BASE_DIR/batch_output"
    fi
    
    # Validate number of jobs
    if ! [[ "$NUM_JOBS" =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: Number of jobs must be a positive integer" >&2
        exit 1
    fi
    
    # Check if apptainer is available
    if ! command -v apptainer > /dev/null 2>&1; then
        echo "Error: apptainer command not found. Please install Apptainer/Singularity" >&2
        exit 1
    fi
    
    # Run batch processing
    run_batch "$BASE_DIR" "$NUM_JOBS" "$CONTAINER_IMAGE" "$OUTPUT_BASE"
}

# Run main function with all arguments
main "$@"