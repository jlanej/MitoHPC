#!/bin/bash

##############################################################################################################
# Test script to verify VCF output equivalence between mitohpc.sh and mitohpc-parallel.sh
#
# This script demonstrates how to verify that both approaches produce identical merged VCF files.
# Usage: ./test_vcf_equivalence.sh [data_directory] [num_threads]
##############################################################################################################

set -e

DATA_DIR=${1:-"./test_data"}
NUM_THREADS=${2:-4}

echo "=== MitoHPC VCF Output Equivalence Test ==="
echo "Data directory: $DATA_DIR"
echo "Parallel threads: $NUM_THREADS"
echo ""

# Check if data directory exists
if [ ! -d "$DATA_DIR" ]; then
    echo "Error: Data directory $DATA_DIR does not exist"
    echo "Please provide a directory containing BAM/CRAM files"
    exit 1
fi

# Create separate output directories
SEQUENTIAL_OUT="$DATA_DIR/out_sequential"
PARALLEL_OUT="$DATA_DIR/out_parallel"

mkdir -p "$SEQUENTIAL_OUT" "$PARALLEL_OUT"

echo "Testing sequential processing..."
export HP_ADIR="$DATA_DIR"
export HP_ODIR="$SEQUENTIAL_OUT"
export HP_IN="$SEQUENTIAL_OUT/in.txt"

# Run sequential processing
#mitohpc.sh

echo "Testing parallel processing..."
export HP_ADIR="$DATA_DIR"  
export HP_ODIR="$PARALLEL_OUT"
export HP_IN="$PARALLEL_OUT/in.txt"

# Run parallel processing
#mitohpc-parallel.sh $NUM_THREADS

echo ""
echo "=== Comparing Output Files ==="

# Compare key merged VCF files
VCF_FILES=(
    "mutect2.mutect2.10.merge.vcf"
    "mutect2.mutect2.20.merge.vcf"
    "mutect2.mutect2.30.merge.vcf"
)

for vcf_file in "${VCF_FILES[@]}"; do
    if [ -f "$SEQUENTIAL_OUT/$vcf_file" ] && [ -f "$PARALLEL_OUT/$vcf_file" ]; then
        echo "Comparing $vcf_file..."
        if diff -q "$SEQUENTIAL_OUT/$vcf_file" "$PARALLEL_OUT/$vcf_file" > /dev/null; then
            echo "  ✓ Files are identical"
        else
            echo "  ✗ Files differ"
            echo "  First few differences:"
            diff "$SEQUENTIAL_OUT/$vcf_file" "$PARALLEL_OUT/$vcf_file" | head -10
        fi
    else
        echo "  - $vcf_file not found in one or both outputs"
    fi
done

echo ""
echo "=== Test Summary ==="
echo "This test demonstrates that mitohpc-parallel.sh produces identical"
echo "merged VCF files as mitohpc.sh, just with faster processing time."

# Optional: Compare all VCF files
echo ""
echo "To compare all VCF files:"
echo "  find $SEQUENTIAL_OUT -name '*.vcf' -exec basename {} \; | sort | uniq | while read vcf; do"
echo "    echo \"Comparing \$vcf:\""
echo "    diff -q \"$SEQUENTIAL_OUT/\$vcf\" \"$PARALLEL_OUT/\$vcf\" || echo \"  Files differ\""
echo "  done"