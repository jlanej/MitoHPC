#!/bin/bash

# Helper script to check memory allocation settings
# Usage: ./check_memory_allocation.sh [threads]

if [ "$1" ]; then
    export HP_P=$1
    echo "Checking memory allocation with HP_P=$HP_P threads (user-specified)"
else
    echo "Checking memory allocation with auto-detected threads"
fi

# Set HP_SDIR if not already set
if [ -z "$HP_SDIR" ]; then
    export HP_SDIR="$(dirname "$(readlink -f "$0")")"
fi

# Source the init script
source "$HP_SDIR/init.sh"

echo "=========================================="
echo "Memory Allocation Summary:"
echo "=========================================="
echo "Threads (HP_P): $HP_P"
echo "Total Memory (HP_MM_TOTAL): $HP_MM_TOTAL"
echo "Per-thread Memory (HP_MM): $HP_MM"
echo ""
echo "Usage breakdown:"
echo "  • Job schedulers (SLURM/SGE): $HP_MM_TOTAL"
echo "  • samtools sort per thread: $HP_MM"
echo "  • samtools sort total usage: $HP_P × $HP_MM = $((HP_P * 2))G"
echo "  • Java heap size: $HP_MM_TOTAL"
echo ""
echo "Example samtools sort command:"
echo "  samtools sort -m $HP_MM -@ $HP_P input.bam"
echo ""
echo "Memory safety check: $([ $((HP_P * 2)) -eq ${HP_MM_TOTAL%G} ] && echo "✅ SAFE" || echo "❌ ISSUE")"