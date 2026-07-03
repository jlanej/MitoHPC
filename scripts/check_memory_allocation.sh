#!/bin/bash

# Helper to inspect MitoHPC memory/CPU allocation.
#
# Usage:
#   ./check_memory_allocation.sh [threads] [concurrent_samples]
#     threads             optional HP_P override (per-sample threads); else init.sh auto-detects
#     concurrent_samples  optional N — if given, also show what resource_budget.sh would assign
#                         when N samples run in parallel (the budgeted, machine-aware allocation
#                         used by mitohpc-parallel.sh / mitohpc-batch-container.sh)

if [ "$1" ]; then
    export HP_P=$1
    echo "Checking memory allocation with HP_P=$HP_P threads (user-specified)"
else
    echo "Checking memory allocation with auto-detected threads"
fi
N_SAMPLES="${2:-}"

# Set HP_SDIR if not already set
if [ -z "$HP_SDIR" ]; then
    export HP_SDIR="$(dirname "$(readlink -f "$0")")"
fi

# Source the init script (per-sample defaults)
source "$HP_SDIR/init.sh"

# Parse the per-thread sort buffer (e.g. "2G" -> 2) so the breakdown reflects the actual HP_MM,
# not a hard-coded 2. samtools `sort -m` is PER THREAD, so the per-sample sort total is HP_P*HP_MM.
MM_G=${HP_MM%[Gg]}
SORT_TOTAL=$(( HP_P * MM_G ))

echo "=========================================="
echo "Per-sample allocation (init.sh defaults):"
echo "=========================================="
echo "Threads (HP_P): $HP_P"
echo "Total Memory (HP_MM_TOTAL): $HP_MM_TOTAL"
echo "Per-thread Memory (HP_MM): $HP_MM"
echo ""
echo "Usage breakdown (one sample):"
echo "  • Job schedulers (SLURM/SGE): $HP_MM_TOTAL"
echo "  • samtools sort per thread: $HP_MM"
echo "  • samtools sort total usage: $HP_P × ${MM_G}G = ${SORT_TOTAL}G"
echo "  • Java heap size (HP_JOPT): $HP_JOPT"
echo ""
echo "Example samtools sort command:"
echo "  samtools sort -m $HP_MM -@ $HP_P input.bam"
echo ""

# Real per-sample consistency note: the sort buffer total should fit within the per-sample
# memory budget HP_MM_TOTAL (which also bounds the java heap). This is a genuine check — unlike
# a tautology, it can FAIL if HP_MM/HP_P/HP_MM_TOTAL are ever set inconsistently.
if [ "$SORT_TOTAL" -le "${HP_MM_TOTAL%[Gg]}" ]; then
    echo "Per-sample sort buffer (${SORT_TOTAL}G) fits within HP_MM_TOTAL ($HP_MM_TOTAL): ✅ OK"
else
    echo "Per-sample sort buffer (${SORT_TOTAL}G) EXCEEDS HP_MM_TOTAL ($HP_MM_TOTAL): ❌ ISSUE"
fi

# If a concurrency level N was given, show the machine-aware budget that the parallel wrappers
# actually apply (resource_budget divides the detected CPU/RAM across N concurrent samples).
if [ -n "$N_SAMPLES" ]; then
    echo ""
    echo "=========================================="
    echo "Budgeted allocation for N=$N_SAMPLES concurrent sample(s):"
    echo "=========================================="
    if [ -r "$HP_SDIR/resource_budget.sh" ]; then
        . "$HP_SDIR/resource_budget.sh"
        resource_budget "$N_SAMPLES"
        echo ""
        echo "Resulting per-sample knobs: HP_P=$HP_P  HP_MM=$HP_MM  HP_MM_TOTAL=$HP_MM_TOTAL"
        echo "  HP_JOPT=$HP_JOPT"
    else
        echo "resource_budget.sh not found next to this script — cannot preview the budget." >&2
    fi
fi
