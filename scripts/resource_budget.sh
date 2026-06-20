#!/usr/bin/env bash
##############################################################################################################
# resource_budget.sh — cgroup-aware CPU/RAM budgeting for PARALLEL MitoHPC runs.
#
# THE PROBLEM IT SOLVES
#   MitoHPC's per-sample resource knobs are sized for ONE sample:
#       HP_P        threads per sample          (init.sh defaults to `nproc` => ALL cores)
#       HP_MM       samtools `sort -m` buffer   (init.sh defaults to 2G — and it is PER THREAD)
#       HP_JOPT     java -Xms/-Xmx for mutect2  (init.sh defaults to (HP_P*2)G)
#   The wrapper scripts (mitohpc-parallel.sh / mitohpc-batch-container.sh) then fan out N
#   samples concurrently with `parallel -j N`. Nothing bounds the AGGREGATE: N concurrent
#   samples can demand  N*HP_P threads,  N*(HP_P*HP_MM)G of sort buffers, and N*Xmx of java
#   heap — none of it checked against the machine. On a box where init.sh sets HP_P=nproc,
#   a single parallel run can request N*nproc threads and tens of GB of heap and thrash or OOM.
#
# THE FIX
#   `resource_budget N` re-derives the per-sample knobs from the *shared* machine budget so
#   that N concurrent samples fit. It:
#     1. detects usable cores + RAM (cgroup v2/v1 aware — a container under a cgroup CPU/mem
#        limit still reports the *host* via nproc and /proc/meminfo, so naive detection
#        over-commits; we read the cgroup limits when present),
#     2. divides that shared budget across the N concurrent samples, and
#     3. exports HP_P / HP_MM / HP_MM_TOTAL / HP_JOPT accordingly, warning if N is so large
#        that even the budgeted share oversubscribes CPU or RAM.
#
# WHY A SEPARATE, SOURCED-AFTER-init.sh FILE
#   init.sh is a frozen upstream file (do not edit). Sourcing this file AFTER `. init.sh`
#   lets us override the per-sample defaults for the parallel wrappers without touching any
#   pre-existing deliverable. It is purely additive and only ever LOWERS resource requests
#   vs. a naive run. An explicitly user-set HP_P is respected, but capped to a safe share.
#
# USAGE
#   . "$HP_SDIR/resource_budget.sh"      # source it (defines the function)
#   resource_budget "$NUM_SAMPLES"       # N = number of samples run concurrently (parallel -j N)
#   # ... then run.sh / filter.sh inherit the budgeted HP_P / HP_MM / HP_MM_TOTAL / HP_JOPT.
#
#   bash resource_budget.sh 25           # run directly to PREVIEW the budget for N=25 (no side effects elsewhere)
##############################################################################################################

# --- detection helpers (all echo a positive integer and always return 0, so they are safe
# --- to use in `x=$(...)` under `set -e`) -------------------------------------------------

# Usable CPU count: prefer the cgroup CPU quota (containers), else nproc, else 4.
_rb_cores() {
    local q p c
    # cgroup v2: "<quota> <period>" or "max <period>"
    if [ -r /sys/fs/cgroup/cpu.max ]; then
        read -r q p < /sys/fs/cgroup/cpu.max 2>/dev/null || true
        if [ "$q" != "max" ] && [ -n "$q" ] && [ -n "$p" ]; then
            if [ "$q" -gt 0 ] 2>/dev/null && [ "$p" -gt 0 ] 2>/dev/null; then
                c=$(( (q + p - 1) / p ))          # ceil(quota/period)
                if [ "$c" -ge 1 ]; then echo "$c"; return 0; fi
            fi
        fi
    fi
    # cgroup v1
    if [ -r /sys/fs/cgroup/cpu/cpu.cfs_quota_us ] && [ -r /sys/fs/cgroup/cpu/cpu.cfs_period_us ]; then
        q=$(cat /sys/fs/cgroup/cpu/cpu.cfs_quota_us 2>/dev/null || echo -1)
        p=$(cat /sys/fs/cgroup/cpu/cpu.cfs_period_us 2>/dev/null || echo 0)
        if [ "$q" -gt 0 ] 2>/dev/null && [ "$p" -gt 0 ] 2>/dev/null; then
            c=$(( (q + p - 1) / p ))
            if [ "$c" -ge 1 ]; then echo "$c"; return 0; fi
        fi
    fi
    nproc 2>/dev/null || echo 4
}

# Usable RAM in whole GB: min(physical MemTotal, cgroup memory limit), else 8.
_rb_ram_g() {
    local phys=0 cg=0 g=0 kb v bytes=""
    if [ -r /proc/meminfo ]; then
        kb=$(awk '/^MemTotal:/{print $2; exit}' /proc/meminfo 2>/dev/null || echo "")
        if [ -n "$kb" ] && [ "$kb" -gt 0 ] 2>/dev/null; then phys=$(( kb / 1024 / 1024 )); fi
    fi
    # cgroup v2 then v1 (v1 "unlimited" is a ~9.2e18 sentinel — treat as unset)
    if [ -r /sys/fs/cgroup/memory.max ]; then
        v=$(cat /sys/fs/cgroup/memory.max 2>/dev/null || echo max)
        if [ "$v" != "max" ] && [ -n "$v" ]; then bytes=$v; fi
    fi
    if [ -z "$bytes" ] && [ -r /sys/fs/cgroup/memory/memory.limit_in_bytes ]; then
        v=$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes 2>/dev/null || echo "")
        if [ -n "$v" ] && [ "$v" -gt 0 ] 2>/dev/null && [ "$v" -lt 9000000000000000000 ] 2>/dev/null; then
            bytes=$v
        fi
    fi
    if [ -n "$bytes" ] && [ "$bytes" -gt 0 ] 2>/dev/null; then cg=$(( bytes / 1024 / 1024 / 1024 )); fi
    # smallest positive estimate wins; default 8G
    if [ "$phys" -gt 0 ]; then g=$phys; fi
    if [ "$cg" -gt 0 ]; then
        if [ "$g" -eq 0 ] || [ "$cg" -lt "$g" ]; then g=$cg; fi
    fi
    if [ "$g" -lt 1 ]; then g=8; fi
    echo "$g"
}

# --- the budgeter -------------------------------------------------------------------------
# resource_budget N   (N = number of samples processed concurrently; defaults to 1)
resource_budget() {
    local n="${1:-1}"
    case "$n" in ''|*[!0-9]*) n=1;; esac
    if [ "$n" -lt 1 ]; then n=1; fi

    local cores ram_g usable fair_p p perjob mm
    cores=$(_rb_cores)
    ram_g=$(_rb_ram_g)

    # keep ~15% RAM headroom for the OS, page cache, and bwa/bcftools/bedtools overhead
    usable=$(( ram_g * 85 / 100 ))
    if [ "$usable" -lt 1 ]; then usable=1; fi

    # threads per sample = fair share of cores across the N concurrent samples.
    fair_p=$(( cores / n ))
    if [ "$fair_p" -lt 1 ]; then fair_p=1; fi
    p="$fair_p"
    # respect an explicit user HP_P, but never above the fair share (that would oversubscribe)
    if [ -n "${HP_P:-}" ] && [ "${HP_P:-0}" -ge 1 ] 2>/dev/null; then
        if [ "$HP_P" -lt "$fair_p" ]; then p="$HP_P"; fi
    fi

    # per-sample RAM budget (GB) = shared usable pool / N concurrent samples
    perjob=$(( usable / n ))
    if [ "$perjob" -lt 2 ]; then perjob=2; fi      # floor so java + sort have headroom

    # samtools `sort -m` is PER THREAD across HP_P threads, so the per-thread buffer must be
    # perjob/HP_P to keep one sample's sort total within its RAM budget.
    mm=$(( perjob / p ))
    if [ "$mm" -lt 1 ]; then
        mm=1
        # RAM too tight for this many threads: cap threads so sort total (~p*1G) ~ perjob
        if [ "$p" -gt "$perjob" ]; then p="$perjob"; fi
    fi

    export HP_P="$p"
    export HP_MM="${mm}G"
    export HP_MM_TOTAL="${perjob}G"
    export HP_JOPT="-Xms512m -Xmx${perjob}G -XX:ParallelGCThreads=$p"

    echo "[resource_budget] detected ${cores} usable cores, ${ram_g}G RAM; budgeting for N=${n} concurrent sample(s)" >&2
    echo "[resource_budget]   per sample: HP_P=$HP_P threads, HP_MM=$HP_MM/thread sort buffer, java $HP_JOPT" >&2
    echo "[resource_budget]   aggregate target: $(( p * n )) threads, ~$(( perjob * n ))G RAM (of ${usable}G usable)" >&2
    if [ $(( p * n )) -gt "$cores" ]; then
        echo "[resource_budget]   WARNING: $(( p * n )) threads > ${cores} cores — CPU oversubscribed; reduce N." >&2
    fi
    if [ $(( perjob * n )) -gt "$usable" ]; then
        echo "[resource_budget]   WARNING: ~$(( perjob * n ))G > ${usable}G usable — RAM may be tight; reduce N." >&2
    fi
}

# Run directly (not sourced) to preview a budget: `bash resource_budget.sh [N]`
if [ "${BASH_SOURCE[0]:-$0}" = "$0" ]; then
    resource_budget "${1:-1}"
fi
