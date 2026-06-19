#!/usr/bin/env bash
set -eu
#########################################################################################
# Aggregate per-sample structural-variant calls into cohort files. SEPARATE from
# getSummary.sh (which is never modified). Invoked only when HP_SV is set.
#
#   $ODIR/sv.concat.vcf   all per-sample $O.sv.vcf records, one header
#   $ODIR/sv.tab          all per-sample $O.sv.tab rows, one header
#
# Usage: getSVSummary.sh [ODIR]   (defaults to $HP_ODIR; uses $HP_IN for sample prefixes)
#########################################################################################

if [ "$#" -lt 1 ]; then ODIR=$HP_ODIR; else ODIR=$1; fi
SDIR=${HP_SDIR:-$(cd "$(dirname "$0")" && pwd)}
test -s "$HP_IN"

# VCF: concatenate, drop per-sample ##sample= lines, dedup, keep header (mirrors getSummary.sh)
awk '{print $3}' "$HP_IN" | sed 's|$|.sv.vcf|' | xargs cat 2>/dev/null | \
  grep -v '^##sample=' | "$SDIR/uniq.pl" | bedtools sort -header > "$ODIR/sv.concat.vcf"

# TAB: one header + all data rows
first=$(awk 'NR==1{print $3; exit}' "$HP_IN").sv.tab
{ [ -s "$first" ] && head -1 "$first"
  awk '{print $3}' "$HP_IN" | sed 's|$|.sv.tab|' | xargs grep -hv '^#' 2>/dev/null || true
} > "$ODIR/sv.tab"

echo "[getSVSummary] wrote $ODIR/sv.concat.vcf and $ODIR/sv.tab" >&2
