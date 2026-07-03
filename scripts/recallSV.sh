#!/usr/bin/env bash
set -eu
#########################################################################################
# Re-run ONLY the mtDNA SV caller (+ optional samplot visualization + cohort report) on the
# PERSISTED circular-aware BAMs ($O.sv.bam, kept during a full run via HP_SV_KEEPBAM),
# skipping the expensive subsample -> realign -> SNV-calling front end of filter.sh.
#
# Use this after tuning the SV caller / SVCONF-SVIMPACT weights / samplot filters to refresh
# the sv.* outputs in MINUTES instead of a multi-hour full re-run. It never touches the frozen
# SNV/CN/haplogroup deliverables.
#
#   recallSV.sh [ODIR] [JOBS]     ODIR defaults to $HP_ODIR (or ./out); JOBS defaults to nproc
#
# Needs the persisted *.sv.bam (run the pipeline once with HP_SV_KEEPBAM=1 first), python3+pysam,
# bcftools/bgzip/tabix (for the cohort report), and samplot (optional, for HP_SV_PLOT).
#########################################################################################

SDIR=${HP_SDIR:-$(cd "$(dirname "$0")" && pwd)}
ODIR=${1:-${HP_ODIR:-out}}
JOBS=${2:-$(nproc 2>/dev/null || echo 4)}
test -d "$ODIR"

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
in=$tmp/in.txt; cmds=$tmp/cmds.sh; : > "$in"; : > "$cmds"

# enumerate persisted BAMs -> (sample, bam, output-prefix); rebuild the in.txt getSVSummary needs
while IFS= read -r bam; do
  pref=${bam%.sv.bam}; s=$(basename "$pref")
  printf '%s\t%s\t%s\n' "$s" "$bam" "$pref" >> "$in"
  printf 'bash %q %q %q %q\n' "$SDIR/callSV.sh" "$s" "$bam" "$pref" >> "$cmds"
done < <(find "$ODIR" -name '*.sv.bam' | sort)

if [ ! -s "$in" ]; then
  echo "[recallSV] no persisted *.sv.bam found under $ODIR." >&2
  echo "[recallSV] re-run the pipeline once with HP_SV_KEEPBAM=1 to preserve the alignments first." >&2
  exit 1
fi
n=$(wc -l < "$in")
echo "[recallSV] re-calling SV on $n persisted BAM(s) under $ODIR with $JOBS job(s)" >&2

if command -v parallel >/dev/null 2>&1; then
  parallel -j "$JOBS" --joblog "$tmp/par.log" < "$cmds"
else
  xargs -P "$JOBS" -I {} bash -c '{}' < "$cmds"
fi

# cohort aggregation + interactive report (reuses the standard summary path)
HP_IN=$in HP_ODIR=$ODIR bash "$SDIR/getSVSummary.sh" "$ODIR"
echo "[recallSV] done -> $ODIR/sv.{tab,merged.vcf.gz,sites.vcf.gz,report.html}" >&2
