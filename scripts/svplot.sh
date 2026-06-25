#!/usr/bin/env bash
set -eu
#########################################################################################
# Generate samplot PNGs for the VISUALIZABLE subset of a sample's SV calls, on the LIVE
# $O.bam. Invoked by callSV.sh while the BAM is still alive (filter.sh deletes it right
# after); writes ${O}.sv.<bp5>_<end>.png per selected call and appends the manifest
# ${O}.sv.plots.tsv (sample, bp5, end, afc, svlen, svconf, png, filter, flags) that getSVSummary.sh
# collects and svReport.py embeds (png is column 7 for back-compat; filter/flags are appended).
#
# DEFAULT-OFF: only runs when HP_SV_PLOT is non-empty (mitohpc-batch-container.sh turns it
# on; the bare pipeline stays unchanged). Degrades gracefully if samplot is unavailable.
#
# Configurable filter (env, all optional):
#   HP_SV_PLOT_PASS=1        only FILTER=PASS calls (set 0 to include all)
#   HP_SV_PLOT_MINAF=0.03    minimum coverage-dosage heteroplasmy (AFC)
#   HP_SV_PLOT_SKIP=HP,DLOOP,NUMT   skip calls whose FP-flags contain ANY of these
#                                   (the dominant control-region/homopolymer/NUMT artifacts;
#                                    REPEAT is NOT skipped — the real common deletion carries it)
#   HP_SV_PLOT_MINSVCONF=    optional extra floor on the confidence score (empty = off)
#   HP_SV_PLOT_MAX=200       safety cap on plots per sample
#   HP_SV_PLOT_ALL=1         escape hatch: plot EVERY call (flips the PASS/MINAF/SKIP DEFAULTS to
#                            permissive so nothing is filtered out; an explicitly-set PASS/MINAF/SKIP
#                            still wins). Pairs with getSVSummary.sh defaulting --plot-dedup to 0, so
#                            the report embeds them all. Rarely needed; used for the example report.
#
# Args:  1: sample   2: BAM (live $O.bam)   3: output prefix $O
#########################################################################################

S=$1; BAM=$2; O=$3
MT=${HP_MT:-chrM}
SDIR=${HP_SDIR:-$(cd "$(dirname "$0")" && pwd)}
RDIR=${HP_RDIR:-$(cd "$SDIR/../RefSeq" && pwd)}
SAMPLOT=${HP_SAMPLOT:-samplot}
TAB="$O.sv.tab"
MAN="$O.sv.plots.tsv"

# Annotation track(s) drawn under each plot (samplot -A; each must be a bgzipped+tabixed BED under
# $RDIR). Default: genes.bed.gz — labels the gene/tRNA/rRNA each deletion removes (the OXPHOS complex
# is implicit in the gene name, e.g. ND*/COX*/CYTB/ATP*). Configurable comma-separated list; e.g.
# HP_SV_PLOT_ANNOT=genes.bed.gz,CDS.bed.gz . Empty disables the track.
ANNOT=${HP_SV_PLOT_ANNOT-genes.bed.gz}
aopt=""
if [ -n "$ANNOT" ]; then
  oldifs=$IFS; IFS=','
  for a in $ANNOT; do
    af="$RDIR/$a"
    if [ -s "$af" ] && { [ -s "$af.tbi" ] || [ -s "$af.csi" ]; }; then aopt="$aopt -A $af"; fi
  done
  IFS=$oldifs
fi

[ -n "${HP_SV_PLOT:-}" ] || exit 0                 # default-off
[ -s "$TAB" ] || exit 0                            # nothing called
test -s "$BAM"
if ! command -v "$SAMPLOT" >/dev/null 2>&1; then
  echo "[svplot] WARNING: '$SAMPLOT' not on PATH — skipping SV visualization for $S" >&2
  exit 0
fi

# HP_SV_PLOT_ALL=1 plots EVERY call: it only changes the DEFAULTS of PASS/MINAF/SKIP to permissive,
# so an explicitly-set HP_SV_PLOT_PASS/MINAF/SKIP still overrides (e.g. ALL=1 + HP_SV_PLOT_SKIP=NUMT
# plots everything except NUMT). The `-` (not `:-`) on SKIP honors an explicit empty value.
if [ -n "${HP_SV_PLOT_ALL:-}" ]; then
  PASS=${HP_SV_PLOT_PASS:-0}
  MINAF=${HP_SV_PLOT_MINAF:-0}
  SKIP=${HP_SV_PLOT_SKIP-}
else
  PASS=${HP_SV_PLOT_PASS:-1}
  MINAF=${HP_SV_PLOT_MINAF:-0.03}
  SKIP=${HP_SV_PLOT_SKIP:-HP,DLOOP,NUMT}
fi
MINSVCONF=${HP_SV_PLOT_MINSVCONF:-}
MAX=${HP_SV_PLOT_MAX:-200}

: > "$MAN"
# select visualizable rows from the tidy tab (parse by column NAME); emit: bp5 end afc svlen svconf filter flags.
# filter/flags are carried so the report can label each plot (esp. in HP_SV_PLOT_ALL mode, where non-PASS /
# WRAP / artifact calls are shown and a viewer needs to know which is which — a WRAP/origin call renders as a
# misleading genome-spanning event under samplot's linear view).
awk -F'\t' -v pass="$PASS" -v minaf="$MINAF" -v skip="$SKIP" -v minsc="$MINSVCONF" -v mx="$MAX" '
  NR==1{ for(i=1;i<=NF;i++) h[$i]=i; next }
  { n++
    if (pass==1 && $h["filter"]!="PASS") next
    if ($h["af_coverage"]+0 < minaf+0) next
    if (minsc!="" && $h["svconf"]!="." && $h["svconf"]+0 < minsc+0) next
    fl="," $h["flags"] ","
    m=split(skip,sk,","); drop=0
    for (j=1;j<=m;j++) if (sk[j]!="" && index(fl, "," sk[j] ",")) { drop=1; break }
    if (drop) next
    if (++k > mx) next
    print $h["pos_bp5"] "\t" $h["end_bp3"] "\t" $h["af_coverage"] "\t" $h["svlen"] "\t" $h["svconf"] "\t" $h["filter"] "\t" $h["flags"]
  }' "$TAB" | while IFS=$'\t' read -r bp5 end afc svlen svconf filter flags; do
  png="$O.sv.${bp5}_${end}.png"
  title="$S  m.$((bp5+1))_${end}del  VAF=$afc"
  if "$SAMPLOT" plot -n "$title" -b "$BAM" -o "$png" -c "$MT" -s "$bp5" -e "$end" -t DEL $aopt >/dev/null 2>&1 \
     && [ -s "$png" ]; then
    # png stays column 7 (CI/back-compat); filter/flags appended after it
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$S" "$bp5" "$end" "$afc" "$svlen" "$svconf" "$png" "${filter:-.}" "${flags:-.}" >> "$MAN"
  else
    echo "[svplot] WARNING: samplot failed for $S m.${bp5}_${end}" >&2
  fi
done

echo "[svplot] $S -> $(grep -c . "$MAN" 2>/dev/null || echo 0) plot(s)" >&2
