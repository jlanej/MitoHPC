#!/usr/bin/env bash
set -eu
#########################################################################################
# Generate samplot PNGs for the VISUALIZABLE subset of a sample's SV calls, on the LIVE
# $O.bam. Invoked by callSV.sh while the BAM is still alive (filter.sh deletes it right
# after); writes ${O}.sv.<bp5>_<end>.png per selected call and appends the manifest
# ${O}.sv.plots.tsv (sample, bp5, end, afc, svlen, svconf, png) that getSVSummary.sh
# collects and svReport.py embeds.
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
#
# Args:  1: sample   2: BAM (live $O.bam)   3: output prefix $O
#########################################################################################

S=$1; BAM=$2; O=$3
MT=${HP_MT:-chrM}
SAMPLOT=${HP_SAMPLOT:-samplot}
TAB="$O.sv.tab"
MAN="$O.sv.plots.tsv"

[ -n "${HP_SV_PLOT:-}" ] || exit 0                 # default-off
[ -s "$TAB" ] || exit 0                            # nothing called
test -s "$BAM"
if ! command -v "$SAMPLOT" >/dev/null 2>&1; then
  echo "[svplot] WARNING: '$SAMPLOT' not on PATH — skipping SV visualization for $S" >&2
  exit 0
fi

PASS=${HP_SV_PLOT_PASS:-1}
MINAF=${HP_SV_PLOT_MINAF:-0.03}
SKIP=${HP_SV_PLOT_SKIP:-HP,DLOOP,NUMT}
MINSVCONF=${HP_SV_PLOT_MINSVCONF:-}
MAX=${HP_SV_PLOT_MAX:-200}

: > "$MAN"
# select visualizable rows from the tidy tab (parse by column NAME); emit: bp5 end afc svlen svconf
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
    print $h["pos_bp5"] "\t" $h["end_bp3"] "\t" $h["af_coverage"] "\t" $h["svlen"] "\t" $h["svconf"]
  }' "$TAB" | while IFS=$'\t' read -r bp5 end afc svlen svconf; do
  png="$O.sv.${bp5}_${end}.png"
  title="$S  m.$((bp5+1))_${end}del  VAF=$afc"
  if "$SAMPLOT" plot -n "$title" -b "$BAM" -o "$png" -c "$MT" -s "$bp5" -e "$end" -t DEL >/dev/null 2>&1 \
     && [ -s "$png" ]; then
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$S" "$bp5" "$end" "$afc" "$svlen" "$svconf" "$png" >> "$MAN"
  else
    echo "[svplot] WARNING: samplot failed for $S m.${bp5}_${end}" >&2
  fi
done

echo "[svplot] $S -> $(grep -c . "$MAN" 2>/dev/null || echo 0) plot(s)" >&2
