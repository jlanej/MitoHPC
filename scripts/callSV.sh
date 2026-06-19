#!/usr/bin/env bash
set -eu
#########################################################################################
# Standalone mitochondrial structural-variant (deletion) caller for MitoHPC.
#
# Thin driver around scripts/callsv.py (Python3 + pysam), which reads the circular-aware
# chrM alignment ($O.bam) in-process and writes ONLY new files:
#     $O.sv.vcf   per-sample deletion calls (VCFv4.2, SVTYPE=DEL)
#     $O.sv.tab   flat table
# It never reads-for-write or deletes any existing deliverable. Default-off: it is only
# invoked when HP_SV=callsv (see filter.sh). See docs/SV_METHODS.md.
#
# Method (v1): split-read junctions (SA:Z tags) corroborated by a coverage drop. A
# deletion is PASS only when both signals agree. Two heteroplasmy estimates (junction
# fraction AFJ + coverage ratio AFC) plus a disagreement QC field (AFDIFF).
#
# Args:  1: sample name   2: BAM (circular-aware chrM)   3: output prefix
#########################################################################################

S=$1
BAM=$2
O=$3

SDIR=${HP_SDIR:-$(cd "$(dirname "$0")" && pwd)}
RDIR=${HP_RDIR:-$(cd "$SDIR/../RefSeq" && pwd)}
MT=${HP_MT:-chrM}
MTLEN=${HP_MTLEN:-16569}
PY=${HP_PYTHON:-python3}            # override to point at a python that has pysam

# tunables (defaults chosen for ~2000x subsampled chrM; override via init.sh)
MINMAPQ=${HP_SV_MINMAPQ:-20}    # min MAPQ for split reads (drops NUMT multi-mappers)
MINJR=${HP_SV_MINJR:-3}         # min junction-supporting reads for PASS
MINSIZE=${HP_SV_MINSIZE:-50}    # min deletion size (bp)
MAXSIZE=${HP_SV_MAXSIZE:-0}     # max deletion size (0 => MTLEN-1)
PAD=${HP_SV_PAD:-25}            # breakpoint clustering / repeat tolerance (bp)
DROP=${HP_SV_DROP:-0.9}         # max medInside/medFlank to confirm a deletion (<=0.9 => >=10% drop)
FLANK=${HP_SV_FLANK:-200}       # flanking window for the coverage ratio (bp)
MINDP=${HP_SV_MINDP:-0}         # min flank depth for PASS (0 => disabled)

test -s "$BAM"
test -s "$RDIR/$MT.fa"

# optional false-positive masks (flags only)
maskopt=""
[ -s "$RDIR/HP.bed.gz" ]    && maskopt="$maskopt --hp $RDIR/HP.bed.gz"
[ -s "$RDIR/NUMT.vcf.gz" ]  && maskopt="$maskopt --numt $RDIR/NUMT.vcf.gz"
[ -s "$RDIR/DLOOP.bed.gz" ] && maskopt="$maskopt --dloop $RDIR/DLOOP.bed.gz"

"$PY" "$SDIR/callsv.py" \
  --bam "$BAM" --ref "$RDIR/$MT.fa" --header "$SDIR/sv.vcf" --sample "$S" \
  --out "$O.sv.vcf" --tab "$O.sv.tab" \
  --chrom "$MT" --mtlen "$MTLEN" --minmapq "$MINMAPQ" --minjr "$MINJR" \
  --minsize "$MINSIZE" --maxsize "$MAXSIZE" --pad "$PAD" --drop "$DROP" \
  --flank "$FLANK" --mindepth "$MINDP" $maskopt

echo "[callSV] $S -> $O.sv.vcf ($(grep -vc '^#' "$O.sv.vcf") records)" >&2
