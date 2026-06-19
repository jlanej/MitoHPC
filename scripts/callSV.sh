#!/usr/bin/env bash
set -eu
#########################################################################################
# Standalone mitochondrial structural-variant (deletion) caller for MitoHPC.
#
# Consumes the circular-aware chrM alignment ($O.bam) and writes ONLY new files:
#     $O.sv.vcf   per-sample deletion calls (VCFv4.2, SVTYPE=DEL)
#     $O.sv.tab   flat table
# It never reads-for-write or deletes any existing deliverable. Default-off: it is only
# invoked when HP_SV is set (see filter.sh). See docs/SV_CALLING.md.
#
# Method (v1): split-read junctions (sa2del.pl) corroborated by a coverage drop
# (svCall.pl). A deletion is PASS only when both signals agree. Two heteroplasmy
# estimates (junction fraction + coverage ratio) plus a disagreement QC field.
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
P=${HP_P:-1}

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

# masks (optional; flags only)
HP_BED=$RDIR/HP.bed.gz
NUMT_VCF=$RDIR/NUMT.vcf.gz
DLOOP_BED=$RDIR/DLOOP.bed.gz
maskopt=""
[ -s "$HP_BED" ]    && maskopt="$maskopt -hp $HP_BED"
[ -s "$NUMT_VCF" ]  && maskopt="$maskopt -numt $NUMT_VCF"
[ -s "$DLOOP_BED" ] && maskopt="$maskopt -dloop $DLOOP_BED"

DP=$O.sv.dp
# per-base depth over chrM (self-contained; equivalent to $O.cvg)
samtools depth -a -r "$MT" -@ "$P" "$BAM" > "$DP"

# split-read deletion junctions (clustered)
samtools view -h -q "$MINMAPQ" -@ "$P" "$BAM" | \
  "$SDIR/sa2del.pl" -chrM "$MT" -mtlen "$MTLEN" -minsize "$MINSIZE" \
                    -maxsize "$MAXSIZE" -pad "$PAD" -minsupport 2 > "$O.sv.jun"

# corroborate with coverage, flag, format (sa2del.pl already emits clusters sorted by bp5=POS)
cat "$SDIR/sv.vcf" | sed "s|^#CHROM|##sample=$S\n#CHROM|" > "$O.sv.vcf"
"$SDIR/svCall.pl" -sample "$S" -ref "$RDIR/$MT.fa" -cvg "$DP" -mtlen "$MTLEN" \
                  -flank "$FLANK" -minjr "$MINJR" -drop "$DROP" -mindepth "$MINDP" \
                  -pad "$PAD" -tab "$O.sv.tab" $maskopt < "$O.sv.jun" >> "$O.sv.vcf"

rm -f "$DP" "$O.sv.jun"
echo "[callSV] $S -> $O.sv.vcf ($(grep -vc '^#' "$O.sv.vcf") records)" >&2
