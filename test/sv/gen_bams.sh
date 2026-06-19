#!/usr/bin/env bash
set -euo pipefail
#########################################################################################
# Generate mock chrM BAMs for SV-module testing by running simulated reads through the
# SAME circular-aware alignment path the MitoHPC pipeline uses to build $O.bam:
#
#   reads --(minimap2 -ax sr)--> chrMC (circularized, 16869bp)
#         --(samtools view -F 0x90C: drop unmapped/mate-unmapped/secondary/supplementary)-->
#         --(circSam.pl: wrap origin-crossing reads into 1..16569, add SA-tagged splits)-->
#         --(samtools sort/index)--> bams/<sample>.bam   (== a faithful $O.bam)
#
# NOTE: the production pipeline aligns with `bwa mem`; minimap2 is used here only to
# generate mock data (it is commonly preinstalled and emits spec-compliant SA tags).
# The split-read + coverage-drop signal the caller consumes is aligner-agnostic.
#
# Usage: gen_bams.sh [fastq_dir] [bams_dir]
#########################################################################################

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
RDIR=${HP_RDIR:-$ROOT/RefSeq}
SDIR=${HP_SDIR:-$ROOT/scripts}

FQDIR=${1:-$HERE/fastq}
BDIR=${2:-$HERE/bams}
MTC=$RDIR/chrMC.fa          # circularized reference (named chrM, 16869bp)
FAI=$RDIR/chrM.fa.fai       # canonical length (16569) for circSam wrap

command -v minimap2 >/dev/null || { echo "ERROR: minimap2 not found" >&2; exit 1; }
test -s "$MTC"
test -s "$FAI"
mkdir -p "$BDIR"

for r1 in "$FQDIR"/*_1.fq; do
  [ -e "$r1" ] || { echo "No *_1.fq in $FQDIR" >&2; exit 1; }
  s=$(basename "$r1" _1.fq)
  r2=$FQDIR/${s}_2.fq
  echo "[gen_bams] aligning $s" >&2
  minimap2 -ax sr -R "@RG\tID:$s\tSM:$s\tPL:ILLUMINA" "$MTC" "$r1" "$r2" 2>/dev/null | \
    samtools view -h -F 0x90C - | \
    "$SDIR/circSam.pl" -ref_len "$FAI" -offset 0 | \
    samtools sort -o "$BDIR/$s.bam" --write-index -
done

echo "[gen_bams] done -> $BDIR" >&2
ls -la "$BDIR"/*.bam
