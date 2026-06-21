#!/usr/bin/env bash
set -euo pipefail
#########################################################################################
# Regenerate the committed REAL-DATA litmus chrM BAMs from 1000 Genomes high-coverage
# (30x) public CRAMs. These are healthy blood samples (no large mtDNA deletion expected),
# used as a specificity check: the caller must yield ZERO PASS calls on them.
#
#   remote CRAM --(samtools view chrM)--> chrM reads  (GRCh38 chrM == rCRS == RefSeq/chrM.fa)
#               --(subsample to ~2000x)--> the pipeline's working depth
#               --(minimap2 -ax sr chrMC -> -F 0x90C -> circSam.pl -> sort)--> $s.chrM.bam
#                 (== a faithful circular-aware $O.bam, the caller's input)
#
# 1000G high-cov data is open-access (consented for public release); committing chrM-only
# subsets is fine. Source index:
#   ftp.1000genomes.ebi.ac.uk/.../1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index
# Each row: col1 = CRAM URL on ftp.sra.ebi.ac.uk (use https://), col10 = sample.
#
# Needs: samtools, minimap2 (network access). Usage: gen_real.sh [target_depth]
# NOTE: full mtDNA depth is ~8,000-22,000x; the spurious-call artifacts the v2 redesign
# fixes are most visible at FULL depth. Raise target_depth (and size) to stress-test; the
# committed BAMs use ~2000x to stay small (<~10 MB each).
#########################################################################################

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../.." && pwd)
RDIR=${HP_RDIR:-$ROOT/RefSeq}; SDIR=${HP_SDIR:-$ROOT/scripts}
DEPTH=${1:-2000}
ref=$RDIR/chrM.fa; mtc=$RDIR/chrMC.fa; fai=$RDIR/chrM.fa.fai
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

# sample -> 1000G high-cov CRAM (https). Add rows to widen the panel (more candidates commented).
SAMPLES="
NA12718	https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram
NA12748	https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239481/NA12748.final.cram
NA12775	https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239482/NA12775.final.cram
"
# more candidates (uncomment / add to widen): from the high-cov sequence.index, col1=CRAM url, col10=sample
#   NA12777  .../ERR3239483/NA12777.final.cram      NA12778  .../ERR3239484/NA12778.final.cram
#   NA12827  .../ERR3239485/NA12827.final.cram      NA18488  .../ERR3239491/NA18488.final.cram

printf '%s\n' "$SAMPLES" | while IFS=$'\t' read -r s url; do
  [ -n "${s:-}" ] || continue
  echo "[gen_real] $s" >&2
  samtools view --reference "$ref" -b -o "$tmp/$s.chrM.bam" "$url" chrM
  dp=$(samtools depth -a "$tmp/$s.chrM.bam" | awk '{x+=$3;n++}END{printf "%.0f",x/n}')
  sf=$(awk -v d="$dp" -v t="$DEPTH" 'BEGIN{f=t/d; if(f>=1)f=0.99; printf "%.4f", f}'); sf=${sf#0.}
  samtools view -s "7.$sf" -b "$tmp/$s.chrM.bam" | samtools sort -o "$tmp/$s.sub.bam" -
  samtools collate -u -O "$tmp/$s.sub.bam" | \
    samtools fastq -1 "$tmp/r1.fq" -2 "$tmp/r2.fq" -0 /dev/null -s /dev/null -n - 2>/dev/null
  minimap2 -ax sr "$mtc" "$tmp/r1.fq" "$tmp/r2.fq" 2>/dev/null | samtools view -h -F 0x90C - | \
    "$SDIR/circSam.pl" -ref_len "$fai" -offset 0 | samtools sort -o "$HERE/$s.chrM.bam" --write-index -
  echo "[gen_real] wrote $HERE/$s.chrM.bam ($(samtools view -c "$HERE/$s.chrM.bam") reads)" >&2
done
ls -la "$HERE"/*.chrM.bam
