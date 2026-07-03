#!/usr/bin/env bash
set -euo pipefail
#########################################################################################
# Realistic POSITIVE control: spike a known del4977 into a REAL chrM wild-type background
# (1000G high-cov, e.g. NA12718) at a target heteroplasmy, realign through the pipeline's
# circular path, and emit the caller's $O.bam. This is a "semi-real" truth set — real WT
# error/coverage/NUMT structure + a deletion of known breakpoints & fraction — used to vet
# that the method RECOVERS a real-background deletion (PASS + COMMON + AFC ~ het) with no
# off-target false positives, which neither the pure-simulation mocks nor the healthy real
# samples (0 PASS) can show on their own.
#
#   real chrM BAM --(samtools fastq)--> WT background reads (W pairs)
#   + del4977 reads (D pairs, D = het/(1-het) * W * len(delGenome)/16569)  [make_testdata model]
#   --> merge --(minimap2 -ax sr chrMC -> -F 0x90C -> circSam.pl -> sort)--> $O.bam (committed)
#
# Usage: gen_spike.sh [het=0.20] [real_bam=NA12718.chrM.bam]   Needs samtools, minimap2, python3.
#########################################################################################

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../.." && pwd)
RDIR=${HP_RDIR:-$ROOT/RefSeq}; SDIR=${HP_SDIR:-$ROOT/scripts}
HET=${1:-0.20}
REALBAM=${2:-$HERE/NA12718.chrM.bam}
BP5=8469; BP3=13447; DELGLEN=$((16569 - (BP3 - BP5 - 1)))   # retained length of the del genome
hh=$(awk -v h="$HET" 'BEGIN{printf "%02.0f", h*100}')
OUT="$HERE/spike_del4977_h${hh}.chrM.bam"
mtc=$RDIR/chrMC.fa; fai=$RDIR/chrM.fa.fai
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
test -s "$REALBAM"

echo "[gen_spike] real WT background: $REALBAM" >&2
samtools collate -u -O "$REALBAM" 2>/dev/null | \
  samtools fastq -1 "$tmp/wt1.fq" -2 "$tmp/wt2.fq" -0 /dev/null -s /dev/null -n - 2>/dev/null
W=$(( $(wc -l < "$tmp/wt1.fq") / 4 ))
D=$(awk -v h="$HET" -v w="$W" -v g="$DELGLEN" 'BEGIN{printf "%d", (h/(1-h))*w*g/16569}')
echo "[gen_spike] WT pairs=$W ; injecting D=$D del4977 pairs for het~$HET" >&2

# generate D del4977 read-pairs from the WT+DEL simulator's emitter (deterministic seed)
MT_PATH="$ROOT/test/sv/make_testdata.py" HP_RDIR="$RDIR" \
  "${HP_PYTHON:-python3}" - "$D" "$BP5" "$BP3" "$tmp" <<'PY'
import os, random, sys, importlib.util
spec = importlib.util.spec_from_file_location("mt", os.environ["MT_PATH"])
mt = importlib.util.module_from_spec(spec); spec.loader.exec_module(mt)
D, bp5, bp3, tmp = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
seq = mt.read_fasta_single(os.path.join(os.environ["HP_RDIR"], "chrM.fa"))
eg = mt.make_deletion(seq, bp5, bp3)
rng = random.Random(7)
with open(os.path.join(tmp, "del1.fq"), "w") as f1, open(os.path.join(tmp, "del2.fq"), "w") as f2:
    mt.emit_reads(f1, f2, eg + eg, D, 150, 300, 450, 0.001, rng, "del4977spike")
PY

cat "$tmp/wt1.fq" "$tmp/del1.fq" > "$tmp/r1.fq"
cat "$tmp/wt2.fq" "$tmp/del2.fq" > "$tmp/r2.fq"
minimap2 -ax sr "$mtc" "$tmp/r1.fq" "$tmp/r2.fq" 2>/dev/null | samtools view -h -F 0x90C - 2>/dev/null | \
  "$SDIR/circSam.pl" -ref_len "$fai" -offset 0 2>/dev/null | samtools sort -o "$OUT" --write-index - 2>/dev/null
echo "[gen_spike] wrote $OUT ($(samtools view -c "$OUT") reads, ~$(samtools depth -a "$OUT"|awk '{s+=$3;n++}END{printf "%.0f",s/n}')x)" >&2
ls -la "$OUT"
