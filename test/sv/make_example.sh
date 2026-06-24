#!/usr/bin/env bash
set -eu
#########################################################################################
# Regenerate the committed example SV outputs in test/sv/example/ from the mock BAMs, so
# the repo carries a representative "taste" of what the SV module emits:
#
#   example/<sample>.sv.vcf / .sv.tab   a per-sample call set (rich VCF + tidy table)
#   example/sv.tab                      cohort tidy/long table
#   example/sv.merged.vcf.gz (+ .tbi)   cohort genotype matrix (NS recurrence)
#   example/sv.sites.vcf.gz  (+ .tbi)   sites-only union (annotation substrate)
#   example/sv.report.html              the interactive cohort report
#
# Machine-specific absolute paths in VCF provenance headers are rewritten to repo-relative
# so the committed files are portable. Needs python3+pysam, bcftools, bgzip, tabix.
#########################################################################################

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
export HP_SDIR=$ROOT/scripts HP_RDIR=$ROOT/RefSeq
EX=$HERE/example
SAMPLE=sv_del4977_h30          # representative per-sample example (canonical common deletion)

# Include the samplot gallery in the example report so it demonstrates the visualization. Needs
# samplot on PATH; svplot.sh degrades gracefully (no gallery) if it is absent. The PNGs are embedded
# as base64 in the committed sv.report.html (portable); the per-sample .png/manifest are NOT committed.
export HP_SV_PLOT=${HP_SV_PLOT:-1}

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/out" "$EX"; : > "$tmp/in.txt"
for b in "$HERE"/bams/*.bam; do
  s=$(basename "$b" .bam); mkdir -p "$tmp/out/$s"
  bash "$HP_SDIR/callSV.sh" "$s" "$b" "$tmp/out/$s/$s" >/dev/null
  printf "%s\t%s\t%s\n" "$s" "$b" "$tmp/out/$s/$s" >> "$tmp/in.txt"
done
export HP_IN=$tmp/in.txt HP_ODIR=$tmp/out
bash "$HP_SDIR/getSVSummary.sh" "$tmp/out" >/dev/null

# strip machine-specific absolute paths -> repo-relative, and normalize the version-stamped
# ##source line so the committed example is stable across regenerations
SED="s|$ROOT/|./|g; s|$tmp/|./out/|g; s|^##source=MitoHPC_callsv .*|##source=MitoHPC_callsv (example output)|; s| --version [^ ]*| --version example|g"
sani(){ case "$1" in
  *.gz) bgzip -dc "$1" | sed "$SED" | bgzip > "$1.s" && mv "$1.s" "$1" && tabix -f -p vcf "$1";;
  *)    sed "$SED" "$1" > "$1.s" && mv "$1.s" "$1";;
esac; }

cp "$tmp/out/$SAMPLE/$SAMPLE.sv.vcf" "$EX/" && sani "$EX/$SAMPLE.sv.vcf"
cp "$tmp/out/$SAMPLE/$SAMPLE.sv.tab" "$EX/"
cp "$tmp/out/sv.tab" "$EX/"
cp "$tmp/out/sv.merged.vcf.gz" "$EX/" && sani "$EX/sv.merged.vcf.gz"
cp "$tmp/out/sv.sites.vcf.gz" "$EX/"  && sani "$EX/sv.sites.vcf.gz"
cp "$tmp/out/sv.report.html" "$EX/"

echo "[make_example] wrote example outputs to $EX" >&2
ls -la "$EX"
