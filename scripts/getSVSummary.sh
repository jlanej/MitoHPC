#!/usr/bin/env bash
set -eu
#########################################################################################
# Aggregate per-sample structural-variant calls into cohort files. SEPARATE from
# getSummary.sh (which is never modified). Invoked only when HP_SV is set.
#
#   $ODIR/sv.tab            tidy/long table: one row per (sample, deletion) for R/pandas
#   $ODIR/sv.merged.vcf.gz  cohort genotype matrix: one row per site, one column per
#                           sample, NS = number of samples carrying it (recurrence substrate)
#   $ODIR/sv.sites.vcf.gz   sites-only union (no genotypes) for annotation (AnnotSV/VEP)
#
# Uses bcftools/bgzip/tabix (already pipeline deps; the per-sample VCFs carry the ##contig
# line that bcftools needs). Each per-sample VCF is bgzip+tabix-indexed first. (Note: a
# single concatenated VCF is not produced — different samples => use the merged matrix or
# the long sv.tab. Exact-match merge can over-split imprecise breakpoints; positional/fuzzy
# merging across a cohort is a future refinement.)
#
# Usage: getSVSummary.sh [ODIR]   (defaults to $HP_ODIR; uses $HP_IN for sample prefixes)
#########################################################################################

if [ "$#" -lt 1 ]; then ODIR=$HP_ODIR; else ODIR=$1; fi
SDIR=${HP_SDIR:-$(cd "$(dirname "$0")" && pwd)}
RDIR=${HP_RDIR:-$(cd "$SDIR/../RefSeq" && pwd)}
MB="$RDIR/mitobreak.tsv.gz"                        # MitoBreak known-breakpoint DB (cohort annotation)
MBTOL="${HP_SV_MITOBREAK_TOL:-20}"                 # per-breakpoint match tolerance (bp)
test -s "$HP_IN"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
list=$tmp/list.txt; : > "$list"

# (1) tidy/long table: header from the first sample + every sample's data rows
first=$(awk 'NR==1{print $3; exit}' "$HP_IN").sv.tab
{ [ -s "$first" ] && head -1 "$first"
  for p in $(awk '{print $3}' "$HP_IN"); do
    [ -s "$p.sv.tab" ] && tail -n +2 "$p.sv.tab"
  done
} > "$ODIR/sv.tab"

# (1b) MitoBreak annotation (cohort-only, additive): append a `mitobreak` column flagging breakpoints
#      previously reported in MitoBreak. Per-sample $O.sv.tab is NOT touched — only this consolidated
#      table. The data dictionary (which documents the `mitobreak` column) is copied alongside ONLY
#      when the column is actually added, so the two can never disagree.
if [ ! -s "$MB" ]; then
  echo "[getSVSummary] WARNING: $MB absent — MitoBreak annotation DISABLED" >&2
elif [ -s "$ODIR/sv.tab" ]; then
  "${HP_PYTHON:-python3}" "$SDIR/svMitoBreak.py" --db "$MB" --tol "$MBTOL" \
    --tab "$ODIR/sv.tab" --tab-out "$tmp/sv.tab.mb" \
    && test -s "$tmp/sv.tab.mb" && mv "$tmp/sv.tab.mb" "$ODIR/sv.tab"
  [ -s "$SDIR/sv.tab.dict.tsv" ] && cp "$SDIR/sv.tab.dict.tsv" "$ODIR/sv.tab.dict.tsv"
fi

# bgzip+tabix each existing per-sample VCF
for p in $(awk '{print $3}' "$HP_IN"); do
  v=$p.sv.vcf
  if [ ! -s "$v" ]; then echo "[getSVSummary] WARN: missing $v" >&2; continue; fi
  g=$tmp/$(basename "$p").vcf.gz
  bgzip -c "$v" > "$g"
  tabix -p vcf -f "$g"
  echo "$g" >> "$list"
done
if [ ! -s "$list" ]; then echo "[getSVSummary] no SV VCFs found" >&2; exit 0; fi

# (2) cohort genotype matrix: one row per site, one column per sample, NS=#samples called
if [ "$(wc -l < "$list")" -gt 1 ]; then
  bcftools merge -m none -l "$list" -Ou | \
    bcftools +fill-tags -Oz -o "$ODIR/sv.merged.vcf.gz" -- -t AN,AC,AF,NS
else
  bcftools view -Oz -o "$ODIR/sv.merged.vcf.gz" "$(head -1 "$list")"
fi
tabix -p vcf -f "$ODIR/sv.merged.vcf.gz"

# (2b) MitoBreak annotation of the merged VCF (adds the MITOBREAK INFO field), done BEFORE the sites
#      union below so sv.sites.vcf.gz inherits it for free. Cohort-only; per-sample VCFs untouched.
if [ -s "$MB" ]; then
  "${HP_PYTHON:-python3}" "$SDIR/svMitoBreak.py" --db "$MB" --tol "$MBTOL" \
    --vcf "$ODIR/sv.merged.vcf.gz" --vcf-out "$tmp/sv.merged.mb.vcf" \
    && bgzip -c "$tmp/sv.merged.mb.vcf" > "$ODIR/sv.merged.vcf.gz" \
    && tabix -p vcf -f "$ODIR/sv.merged.vcf.gz"
fi

# (3) sites-only union (drop genotypes) for annotation tools
bcftools view -G -Oz -o "$ODIR/sv.sites.vcf.gz" "$ODIR/sv.merged.vcf.gz"
tabix -p vcf -f "$ODIR/sv.sites.vcf.gz"

# (3b) collect the per-sample samplot manifests (if HP_SV_PLOT produced any) into one cohort manifest
popt=""
: > "$ODIR/sv.plots.tsv"
for p in $(awk '{print $3}' "$HP_IN"); do
  [ -s "$p.sv.plots.tsv" ] && cat "$p.sv.plots.tsv" >> "$ODIR/sv.plots.tsv"
done
[ -s "$ODIR/sv.plots.tsv" ] && popt="--plots $ODIR/sv.plots.tsv"

# (4) interactive, self-contained HTML report (svReport.py is stdlib-only — no pysam needed;
#     embeds any samplot PNGs as base64 so the report stays single-file/offline).
#     HP_SV_PLOT_ALL=1 (the "show everything" mode) embeds EVERY plot: default --plot-dedup to 0 (no
#     representative collapsing) and tell the report its gallery is unfiltered/unsubsampled (--plot-all).
N=$(grep -vc '^#' "$HP_IN")
gopt=""; [ -s "$RDIR/genes.bed.gz" ] && gopt="--genes $RDIR/genes.bed.gz"
if [ -n "${HP_SV_PLOT_ALL:-}" ]; then DEDUP=${HP_SV_PLOT_DEDUP:-0}; allopt="--plot-all"
else                                  DEDUP=${HP_SV_PLOT_DEDUP:-25}; allopt=""; fi
"${HP_PYTHON:-python3}" "$SDIR/svReport.py" --tab "$ODIR/sv.tab" --nsamples "$N" \
  --mtlen "${HP_MTLEN:-16569}" --chrom "${HP_MT:-chrM}" --plot-dedup "$DEDUP" $allopt \
  $gopt $popt --out "$ODIR/sv.report.html" \
  || echo "[getSVSummary] WARN: HTML report generation failed" >&2

echo "[getSVSummary] wrote $ODIR/{sv.tab, sv.merged.vcf.gz, sv.sites.vcf.gz, sv.report.html}" >&2
