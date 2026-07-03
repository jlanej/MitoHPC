# Example SV output

A committed "taste" of what the MitoHPC structural-variant module emits, generated from the mock
BAMs in [`../bams/`](../bams) (a 10-sample toy cohort). Open `sv.report.html` in a browser for the
interactive view; the rest are the machine-readable deliverables.

| File | What it is |
|---|---|
| **`sv.report.html`** | Self-contained interactive cohort report — circular mtDNA map + linear genome browser (gene/OXPHOS-complex annotation, per-position deletion-frequency track, VAF-coloured calls), live filtering, summary stats, recurrence table. Open offline in any browser. |
| `sv_del4977_h30.sv.vcf` | A representative **per-sample VCF** (the canonical common deletion): spec-correct headers (`##contig`/`##reference`/provenance), sample-named genotype column, rich `INFO` (`HOMLEN`/`HOMSEQ`/`DELCLASS`/`CIPOS`/`SVCLAIM`/`GENE`/`COMMON`/`HGVS`), `FORMAT GT:DP:AD:AF:SR`. |
| `sv_del4977_h30.sv.tab` | The same call as a tidy one-row-per-deletion table. |
| `sv.tab` | **Cohort** tidy/long table (one row per sample-deletion) for R/pandas. |
| `sv.merged.vcf.gz` (+`.tbi`) | **Cohort genotype matrix** (`bcftools merge`): one row per site, one column per sample, `NS` = number of samples carrying it (recurrence). |
| `sv.sites.vcf.gz` (+`.tbi`) | Sites-only union for annotation (AnnotSV/VEP). |

Machine-specific paths in the VCF provenance headers are rewritten to repo-relative, and the
version-stamped `##source` line is normalised, so these files are portable and stable.

## Regenerate

```bash
bash test/sv/make_example.sh        # needs python3+pysam, bcftools, bgzip, tabix
```

`test/sv/run_test.py` checks that these committed examples match the **current** output schema
(`example_report` / `example_tab_schema` / `example_vcf_schema`); if you change the VCF/tab/report
format, regenerate them or the test will flag the drift.
