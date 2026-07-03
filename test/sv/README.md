# SV module — mock data & test harness

Self-contained, real-time evaluation of the MitoHPC structural-variant caller
(`scripts/callSV.sh` and friends). No full pipeline run is required.

## Layout

| Path | Committed | Purpose |
|---|---|---|
| `bams/*.bam(.csi)` | ✅ (21 BAMs, ~31 MB total) | Mock `$O.bam`-equivalents: circular-aware chrM alignments carrying split-read junctions + coverage signal |
| `truth.tsv` | ✅ | Ground truth per (sample, event): `sample kind bp5 bp3 svlen het depth expect` (`expect` = today's caller behavior) |
| `make_testdata.py` | ✅ | Read simulator: WT + per-event circular genomes (deletion / duplication / inversion / inverted-dup / dup-del / origin-crossing) → FASTQ |
| `gen_bams.sh` | ✅ | Aligns FASTQ → BAM through the pipeline's circular path (minimap2 → `circSam.pl` → sort) |
| `run_test.py` | ✅ | The harness (python3 + pysam): scenarios vs `truth.tsv`, degenerate inputs, cohort, VCF-spec gate, example-schema check |
| `run_test.sh` | ✅ | Thin wrapper → `run_test.py` |
| `example/` | ✅ | Committed example outputs (per-sample VCF/tab, cohort VCFs, `sv.report.html`) — a "taste" of what the module emits. See `example/README.md` |
| `make_example.sh` | ✅ | Regenerates `example/` from the mock BAMs (paths sanitized to repo-relative) |
| `fastq/`, `out/` | ❌ (gitignored) | Regenerable intermediates |

## Samples (21) + robustness checks

> **[`TEST_BAMS.md`](TEST_BAMS.md) is the authoritative catalog** of all 21 mock BAMs — what each one
> is, *why* that scenario exists (the caller behavior it pins down), and what the harness asserts.
> Below is the at-a-glance grouping.

- **Deletions (PASS-able):** `sv_del4977_h30/h05` (common, 30/5%), `sv_del6000_h50` (non-repeat),
  `sv_multidel` (two), `sv_homoplasmy` (95%), `sv_dloop` (D-loop flag), `sv_lowcov` (40×), `sv_del_500`
  (small detectable), `sv_del_45` (< minsize → 0 records), `sv_del_13kb` (majority-arc w/ drop → PASS).
- **Origin-crossing (WRAP-withheld):** `sv_origin` (clips OriH), `sv_del_origin_spares` (spares) — the
  [origin-resolution](../docs/SV_EVENT_TYPES.md) regression pair.
- **Forward-looking (not yet callable; design in [`../docs/SV_EVENT_TYPES.md`](../docs/SV_EVENT_TYPES.md)):**
  duplications `sv_dup` (1 kb) / `sv_dup_large` (5 kb); complex `sv_dupdel` (a documented spurious-PASS
  gap) / `sv_invdup`; inversions `sv_inv_small/large/origin/lowhet` (0 records — the strand-filter blind
  spot). Each carries an `expect` in `truth.tsv`; a signature pre-assertion confirms the BAM really
  holds its signal first.
- **Control:** `sv_wt` (wild-type specificity, also the degenerate-input substrate).

Plus: degenerate inputs (empty / unindexed / wrong-contig / wrong-`mtlen` BAM → clean error, never a
traceback), cohort aggregation (`getSVSummary.sh` merge matrix + sites + recurrence), and a
`bcftools` VCF-spec gate. The del4977 breakpoint lands inside the 13 bp direct repeat (called
~`8482/13447`), matched within tolerance.

## Run the test

Needs `python3` with `pysam` (the caller does everything in-process). The cohort + VCF-spec gates
also use `bcftools`/`bgzip`/`tabix` when present (skipped otherwise). Point `HP_PYTHON` at an
interpreter that has `pysam` if your default `python3` doesn't:

```bash
bash test/sv/run_test.sh                                   # -> "ALL TESTS PASSED"
HP_PYTHON=/path/to/venv/bin/python bash test/sv/run_test.sh   # if pysam is in a venv
```

## Regenerate the BAMs (optional)

Requires `minimap2` + `samtools` + `perl` (and `RefSeq/chrMC.*`):

```bash
python3 test/sv/make_testdata.py -ref RefSeq/chrM.fa -out test/sv/fastq
bash    test/sv/gen_bams.sh        test/sv/fastq test/sv/bams
cp      test/sv/fastq/truth.tsv    test/sv/truth.tsv
```

The simulator is deterministic (fixed seed), so regenerated BAMs are reproducible.
The production pipeline aligns with `bwa mem`; minimap2 is used here only to mint mock
data — the split-read + coverage-drop signal the caller consumes is aligner-agnostic.
