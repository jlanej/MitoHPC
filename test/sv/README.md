# SV module — mock data & test harness

Self-contained, real-time evaluation of the MitoHPC structural-variant caller
(`scripts/callSV.sh` and friends). No full pipeline run is required.

## Layout

| Path | Committed | Purpose |
|---|---|---|
| `bams/*.bam(.csi)` | ✅ (~13 MB total) | Mock `$O.bam`-equivalents: circular-aware chrM alignments carrying split-read junctions + coverage signal |
| `truth.tsv` | ✅ | Ground truth per (sample, event): `sample kind bp5 bp3 svlen het depth` |
| `make_testdata.py` | ✅ | Read simulator: WT + per-event circular genomes (deletion / duplication / origin-crossing) → FASTQ |
| `gen_bams.sh` | ✅ | Aligns FASTQ → BAM through the pipeline's circular path (minimap2 → `circSam.pl` → sort) |
| `run_test.py` | ✅ | The harness (python3 + pysam): scenarios vs `truth.tsv`, degenerate inputs, cohort, VCF-spec gate, example-schema check |
| `run_test.sh` | ✅ | Thin wrapper → `run_test.py` |
| `example/` | ✅ | Committed example outputs (per-sample VCF/tab, cohort VCFs, `sv.report.html`) — a "taste" of what the module emits. See `example/README.md` |
| `make_example.sh` | ✅ | Regenerates `example/` from the mock BAMs (paths sanitized to repo-relative) |
| `fastq/`, `out/` | ❌ (gitignored) | Regenerable intermediates |

## Samples (10) + robustness checks

| Sample | Construction | Checks |
|---|---|---|
| `sv_del4977_h30` / `_h05` | common deletion @30% / @5% | PASS + `REPEAT`/`COMMON`/`HOMLEN=13`/`DELCLASS=I`/genes; 5% → `no_cvg_drop` tier |
| `sv_del6000_h50` | non-repeat deletion @50% | PASS, no `REPEAT`/`COMMON` |
| `sv_multidel` | **two** deletions (del4977 + del6000) | both detected as separate PASS records |
| `sv_homoplasmy` | common deletion @95% | PASS, `AFJ→1.0`, no divide-by-zero |
| `sv_dup` | tandem duplication | **zero PASS** (`CVGR>1` → `no_cvg_drop`) |
| `sv_origin` | origin-crossing deletion | **zero PASS**, all coords ≤ contig (valid VCF) |
| `sv_dloop` | 5′ breakpoint in the D-loop | PASS + `DLOOP` flag |
| `sv_lowcov` | common deletion @50%, 40× depth | still detected |
| `sv_wt` | wild-type | 0 PASS (specificity) |

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
