# SV module — mock data & test harness

Self-contained, real-time evaluation of the MitoHPC structural-variant caller
(`scripts/callSV.sh` and friends). No full pipeline run is required.

## Layout

| Path | Committed | Purpose |
|---|---|---|
| `bams/*.bam(.bai)` | ✅ (~1–2 MB each) | Mock `$O.bam`-equivalents: circular-aware chrM alignments carrying split-read junctions + coverage drops |
| `truth.tsv` | ✅ | Ground truth (breakpoints, svlen, heteroplasmy) per sample |
| `make_testdata.py` | ✅ | Read simulator (WT + deleted circular genomes → FASTQ) |
| `gen_bams.sh` | ✅ | Aligns FASTQ → BAM through the pipeline's circular path (minimap2 → `circSam.pl` → sort) |
| `run_test.sh` | ✅ | Runs the caller on each BAM and compares to `truth.tsv` |
| `fastq/`, `out/` | ❌ (gitignored) | Regenerable intermediates |

## Samples

| Sample | Deletion | Heteroplasmy | Tests |
|---|---|---|---|
| `sv_del4977_h30` | common deletion m.8470–13446 (~4977 bp, 13 bp direct repeat) | 30% | PASS call, `REPEAT` flag, heteroplasmy |
| `sv_del4977_h05` | common deletion | 5% | low-heteroplasmy detection (split-only `no_cvg_drop` tier) |
| `sv_del6000_h50` | non-repeat deletion m.6000–10998 (~4999 bp) | 50% | generality, **no** `REPEAT` flag |
| `sv_wt` | none | 0% | specificity (zero PASS calls) |

The del4977 breakpoint lands inside the 13 bp direct repeat (called ~`8482/13447`,
svlen 4964) — the real-world microhomology ambiguity — which `run_test.sh` matches
within tolerance and which trips the `REPEAT` flag.

## Run the test

```bash
bash test/sv/run_test.sh          # -> "ALL TESTS PASSED"
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
