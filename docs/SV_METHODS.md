# MitoHPC mtDNA Deletion Calling — Method (as implemented)

**Status:** v1 (single large-scale deletions) · default **off** (`HP_SV`) · purely additive
**Applies to code:** `scripts/{callsv.py, callSV.sh, sv.vcf, getSVSummary.sh}` (core in Python 3 + `pysam`)
**Companion docs:** design/literature → [`SV_CALLING.md`](SV_CALLING.md) · guardrails → [`../CLAUDE.md`](../CLAUDE.md)

> ⚠️ **Maintenance contract:** this document describes the code *as it actually runs*. If you
> change the algorithm, thresholds, output schema, or wiring, update this file in the **same
> change** and add a line to the [Changelog](#changelog). Numbers and formulas here are quoted
> from the scripts and must stay in sync with them.

---

## 1. Intuition (the 30-second version)

A mitochondrial **deletion** leaves two fingerprints in aligned short reads:

```
                      deletion (absent from mutant molecules)
                  ┌───────────────────────────────────┐
  reference  ──────●bp5                             bp3●──────────────
                    \                                 /
  a read that        \____________  ___________ ____/        ← SPLIT READ:
  spans the junction              \/                          left part maps to bp5,
                                  (read)                       right part maps to bp3
                                                               (carried as an SA:Z tag)

  depth  ▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▁▁▁▁▁▁▁▁▁▁▁▁▁▁▆▆▆▆▆▆▆▆▆▆▆▆▆        ← COVERAGE DROP:
                          └ only wild-type molecules cover the deleted span ┘
```

The caller finds **both** signals and requires them to agree:

1. **Split reads** (`SA:Z:` tags) pin the breakpoints `bp5 → bp3` to base-pair resolution.
2. The **coverage drop** between `bp5` and `bp3` confirms a real dosage loss (and is what a
   nuclear-segment/chimera artifact will *not* produce).

Heteroplasmy (the mutant fraction) is estimated **two independent ways** — from the junction-read
fraction and from the coverage ratio — and their disagreement is reported as a QC flag. A
deletion is **PASS** only when a clustered split junction *and* a coverage drop coincide.

This is the standard split-read + read-depth corroboration principle (cf. DELLY/LUMPY), specialized
to the small circular mitochondrial genome and implemented in **Python 3 with `pysam`** (the
de-facto htslib binding — the same C library behind `samtools`), so the BAM is parsed in-process
with no Perl and no subprocess shelling.

---

## 2. Where it runs, what it reads, what it writes

The module is a **standalone, additive** stage. It is invoked per sample from `filter.sh` only
when `HP_SV=callsv`, in a gated block placed **after** the GRIDSS block and **before** the BAM
cleanup (`rm -f $O.bam*`), so it sees the live alignment:

```bash
# scripts/filter.sh (gated; default off)
if [ $HP_SV ] && [ "$HP_SV" == "callsv" ] ; then
  if [ ! -s $O.sv.vcf ] ; then
    callSV.sh $S $O.bam $O
  fi
fi
```

| Reads (never modified) | Writes (new `*.sv.*` only) |
|---|---|
| `$O.bam` — MitoHPC's circular-aware, subsampled, deduplicated chrM alignment (1..16569) | `$O.sv.vcf` — per-sample deletion calls (VCFv4.2, `SVTYPE=DEL`) |
| `RefSeq/$HP_MT.fa` — reference (for the VCF REF base) | `$O.sv.tab` — flat table (same calls) |
| `RefSeq/{HP,DLOOP}.bed.gz`, `RefSeq/NUMT.vcf.gz` — FP masks (flags only) | (cohort) `$ODIR/sv.concat.vcf`, `$ODIR/sv.tab` via `getSVSummary.sh` |

With `HP_SV` empty the block is skipped and **no existing deliverable changes** (verified by the
diff in `CLAUDE.md` §0 / the additive wiring). `getSummary.sh` is never edited; cohort SV
aggregation lives in a separate `getSVSummary.sh`, gated on `HP_SV`.

---

## 3. Pipeline (`callSV.sh` → `callsv.py`)

`callSV.sh $S $BAM $O` is a thin bash driver: it resolves the `HP_SV_*` thresholds (defaults in
brackets, set in `init.sh`) and the `RefSeq` mask paths, then runs the Python caller once:

```
callSV.sh ──▶ $HP_PYTHON(=python3) scripts/callsv.py --bam $O.bam --ref chrM.fa
                                                      --header scripts/sv.vcf --sample $S
                                                      --out $O.sv.vcf --tab $O.sv.tab [knobs/masks]
                         │
   ┌─────────────────────┴───────────────── all in one pysam pass ───────────────────────┐
   │ extract_junctions(): iterate $O.bam (MAPQ≥20), SA:Z split reads → clustered junctions │
   │ per_base_depth():    pysam count_coverage(quality_threshold=0) ≡ `samtools depth -a`  │
   │ call():              corroborate + heteroplasmy + flags → $O.sv.vcf + $O.sv.tab       │
   └──────────────────────────────────────────────────────────────────────────────────────┘
```

There are **no temp files and no subprocesses** — `pysam` reads the BAM, computes depth, and reads
the reference/masks in-process. The single dependency is `pysam` (installed in the image via a
pinned pip manylinux wheel that bundles htslib; see §12). `callsv.py` runs on Python 3.8+ and is
parity-tested against the original Perl implementation (field-for-field identical on the mock BAMs).

---

## 4. Stage A — split-read junctions (`callsv.py: extract_junctions`)

Iterates `$O.bam` via `pysam.AlignmentFile.fetch(chrom)` and returns clustered junctions
`(bp5, bp3, SVLEN, JR, strand)`.

### 4.1 Which reads are used
For each alignment record, keep it only if it is a **primary** alignment carrying an `SA:Z:` tag:
- skip `is_unmapped` (`0x4`), `is_secondary` (`0x100`), `is_supplementary` (`0x800`);
- require `read.mapping_quality ≥ HP_SV_MINMAPQ` (default **20**) — the first NUMT/multimapper guard;
- require `read.has_tag("SA")` (the supplementary mate was dropped upstream by `-F 0x90C`, but the
  tag remains on the primary — this is exactly the `$O.bam` representation).

### 4.2 Junction definition
Two segments are reconstructed: the **primary** (from its `POS`+`CIGAR`) and the **first `SA`
entry** (`rname,pos,strand,CIGAR,…`). Reference span = sum of CIGAR ops that consume reference
(`M/D/N/=/X`), 1-based inclusive. A junction is kept only when both segments are:
- on the **same contig**, and
- the **same strand** (same orientation ⇒ deletion; opposite ⇒ inversion, **not** called in v1).

With the upstream segment = smaller start:

```
bp5   = reference END   of the upstream segment   (last retained base before the deletion)
bp3   = reference START of the downstream segment  (first retained base after the deletion)
SVLEN = bp3 - bp5 - 1                               (number of deleted bases)
```

Junctions with `SVLEN < HP_SV_MINSIZE` (default **50**) or `> mtlen-1` are discarded (small indels
are the SNV caller's territory; the upper bound excludes the whole-genome/origin artifact).

> **Direct-repeat microhomology.** When a deletion is flanked by a direct repeat (e.g. the 13 bp
> repeat of the common deletion), the aligner places the breakpoint *anywhere within the repeat*,
> so the read shows a few bases of overlap between the two arms and `bp5/bp3` shift by up to the
> repeat length. This is expected and absorbed by clustering (`HP_SV_PAD`) and surfaced by the
> `REPEAT` flag — not an error.

### 4.3 Clustering
Each supporting read contributes `(bp5, bp3, read-id)` where read-id is mate-aware (`name/1`,
`name/2`). Points are sorted by `(bp5, bp3)` and grouped by **single-linkage greedy clustering to
the cluster seed**: a point joins the current cluster if both `|bp5 − seed.bp5| ≤ HP_SV_PAD` and
`|bp3 − seed.bp3| ≤ HP_SV_PAD` (default pad **25**); otherwise it starts a new cluster (its own
seed). Per cluster:
- `bp5`, `bp3`, `strand` = the **mode** (most frequent value) across members;
- `JR` = number of **distinct** supporting reads.

Clusters with `JR < 2` are dropped here (`-minsupport 2`); the PASS threshold `HP_SV_MINJR`
(default 3) is applied later, so the 2-read tier is still visible as a non-PASS call.

---

## 5. Stage B — coverage corroboration, heteroplasmy, flags (`callsv.py: call`)

**Input:** clustered junctions + per-base depth + reference (`--ref`) + masks. Depth comes from
`pysam.AlignmentFile.count_coverage(chrom, 0, mtlen, quality_threshold=0)` summed over A/C/G/T into
a 1-based array `dep[1..mtlen]`; with the default `read_callback="all"` (skips unmapped/secondary/
QC-fail/dup) this matches `samtools depth -a` (verified field-for-field on the mock BAMs). All
position windows **wrap modulo mtlen** so the circular origin is handled.

### 5.1 Coverage statistics (per junction)
```
medInside = median dep[ bp5+1 .. bp3-1 ]                                  (the deleted span)
medFlank  = median dep[ bp5-FLANK+1 .. bp5 ]  ∪  dep[ bp3 .. bp3+FLANK-1 ] (FLANK = HP_SV_FLANK[200])
CVGR      = medInside / medFlank                                          (1 if medFlank = 0)
```

### 5.2 Heteroplasmy (two estimates + disagreement)
```
span = round( (dep[bp5] + dep[bp3]) / 2 )      # local wild-type read depth at the breakpoints
SR   = max(span - JR, 0)                        # spanning (intact) reads
AFJ  = JR / (JR + SR)                            # (1) junction fraction      → primary AF
AFC  = clamp(1 - CVGR, 0, 1)                     # (2) coverage ratio         → orthogonal estimate
AFDIFF = |AFJ - AFC|                             # QC: large ⇒ amplification bias or DUP-as-DEL
```
`AFJ` is the reported `AF`; `AFC` and `AFDIFF` are reported in `INFO` so every number is auditable
and the two methods can be cross-checked.

### 5.3 PASS / FILTER logic
A call is **PASS** only if **none** of these fire:

| FILTER reason | Condition |
|---|---|
| `lowJR` | `JR < HP_SV_MINJR` (default 3) |
| `no_cvg_drop` | `CVGR > HP_SV_DROP` (default **0.9** ⇒ requires ≥10% coverage drop) |
| `WRAP` | a breakpoint within `originpad` (20 bp) of the origin (1 or mtlen) |
| `lowDP` | `HP_SV_MINDP > 0` and `medFlank < HP_SV_MINDP` (default 0 ⇒ disabled) |

The `no_cvg_drop` tier is where genuine **low-heteroplasmy** deletions land (below ~10% the
coverage dip is within noise) — they are still reported with full junction evidence, just not PASS.

### 5.4 False-positive / annotation flags (`INFO`)
| Flag | Meaning (fires if either breakpoint matches) |
|---|---|
| `REPEAT` | within `HP_SV_PAD` of the del4977 13 bp direct repeat (m.8470–8482 or m.13447–13459) |
| `WRAP` | within 20 bp of the artificial origin |
| `HP` | inside a homopolymer run (`RefSeq/HP.bed.gz`) |
| `DLOOP` | inside the control region (`RefSeq/DLOOP.bed.gz`) |
| `NUMT` | exact position in `RefSeq/NUMT.vcf.gz` |

Flags annotate; they do not by themselves reject a call (except `WRAP`, which also blocks PASS).
The primary NUMT defense is upstream: `$O.bam` was already built by competing reads against the
NUMT reference (`filter.sh`), and the mandatory coverage-drop gate removes chimeras that lack a
real dosage loss.

---

## 6. Output schema

### 6.1 `$O.sv.vcf` (and cohort `sv.concat.vcf`)
Header from `scripts/sv.vcf` (with a `##sample=` line injected). One record per cluster:

```
#CHROM POS   ID REF ALT    QUAL FILTER       INFO
chrM   8482  .  A   <DEL>  .    PASS         SM=sv_del4977_h30;SVTYPE=DEL;END=13446;SVLEN=-4964;
                                             JR=133;SR=368;AFJ=0.265;AFC=0.254;AFDIFF=0.011;
                                             CVGR=0.746;REPEAT   GT:DP:AF  0/1:566:0.265
```
- `POS = bp5`, `END = bp3-1` (last deleted base), `SVLEN = -(deleted bases)`, `REF` = reference
  base at `bp5`, `ALT = <DEL>`.
- `FORMAT = GT:DP:AF` → `0/1 : round(medFlank) : AFJ`.

### 6.2 `$O.sv.tab` (and cohort `sv.tab`)
`#sample  chrom  bp5  bp3  svlen  JR  SR  AFJ  AFC  AFDIFF  CVGR  FLANKDP  FILTER  flags`

### 6.3 Cohort aggregation (`getSVSummary.sh`)
Concatenates per-sample VCFs (drops `##sample=`, dedups, `bedtools sort -header`) into
`$ODIR/sv.concat.vcf`, and per-sample tabs (one header) into `$ODIR/sv.tab`. Separate from
`getSummary.sh`; emitted by `run.sh` only when `HP_SV` is set.

---

## 7. Parameters (all `HP_SV_*`, set in `init.sh`)

| Variable | Default | Meaning / effect |
|---|---|---|
| `HP_SV` | *(empty)* | `callsv` enables the module; empty = off |
| `HP_SV_MINMAPQ` | 20 | min MAPQ for split reads (NUMT multimapper guard) |
| `HP_SV_MINJR` | 3 | min distinct junction reads for a PASS deletion |
| `HP_SV_MINSIZE` | 50 | min deletion size (bp); separates from small indels |
| `HP_SV_MAXSIZE` | 0 | max deletion size (bp); `0` ⇒ `mtlen-1` |
| `HP_SV_PAD` | 25 | breakpoint clustering + direct-repeat tolerance (bp) |
| `HP_SV_DROP` | 0.9 | max `medInside/medFlank` for PASS (≤0.9 ⇒ ≥10% drop) |
| `HP_SV_FLANK` | 200 | flank window (bp) for the coverage ratio |
| `HP_SV_MINDP` | 0 | min flank depth for PASS (0 = disabled) |

`HP_SV_DROP` is the key sensitivity/specificity knob and is an open tuning question (see
`SV_CALLING.md` §11); it should be calibrated against a spiked dilution series.

---

## 8. Circular-genome handling

The module inherits circular correctness from the existing pipeline rather than re-implementing it:
- `$O.bam` was produced by `circSam.pl` from reads aligned to the **circularized** reference
  `chrMC` (`HP_E=300` bp appended), so a read crossing the artificial origin already has its parts
  wrapped into 1..16569.
- In `callsv.py`, every coverage window wraps modulo `mtlen`, so flanks straddling 16569/1 are
  computed correctly.
- v1 **does not** disambiguate a circular deletion from its complementary-arc duplication; junctions
  at the origin are flagged `WRAP` and kept out of PASS (deferred to a future tier — see roadmap in
  `SV_CALLING.md` §10).

---

## 9. Worked example (common deletion @ 30% heteroplasmy)

A read drawn across the deletion junction in a mutant molecule aligns to wild-type chrM as a split:
```
primary:  POS 8345  CIGAR 138M12S       (left arm, ends at 8482)
SA tag:   chrM,13447,+,125S25M          (right arm, starts at 13447)   ← 13 bp arm overlap = repeat
```
→ `extract_junctions`: `bp5=8482, bp3=13447, SVLEN=4964`. With 133 such reads clustered: `JR=133`.
→ `call`: `medFlank≈566`, `medInside≈422` ⇒ `CVGR=0.746`; `AFC=0.254`; `span≈501`,
   `SR≈368`, `AFJ=0.265`. `CVGR 0.746 ≤ 0.9` and `JR 133 ≥ 3` ⇒ **PASS**; breakpoints in the
   repeat ⇒ `REPEAT`. The two heteroplasmy estimates agree (`AFDIFF=0.011`) and bracket the true
   0.30. (Conventionally the common deletion is "4977 bp"; split reads report 4964 because the
   shared 13 bp repeat is counted once — both describe the same event.)

---

## 10. Validation & mock data (`test/sv/`)

Real-time, self-contained evaluation — no full pipeline run needed:

```bash
bash test/sv/run_test.sh        # -> ALL TESTS PASSED
```

`make_testdata.py` simulates paired-end reads from **wild-type + deleted circular genomes** at a
known heteroplasmy (so the coverage ratio outside vs inside the deletion equals the spiked
fraction by construction); `gen_bams.sh` aligns them through the pipeline's own circular path
(`minimap2 -ax sr → samtools view -F 0x90C → circSam.pl → sort`) to produce faithful `$O.bam`
files (committed, 1–2 MB each). `run_test.sh` runs `callSV.sh` and matches calls to `truth.tsv`
within tolerance.

| Sample | Truth | Result |
|---|---|---|
| common deletion @ 30% | DEL m.~8470–13446, 0.30 | **PASS**, `REPEAT`, AFJ 0.265 |
| common deletion @ 5% | same, 0.05 | detected, `no_cvg_drop` tier, AFJ 0.041 |
| non-repeat deletion @ 50% | DEL m.~6000–10998, 0.50 | **PASS**, no `REPEAT`, AFJ 0.404 |
| wild-type | none | 0 calls (specificity) |

See [`../test/sv/README.md`](../test/sv/README.md) for layout and regeneration.

---

## 11. Known limitations (v1)

- **Recall floor.** Uses `SA:Z:` split reads only (no soft-clip-only clustering, no local
  assembly); very low-heteroplasmy junctions with few split reads can be missed. The coverage-drop
  PASS gate intentionally relegates sub-~10% events to the `no_cvg_drop` tier.
- **DEL vs DUP.** No origin-of-replication logic; a complementary-arc duplication could be reported
  as a deletion. Mitigated by the coverage-drop gate (biases to true deletions) and `WRAP`/`DLOOP`.
- **Heteroplasmy is approximate.** `SR` is a coverage proxy (depth at the breakpoints − `JR`), not
  an exact intact-spanning-pair count; `AFC` is mildly biased near breakpoints by the coverage
  transition and near the D-loop. Reporting both estimates + `AFDIFF` exposes this.
- **Subsampling.** `HP_L` (~2000×) caps the lowest detectable heteroplasmy vs deep dedicated assays.
- **Deletions only.** No duplications, insertions, inversions, or multiple/complex rearrangements.

Planned v2/v3 work (soft-clip clustering, exact spanning counts, competitive NUMT re-scoring,
DEL/DUP disambiguation, optional eKLIPse/long-read engines) is in `SV_CALLING.md` §10.

---

## 12. How to run

**In the pipeline** (per sample, then cohort): set in `init.sh`
```bash
export HP_SV=callsv          # (optionally override HP_SV_* thresholds)
```
then run normally (`run.sh > run.all.sh; bash run.all.sh`). Produces `$O.sv.vcf`/`$O.sv.tab` per
sample and `$ODIR/sv.concat.vcf`/`$ODIR/sv.tab` for the cohort.

**Standalone** on any chrM BAM (needs `python3` with `pysam`):
```bash
HP_SDIR=scripts scripts/callSV.sh SAMPLE path/to.bam out/SAMPLE
# point at a specific interpreter if needed:
HP_SDIR=scripts HP_PYTHON=/path/to/venv/bin/python scripts/callSV.sh SAMPLE path/to.bam out/SAMPLE
```

**Dependency / Docker:** the only new dependency is `pysam` (Python), installed in the image by
`install_sysprerequisites.sh` (`pip install pysam==0.24.0`, a manylinux wheel that bundles htslib —
no compiler needed) and checked by `checkInstall.sh`. CI installs it via `actions/setup-python` +
pip (`.github/workflows/sv-test.yml`) and also exercises it inside the built image
(`docker-publish.yml`). `samtools`/`bedtools` remain installed for the rest of the pipeline but the
SV caller no longer shells out to them.

---

## Changelog

- **v1.1 (Python/pysam port):** reimplemented the two Perl cores (`sa2del.pl`, `svCall.pl`) as a
  single Python 3 + `pysam` module `scripts/callsv.py`; `callSV.sh` is now a thin driver
  (`HP_PYTHON` override). BAM iteration, SA/CIGAR parsing, and per-base depth (`count_coverage`,
  `quality_threshold=0`) run in-process — no Perl, no `samtools`/temp-file shelling in the SV path.
  **Field-for-field parity** with the Perl v1 verified on the mock BAMs; algorithm, thresholds, and
  output schema unchanged. `pysam` added to the Docker image + CI.
- **v1 (initial):** single large-scale deletion caller — `SA:Z:` split-read clustering +
  coverage-drop corroboration and dual heteroplasmy estimates (originally `sa2del.pl`/`svCall.pl`),
  additive `HP_SV` wiring, cohort `getSVSummary.sh`, and the `test/sv/` mock-data harness.
