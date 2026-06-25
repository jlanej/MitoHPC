# Mitochondrial Structural-Variant (SV) Calling — State of the Art & MitoHPC Design

> Companion to the guardrails in [`/CLAUDE.md`](../CLAUDE.md). This document collates the
> state-of-the-art for mtDNA structural-variant detection (methods, tools, call types, pros/cons,
> verified facts) and proposes a **simple, robust, defensible** way to add a **standalone,
> additive** SV module to MitoHPC.
>
> **Hard constraints (from `CLAUDE.md`):** the SV module is purely additive and must never alter
> any existing deliverable; it must reuse MitoHPC's circular-genome machinery; it is gated behind
> a default-off flag. Reusing the existing latent split-read BED (`$O.sa.bed`) is **optional** —
> we are free to compute a better signal if it improves results.
>
> _Last updated 2026-06-18. Built from a multi-source literature review with adversarial
> fact-checking; key claims are tagged with a confidence level and primary sources are listed in
> §12._

---

## 1. Scope & goal

Add the ability to **detect, quantify, and report mitochondrial structural variants** — primarily
**single large-scale deletions (SLSMDs)** for v1 — from the short-read alignments MitoHPC already
produces, without touching the SNV/CN/haplogroup deliverables. The north star is a caller that is
**robust and defensible first, complete later**: narrow, correct, auditable v1 → richer v2/v3.

The dominant clinically-relevant mtDNA SV is the **single large deletion**, so v1 targets exactly
that one call class with high confidence. **Tandem duplications and inversions are now implemented as
opt-in** (`--call-dup` / `--call-inv`, default off; see `SV_EVENT_TYPES.md`); multiple/complex
rearrangements, dispersed duplications, insertions, and the origin DEL-vs-DUP resolution remain
deferred (see roadmap §10).

---

## 2. Biology & clinical context (why deletions are the priority)

- **mtDNA is a ~16,569 bp circular, high-copy, heteroplasmic genome.** Variants exist as a mixture
  ("heteroplasmy"); SVs have a heteroplasmy fraction just like SNVs.
- **Single large-scale deletions (SLSMDs)** range ~1.1–10 kb and cause three overlapping,
  age/tissue-dependent phenotypes: **Pearson syndrome** (infancy, marrow failure, multisystem),
  **Kearns–Sayre syndrome** (onset <20 y; retinopathy + CPEO + cardiac conduction block), and
  **isolated CPEO** (deletion confined to skeletal muscle). _[high confidence]_
- **The "common deletion" (ΔmtDNA4977 / m.8470_13447):** removes exactly **4977 bp** and is the
  single most frequent SLSMD (~30–37% of patients). It is flanked by a **perfect 13-bp direct
  repeat `5'-ACCTCCCTCACCA-3'`** at **m.8470–8482** (ATP8/ATP6 region) and **m.13447–13459**
  (ND5). One repeat copy is lost, one retained, so the actually-removed span is m.8470–13446 (≡
  m.8482–13459); literature quotes it either way. _[high confidence — Schon et al., Science 1989;
  Yusoff et al. 2019]_. **MitoHPC already ships synthetic junction references for this event:**
  `RefSeq/chrM.8455-8482_13460-13474.fa` (`del4977_43mer`, plus `_33mer`/`_23mer`) — ready-made
  positive controls. The 43-mer `CACAAACTACCACCTACCTCCCTCACCATTGGCAGCCTAGCAT` embeds the repeat.
- **Breakpoint classes:** ~70% are **Class I** (flanked by perfect direct repeats; slipped-strand
  mispairing / homologous-recombination mechanism), ~30% **Class II** (imperfect or no repeats).
  Breakpoint regions are significantly enriched for sequence homology. _[high confidence]_
- **Heteroplasmy thresholds:** biochemical dysfunction is classically cited at **>~60%** mutant
  mtDNA in the relevant tissue, though dysfunction is documented far lower (symptomatic at
  25–33%, occasionally <10%). For a screening caller, the useful target is to **reliably detect
  and quantify down to the low single-digit percent and below**. _[high confidence]_

**Sequencing signals a deletion produces (the substrate for any caller):** _[high confidence]_
1. **Read-depth drop** — a sharp, sustained coverage dip between the two breakpoints.
2. **Split / soft-clipped reads spanning the junction** — a read 3'-clipped marks the 5'
   breakpoint; 5'-clipped marks the 3' breakpoint. These give base-pair breakpoints.
3. **Discordant read pairs** — abnormally large apparent insert size flanking the deletion.

---

## 3. Algorithmic primitives (the building blocks)

Every SV caller combines a subset of these. For a small, high-copy, circular genome the trade-offs
differ from nuclear WGS.

| Primitive | What it detects | Strength on chrM | Weakness on chrM |
|---|---|---|---|
| **Split-read / soft-clip (SR)** | Base-pair breakpoints of the junction | Reaches **sub-1%** heteroplasmy at high depth; precise breakpoints | NUMT/chimeric reads create false junctions; needs reads that actually cross the junction |
| **Read-depth / coverage (RD)** | Dosage drop across the deleted span | Simple, robust, already computed (`$O.cvg`); good heteroplasmy at ≥~5% | **Loses power below ~5%**; confounded by D-loop/coverage non-uniformity and PCR/primer artifacts |
| **Discordant read-pair (PE)** | Approximate breakpoint region | Cheap | Poor breakpoint resolution; chrM's own insert-size distribution distorts the model |
| **Local/breakend assembly** | Precise complex junctions | Highest precision (GRIDSS/Manta) | Heavy for a 16.5 kb target; tools assume linear genome |

**Heteroplasmy estimation — two standard families** _[high confidence; both standard, not always
concordant]_:
- **Junction (split-read) fraction:** `AF_junction = J / (J + S)` where `J` = distinct
  junction-supporting reads, `S` = wild-type reads spanning the breakpoints intact. Used by
  MitoSAlt, eKLIPse, MitoDel. Reaches sub-1%.
- **Coverage ratio:** `AF_cov = 1 − (mean depth inside deletion / mean depth in flanking WT)`.
  Used by long-read/IGV workflows. Reliable ≥~5%.
- **Caveat:** PCR/long-range-amplification bias toward shorter (deleted) molecules can inflate the
  junction fraction (~4× vs qPCR in Splice-Break). Reporting **both** estimates plus their
  disagreement is the defensible move — it surfaces amplification bias and
  deletion-vs-duplication ambiguity rather than hiding them.
- MitoDel's binomial heteroplasmy with analytic 95% CI: `q = (n/r) × (16569/l)` (n=supporting
  reads, r=total reads, l=deletion length) is a clean, citable formulation.

**Circular-genome handling — two established strategies** _[high confidence; standard practice]_:
1. **Rotate/shift the reference** by an offset, re-call, map coordinates back. GATK Mutect2
   mito-mode and gnomAD v3.1 realign chrM to a reference **shifted by 8000 bp** to recover the
   control region; mitoCaller breaks at position 8000.
2. **Duplicate/concatenate (extend)** the reference so an origin-spanning read aligns across the
   seam, then remap back. CircularMapper appends N bp; some aligners (GSNAP `--circular`) do it
   natively.
   **MitoHPC already does both:** `chrMC` (concatenate, `HP_E=300`) for alignment + `chrMR`
   (rotate) for origin-region calling, with `circSam.pl` wrapping spanning reads into 1..16569 and
   the `(POS−HP_E) % HP_MTLEN` coordinate fix. **An SV module should reuse this, not reinvent it.**

---

## 4. Tool landscape

### 4a. Dedicated mtDNA deletion/SV tools (the relevant prior art)

| Tool | Signal | Calls | Detection limit | Circular handling | Maintenance | Notes |
|---|---|---|---|---|---|---|
| **eKLIPse** | Soft-clip + BLASTN realignment | Deletions (dups only in unpublished v2.1) | **<0.5%** mutant | Soft-clip junction recognition | Peer-reviewed (Genet Med 2019), **Python 2.7 (EOL)** | BAM/SAM in; per-gene cumulative deletion %; Circos plots; clinically used |
| **MitoSAlt** | Split-read via LAST (HISAT2 prefilter) | **Deletions AND duplications**, single-bp | **~0.5%** | **OriH/OriL logic disambiguates DEL vs complementary-arc DUP** | Peer-reviewed (PLOS Genet 2020), research pipeline | The most principled circular handling; multi-step (HISAT2+LAST) |
| **MitoDel** | Soft-clip/unmapped → BLAT | Deletions only | **~0.006% theoretical; 100% sens/0 FP at 0.1%** | 1..16569 only, no explicit origin logic | BMC Bioinf 2017; **lab-page only, no repo** | Binomial heteroplasmy + 95% CI; ≥10 supporting reads default |
| **Splice-Break2** | RNA-seq splice-junction aligner (MapSplice2) repurposed | Deletions (breakpoints) | floor ~0.1% | rCRS; masks first/last 100 bp | **Actively maintained** (USC/Hjelm), GitHub v3.0.1 | DNA-seq + RNA-seq; "deletion read %" is biased vs true heteroplasmy; needs ≥5000× |
| **MitoSV-seq** | RCA + Pindel/VarScan2 | Full SV spectrum (DEL/INS/DUP/INV) + SNV | <1% | RCA favors intact circles (~20× fewer NUMTs) | Research method (2020) | Requires wet-lab RCA + restriction digest |
| **Cas9/baldur** (ONT) | Cas9 linearize+barcode, ML caller | Deletions (incl. multiple), phasing | 88% sens @ ≥0.5% VAF, 95% precision | Cut-site removes mapping ambiguity | Research (Nat Commun 2022) | Long-read; bespoke ML caller |

**Not SV callers (avoid mis-attribution — verified):** `mity` is **FreeBayes-based, SNV/indel
only — it does NOT call SVs and is NOT GRIDSS-based** (the GRIDSS↔mito confusion is actually
**svaNUMT**, a NUMT-annotation tool). `MToolBox`, `mtDNA-Server 2 / Mutserve2`, and `mitoCaller`
are SNV/heteroplasmy callers; rare-disease pipelines pair them with **MitoSAlt** for deletions.
`HelixMTdb` is a population frequency **database**, not a caller. _[high confidence]_

### 4b. Generic short-read SV callers — and why they struggle on chrM

| Tool | Signals | Key chrM problem | _[conf]_ |
|---|---|---|---|
| **GRIDSS/GRIDSS2** | Break-end assembly + SR + PE | No circular/origin handling → can't separate circular DEL from complementary-arc DUP; heavyweight (JVM+assembly) for 16.5 kb; needs variant-supporting reads (misses very low VAF). *(This is the engine behind MitoHPC's existing `HP_V=gridss` stub.)* | high |
| **Manta** | PR + SR + local assembly (no RD) | **MaxDepth filter removes breakends where depth >3× median** → kills calls on uniform ultra-high-depth chrM; WGS guidance excludes chrM from `--callRegions` | high |
| **DELLY2** | PE + SR + RD | **chrM/MT is in DELLY's default exclude templates** (`human.hg38.excl.tsv`) → out-of-the-box DELLY does not call chrM; needs a genome-wide insert-size model chrM distorts | high |
| **LUMPY/smoove** | PE + SR + RD jointly | PE/RD assume linear, near-uniform depth — both violated on circular high-copy chrM | high |
| **BreakDancer** | PE only | Poor breakpoint resolution; depends on well-behaved insert size; superseded | high |
| **Pindel** | SR (pattern growth) | Good precision but no circular logic; loses power at low VAF and on the origin junction | high |

**Takeaway:** general SV callers are built for **large linear nuclear genomes** and assume
roughly uniform depth and a linear reference. On chrM they are confounded by (a) the circular
origin seam, (b) extreme/uniform depth (Manta's MaxDepth, DELLY's exclusion), and (c) low-VAF
heteroplasmy. Several **exclude chrM by default**. This is the core reason a small, purpose-built
caller that reuses MitoHPC's circular-aware alignment is more defensible than bolting on a generic
tool.

### 4c. Long-read callers (future direction)

`Sniffles2` is the **documented** caller in published mtDNA long-read (ONT/PacBio) deletion
workflows (with minimap2/NGMLR; mosaic mode for ~5–20% VAF). `cuteSV` is a strong generic
long-read caller but **has no documented mtDNA-specific use** (the claim is plausible
extrapolation only); generic long-read benchmarks **explicitly filter out chrM**. `SVIM` and
`SVision-pro` are notable for distinguishing tandem vs interspersed duplications and complex
rearrangements; `dysgu` spans short+long reads. Long reads resolve complex/duplicated
rearrangements short reads cannot, with deletion heteroplasmy correlating with ddPCR at r²=0.98
(>~10% reliable), but ONT's ~4–6% per-base error sets a ~5% floor for point variants.
_[partially-supported; Sniffles for mtDNA = supported, cuteSV for mtDNA = unsubstantiated]_

---

## 5. NUMTs & false positives _[high confidence]_

**NUMTs (nuclear-mitochondrial segments)** are a recognized primary source of false-positive mtDNA
calls: NUMT-derived reads mis-align to chrM, worst at **low VAF and low mtDNA copy number** (a
mis-mapped heterozygous NUMT yields VAF ≈ 1/(1+mtCN)). Mitigations used by leading pipelines:
- **Cross/competitive alignment** between mt and nuclear references (MToolBox realigns to nuclear;
  **MitoHPC already extracts NUMT-locus reads via `$HP_RNUMT` and competes them against chrM/NUMT
  references** — `filter.sh:94-101`). This is MitoHPC's strongest built-in NUMT defense.
- **Curated NUMT databases** (`RefSeq/NUMT.vcf.gz`; MitoHPC flags ~382 NUMT SNVs / 88 regions).
- **VAF and copy-number thresholds** (gnomAD drops mtCN<50; removes <1% heteroplasmy).
- **Mapping-quality** filtering (auxiliary, not the dominant mechanism).

For an SV module, the practical implication: **most NUMT reads are already removed upstream**, and
the **dual-signal requirement (a real coverage dropout must accompany the junction)** removes most
remaining NUMT/chimeric artifacts, because a NUMT chimera does not produce a true dosage drop.

---

## 6. Detection limits — what's realistically achievable _[high confidence, corrected]_

- **Split-read methods reach sub-1%** heteroplasmy: MitoDel 100% sens / 0 FP at 0.1% (theoretical
  ~0.006%); eKLIPse <0.5%; MitoSAlt ~0.5%. These are **read-depth dependent** (often >3000× for
  1.5–3% MAF).
- **Coverage-drop alone loses power below ~5%.**
- **Both NGS approaches beat Southern blot** (~≥10% only) and even ddPCR ("cannot reliably detect
  <10% heteroplasmy" for an SLSMD).
- **MitoHPC's subsampling** (`HP_L=222000` ≈ 2000× chrM) caps the lowest detectable heteroplasmy
  vs. deep dedicated assays — a deliberate, documentable trade-off. (v2 can offer an un-subsampled
  SV path.)

---

## 7. Validation & truth sets _[high confidence]_

- **MITOMAP "Reported mtDNA Deletions" (DeletionsSingle)** — curated catalog of literature
  deletions with sequenced breakpoints; a **reference**, not a validated ground-truth benchmark
  (incomplete; breakpoints approximate near repeats).
- **MitoBreak** — breakpoint-precise mtDNA-rearrangement database (PMID 24170808); use alongside
  MITOMAP for known-vs-novel flagging.
- **Simulators** — MitoSAlt ships a simulation framework (e.g. 200 events at 0.5% heteroplasmy,
  deletion sizes 50/500/2000 bp at ~6000×). The right way to tune thresholds.
- **Built-in positive control for MitoHPC:** the `del4977` junction FASTAs already in `RefSeq/`.
  A spiked dilution series of the common deletion is the natural acceptance test.

---

## 8. Design options for MitoHPC (three approaches compared)

All three respect the additive constraint identically (mirror the `HP_V=gridss` precedent: a
default-off flag, a single gated block in `filter.sh` placed **before the BAM cleanup**, a
separate aggregator, all outputs under an `*.sv.*` namespace, `getSummary.sh` untouched). They
differ in **where the split-read signal comes from** and **how much they depend on external tools**.

| | **A. `sa2sv`** (reuse `$O.sa.bed`) | **B. `CoSplitDel`** (recompute from live BAM) | **C. eKLIPse wrapper** (+ in-house fallback) |
|---|---|---|---|
| Split signal | Existing `$O.sa.bed` (SA tags only) | **MAPQ-filtered, recomputed from live `$O.bam`** (+ soft-clip capable) | eKLIPse soft-clip; `$O.sa.bed` fallback |
| Coverage corroboration | Yes (optional → promote to mandatory) | Yes (mandatory AND-gate) | Yes (corroborating) |
| Heteroplasmy | Junction + coverage | **Junction + coverage + disagreement QC** | Junction + coverage |
| Dependencies | One (`pysam`; pip manylinux wheel) | One (`pysam`) | **eKLIPse (Python 2.7 EOL) — real packaging debt** |
| Recall ceiling | **Bounded by SA tags** (no soft-clip-only junctions) | Higher (MAPQ filter + soft-clip option) | High (eKLIPse) but fragile |
| Simplicity (panel score) | **9** | 7–8 | 6 |
| Defensibility (panel score) | 7 | **8** | 7 |
| Re-runnable post-hoc over `out/`? | Yes (sa.bed+cvg persist) | Core path needs the live BAM | Partly |

**Why not the generic engines for v1:** §4b — they exclude chrM, choke on uniform high depth, and
can't handle the circular origin. **Why not eKLIPse for v1:** it's unmaintained Python 2.7 (EOL),
adds Docker/packaging debt, and detects deletions only in its published form — more fragile than an
in-house split+coverage caller whose inputs MitoHPC already computes. eKLIPse is a good **optional
v2 engine** behind the same flag (its peer-reviewed validation adds value without becoming a hard
dependency).

---

## 9. Recommended v1 — `callSV` (blended, deletion-only, dual-signal)

> **Status: implemented.** Core caller `callsv.py` (Python 3 + `pysam`) with a thin `callSV.sh`
> driver, plus `sv.vcf` header and `getSVSummary.sh`; wired via `HP_SV` (default off). `pysam` is
> installed in the Docker image and CI. The as-built method (formulas, schema, parameters) is
> documented in **[`SV_METHODS.md`](SV_METHODS.md)**. The output is general VCF/SV best practice
> (`##contig`/`##reference`/provenance, sample-named column, `HOMLEN`/`DELCLASS`/`CIPOS`, `SVCLAIM`,
> affected `GENE`s, `COMMON`, `HGVS`, cohort `bcftools merge` recurrence matrix + sites union).
> Evaluated by `bash test/sv/run_test.sh` over **21 committed mock BAMs** (~40 checks): deletions
> across a size range (45 bp→13 kb, repeat/non-repeat, multi-deletion, near-homoplasmy, D-loop,
> low-coverage), **tandem duplications and origin-crossing both correctly not PASSed**, plus
> forward-looking **duplication / inversion / complex** fixtures (design in `SV_EVENT_TYPES.md`;
> mostly 0-record / 0-PASS until those paths land), wild-type specificity, degenerate-input
> robustness (no tracebacks), cohort recurrence, a `bcftools` VCF-spec gate, and real-data vetting.

A blend that takes **B's signal source** (recompute a clean split-read signal from the **live
`$O.bam`**, since we must hook in before BAM cleanup anyway — this is the "build our own, better
signal" latitude granted in `CLAUDE.md` §3) with **A's simplicity** and the **mandatory
dual-signal gate** that makes it defensible. Deletion-only, high-specificity, fully auditable.

> Rationale for not simply reusing `$O.sa.bed`: `sam2bedSA.pl` uses **only `SA:Z:` tags** (no
> soft-clip-only reads, no MAPQ filter, no breakpoint refinement) — it's a coarse signal with a
> real recall ceiling. Because the SV step runs **while `$O.bam` is still alive**, recomputing a
> MAPQ-filtered (and optionally soft-clip-augmented) breakpoint set from the BAM costs almost
> nothing and is strictly better. `$O.sa.bed` remains a zero-dependency fallback / post-hoc option.

### 9.1 Signals & call logic
1. **Split junctions:** `samtools view -h -q $HP_SV_MINMAPQ $O.bam | <breakpoint extractor> →
   candidate junctions`. Keep same-`chrM`, same-strand segment pairs whose inner gap is a
   deletion; cluster breakpoints within `±HP_SV_PAD` (absorbs the 13-bp direct-repeat ambiguity);
   count distinct supporting read names; drop clusters with `< HP_SV_MINJR` reads.
2. **Coverage corroboration (mandatory for PASS):** from `$O.cvg`, compute median depth **inside**
   `(up,dn)` vs **flanking WT windows** (`HP_SV_FLANK`, wrapped mod `HP_MTLEN`). A cluster is
   **PASS** only if it has **both** `≥HP_SV_MINJR` junction reads **AND**
   `medInside/medFlank ≤ HP_SV_DROP`. Split-only or coverage-only candidates are emitted with a
   **non-PASS FILTER tag** (never silently dropped — auditability), recovering the sub-5% range
   where coverage loses power as a lower-confidence tier.
3. This is the classic **DELLY/LUMPY corroboration principle** — the single most defensible thing
   a short-read v1 deletion caller can do on chrM, and it neutralizes most NUMT/chimera artifacts.

### 9.2 Heteroplasmy (report both + disagreement)
- `AF_junction = J / (J + S)`, with `S = round(0.5·(depth[up]+depth[dn])) − J` from `$O.cvg`
  (v2: exact intact-spanning count from the BAM).
- `AF_cov = 1 − medInside/medFlank`.
- **Primary `AF` = `AF_junction`**; `AF_cov` and `|AF_junction − AF_cov|` (a QC/defensibility
  field) go in `INFO`. The formula is documented in the VCF header.

### 9.3 Circularity
Reuse only. `$O.bam` is built via `circSam.pl` from the circularized `chrMC` (`HP_E=300`), so
origin-crossing junctions already appear wrapped into 1..16569. Coverage windows wrap mod
`HP_MTLEN`. Origin-spanning junctions are **flagged `WRAP`** and kept out of PASS in v1 (no
deletion-vs-complementary-arc-duplication disambiguation yet — documented limitation; the
coverage-dropout gate already biases toward true deletions).

### 9.4 False-positive controls (all over existing `RefSeq/` assets)
- **Min supporting reads** `HP_SV_MINJR` (default 3); **min/max size** `HP_SV_MINSIZE` (default
  ~50–500 bp — see open questions) to `HP_MTLEN−1`; **coverage-drop** `HP_SV_DROP`.
- **MAPQ prefilter** `HP_SV_MINMAPQ` (default 20) on split reads → drops NUMT multi-mappers.
- **NUMT:** flag breakpoints overlapping `NUMT.vcf.gz`; rely on the upstream competitive NUMT
  realignment (`filter.sh:94-101`); v2 can re-score junction reads against `NUMT.fa` vs `chrMC`.
- **Homopolymers / control region:** flag/drop breakpoints in `HP.bed.gz`, flag `DLOOP.bed.gz`.
- **Direct-repeat awareness:** left-align and **`REPEAT`-flag** the canonical `del4977` using the
  `RefSeq/chrM.8455-8482_13460-13474.fa` coordinates (annotate, never drop).

### 9.5 Outputs (purely additive)
- Per sample: **`$O.sv.vcf`** (VCFv4.2; `SVTYPE=DEL`, `END`, `SVLEN`, `JR`, `SR`, `AFJ`, `AFC`,
  `AFDIFF`, `REPEAT/NUMT/HP/DLOOP/WRAP` flags; `FORMAT GT:DP:AF`) and **`$O.sv.tab`** (flat TSV).
  Run through `filterVcf.pl -sample $S -source $HP_SV -depth $HP_DP` for schema consistency.
- Cohort: a **separate `getSVSummary.sh`** (never edit `getSummary.sh`) cats per-sample files into
  `$ODIR/sv.concat.vcf` + `$ODIR/sv.tab` using the existing `awk|sed|xargs cat|uniq.pl|bedtools
  sort -header` idiom, invoked only when the flag is set.

### 9.6 Wiring (mirror `HP_V` exactly)
- **`init.sh`:** new `HP_SV` (default empty; v1 accepts one value, e.g. `callsv`) + tunables
  `HP_SV_MINJR=3`, `HP_SV_MINSIZE`, `HP_SV_PAD=5`, `HP_SV_DROP`, `HP_SV_FLANK=200`,
  `HP_SV_MINMAPQ=20`, printed/exported like the `HP_V` block.
- **`run.sh`:** add a validation line beside the `HP_V` check; emit a flag-gated `getSVSummary.sh`
  line **after** the existing `getSummary.sh` line.
- **`filter.sh`:** **one** gated block, placed immediately after the GRIDSS block and **before
  `rm -f $O.bam*` (line 232)**: `if [ $HP_SV ] ; then callSV.sh $S $O.bam $O ; fi`. No other edit.
- **New files (as built):** `scripts/callsv.py` (Python 3 + `pysam` — junction extraction,
  in-process depth, corroboration, VCF/tab), `scripts/callSV.sh` (thin driver), `scripts/sv.vcf`
  (header, mirrors `gridss.vcf`), `scripts/getSVSummary.sh`.
- **Acceptance:** (1) with `HP_SV` empty, `out/` is **byte-for-byte identical** (diff
  `examples1/`, `examples2/` before/after — wire into CI). (2) the `del4977` junction smoke test
  produces a PASS `DEL` at ~m.8470_13447 with the `REPEAT` flag set.

---

## 10. Phased roadmap

- **v1 — `callSV` (minimal, deletion-only):** live-BAM split junctions + mandatory coverage-drop
  AND-gate; two heteroplasmy estimates + disagreement QC; FP controls via existing
  `HP.bed.gz`/`NUMT.vcf.gz`/`DLOOP.bed.gz`; `del4977` `REPEAT` annotation; `WRAP` out of PASS;
  `$O.sv.vcf`/`$O.sv.tab` + separate `getSVSummary.sh`. **Bar:** byte-identical `out/` when off;
  `del4977` smoke test passes.
- **v2 — sensitivity, precision, exact heteroplasmy (still additive):** (a) soft-clip-cluster
  support (reads clipped without an `SA` partner) → recall toward eKLIPse/MitoSAlt; (b) split-only
  lower-confidence tier for sub-5% events; (c) **exact intact-spanning read count** from the BAM
  for the heteroplasmy denominator; (d) competitive NUMT re-scoring of junction reads
  (`filter.sh:94-101` trick); (e) optional **`HP_SV=eklipse`** engine behind the same flag with a
  `fixeklipseVcf.pl` shim into the identical `$O.sv` namespace; (f) un-subsampled SV path for very
  low heteroplasmy.
- **v3 — event richness & cross-platform** — full design in **[`SV_EVENT_TYPES.md`](SV_EVENT_TYPES.md)**
  (DUP/INV/complex + origin-resolution, harmonized): (a) **DEL-vs-complementary-arc-DUP disambiguation**
  (MitoSAlt OriH/OriL logic); (b) multiple/concurrent deletions + basic
  junction phasing; (c) inversions/insertions where `SA` orientation supports them; (d) optional
  **long-read engine** (`HP_SV=sniffles` over NGMLR/minimap2) for complex/duplicated
  rearrangements; (e) annotate against **MITOMAP/MitoBreak** for known-vs-novel flagging. Each is a
  new `HP_SV` value with its own shim into the frozen `*.sv.*` schema (append-only INFO/FORMAT).

---

## 11. Risks & open questions

**Risks**
- **Recall ceiling at low heteroplasmy** — the coverage-drop AND-gate suppresses sub-5% deletions
  whose dip is invisible at 2000×; deliberate specificity-over-sensitivity for v1. _Mitigation:_
  keep split-only candidates behind a FILTER tag; v2 adds the split-only tier + un-subsampled path.
- **DEL-vs-DUP ambiguity** on the circular genome — v1 reports intra-chrM junctions as `DEL`;
  deferred to v3. _Mitigation:_ `WRAP`/`DLOOP` flags + coverage gate bias toward true deletions.
- **Heteroplasmy denominator is a coverage proxy** in v1 → approximate AF, biased in the
  D-loop/uneven coverage and by `HP_L` subsampling. _Mitigation:_ report `AF_cov` + `AFDIFF`; v2
  computes an exact spanning count.
- **Residual NUMT chimeras** — v1 only flags `NUMT.vcf.gz` overlaps. _Mitigation:_ dual-signal gate
  removes most (no real dropout); v2 adds competitive re-scoring.
- **Integration-slot fragility** — the block must sit before `rm -f $O.bam*`. _Mitigation:_ a
  misordered block degrades gracefully (sa.bed/cvg persist); CI `out/` diff catches regressions.
- **Schema lock-in** — once `$O.sv.*` columns ship, evolve **append-only** (never rename).

**Open questions (tune before locking defaults)**
- Best `HP_SV_DROP` at ~2000× subsampling? (panel proposed 0.5 vs 0.95 — a 2× spread; tune against
  a spiked `del4977` dilution series.)
- Emit split-only candidates as a non-PASS tier in v1, or withhold until v2?
- `HP_SV_MINSIZE` lower bound — 50 bp vs 500 bp (cleaner separation from the SNV/indel caller)?
- Enforce a minimum flank depth (`HP_DP`) on WT windows to avoid spurious ratios near the D-loop?
- Does `WRAP` exclusion risk dropping real near-origin deletions? Is `$OR.bam` corroboration needed
  even in v1?
- Add a second, **non-repeat-flanked** synthetic deletion as a positive control so the caller
  isn't overfit to the 13-bp direct repeat?
- Cohort `sv.concat.vcf`: position-merge across samples (recurrent-deletion detection) in v1, or
  per-sample concat until v2?

---

## 12. Key references (primary sources)

**mtDNA deletion biology**
- Schon et al., *Science* 1989 — direct repeat is a hotspot for the common deletion.
- Yusoff et al. 2019 — ΔmtDNA4977 overview: https://pmc.ncbi.nlm.nih.gov/articles/PMC6478002/
- GeneReviews, *Single Large-Scale mtDNA Deletion Syndromes*: https://www.ncbi.nlm.nih.gov/books/NBK1203/

**Dedicated mtDNA SV tools**
- eKLIPse — Goudenège et al., *Genet Med* 2019: https://www.nature.com/articles/s41436-018-0350-8 · repo https://github.com/dooguypapua/eKLIPse
- MitoSAlt — Basu et al., *PLOS Genet* 2020: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009242 · https://pmc.ncbi.nlm.nih.gov/articles/PMC7769605/
- MitoDel — Marshall et al., *BMC Bioinformatics* 2017: https://pmc.ncbi.nlm.nih.gov/articles/PMC5657046/
- Splice-Break — *NAR* 2019: https://academic.oup.com/nar/article/47/10/e59/5380497 · Splice-Break2 https://github.com/brookehjelm/Splice-Break2
- MitoSV-seq 2020: https://pmc.ncbi.nlm.nih.gov/articles/PMC7334819/ · Cas9/baldur (ONT) 2022: https://pmc.ncbi.nlm.nih.gov/articles/PMC9537161/

**Generic SV callers**
- GRIDSS — Cameron et al., *Genome Res* 2017 · https://github.com/PapenfussLab/gridss
- Manta — Chen et al., *Bioinformatics* 2016 · https://github.com/Illumina/manta
- DELLY — Rausch et al., *Bioinformatics* 2012: https://pmc.ncbi.nlm.nih.gov/articles/PMC3436805/ · https://github.com/dellytools/delly
- LUMPY — *Genome Biol* 2014 · smoove https://github.com/brentp/smoove

**Long-read**
- Sniffles2 — *Nat Biotechnol* 2024: https://www.nature.com/articles/s41587-023-02024-y
- Nanopore mtDNA deletions — *Front Genet* 2023: https://pmc.ncbi.nlm.nih.gov/articles/PMC10344361/
- PacBio vs short-read mtDNA — *IJMS* 2025: https://www.mdpi.com/1422-0067/27/8/3562
- cuteSV https://github.com/tjiangHIT/cuteSV · SVIM https://github.com/eldariont/svim · dysgu https://github.com/kcleal/dysgu

**NUMTs, circular handling, population/SNV pipelines**
- NUMT confounding — Wei et al., *Front Cell Dev Biol* 2019: https://www.frontiersin.org/articles/10.3389/fcell.2019.00201/full
- gnomAD v3.1 mtDNA: https://gnomad.broadinstitute.org/news/2020-11-gnomad-v3-1-mitochondrial-dna-variants/ · Genome Research 2022 (Laricchia/Lake et al.)
- GATK Mutect2 mitochondria mode: https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels
- CircularMapper — Peltzer et al. 2016: https://github.com/apeltzer/CircularMapper
- mity (SNV/indel only, **not** SV): https://github.com/KCCG/mity

**Truth sets**
- MITOMAP Reported Deletions: https://www.mitomap.org/foswiki/bin/view/MITOMAP/DeletionsSingle
- MitoBreak (PMID 24170808)

**MitoHPC**
- Battle et al., *NAR Genomics & Bioinformatics* 2022: https://academic.oup.com/nargab/article/4/2/lqac034/6586827 · https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9112767/
