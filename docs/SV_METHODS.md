# MitoHPC mtDNA Deletion Calling — Method (as implemented)

**Status:** v1 (single large-scale deletions) · default **off** (`HP_SV`) · purely additive
**Applies to code:** `scripts/{callsv.py, callSV.sh, sv.vcf, getSVSummary.sh}` (core in Python 3 + `pysam`)
**Companion docs:** design/literature → [`SV_CALLING.md`](SV_CALLING.md) · guardrails → [`../CLAUDE.md`](../CLAUDE.md)

> ⚠️ **Maintenance contract:** this document describes the code *as it actually runs*. If you
> change the algorithm, thresholds, output schema, or wiring, update this file in the **same
> change** and add a line to the [Changelog](#changelog). Numbers and formulas here are quoted
> from the scripts and must stay in sync with them.

---

## 0. Methods abstract — two levels (for reviewers / grant text)

Two self-contained descriptions of the method, written to be lifted directly into a manuscript or
grant. **Level 1** gives a genomics reader the idea in a few sentences; **Level 2** is a dense,
parameterized methods paragraph for a preliminary-data / methods section. Everything below is
expanded, with code references, in §§1–10.

### Level 1 — intuitive (the idea in a paragraph)

We detect **large mitochondrial deletions** directly from the circular-aware chrM alignment that
MitoHPC already builds for SNV/heteroplasmy and copy-number estimation, so no additional alignment
is required. A true deletion leaves **two complementary signatures** in short-read data, and we
require **both**: (i) *split reads* — a read crossing the deletion junction aligns to the wild-type
reference in two pieces, the second recorded as a supplementary-alignment (`SA:Z`) tag, pinning the
two breakpoints (bp5 → bp3) to base-pair resolution; and (ii) a *coverage drop* — only wild-type
molecules cover the deleted span, so read depth falls across it in proportion to the mutant
fraction. We extract split-read junctions, cluster the reads that support the same breakpoint pair,
and call a deletion **PASS by either path**: a clustered junction corroborated by a coverage drop, *or*
strong, clean split-read evidence alone (since mtDNA read depth is finicky) — the split-read +
read-depth corroboration principle of general-purpose SV callers (DELLY/LUMPY/Manta), specialized to
the small, high-copy, **circular** mitochondrial genome (see §0.5 for the precedent map).
Heteroplasmy (mutant fraction) is reported as the **coverage-dosage** fraction (the depth deficit
across the deletion — the standard measure for large mtDNA deletions), with the **junction-read
fraction** as corroborating evidence; the dual-path PASS is what makes a call robust at any
sequencing depth.
Because the underlying alignment is made against a *circularized* reference (the first 300 bp
re-appended), deletions spanning the artificial linear origin are represented as split alignments
and handled natively, without ad-hoc linearization. v1 targets large single deletions; duplications,
inversions, insertions, and true small origin-crossing deletions are out of scope (the former two
appear only as negative controls — see Level 2).

### Level 2 — grant-ready methods (preliminary data)

*Mitochondrial large-deletion calling.* Large mtDNA deletions are called by a purpose-built
module (`callsv.py`; Python 3 / `pysam`) operating on the circularized, deduplicated, subsampled
chrM alignment produced by MitoHPC's existing realignment (reads mapped to a chrM reference extended
by 300 bp so origin-spanning reads map contiguously; the revised Cambridge Reference Sequence — rCRS,
NC_012920.1 — 16,569 bp). The caller integrates two complementary lines of evidence. *(i) Split-read junctions:* for each primary alignment
with mapping quality ≥ 20 (secondary/supplementary records excluded) carrying a supplementary-
alignment (`SA:Z`) tag, the primary and first supplementary segments are reconstructed from their
CIGAR strings; same-contig, same-strand segment pairs define a candidate deletion junction (bp5 =
last retained base upstream, bp3 = first retained base downstream). Reads supporting concordant
breakpoints are grouped by greedy single-linkage clustering to the cluster seed (the first member,
±25 bp), with cluster breakpoints taken as the per-coordinate mode and junction support (JR) as the
number of *distinct* supporting reads. Only deletions the aligner represents as a supplementary/split
alignment are visible, so the practical lower size limit is aligner-set (≈ hundreds of bp at 150 bp
reads); the `HP_SV_MINSIZE` = 50 bp parameter is a lower guard, not the detection limit. *(ii)
Read-depth corroboration:* per-base depth is computed in-process (`pysam` count_coverage) and a
deletion must exhibit a coverage drop — median depth across the deleted span divided by median depth
in the flanking 200 bp windows ≤ 0.9 (≥ 10% dosage loss). A deletion is reported as **PASS** only
when a clustered junction (≥ 3 distinct split reads) and a coverage drop co-occur; junction-only
events (e.g. sub-~10% heteroplasmy, where the depth dip is within noise) are retained but flagged
rather than PASS. All genomic windows wrap modulo 16,569 bp, so breakpoints and flanks straddling
the origin are computed correctly. The thresholds are configurable defaults intended for calibration
against a heteroplasmy dilution series; the coverage-drop gate (not the size/support minima) sets the
practical lower heteroplasmy limit.

Heteroplasmy is reported as the **coverage-dosage** fraction — the standard for large mtDNA
deletions (eKLIPse, MitoSAlt, Damas et al.): AFC = 1 − trimmed-median(depth inside the deletion) /
trimmed-median(depth in flanking windows), computed over windows that exclude a transition pad at
each breakpoint and mask the control region (D-loop), origin, homopolymers, and NUMT-like sites so
those fragile regions cannot fake a dosage loss. A second, complementary estimate — the junction VAF
AFJ = JR/(JR + SR), where SR is the count of wild-type reads aligned **reference-contiguously across
the breakpoint** (the true spanning population, not the whole pileup) — both corroborates the event
and serves as a depth-robust gate: a call PASSes only when the dosage drop is matched by a
proportional junction (AFJ ≥ a fraction of AFC and above an absolute floor), which rejects coverage
"bowls" with no real junction at any sequencing depth. False-positive control is layered: NUMT
(nuclear mitochondrial DNA segment) paralog reads are suppressed upstream by competitive alignment
against a NUMT reference during
MitoHPC realignment (inherited from that step); the mapping-quality filter, minimum junction
support, and the junction-vs-dosage consistency gates add specificity; breakpoints in the D-loop,
NUMT-like sites, or near the origin require strong junction support to PASS; and very large
deletions with weak junction support are rejected. v1.x reports deletions only: a tandem duplication
yields a coverage *gain* and is rejected, and true small origin-crossing deletions are flagged out of
PASS (deletion-versus-duplication is not yet disambiguated there).

Output follows general SV/VCF best practice for interoperability and reproducibility: a
spec-correct VCFv4.2 per sample with full provenance (tool and `pysam` versions, reference path,
`##contig` with sequence MD5, the exact command line, and one header line per parameter), the sample
as the genotype column, and rich annotation — breakpoint microhomology/direct-repeat length and
sequence (HOMLEN/HOMSEQ, with IMPRECISE/CIPOS/CIEND), a homology class (DELCLASS), an evidence claim
(SVCLAIM = junction and/or depth), recognition of the canonical MITOMAP common deletion (del4977;
COMMON), affected mtDNA features (GENE/NGENE), HGVS notation, and FORMAT GT:DP:AD:AF:SR. A tidy long
table accompanies each VCF. Across a cohort, calls are aggregated (`bcftools`) into a genotype matrix
with per-site recurrence (NS), a sites-only union for downstream annotation (AnnotSV/VEP), and a
self-contained, dependency-free interactive HTML report (circular and linear genome browsers with
gene / OXPHOS (oxidative phosphorylation) complex annotation, a per-position deletion-frequency track, and recurrence/summary
views). The module adds a single dependency (`pysam`, a pinned manylinux wheel bundling htslib;
containerized) and is purely additive and **default-off**: existing SNV, heteroplasmy, copy-number,
and haplogroup deliverables are byte-for-byte unchanged.

The implementation has been validated **in silico** as a proof-of-concept under idealized
conditions: paired-end reads are generated from defined mixtures of wild-type and deletion-bearing
circular genomes and aligned through the pipeline's own circular path, so the inside-vs-flank
coverage ratio equals the spiked heteroplasmy by construction. A committed regression suite (24
checks) confirms recovery of spiked deletions — including the common deletion at 30% and 5%
heteroplasmy, multiple simultaneous deletions, near-homoplasmy (95%), a control-region deletion, and
a low-coverage (40×) case — and the expected negatives: wild-type samples, tandem duplications
(coverage gain), and origin-crossing artifacts all yield zero PASS calls, degenerate inputs fail
cleanly without tracebacks, and every emitted VCF passes a `bcftools` specification gate. It also
includes real-data vetting: three healthy 1000G high-coverage chrM samples (0 PASS, specificity) and
a del4977 spiked into a real wild-type background (recovered PASS + COMMON, positive control). This
establishes algorithmic correctness, not real-world performance; quantitative benchmarking on real
data is the planned next step — a heteroplasmy × depth titration to establish sensitivity and the
limit of detection (LoD), an empirical per-genome false-positive rate (including NUMT-stress samples),
breakpoint-accuracy statistics, at least one orthogonally confirmed positive control (e.g. a single
large-scale deletion validated by long-range PCR/ddPCR/Southern blot), and a head-to-head
concordance against an established mtDNA deletion caller (e.g. eKLIPse, MitoSAlt).

> *Citations to add when this text is placed in a grant/manuscript:* common deletion (del4977) —
> Schon et al., *Science* 1989; MITOMAP — Lott et al. 2013 / mitomap.org; rCRS (NC_012920.1) —
> Andrews et al., *Nat Genet* 1999; split-read + read-depth SV precedent — DELLY (Rausch 2012),
> LUMPY (Layer 2014), Manta (Chen 2016); mtDNA-specific deletion callers for the comparison —
> eKLIPse (Goudenège 2019), MitoSAlt (Basu 2020). (Verify each before submission.)

---

## 0.5 Relationship to established methods

This caller is a deliberate recombination of two mature lineages — **mtDNA-specific deletion
callers** (eKLIPse, MitoSAlt, Damas et al. / MitoBreak) and **general short-read SV callers**
(DELLY, LUMPY, Manta, GRIDSS) — specialized to MitoHPC's small, high-copy, *circular* chrM
alignment. Nothing here is novel in isolation; the contribution is the fit and a handful of
mtDNA-specific corrections. Verify every citation before placing this in a manuscript.

| Design choice (this caller) | Closest precedent | What we take | Where we deviate / why |
|---|---|---|---|
| Breakpoints from `SA:Z` split reads + soft-clip harvesting (§4.1, §4.4) | eKLIPse (soft-clipping), MitoSAlt (split/gapped-read realignment); LUMPY/Manta split-read (`SR`) | Split/clipped reads pin breakpoints to base-pair resolution | We read the aligner's `SA` tags off the realigned `$O.bam` and *harvest* short clips only to **reinforce** an existing SA junction, never to seed one. No local re-alignment of clipped reads (eKLIPse) and no assembly (Manta/GRIDSS). |
| **Primary** heteroplasmy = masked, trimmed coverage-dosage `AFC` (§5.2) | Damas et al. 2014 / MitoBreak; eKLIPse & MitoSAlt also report a dosage fraction | `1 − inside/flank` depth ratio as the field-standard large-deletion heteroplasmy | We mask D-loop/origin/homopolymer/NUMT and exclude a transition pad so fragile regions can't fake a drop. eKLIPse (soft-clip-derived) is reported to **underestimate** heteroplasmy; our masked/trimmed `AFC` tracks the spiked fraction in-silico (del4977 @30% → 0.27) — real-data concordance is future work (§0). |
| Corrected junction VAF `AFJ = JR/(JR+SR)`, `SR` = true reference-contiguous spanning templates (§5.2) | General-caller VAF intuition (junction reads / molecules at the locus) | Junction-read fraction as a depth-robust second estimate | The fix vs a naive `JR/depth`: at mtDNA depth (8–22k×) `JR/depth ≈ 0` even for a real 30% deletion. We count only wild-type reads aligned **contiguously across the breakpoint** (`count_spanning_boundaries`), so `AFJ` tracks heteroplasmy at any depth. |
| Evidence-type separation: `SVCLAIM=J`/`DJ` PASS paths + the `JSUP`/`SRCONS`/`SRSB` lens (§4.5, §5.3) | LUMPY (exposes `SU`/`PE`/`SR`), Manta (`SR` vs `PR`); general-caller per-call quality scoring | Keep split-read support as its own, depth-independent dimension and surface it | We record *which* evidence drove a call and a depth-independent junction-quality lens. Because a direct repeat slides the breakpoint but conserves deletion size, `SRCONS` is measured on **svlen**, not position (microhomology-invariant) — an mtDNA-specific twist. We have **no** paired-end / discordant (`PE`/`PR`) channel. |
| Circular handling + del/dup ambiguity (§8, §5.4) | MitoSAlt (on a *circular* genome every split read is a deletion **or** a duplication; Basu 2020) | Treat chrM as circular; origin-spanning junctions are first-class | We **inherit** circularity from MitoHPC's `chrMC` realignment (first 300 bp re-appended) rather than re-linearizing. We do **not** yet resolve the del-vs-dup ambiguity MitoSAlt flags: origin junctions are `WRAP`-flagged out of PASS, and a duplication's coverage *gain* is rejected, not called. |
| NUMT (nuclear-mtDNA paralog) control (§5.4) | Standard mtDNA-pipeline concern; general callers rely on MAPQ/mappability | Suppress paralog reads before calling | Load-bearing defense is **upstream**: `$O.bam` was built by competitive alignment against a NUMT reference during MitoHPC realignment. In the caller, MAPQ≥20, a `NUMT` flag, and the junction↔dosage consistency gate are secondary guards. |
| Scope: single large-scale **deletions** only (§11) | eKLIPse / MitoSAlt target the same clinical class (Kearns-Sayre / CPEO / Pearson; del4977) | The clinically dominant mtDNA SV class | Narrower than MitoSAlt (also calls duplications) and the general callers (inversions/insertions/complex). |

**What we deliberately do *not* do** (honest scope, not oversight):
- **No local/breakend assembly** — unlike Manta (assembly of `SR`/`PR`) or GRIDSS (breakend assembly, ~90% precision). Breakpoints are read off the aligner's split/clip records directly.
- **No paired-end / discordant-pair signal** — unlike DELLY/LUMPY (`PE`), Manta (`PR`). Split-read + read-depth only; on small high-copy mtDNA that pair is sufficient and simpler.
- **No duplication / inversion / insertion / complex calling** — unlike MitoSAlt (duplications) and general callers. Opposite-strand `SA` (inversion signature) and coverage *gains* (tandem-dup signature) are recognized and **rejected** so they can't masquerade as deletions, but are not reported as their own events.
- **No probabilistic multi-evidence model / breakpoint likelihood** — unlike LUMPY's probabilistic integration or GRIDSS's quality model; PASS is a transparent set of thresholded gates (§5.3), chosen for auditability over a single opaque score.

**Net positioning.** Closest in spirit to **MitoSAlt** (split-read-driven, heteroplasmy that tracks
truth, circular-aware) but deletion-only and assembly-free; takes **dosage-as-primary-heteroplasmy**
from **Damas/MitoBreak** with explicit fragile-region masking; takes **split-read-support-as-its-own-
dimension** (`SVCLAIM`, `JSUP`) from **LUMPY/Manta**; and corrects the naive junction-VAF denominator
that breaks at mitochondrial depth. Head-to-head concordance against eKLIPse / MitoSAlt on real data
is named as required future work (§0) and is not yet done.

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

Heteroplasmy (the mutant fraction) is estimated **two complementary ways** — from the junction-read
fraction and from the coverage ratio (they share the breakpoint depth, so their disagreement is a QC
flag, not an independent confirmation). A deletion is **PASS** by **either** path: a clustered split
junction corroborated by a coverage drop (`DJ`), *or* strong, clean split-read evidence alone (`J` —
because mtDNA read depth is finicky); see §5.3.

This is the split-read + read-depth corroboration principle of general SV callers (read-depth is
integrated explicitly by LUMPY; DELLY/Manta establish the split-read breakpoint precedent), with the
coverage-dosage heteroplasmy of the mtDNA-specific callers (eKLIPse/MitoSAlt/Damas) — see §0.5 —
specialized to the small circular mitochondrial genome and implemented in **Python 3 with `pysam`**
(the de-facto htslib binding — the same C library behind `samtools`), so the BAM is parsed in-process
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
| `RefSeq/{HP,DLOOP}.bed.gz`, `RefSeq/NUMT.vcf.gz` — FP masks (flags only) | (cohort, via `getSVSummary.sh`) `$ODIR/{sv.tab, sv.merged.vcf.gz,` |
| `RefSeq/genes.bed.gz` — mtDNA features for `GENE`/`NGENE` annotation | `sv.sites.vcf.gz, sv.report.html}` |

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
   ┌──────────────────────┴──────── in-process (pysam); no temp files, no subprocess ──────────────┐
   │ extract_junctions():          $O.bam SA:Z split reads → clusters (+ soft-clip harvest rescan)  │
   │ per_base_depth():             pysam count_coverage(quality_threshold=0) ≡ `samtools depth -a`  │
   │ count_spanning_boundaries():  one pass → wild-type spanning reads per breakpoint (the AFJ den.) │
   │ call():                       corroborate + heteroplasmy + lens + flags → $O.sv.vcf + $O.sv.tab │
   └──────────────────────────────────────────────────────────────────────────────────────┘
```

There are **no temp files and no subprocesses** — `pysam` reads the BAM, computes depth, and reads
the reference/masks in-process. The single dependency is `pysam` (installed in the image via a
pinned pip manylinux wheel that bundles htslib; see §12). `callsv.py` runs on Python 3.8+ and is
parity-tested against the original Perl implementation (field-for-field identical on the mock BAMs).

---

## 4. Stage A — split-read junctions (`callsv.py: extract_junctions`)

Iterates `$O.bam` via `pysam.AlignmentFile.fetch(chrom)` and returns clustered junctions
`(bp5, bp3, SVLEN, JR, strand, SRCONS, SRSB)` (the last two feed the evidence lens, §4.5).

### 4.1 Which reads are used
For each alignment record, keep it only if it is a **primary** alignment carrying an `SA:Z:` tag:
- skip `is_unmapped` (`0x4`), `is_secondary` (`0x100`), `is_supplementary` (`0x800`);
- require `read.mapping_quality ≥ HP_SV_MINMAPQ` (mapping quality, MAPQ; default **20**) — a first
  NUMT/multimapper guard, near-inert on the already NUMT-competitively-realigned chrM `$O.bam` (the
  upstream competitive alignment is the load-bearing NUMT defense, §5.4);
- require `read.has_tag("SA")` (the supplementary mate was dropped upstream by `-F 0x90C`, but the
  tag remains on the primary — this is exactly the `$O.bam` representation).

> **Soft-clip harvesting recovers the SA-only `JR` bias (§4.4).** The aligner emits an `SA:Z:` tag
> only when the clipped arm is long enough; a deletion read whose clipped arm is too short (~11–15% of
> breakpoint-clipped reads on the test data, clips ≤16 bp) carries the junction as a **soft-clip with
> no SA tag**, and was dropped from `JR` (and, being clipped, also from `SR`) — a one-sided
> `AFJ = JR/(JR+SR)` deflation. The caller now **adds these clipped reads back** onto an existing
> SA-supported breakpoint (`HP_SV_MINCLIP`, §4.4), raising `JR`/`AFJ` by ~10–16% and the PASS LoD from
> ~8% to ~7% with no loss of specificity (real healthy samples stay 0 PASS).

### 4.2 Junction definition
Two segments are reconstructed: the **primary** (from its `POS`+`CIGAR`) and the **first `SA`
entry** (`rname,pos,strand,CIGAR,…`; a read with multiple `SA` segments — a complex/multi-junction
read — is truncated to its first SA segment in v1, sufficient for the single-large-deletion scope).
Reference span = sum of CIGAR ops that consume reference (`M/D/N/=/X`), 1-based inclusive. A junction
is kept only when both segments are:
- on the **same contig**, and
- the **same strand** (same orientation ⇒ deletion; opposite ⇒ inversion, **not** called in v1).

With the upstream segment = smaller start:

```
bp5   = reference END   of the upstream segment   (last retained base before the deletion)
bp3   = reference START of the downstream segment  (first retained base after the deletion)
SVLEN = bp3 - bp5 - 1                               (number of deleted bases)
```

Junctions with `SVLEN < HP_SV_MINSIZE` (default **50**) or `> mtlen-1` are discarded (small indels
are the SNV caller's territory; the upper bound excludes the whole-genome/origin artifact). The
**deleted bases are `bp5+1 .. bp3-1`** (`POS=bp5` is the last *retained* base, per the VCF `DEL`
convention; `END=bp3-1`); `GENE`/`NGENE` and `HGVS` are computed over that span.

> **Direct-repeat microhomology.** When a deletion is flanked by a direct repeat (e.g. the 13 bp
> repeat of the common deletion), the aligner places the breakpoint *anywhere within the repeat*,
> so the read shows a few bases of overlap between the two arms and `bp5/bp3` shift by up to the
> repeat length. This is expected and absorbed by clustering (`HP_SV_PAD`) and surfaced by the
> `REPEAT` flag — not an error.

### 4.3 Clustering
Each supporting read contributes `(bp5, bp3, read-id)` where read-id is the **template** name
(`query_name`) — the same unit as the `SR` spanning count, so `AFJ = JR/(JR+SR)` compares like with
like (a fragment is one molecule regardless of which mate split). Points are sorted by `(bp5, bp3)`
and grouped by **single-linkage greedy clustering**: a point joins the current cluster if it is
within `HP_SV_PAD` (default **25** bp) of that cluster's **last-added member** in both coordinates;
otherwise it starts a new cluster. (Linking to the last member, not a fixed seed, lets a breakpoint
smear spread *transitively* up to PAD per step, so a real event spread by microhomology/error is not
silently fragmented into sub-`minsupport` pieces.) Per cluster:
- `bp5`, `bp3`, `strand` = the **mode** (most frequent value) across members;
- `JR` = number of **distinct** supporting templates.

> **Caveat.** Clustering is single-linkage with a fixed `HP_SV_PAD`; a real breakpoint smear *wider*
> than PAD across a large direct repeat could still fragment one event, and conversely two genuinely
> distinct junctions closer than PAD could merge. PAD=25 comfortably covers the del4977 13 bp repeat.

Clusters with `JR < 2` are dropped here (`-minsupport 2`); the PASS threshold `HP_SV_MINJR`
(default 3) is applied later, so the 2-read tier is still visible as a non-PASS call.

### 4.4 Soft-clip harvesting (reinforce-only)
SA-only `JR` undercounts the junction (§4.1). After the SA clusters are fixed, a second pass adds
**clipped reads (CIGAR `S` soft- or `H` hard-clip) back onto an existing SA-supported breakpoint**: a
read with a trailing clip (≥ `HP_SV_MINCLIP`, default **10** bp) ending within `HP_SV_PAD` of `bp5`, or
a leading clip starting within `HP_SV_PAD` of `bp3`, is a read crossing that junction; its template is added to the
cluster's support. This **only reinforces a junction already evidenced by ≥`minsupport` split reads —
it can never create a call** (so the false-positive surface is unchanged: real healthy samples stay 0
PASS), while recovering the ~11–15% of breakpoint-clipped reads the SA tag omits. Effect: `JR`/`AFJ`
rise ~10–16% and the PASS limit of detection improves from ~8% to ~7% (del4977 @7% now PASSes). Set
`HP_SV_MINCLIP=0` to disable. These harvested reads are clipped at the boundary, so they are still
correctly *excluded* from `SR`, keeping `AFJ = JR/(JR+SR)` consistent.

> **Why harvesting makes the caller robust to the *aligner settings*.** Whether a junction read's
> short arm becomes an `SA` tag or a bare soft-clip is governed by the aligner — for `bwa mem`, the
> supplementary-alignment score threshold `-T` (default **30** → the shorter arm needs ~30 bp) and the
> clipping penalty `-L`. Production MitoHPC uses `bwa mem -Y` with both at default (`filter.sh`); the
> mocks use `minimap2 -ax sr`. A dedicated harness (`test/sv/aligntest.py`, committed result
> `test/sv/real/aligntest.tsv`) aligns a simulated del4977 under `bwa` (default, `-T20`, `-T15`,
> `-L2`) and `minimap2 -ax sr` and runs the caller. Result: **with harvesting on, every configuration
> recovers the *identical* `JR`/`AFJ`, PASS outcome, and breakpoint** (e.g. @10% het all give
> `JR=473`, PASS). With harvesting *off*, the aligner choice matters by ~5% (`bwa` default 416 vs
> `bwa -T15` 435 vs `minimap2` 414 junction reads) — exactly the short-arm reads `-T 30` leaves as
> soft-clips. **So the production `bwa` defaults are not limiting deletion detection**, and lowering
> `-T` is unnecessary; harvesting closes the gap. Re-run the harness if your read length or aligner
> differs.

### 4.5 Split-read evidence lens (`JSUP` / `SRCONS` / `SRSB`)
Split reads are **positional** evidence — they pin a breakpoint to a single base — so a handful of
reads all clipping at the *same* base is strong evidence of a real junction even when read depth shows
nothing (this is why general SV callers expose split-read support as its own dimension, e.g. LUMPY's
`SU/SR` and Manta's `SR` vs `PR`). To make low-level events **curatable apart from depth**, the caller
emits a depth-independent junction-quality lens per call:

```
SRCONS = fraction of supporting reads agreeing on the deletion SIZE (svlen) within HP_SV_SRTOL[5] bp
         — measured on svlen, NOT the absolute breakpoint, because a direct repeat slides bp5/bp3
         together but CONSERVES svlen, so a clean del4977 scores ~1.0 despite its 13 bp microhomology
SRSB   = strand balance min(+,-)/total of the supporting reads (0 = one-strand-only; 0.5 = balanced)
JSUP   = HIGH  if SRCONS ≥ HP_SV_SRMINCONS[0.7] and JR ≥ HP_SV_MINJR and SRSB ≥ HP_SV_SRMINSB[0.1]
         MOD   if SRCONS ≥ HP_SV_SRMINCONS  (clean, consistent junction but low-count or one-strand —
                                             a CREDIBLE low-level event worth manual review)
         LOW   otherwise (scattered breakpoint sizes — a likely mapping / NUMT artifact)
```

`JSUP` is a **lens, not a gate** — it does not change PASS. Use it to triage: filter to `JSUP ∈
{HIGH, MOD}` (with PASS off) to surface clean, low-heteroplasmy junctions that lack a depth signal
(the tumor del4977-at-5% regime), and to down-rank `LOW`/strand-biased clusters as artifacts. The
`sv.report.html` exposes it directly (a `JSUP` filter + a "min split reads (JR)" slider; §6.3).
Computed from the SA reads (the well-characterised evidence); soft-clip-harvested reads count toward
`JR` but not `SRCONS`/`SRSB`.

---

## 5. Stage B — coverage corroboration, heteroplasmy, flags (`callsv.py: call`)

**Input:** clustered junctions + per-base depth + reference (`--ref`) + masks. Depth comes from
`pysam.AlignmentFile.count_coverage(chrom, 0, mtlen, quality_threshold=0)` summed over A/C/G/T into
a 1-based array `dep[1..mtlen]`; with the default `read_callback="all"` (skips unmapped/secondary/
QC-fail/dup) this matches `samtools depth -a` (verified field-for-field on the mock BAMs). All
position windows **wrap modulo mtlen** so the circular origin is handled.

### 5.1 Coverage statistics (per junction) — masked, transition-excluded, trimmed
```
inside = [ bp5+TRANS+1 .. bp3-TRANS-1 ]           (deleted span, minus a TRANS pad at each end)
flank  = [ bp5-TRANS-FLANK+1 .. bp5-TRANS ]  ∪  [ bp3+TRANS .. bp3+TRANS+FLANK-1 ]
         TRANS = HP_SV_TRANS[150] ;  FLANK = HP_SV_FLANK[200]
         positions in the D-loop / origin / homopolymers / NUMT are EXCLUDED from both windows
CVGR   = trimmedMedian(inside) / trimmedMedian(flank)        (trimmed = drop 15% of each tail; 1 if 0)
```
(Full derivation + the dosage AF and fallback in §5.2. The transition pad keeps the breakpoint smear
out of the dosage; the masks remove the control-region/NUMT "bowls" that would fake a drop.)

### 5.2 Heteroplasmy (coverage-dosage primary + corrected junction VAF)
The **primary** heteroplasmy is the **coverage-dosage** estimate — the field standard for large
mtDNA deletions (eKLIPse, MitoSAlt, Damas et al.): the fractional depth loss across the deleted
span, over **masked, transition-excluded** windows so fragile regions can't fake a drop.
```
AFC = clamp( 1 − trimmedMedian(inside) / trimmedMedian(flank), 0, 1 )   # PRIMARY → reported AF
  inside = [bp5+TRANS+1 .. bp3−TRANS−1] ;  flank = FLANK bp beyond a TRANS pad at each breakpoint
  TRANS = HP_SV_TRANS[150] ;  FLANK = HP_SV_FLANK[200] ;  trimmed median drops 15% of each tail
  positions in the D-loop / origin / homopolymers (HP.bed.gz) / NUMT (NUMT.vcf.gz) are EXCLUDED from
  both windows
CVGR = trimmedMedian(inside) / trimmedMedian(flank)                     # 1 = no drop
```
> **Dosage fallback.** Each window must have ≥ `MINBASE` (fixed **50**) usable (unmasked) bases **and**
> a positive flank median for the dosage to be trustworthy. Otherwise (a deletion smaller than ~2×TRANS,
> or an over-masked interior) the call is marked *dosage-not-estimable*: `CVGR=1`, `AFC` falls back to
> `AFJ`, and `SVCLAIM=J`. So sub-~300 bp deletions are quantified by the junction, not dosage.

The **junction VAF** corroborates and is the depth-robust gate input:
```
SR  = wild-type reads aligned reference-CONTIGUOUSLY across a breakpoint (true spanning count,
      min over bp5/bp3) — NOT the v1 `depth − JR`, which used the whole pileup
AFJ = JR / (JR + SR)                                                    # corrected junction fraction
AFDIFF = |AFJ − AFC|                                                    # junction-vs-dosage QC
```
*How `SR` is counted (`count_spanning_boundaries`, ONE BAM pass for all breakpoints):* a read template
spans a boundary `B` if some aligned block covers **both** `B` and `B+1` (continuously matched across
it); a read soft-clipped or SA-split *at* `B` has a block that ends at `B` and is correctly excluded
(it carries the deletion, not the wild type). A template's primary + supplementary blocks are
**unioned** by `query_name`, so an origin-crossing wild-type read (split by `circSam.pl`) still counts;
`SR` for a junction is the **min** of the two breakpoints' counts (conservative). The spanning boundary
lies inside the `bp5..bp3` microhomology zone, so `SR` can wobble by a few reads if clustering picks the
opposite repeat edge — bounded by `HOMLEN` and benign in practice.
> **Why this changed (v2).** v1 divided junction reads by the *entire* pileup, so on real high-copy
> mtDNA (8,000–22,000×) every `AFJ` collapsed below 1% and junction-noise deletions slipped through
> PASS. v2 reports the dosage `AF` and a *correctly normalised* `AFJ`; the two now agree for real
> deletions (e.g. del4977 @30%: AFC 0.27 / AFJ 0.25) and diverge for artifacts (AFJ ≈ 0 vs AFC > 0).
> See the [Changelog](#changelog).

### 5.3 PASS / FILTER logic — two evidence paths
A call is **PASS** if **either** evidence path is satisfied (`SVCLAIM` records which). Both gate on
the **corrected** `AFJ`, so neither re-admits the depth-noise artifacts (whose `AFJ ≈ 0`).

**(DJ) depth + junction** — a coverage drop corroborated by a proportional junction. PASS when
**none** of these fire:

| FILTER reason | Condition |
|---|---|
| `lowJR` | `JR < HP_SV_MINJR` (3) |
| `no_cvg_drop` | `CVGR > HP_SV_DROP` (0.9 ⇒ requires ≥10% dosage drop) |
| `WRAP` | a breakpoint within `originpad` (20 bp) of the origin |
| `lowDP` | `HP_SV_MINDP > 0` and flank depth `< HP_SV_MINDP` (default 0 ⇒ disabled) |
| `low_dosage` | `AFC < HP_SV_MINAF` (0.03) — the dosage drop itself must clear MINAF |
| `lowAFJ` | `AFJ < HP_SV_MINAFJ` (0.02) — depth-robust junction-support floor |
| `unexplained_drop` | `AFJ < HP_SV_AFFRAC × AFC` (0.30) — the coverage drop must be junction-corroborated |
| `fragile_weakJ` | breakpoint in D-loop/NUMT/origin **and** (`AFJ < HP_SV_STRONGAFJ`[0.05] or `JR < HP_SV_STRONGJR`[10]) |
| `bigdel_weakJ` | `SVLEN ≥ HP_SV_BIGDEL` (8000) **and** (`JR < HP_SV_BIGMINJR`[8] or `AFJ < HP_SV_BIGMINAFJ`[0.02]) |

**(J) junction-strong** — clean, well-supported split reads PASS **without** a coverage drop, because
mtDNA read depth is finicky and a high-confidence junction is the highest-quality signal. PASS when
`JR ≥ HP_SV_JMINJR` (8) **and** `AFJ ≥ HP_SV_JMINAFJ` (0.05) **and** `SVLEN < HP_SV_BIGDEL` **and**
not at the origin **and** not a coverage *gain* (`CVGR ≤ 1 + HP_SV_GAINPAD` ⇒ excludes duplications)
**and** (if in a fragile D-loop/NUMT region) the junction is strong (`AFJ ≥ STRONGAFJ`, `JR ≥ STRONGJR`).
A J-path PASS carries `SVCLAIM=J`; a DJ PASS carries `SVCLAIM=DJ`. Filter to `SVCLAIM=DJ` for
depth-corroborated calls only, or keep `J` to surface clean junctions where depth is unreliable.

`unexplained_drop` + `lowAFJ` (DJ) and the `AFJ`/`JR` floors (J) are the load-bearing, depth-invariant
artifact filters: a coverage dip with no proportional junction (a NUMT / mappability / control-region
bowl, `AFJ ≈ 0`) is rejected on **both** paths. Genuine very-low-heteroplasmy deletions (where neither
a clean dosage drop nor a strong junction is present) still land in a non-PASS tier with full evidence.

**Sensitivity (titration + real-background spike).** A heteroplasmy × depth titration of del4977
(`test/sv/titration.py`; committed `test/sv/real/titration_del4977.tsv`) and a del4977 spiked into the
**real** NA12718 background (`gen_spike.sh`) agree: the deletion is **detected 100% down to 2%**
heteroplasmy (real variants are never lost — below the PASS threshold they are surfaced as non-PASS
records with full evidence), and **PASSes at ≥7%** with soft-clip harvesting on (§4.4; the J path
delivers the 7–8% tier, DJ takes over at ≥10%). `AFC` tracks the spiked fraction (0.04 @5%, 0.05 @8%,
0.08 @10%, 0.50 @50%). A **sensitive mode** `HP_SV_JMINAFJ=0.04` reaches a **~5% PASS LoD** (the 5%
spike then PASSes via the J path) and was verified to keep the healthy real samples at **0 PASS** —
the real-data artifacts sit at `AFJ ≈ 0`, far below, so there is headroom. Use it for tissue/tumor
cohorts (see the tissue note below); calibrate per cohort.

> **Tissue matters more than the caller.** The "common deletion" del4977 is a low-abundance biomarker,
> not a high-heteroplasmy clonal event in most samples: ~0.01–0.2% in **blood/buccal**, 0.0001–0.14% in
> aged **muscle**, and 0.0001–**7%** in **solid tumor** (often *lower* than adjacent normal — it is
> selected against in proliferating cells). These are **below** even a ~5% short-read-WGS PASS LoD, so
> **0 PASS large deletions in a blood/buccal cohort is the expected, correct result** — not a caller
> failure (no short-read method, incl. eKLIPse/MitoSAlt, calls a 0.1% deletion). Callable
> heteroplasmy (20–90%) occurs in **single-large-scale-deletion disease tissue** (Kearns-Sayre / CPEO /
> Pearson — muscle; adult blood is often negative). The high end of tumor del4977 (~5–7%) is exactly
> where soft-clip harvesting + the sensitive mode help.

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

### 5.5 Two orthogonal per-call scores: `SVCONF` (confidence) and `SVIMPACT` (biological impact)

Each call carries two independent 0–100 scores answering different questions. **Confidence** = "is this
a real deletion or an artifact?" **Impact** = "how functionally damaging would it be *if real*?" They are
deliberately decoupled: a low-VAF real deletion can be low-confidence yet high-impact, and a strong
control-region artifact can be high-evidence yet low-confidence. Neither feeds the PASS/FILTER decision —
both are additive `INFO`/tab fields for curation and ranking. Every input is a ratio, a saturating count,
or a categorical biology fact, so both are **depth-stable** (the property a LoD sweep needs).

**`SVCONF` — call confidence** (`callsv.py: svconf_score`). Built only from existing per-call signals:

```
SVCONF = clamp( Q + H + DJ − PEN , 0, 100)        # '.' (NA) for WRAP/origin calls
  Q  = 14·SRCONS + 10·min(SRSB/0.40,1) + 8·log1p(min(JR,20))/log1p(20)   # evidence quality, depth-independent
  H  = 40·min(het/0.30, 1),   het = AFC if dosage estimable else AFJ      # heteroplasmy magnitude (monotone)
  DJ = 16·max(0, 1 − |AFJ−AFC|/max(AFJ,AFC))   only when a real coverage drop corroborates (else 0)
  PEN= 16·[nfragile≥1] + 16·[nfragile≥2]        # fragile categories (DLOOP/HP/NUMT/WRAP) at either breakpoint
```

`H` rises monotonically with heteroplasmy (saturating at 30%); `DJ` uses *relative* junction-vs-dosage
agreement (not absolute `AFDIFF`, which grows with het and would invert the score at high VAF); `PEN`
double-weights the both-ends-fragile homopolymer signature of the recurrent control-region artifact.
Bands: **HIGH 70–100**, **MEDIUM 40–69**, **LOW 15–39**, **VERY-LOW 0–14**. Worked: del4977 @30% (DJ,
non-fragile) → **83**; the same del @5% (junction-only, low dosage) → **35**; a `cvg-gain` duplication
mis-as-DEL → **32**; the recurrent m.310_955 control-region artifact → **~0–8**.

**`SVIMPACT` — biological impact if real** (`callsv.py: svimpact_score`). mtDNA-deletion severity is
near-categorical, so the score is the **max of calibrated biology floors** plus a small continuous
tie-breaker (constraint intensity + genome fraction) that orders calls *within* a band:

| Floor | Condition |
|---|---|
| **95** | a replication origin removed (OriH m.110–441 or OriL m.5721–5798 fully deleted) → replication-incompetent |
| 70 | an origin partially overlapped |
| **78** | ≥3 tRNA/rRNA genes fully deleted (e.g. del4977) — massive translation loss |
| 62 | **any** one tRNA or rRNA fully deleted — translation-lethal (a single lost tRNA halts all mt-protein synthesis) |
| 60 | ≥2 distinct OXPHOS complexes disrupted (Complex I/III/IV/V from gene names `ND*`/`CYTB`/`COX*`/`ATP*`) |
| 45 | a full protein-coding gene deleted |
| 25 | a protein-coding gene partially deleted |

`tie = 8·intensity + 4·size_frac`, where `intensity = clamp((mean MLC over the deleted span − 0.10)/0.65, 0, 1)`
(MLC = the Yale per-base **M**itochondrial **L**ocal **C**onstraint score, `RefSeq/MLC.vcf.gz`, averaged over
its ALT rows, missing positions skipped) and `size_frac = min((SVLEN/16569)/0.5, 1)`. Bands: **SEVERE ≥80**,
**HIGH 50–79**, **MODERATE 20–49**, **LOW <20**. Worked: del4977 → **85 (SEVERE)**; an OriL-removing
major-arc deletion → **100**; a partial single-CDS nick → **~32 (MODERATE)**; a control-region-only
deletion sparing the origins → **LOW**. `SVIMPACT` uses no read counts or VAF, so it is identical for the
same deletion at any heteroplasmy or depth.

> **Calibration note.** The `SVCONF` weights and the `SVIMPACT` floor/band cut-points are expert-set
> starting values. They are intended to be tuned against a labelled set — a heteroplasmy × depth LoD
> sweep for `SVCONF` monotonicity/separation, and MITOMAP/KSS-Pearson-CPEO pathogenic single deletions
> vs benign control-region indels for `SVIMPACT` — before being treated as fixed. The structure (which
> signals, which direction) is the robust part; the numbers are the tunable part.

---

## 6. Output schema

The output follows general VCF/SV best practice (not this repo's other VCFs): a spec-correct,
tool-compatible (`bcftools`/IGV/AnnotSV), scientifically rich, reproducible single-sample VCF, plus
a tidy long table and cohort artifacts. **VCF 4.2 with negative `SVLEN` for `DEL`** (the widely
supported convention; do not switch to a positive `SVLEN` unless you also bump `##fileformat` to
4.4 — that is the one genuinely-wrong combination).

### 6.1 `$O.sv.vcf` (per sample)
`callsv.py` injects dynamic provenance/contig headers, then the static field definitions from
`scripts/sv.vcf`, then a `#CHROM` line whose **genotype column is the real sample name** (the
sample is NOT an `INFO` field). Header lines emitted:
`##fileformat`, `##fileDate`, `##source=MitoHPC_callsv <version> (pysam <v>)`, `##reference`,
`##contig=<ID=chrM,length=16569,md5=…>`, `##sample`, `##callsv_command="…"`, one
`##callsv_param_HP_SV_*` per threshold, then the `##ALT/##FILTER/##INFO/##FORMAT` definitions.

Example PASS record (del4977 @30%, from `test/sv/example/`):
```
chrM  8482  .  A  <DEL>  .  PASS  SVTYPE=DEL;END=13446;SVLEN=-4964;SVCLAIM=DJ;IMPRECISE;
   CIPOS=0,13;CIEND=0,13;HOMLEN=13;HOMSEQ=ACCTCCCTCACCA;DELCLASS=I;
   GENE=ATP8:P,ATP6:F,COX3:F,TRNG:F,ND3:F,TRNR:F,ND4L:F,ND4:F,TRNH:F,TRNS2:F,TRNL2:F,ND5:P;
   NGENE=12;COMMON;HGVS=NC_012920.1:m.8483_13446del;JR=154;SR=401;AFJ=0.277;AFC=0.268;
   AFDIFF=0.010;CVGR=0.732;SRCONS=1.000;SRSB=0.496;JSUP=HIGH;REPEAT   GT:DP:AD:AF:SR   0/1:575:401,154:0.268:154
```
(`SR=401` is the wild-type **spanning** count; `AF=0.268` is `AFC`, the coverage-dosage heteroplasmy;
`JR=154` reflects soft-clip harvesting; `JSUP=HIGH` from `SRCONS=1.0`/`SRSB=0.50` — a clean, balanced junction.)

| Field | Meaning |
|---|---|
| `POS / END / SVLEN` | `bp5` / last deleted base (`bp3-1`) / `-(deleted bases)` (negative, 4.2) |
| `SVCLAIM` | `DJ` (junction + coverage agree) or `J` (split-read only) — VCF 4.4 evidence claim. Emitted on **every** record: initialised from the dosage signal (`DJ` when a dosage drop is present and estimable, else `J`), and forced to `J` on a junction-strong PASS. On a non-PASS record it indicates the *evidence type*, not the PASS path. |
| `IMPRECISE`,`CIPOS`,`CIEND` | set when `HOMLEN>0`; CI = `0,HOMLEN` (breakpoint slides within the repeat) |
| `HOMLEN`,`HOMSEQ` | breakpoint microhomology / direct-repeat length + sequence (13 / `ACCTCCCTCACCA` for del4977) |
| `DELCLASS` | `I` (perfect repeat ≥5 bp) / `II` (1–4 bp microhomology) / `III` (none) |
| `GENE`,`NGENE` | mtDNA features deleted, `name:F` (fully) or `name:P` (partial), from `RefSeq/genes.bed.gz` |
| `COMMON` | matches the MITOMAP common deletion del4977 (m.8470_13447, within tolerance) |
| `HGVS` | approximate `NC_012920.1:m.<a>_<b>del` |
| `JR`,`SR` (INFO) | junction (split) reads / wild-type **spanning** reads (`SR` = AFJ denominator). Note: the **FORMAT** `SR` is a *different* quantity (split reads = `JR`, Manta-style); the two share the token by VCF convention |
| `AFC`,`AFJ`,`AFDIFF`,`CVGR` | **primary heteroplasmy** (coverage dosage) / junction fraction (evidence) / disagreement QC / coverage ratio |
| `JSUP`,`SRCONS`,`SRSB` | split-read evidence lens (§4.5): junction-support tier `HIGH`/`MOD`/`LOW` / size-consistency / strand balance — depth-independent, for curating low-level junctions |
| flags | `REPEAT NUMT HP DLOOP WRAP` (advisory breakpoint-region flags) |
| `FORMAT GT:DP:AD:AF:SR` | `0/1 : round(maskedFlankDepth) : SR,JR : AFC : JR` — `AF` carries the coverage-dosage heteroplasmy (`AFC`); `AD` = REF(spanning),ALT(junction); FORMAT `SR` = split reads (`JR`), Manta-style |

### 6.2 `$O.sv.tab` (tidy/long, one row per sample-deletion)
Header (parse by **name**, not position):
`sample chrom pos_bp5 end_bp3 svlen svclaim jr sr af_junction af_coverage afdiff cvgr flank_dp
homlen homseq delclass common ngene gene_list hgvs filter flags srcons srsb jsup`.
Null convention: numeric columns always populated; `homseq` empty when `homlen=0`; `gene_list`/`flags`
comma-joined (or `.` when empty); `common` is `0/1`; `srcons`/`srsb` are 3-decimal floats and `jsup` ∈
{`HIGH`,`MOD`,`LOW`} (the split-read evidence lens, §4.5). `svlen` here is the **positive** deletion
length (the VCF carries the signed `SVLEN`).

### 6.3 Cohort aggregation (`getSVSummary.sh`, gated on `HP_SV`)
Separate from `getSummary.sh` (never touched). bgzip+tabix-indexes each per-sample VCF, then writes:
- **`$ODIR/sv.tab`** — concatenated tidy long table (one header) for R/pandas.
- **`$ODIR/sv.merged.vcf.gz`** — `bcftools merge` cohort genotype matrix: one row per site, one
  column per sample, `NS` = number of samples carrying it (the **recurrence** substrate).
- **`$ODIR/sv.sites.vcf.gz`** — sites-only union (`bcftools view -G`) for annotation (AnnotSV/VEP).
- **`$ODIR/sv.report.html`** — a self-contained, offline **interactive report** (`svReport.py`,
  vanilla SVG/JS, no dependencies): a circular mtDNA map + linear genome browser with gene /
  OXPHOS-complex annotation, a per-position deletion-frequency track, VAF-coloured calls, live
  filtering (PASS / VAF / class / common / sample / **split-read evidence `JSUP` + min-JR**, §4.5),
  summary cards, VAF & size histograms, and a
  recurrence table. Follows the system light/dark theme with a manual toggle. Open it in any browser.

(No single mixed-sample concatenated VCF is produced — different sample columns can't share one VCF;
use the merged matrix or the long table. Exact-match merge can over-split imprecise breakpoints across
a cohort; positional/fuzzy merging is a future refinement.)

---

## 7. Parameters (all `HP_SV_*`)

The original eight (`HP_SV`, `MINMAPQ`, `MINJR`, `MINSIZE`, `MAXSIZE`, `PAD`, `DROP`, `FLANK`, `MINDP`)
are declared in `init.sh`; the v2+ tunables below are resolved via `callSV.sh` `${VAR:-default}`
fallbacks (so they need no edit to the frozen `init.sh` — set them in the environment to override).
Fixed, non-`HP_SV_` constants that also gate calls: cluster floor `minsupport=2` (§4.3), origin guard
`originpad=20` (the `WRAP` width), dosage-window floor `MINBASE=50` (§5.2), and the del4977 13 bp
direct-repeat windows m.8470–8482 / 13447–13459 (§5.4).

| Variable | Default | Meaning / effect |
|---|---|---|
| `HP_SV` | *(empty)* | `callsv` enables the module; empty = off |
| `HP_SV_MINMAPQ` | 20 | min MAPQ for split reads (NUMT multimapper guard) |
| `HP_SV_MINCLIP` | 10 | min soft-clip (bp) to harvest a clipped read onto an existing SA junction (§4.4); `0` disables |
| `HP_SV_MINJR` | 3 | min distinct junction reads for a PASS deletion |
| `HP_SV_MINSIZE` | 50 | min deletion size (bp); separates from small indels |
| `HP_SV_MAXSIZE` | 0 | max deletion size (bp); `0` ⇒ `mtlen-1` |
| `HP_SV_PAD` | 25 | breakpoint clustering + direct-repeat tolerance (bp) |
| `HP_SV_DROP` | 0.9 | max masked `CVGR` for PASS (≤0.9 ⇒ ≥10% drop) |
| `HP_SV_FLANK` | 200 | flank window (bp) for the coverage dosage |
| `HP_SV_MINDP` | 0 | min flank depth for PASS (0 = disabled) |
| `HP_SV_TRANS` | 150 | transition pad excluded from the dosage windows (≥ read length) |
| `HP_SV_MINAF` | 0.03 | min coverage-dosage AF (AFC) for PASS |
| `HP_SV_MINAFJ` | 0.02 | min corrected junction VAF (AFJ) for PASS (depth-robust support floor) |
| `HP_SV_AFFRAC` | 0.30 | `AFJ` must be ≥ `AFFRAC × AFC` (coverage drop must be junction-corroborated) |
| `HP_SV_STRONGAFJ` | 0.05 | junction strength to PASS in a fragile (D-loop/NUMT/origin) region |
| `HP_SV_STRONGJR` | 10 | junction-read count to PASS in a fragile region |
| `HP_SV_BIGDEL` | 8000 | "very large" deletion threshold (bp) |
| `HP_SV_BIGMINJR` | 8 | min `JR` for a very large deletion |
| `HP_SV_BIGMINAFJ` | 0.02 | min corrected `AFJ` for a very large deletion |
| `HP_SV_JMINJR` | 8 | min `JR` for a **junction-strong** (depth-independent) PASS |
| `HP_SV_JMINAFJ` | 0.05 | min corrected `AFJ` for a junction-strong PASS (lower ⇒ more sensitive, less specific) |
| `HP_SV_GAINPAD` | 0.10 | coverage-gain tolerance; `CVGR > 1+GAINPAD` ⇒ duplication, blocks the junction-strong path |
| `HP_SV_SRTOL` | 5 | bp tolerance on per-read deletion size for the split-read consistency `SRCONS` (§4.5) |
| `HP_SV_SRMINCONS` | 0.7 | min `SRCONS` for `JSUP`=MOD/HIGH (a clean, consistent junction) |
| `HP_SV_SRMINSB` | 0.1 | min strand balance `SRSB` for `JSUP`=HIGH |

> **v2 calibration caveat.** The dosage/consistency thresholds above are defaults validated on the
> simulated mocks plus real 1000G high-coverage chrM (healthy → 0 PASS); they are **not** yet locked
> by a heteroplasmy × depth titration. All are `HP_SV_*` so they recalibrate without code change.

`HP_SV_DROP` is the key sensitivity/specificity knob and is an open tuning question (see
`SV_CALLING.md` §11); it should be calibrated against a spiked dilution series.

---

## 8. Circular-genome handling

The module inherits circular-correct **detection input** from the existing pipeline rather than
re-implementing it — but note (third bullet) this is detection input, not full origin resolution:
- `$O.bam` was produced by `circSam.pl` from reads aligned to the **circularized** reference
  `chrMC` (`HP_E=300` bp appended), so a read crossing the artificial origin already has its parts
  wrapped into 1..16569.
- In `callsv.py`, every coverage window wraps modulo `mtlen`, so flanks straddling 16569/1 are
  computed correctly.
- Junction `svlen` is, however, computed **linearly** (`bp3-bp5-1`, no modulo). A deletion whose
  *deleted arc* spans the origin therefore has no positive-length linear form: it is detected-but-
  **suppressed** (negative svlen → dropped at the `minsize` gate, or `WRAP`-flagged out of PASS),
  never mis-called. So origin-spanning deletions are not yet *resolved*, only safely withheld.
- Relatedly, v1 **does not** disambiguate a circular deletion from its complementary-arc duplication;
  junctions at the origin are flagged `WRAP` and kept out of PASS (deferred to a future tier — see
  roadmap in `SV_CALLING.md` §10).

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

`make_testdata.py` simulates paired-end reads from **wild-type + event circular genomes** at a
known heteroplasmy (so the coverage ratio outside vs inside a deletion equals the spiked fraction by
construction); `gen_bams.sh` aligns them through the pipeline's own circular path
(`minimap2 -ax sr → samtools view -F 0x90C → circSam.pl → sort`) to produce faithful `$O.bam`
files (committed, ~13 MB total). `run_test.py` (invoked by `run_test.sh`) runs the caller and checks
calls against `truth.tsv`, then exercises degenerate inputs and (when `bcftools` is present) cohort
aggregation, a VCF-spec gate, and a schema check on the committed `example/` outputs. **24 checks**:

| Scenario | What it verifies |
|---|---|
| del4977 @30% / @5% | PASS + `REPEAT`/`COMMON`/`HOMLEN=13`/`DELCLASS=I`/genes; low-het → `no_cvg_drop` tier |
| non-repeat deletion @50% | PASS, no `REPEAT`, `DELCLASS` from incidental microhomology |
| **multiple deletions** | both deletions detected as separate PASS records (no merge/cross-talk) |
| **near-homoplasmy @95%** | PASS, `AFJ→1.0` (no divide-by-zero) |
| **tandem duplication** | **zero PASS** (coverage *gain*, `CVGR>1` → `no_cvg_drop`) |
| **origin-crossing deletion** | **zero PASS**, all coords ≤ contig length (valid VCF) |
| D-loop breakpoint | PASS + `DLOOP` flag |
| low coverage (40×) | still detected (cohort depth variability) |
| wild-type | 0 PASS (specificity) |
| **degenerate inputs** | empty BAM → 0 records; unindexed/wrong-contig/wrong-`mtlen` → clean one-line error, **never a traceback** |
| **cohort** | `getSVSummary.sh` builds the merge matrix + sites union; recurrence (`NS≥2`) detected |
| **VCF spec** | `bcftools view` accepts every per-sample VCF (no undefined-contig/INFO warnings) |
| **real-data specificity** | committed **1000G high-coverage** chrM (healthy: `test/sv/real/NA*.chrM.bam`) → **0 PASS** — a real-world false-positive guard (real NUMT/D-loop/error structure) the mocks cannot give |
| **real-background positive control** | **del4977 spiked into a real WT background** (`test/sv/real/spike_del4977_h20.chrM.bam`, via `gen_spike.sh`) → recovered **PASS + `COMMON`**, `AFC`≈truth, **no off-target PASS** — real error/coverage + known truth |

24 checks total (20 scenarios + 3 healthy real-data specificity + 1 del4977-into-real-background
positive control). See [`../test/sv/README.md`](../test/sv/README.md) and
[`../test/sv/real/README.md`](../test/sv/real/README.md) for layout and regeneration.

### 10.1 Quantitative LoD &amp; accuracy evaluation (`lod_sweep.py` + `lod_report.py`)

The pass/fail suite above proves the caller *works*; this evaluation **quantifies how well**, on a
heteroplasmy × depth grid, and produces the figures that defend the `SVCONF` calibration. It follows
**CLSI EP17-A2** (limit of detection as a *surface*, not a scalar) and in-silico spike-in benchmark
practice (BAMSurgeon/GIAB philosophy: known truth, realistic noise, negatives, replicate units).

**Two arms, one grid.** `test/sv/lod_sweep.py` runs every replicate through the *real* circular
pipeline path (`minimap2 → -F 0x90C → circSam.pl → callSV.sh`) and writes one tidy row per call to
`real/lod_sweep.tsv`:
- **SIM** — a simulated wild-type + event mixture at a target VAF/depth (cheap → carries the full grid
  and the confidence intervals).
- **REAL** — the same event **spiked into a real 1000G WT chrM background** (NA12718/48/75: extract WT
  FASTQ, inject *D* event read-pairs, realign) → real error/NUMT/coverage structure validates that the
  simulated grid is not a simulator artifact (sim↔real concordance, F8).
- **Variants:** `del4977` (repeat-mediated, COMMON), `NONREP` (~5 kb, no repeat — recovery is not
  repeat-specific), `ORIGIN` (origin-crossing → must be WRAP-suppressed), and `HP_ARTIFACT` =
  `del(305,965)` whose **both breakpoints sit in the D-loop poly-C homopolymer** (`nfragile=2`) — the
  dominant real-cohort false positive, used as the **hard negative** that defends the fragile penalty
  (sized > `MINSIZE` so it is actually emitted and can be scored).
- **Negatives:** the `VAF=0` column of each arm + the three healthy real BAMs (specificity / LoB).
- **Reproducibility:** every replicate's integer seed is a deterministic `crc32(variant,vaf,depth,rep)`,
  injective over the grid (asserted) and `PYTHONHASHSEED`-independent, and is written into the TSV.

**Derived metrics & figures** (`test/sv/lod_report.py` → `real/lod_report/index.html`, self-contained,
offline; matplotlib/numpy/scipy are *dev-harness* deps, never pipeline runtime deps): per-cell
sensitivity with **Wilson** CIs; **LoD50/LoD95 per depth via both probit and logistic** regression
(model-robustness) with a **cluster bootstrap over replicate units** for the CI; ROC + **precision-recall
(AUPRC headline, class-imbalanced)** and **MCC** at the PASS operating point; **calibration** reported
honestly as the raw reliability diagram + ECE/Brier **and** an **isotonic (PAVA) recalibration** on a
held-out split (so the raw hand-weighted `SVCONF` is treated as a *ranking* score and only becomes a
probability through the documented map); heteroplasmy accuracy as per-VAF bias + Bland-Altman (with the
low-VAF `AFC`→0 censoring flagged, and the **empirical** PASS rate — not the near-separable model fit —
as the primary LoD read-out), plus a **false-positive / precision** analysis — a confusion matrix
(precision, recall, FPR, F1, MCC) at the FILTER decision vs FILTER+`SVCONF`, showing the FILTER alone
admits the control-region artifact (precision ≈0.83, FPR ≈0.73 against that hard negative) and the
confidence gate recovers precision (≈0.91, FPR ≈0.34) while genuine wild-type yields **zero** calls.
Ten figures, **numbered in the order presented** and each carrying a plain-language caption (the number
lives only in the caption, so reordering never desyncs the labels): the LoD surface (**detection** rate
*and* PASS rate over heteroplasmy × depth — detection extends below the PASS threshold), the
probit/logistic dose-response + LoD95, Bland-Altman heteroplasmy accuracy, the control-region-artifact
false-positive behaviour vs spike level, `SVCONF`-vs-heteroplasmy monotonicity/depth-overlap, a direct
true-deletion-vs-artifact confidence comparison (with IQR error bars), TP-vs-artifact separation, ROC+PR,
calibration, and sim↔real concordance. Per-cell rates with Wilson CIs are committed in `lod_cells.tsv`.

> The committed `real/lod_sweep.tsv` + `real/lod_report/` come from the **`--quick`** grid (a tractable
> multi-replicate run, regenerable in one command). The publication-grade **`--full`** grid (14 VAF × 5
> depth × ≥30 replicates + 3 spike backgrounds) is documented for a cluster run; both write the same
> schema and the report regenerates identically.

### 10.2 `SVCONF` — what each term is and **why** it is present (reviewer view)

`SVCONF = clamp(Q + H + DJ − PENALTY, 0, 100)` (`.`/NA for WRAP/origin). Each term answers a specific
failure mode, and each is demonstrated by a specific figure above:

| Term | What it is | Why it is present | Shown by |
|---|---|---|---|
| **Q** evidence quality | `14·SRCONS + 10·min(SRSB/0.40,1) + 8·log1p(min(JR,20))/log1p(20)` | A true junction has size-consistent (`SRCONS`≈1), strand-balanced (`SRSB`≈0.5) split reads; homopolymer/mapping artifacts give inconsistent sizes and/or one-strand clips. `JR` is **log-saturated at 20** so ultra-high mtDNA depth cannot inflate confidence (depth-stability). | F4, F5 |
| **H** heteroplasmy magnitude | `40·min(het/0.30, 1)`, `het=AFC` (else `AFJ`) | Confidence must **rise with heteroplasmy** (more mutant molecules ⇒ more believable), expressed as a depth-invariant **ratio** (not a count) and ceilinged at 30% so one term can't dominate. The central claim the LoD sweep validates. | F3, F6 |
| **DJ** junction↔dosage agreement | `16·max(0, 1 − \|AFJ−AFC\|/max(AFJ,AFC))`, only with a corroborating coverage drop | A real deletion makes the junction VAF and the coverage-dosage AF **agree** (two orthogonal estimators of one molecular fraction); an artifact often has a junction with no proportional depth drop. **Relative-normalized** so it doesn't grow with het — the fix that keeps `SVCONF` monotone. | F3, F5 |
| **PENALTY** fragile demotion | `−16 if nfragile≥1, −16 more if nfragile≥2` (DLOOP/HP/NUMT/WRAP) | Targets the **dominant real false positive** — low-VAF control-region homopolymer pseudo-deletions, which trip both DLOOP and HP (`nfragile=2` ⇒ full −32). The `HP_ARTIFACT` hard-negative panel is the evidence it earns its points (origin/WRAP calls are additionally forced to `SVCONF='.'`). | F4, F5 |

`SVCONF` is a **ranking/confidence** score, orthogonal to heteroplasmy (`AF`) and to biological impact
(`SVIMPACT`); its 0–100 value becomes a calibrated probability only through the F6 isotonic map. Weights
are expert-set starting points to be re-tuned on the `--full` grid.

---

## 11. Known limitations (v1)

- **Recall floor.** Uses `SA:Z:` split reads only (no soft-clip-only clustering, no local
  assembly); very low-heteroplasmy junctions with few split reads can be missed. The coverage-drop
  PASS gate intentionally relegates sub-~10% events to the `no_cvg_drop` tier.
- **Minimum size is aligner-bounded.** Only deletions the aligner represents as a *split read*
  (`SA` tag) are seen; smaller deletions that fit inside one gapped alignment (CIGAR `D`) are not
  detected, so the effective floor (~hundreds of bp with 150 bp reads) is set by the aligner, not
  by `HP_SV_MINSIZE`. (v2: also harvest large CIGAR-`D` operations.)
- **DEL vs DUP / origin-crossing.** No origin-of-replication logic. A tandem duplication yields a
  coverage *gain* (`CVGR>1` → `no_cvg_drop`, never PASS), and an origin-crossing deletion is
  reported as its large complementary arc — also `no_cvg_drop` (no coverage drop in the claimed
  span), never PASS. SA coordinates in the chrMC extension are wrapped into `1..mtlen` so VCF
  `POS`/`END` always stay within the contig. The true small origin-crossing deletion is not yet
  resolved (deferred to DEL/DUP disambiguation).
- **Heteroplasmy estimators have known biases.** The primary `AFC` (coverage dosage) is biased for
  **overlapping** deletions (a second event's overlap deepens the dip — see `sv_multidel`, where the
  junction `AFJ` is the accurate per-deletion estimate) and is unavailable for deletions smaller than
  ~2×`HP_SV_TRANS` (the interior window vanishes → junction-only `AFJ`). `AFJ`'s denominator `SR`
  counts SA-representable spanning reads, so very low-heteroplasmy short-arm events can be
  under-supported. Reporting both `AFC` and `AFJ` (+ `AFDIFF`) exposes these; thresholds are not yet
  titration-locked (§7 caveat).
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
sample and `$ODIR/{sv.tab, sv.merged.vcf.gz, sv.sites.vcf.gz, sv.report.html}` for the cohort.

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

- **v2.6 (alignment-settings robustness test):** addressed the question "are MitoHPC's `bwa mem`
  settings stopping us from *seeing* deletions via split reads?" with a dedicated harness
  (`test/sv/aligntest.py`, committed `test/sv/real/aligntest.tsv`) that aligns a simulated del4977
  under `bwa mem` (production default `-T30/-L5`, plus `-T20`/`-T15`/`-L2`) and `minimap2 -ax sr`,
  then calls. Result: **with soft-clip harvesting on, every configuration recovers the identical
  `JR`/`AFJ`, PASS, and breakpoint** (@10% all give `JR=473`, PASS); the aligner's `-T` threshold only
  matters with harvesting *off* (~5%: `bwa` default 416 vs `bwa -T15` 435 vs `minimap2` 414). So the
  production `bwa` defaults are **not** limiting detection and `-T` needs no change — §4.4 documents
  this. Docs/test only; no caller change.
- **v2.5 (doc sync + framing against established callers):** added **§0.5 "Relationship to established
  methods"** — a comparison table positioning the caller against the mtDNA-specific tools (eKLIPse,
  MitoSAlt, Damas/MitoBreak) and general SV callers (DELLY, LUMPY, Manta, GRIDSS), mapping each design
  choice to its precedent and stating, honestly, what we deliberately don't do (no assembly, no
  paired-end channel, no duplication/inversion calling, no probabilistic model). Re-synced the docs to
  the v2.3/v2.4 code via a verification pass: refreshed the §6.1 example record (JR 133→154, +
  `SRCONS`/`SRSB`/`JSUP`) and the §6.2 tab header (+ 3 columns); documented the `count_spanning_boundaries`
  SR algorithm, the AFC `MINBASE=50` fallback, the `SVCLAIM` two-stage semantics, and the fixed
  call-gating constants; corrected "two independent"→"complementary", the single-path PASS framing→two
  paths (DJ/J), and "20 checks"→24; broadened the VCF `PASS` FILTER description. No algorithm change.
- **v2.4 (split-read evidence lens — `JSUP`/`SRCONS`/`SRSB`):** split reads are *positional* evidence
  (single-base breakpoint), so a few reads all clipping at the same base is strong even when depth is
  silent (cf. LUMPY `SR`, Manta `SR` vs `PR`). The caller now emits a depth-independent junction-quality
  lens per call: `SRCONS` (fraction of reads agreeing on the deletion *size* — microhomology-invariant,
  so del4977 scores ~1.0 despite its 13 bp repeat), `SRSB` (strand balance; one-strand-only flags
  artifacts), and a `JSUP` tier (HIGH / MOD = *credible low-level* / LOW = scattered artifact). It is a
  **lens, not a gate** (PASS unchanged). The `sv.report.html` gains a `JSUP` filter + a "min split reads
  (JR)" slider so a handful of consistent low-heteroplasmy junctions (the tumor del4977-at-5% regime)
  can be curated apart from depth-supported calls. New `HP_SV_SRTOL/SRMINCONS/SRMINSB`; new INFO
  `SRCONS/SRSB/JSUP` + 3 tab columns. Suite 24/24; mock del4977 → JSUP=HIGH, SRCONS=1.0.
- **v2.3 (soft-clip harvesting — low-heteroplasmy sensitivity):** prompted by a real cohort
  (blood/buccal/tumor) showing 0 PASS deletions. Research confirmed this is expected for blood/buccal
  (del4977 ~0.01–0.2%) and that tumor del4977 (~0.0001–7%) sits just below the old ~8% LoD. The
  SA-tag-only `JR` was dropping ~11% of breakpoint-clipped reads (clips ≤16 bp); `extract_junctions`
  now **harvests soft-clipped reads onto an existing SA-supported junction** (`HP_SV_MINCLIP`=10,
  reinforce-only — cannot create a call), raising `JR`/`AFJ` ~10–16% and the **PASS LoD from ~8% to
  ~7%** (del4977 @7% now PASSes; @5% via the documented sensitive mode `HP_SV_JMINAFJ=0.04`). Verified:
  mock positives still PASS (del4977 AFJ 0.25→0.28), mock negatives + 3 real healthy samples still
  **0 PASS** (specificity unchanged), suite 24/24. Detection is 100% to 2% throughout. Added a
  tissue-prevalence note (§5.3) so "0 PASS in blood" is read correctly, not as a caller failure.
- **v2.2 (split-read robustness, real-background positive control, abbreviation/doc audit):** an
  adversarial multi-agent audit (abbreviations, doc-vs-code, split-read expert review, real-data
  validity) drove: (code) `JR` now dedupes on the **template** (`query_name`) like `SR`, so
  `AFJ=JR/(JR+SR)` compares like-with-like; clustering links to the **last-added member** (transitive
  single-linkage), not a fixed seed, removing a silent false-negative where a wide breakpoint smear
  fragmented into sub-`minsupport` pieces (mock results byte-identical). (vetting) a **del4977 spiked
  into a real chrM wild-type background** (`test/sv/real/gen_spike.sh`) is a realistic positive control
  — recovered PASS + `COMMON`, `AFC` ≈ truth, no off-target FP — committed + asserted in the suite
  (now **24 checks**). (docs) expanded every reader-facing abbreviation (NUMT, rCRS, OXPHOS, LoD, …);
  fixed doc-vs-code drift (FORMAT `AF`=`AFC`, the §5.1 masked/trimmed `CVGR`, the example record);
  documented the SA-tag-only `JR` one-sided `AFJ` deflation, first-SA-only, and the clustering caveat.
  Algorithm behaviour on the mocks unchanged; default-off unchanged.
- **v2.1 (junction-strong PASS path + sensitivity titration):** adds a second, depth-independent
  PASS path so clean, well-supported split reads PASS **without** a coverage drop (`SVCLAIM=J`) —
  mtDNA read depth is finicky, and a high-confidence junction is the highest-quality signal. Gated on
  the corrected `AFJ` (`HP_SV_JMINJR`[8] / `HP_SV_JMINAFJ`[0.05]), excludes origin-crossing arcs
  (`SVLEN<HP_SV_BIGDEL`), coverage *gains* (`HP_SV_GAINPAD` ⇒ duplications), and fragile-region-weak
  calls, so it does **not** re-admit the real artifacts (`AFJ ≈ 0`) — verified: mock negatives and
  real healthy samples stay 0 PASS. Detection is unchanged (`extract_junctions` untouched); the J
  path only *adds* sensitivity for real deletions whose dosage drop is lost in depth noise. New
  `test/sv/titration.py` (heteroplasmy × depth) measures the PASS LoD (~5–10% at 1–4k×) and confirms
  `AFC` tracks the spiked fraction (≈0.11/0.20/0.51 @10/20/50%). Suite unchanged at 22 checks.
- **v2.0 (heteroplasmy + specificity overhaul — real-data driven):** fixes two defects exposed on a
  real 1384-sample cohort and reproduced on 1000G high-coverage chrM. (1) **Junction VAF was
  structurally wrong:** `SR = total_pileup_depth − JR` put the whole pileup (thousands ×) in the
  denominator, so `AFJ ≈ JR/depth ≈ 0` at mitochondrial depth and every reported VAF read <1%. Now
  `SR = count of wild-type reads aligned reference-contiguously across the breakpoint` (`callsv.py:
  count_spanning_boundaries`, primary+supplementary unioned for origin-crossing reads), so `AFJ` tracks
  heteroplasmy at any depth. (2) **Primary AF is now coverage-dosage** (`AFC`), computed as a
  trimmed median over D-loop/origin/HP/NUMT-masked, transition-excluded windows (eKLIPse/MitoSAlt/
  Damas convention) — `FORMAT/AF` carries it. (3) **PASS now enforces junction↔dosage consistency**
  via depth-robust gates `low_dosage`/`lowAFJ`/`unexplained_drop`/`fragile_weakJ`/`bigdel_weakJ`
  (5 new `HP_SV_*` knobs + 5 new FILTER ids), rejecting coverage "bowls" with no proportional
  junction. Validated: real healthy samples 6/10/1→**0 PASS**; mock positives keep PASS with correct
  VAF (del4977 @30% AFC 0.27/AFJ 0.25; @50%/@95%/D-loop all recovered). New real-data litmus assets
  + specificity test under `test/sv/real/` (suite now 22 checks). `sv.vcf`/`callSV.sh`/`svReport.py`/
  docs updated; `callsv.py` only — no frozen file touched; default-off unchanged.
- **v1.4 (reviewer/grant abstract + provenance):** added §0, a two-level "Methods abstract"
  (Level 1 intuitive + Level 2 grant-ready preliminary-data text) written to be lifted into a
  manuscript/grant, with honest scope and validation framing (in-silico proof-of-concept + an
  explicit real-data validation plan; heteroplasmy estimates described as *complementary* — they
  share the breakpoint depth signal — not "independent"; the aligner-bounded effective size floor
  surfaced; deletions-only scope stated up front). Verified the whole doc against the code via an
  adversarial multi-agent pass and fixed stale references (no `sv.concat.vcf` is produced; §2 and
  §12 now list the real cohort outputs incl. `sv.report.html`). **Provenance fix:** the VCF
  `##callsv_param_*` header now names the real env var (`HP_SV_MINDP`, was the literal-uppercased
  `MINDEPTH`); `callsv.py` only — algorithm, thresholds, and the rest of the schema unchanged;
  committed `example/` regenerated; tests remain 20/20. Default-off behavior unchanged.
- **v1.3 (interactive cohort report):** `scripts/svReport.py` builds a self-contained, offline
  interactive `sv.report.html` (circular mtDNA + linear genome browser, gene/OXPHOS-complex
  annotation, per-position deletion-frequency map, VAF-coloured calls, live filtering, summary
  stats, recurrence table); emitted by `getSVSummary.sh`. Committed example outputs under
  `test/sv/example/` (per-sample VCF/tab, cohort VCFs, report) via `make_example.sh`. Test suite →
  20 checks (adds `html_report` + 3 example-schema gates). No dependency beyond Python stdlib.
- **v1.2 (best-practice output + cohort robustness):** rich, spec-correct VCF — `##contig`/
  `##reference`/provenance headers, **sample-named genotype column** (dropped `INFO/SM`),
  `HOMLEN`/`HOMSEQ`/`DELCLASS`/`IMPRECISE`/`CIPOS`/`CIEND` (breakpoint microhomology), `SVCLAIM`,
  `COMMON` (del4977), `GENE`/`NGENE` (affected mtDNA features), `HGVS`, `FORMAT GT:DP:AD:AF:SR`;
  tidy long `$O.sv.tab`. Cohort `getSVSummary.sh` now builds a `bcftools merge` matrix
  (`sv.merged.vcf.gz`, `NS` recurrence) + sites union (`sv.sites.vcf.gz`) + long `sv.tab`. Caller
  fix: **wrap SA-tag coordinates** into `1..mtlen` (origin-crossing reads no longer emit
  out-of-contig `POS`/`END`). Tests expanded to 16 checks (multiple deletions, near-homoplasmy,
  tandem-dup-not-called, origin-crossing, D-loop, low coverage, degenerate inputs, cohort
  recurrence, bcftools spec gate) via a Python harness `test/sv/run_test.py`. This is a deliberate
  **schema change** from v1.1 (so the "field-for-field parity with perl" claim now applies only to
  the core numeric fields, not the VCF/tab layout). Default-off behavior unchanged.
- **v1.1 (Python/pysam port):** reimplemented the two Perl cores (`sa2del.pl`, `svCall.pl`) as a
  single Python 3 + `pysam` module `scripts/callsv.py`; `callSV.sh` is now a thin driver
  (`HP_PYTHON` override). BAM iteration, SA/CIGAR parsing, and per-base depth (`count_coverage`,
  `quality_threshold=0`) run in-process — no Perl, no `samtools`/temp-file shelling in the SV path.
  **Field-for-field parity** with the Perl v1 verified on the mock BAMs; algorithm, thresholds, and
  output schema unchanged. `pysam` added to the Docker image + CI.
- **v1 (initial):** single large-scale deletion caller — `SA:Z:` split-read clustering +
  coverage-drop corroboration and dual heteroplasmy estimates (originally `sa2del.pl`/`svCall.pl`),
  additive `HP_SV` wiring, cohort `getSVSummary.sh`, and the `test/sv/` mock-data harness.
