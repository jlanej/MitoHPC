# SV event types — unified design for calling deletions, duplications, inversions & complex events

**Status:** the single harmonized spec for extending the MitoHPC SV caller beyond deletions.
**Deletions are implemented and validated.** Duplications, inversions, complex (dup-del / inverted-dup)
events, and circular DEL-vs-DUP *resolution* are **deferred — designed here, not yet built.** Every
addition is opt-in and additive (frozen DEL behavior unchanged). This doc is the place all future
event-type work is harmonized; the as-built deletion method lives in [`SV_METHODS.md`](SV_METHODS.md),
the literature/roadmap in [`SV_CALLING.md`](SV_CALLING.md).

The test fixtures that pin each class (including the forward-looking ones that don't fully resolve yet)
are catalogued in [`../test/sv/TEST_BAMS.md`](../test/sv/TEST_BAMS.md) and §9 below.

---

## 1. The unifying idea — orientation × coverage

Every simple SV is one **breakpoint junction** plus a **coverage signature**, and the *combination*
classifies it. A short-read caller has exactly two orthogonal signals — the **split-read junction
orientation** (from `SA:Z` strand + reference order) and the **coverage change** across the affected
arc — and each event type is a distinct point in that 2-D space:

| Event | Split-read junction (the two segments) | Coverage over the arc | Linear `svlen` | How to classify |
|---|---|---|---|---|
| **Deletion** | same-strand, **forward** ref order (seg A ends `bp5`, seg B resumes at `bp3 > bp5`) | **drop** → `~base·(1−AF)` | `bp3−bp5−1 > 0` | junction **and/or** dosage drop *(implemented)* |
| **Tandem duplication** | same-strand, **reverse** ref order (the "everted"/backward jump: anchor at the higher coord, clip to the lower) | **gain** → `~base·(1+AF)` | negative (linear) | reverse-order junction **co-occurring with a gain** |
| **Inversion (balanced)** | **opposite-strand** (`SA` strand ≠ primary strand; FLAG `0x10` flips) — two reciprocal junctions (`INV3`/`INV5`) | **none** (copy-number-neutral, ratio ≈ 1) | n/a | opposite-strand junction(s), **no** coverage change |
| **Inverted dup (fold-back)** | **opposite-strand** | **gain** | n/a | opposite-strand junction **+ a gain** (≠ balanced INV) |

The single most important consequences:
- **A junction alone never decides DEL vs DUP** — on a circle a junction `(a,b)` is *simultaneously* a
  deletion of arc `a→b` and a duplication of the complementary arc `b→a`. **Coverage direction
  adjudicates** (drop ⇒ DEL, gain ⇒ DUP); near the origin, where the local window can't, the
  origin-preservation rule does (§5).
- **Inversions are a strand problem, not a distance problem.** They are copy-number-neutral, so depth
  gives *zero* signal; only the junction VAF is usable. (And our caller is currently **blind** to them
  — see §2.)

---

## 2. Current caller behavior, per class (the starting point)

| Class | What `callsv.py` does today |
|---|---|
| Deletion | Fully supported (junction + dosage; two PASS paths). The only PASS-able class. |
| Tandem dup | The reverse-order junction *is* found, but the deletion-shaped record is **blocked** from PASS by the coverage-**gain** guard (`cvg_gain = ratio > 1+gainpad` → `no_cvg_drop`; `j_pass` requires `not cvg_gain`). → a **non-PASS DEL-shaped record**. |
| Inversion | **Architecturally invisible.** `extract_junctions` (`if sref != chrom or sstrand != strand: continue`) **discards every opposite-strand `SA` segment**, so inversions emit **zero records** — not even a mis-classified one. |
| Complex (dup-del / inv-dup) | The internal **del** junction may surface as a (likely non-PASS) record; the inverted arm is invisible; net ≈ 0 PASS. |
| Origin-crossing / majority-arc | Reported as the linearized complement, `WRAP`-flagged, non-PASS, `SVCONF='.'` (detected, withheld, never resolved — §5). |

So the forward-looking test BAMs (§9) should mostly be **0 PASS / non-PASS / 0-record today**, by design.

---

## 3. Per-event method (detection + classification)

### 3.1 Deletion — implemented (reference only)
Same-strand forward-order junction; coverage drop. PASS via **DJ** (dosage drop + proportional
junction) or **J** (strong clean junction alone, `svlen < BIGDEL`, not `WRAP`, not a gain).
`AFC = 1 − inside/outside`, `AFJ = JR/(JR+SR)`. See `SV_METHODS.md` §4–5. **No change needed.**

### 3.2 Duplication — `SVTYPE=DUP` (deferred)
**Signature:** the mirror of a deletion — a **reverse-order** junction (segments in *decreasing*
reference order; clipped pieces map **inside** `[bp5,bp3]`) **co-located with a coverage GAIN**.
**Both** are required: a reverse-order junction *alone* is the DEL-of-complement reading (§5); a gain
*alone* without a junction is more likely a NUMT artifact.

**Heteroplasmy must flip sign:** `AFC_dup = inside/outside − 1` (the deletion's `1 − inside/outside`
would be *negative* on a gain). `AFJ = JR/(JR+SR)` is unchanged. The identity `1 + h_dup = 1/(1 − h_del)`
holds (a 50% tandem dup shows `CVGR ≈ 1.5`, **not** 0.5).

**Subtypes** (increasing difficulty): **tandem** (single adjacent junction, gain co-localized) →
**partial / dup-del** (a full monomer + an inserted partial arc; the dup junction sequence is
*identical* to its complementary deletion's — only coverage direction adjudicates; the clinically
dominant mtDNA class — KSS/Pearson/SLSMDS) → **dispersed** (donor gain spatially **decoupled** from the
acceptor-site split reads; NUMT-shaped; hardest).

**Implementation sketch:** branch in the classifier on `reverse-order ∧ gain` → emit `SVTYPE=DUP`,
`SVCLAIM=DJ`, sign-aware `AFC`. Keep tandem first; dispersed/dup-del behind further work.

### 3.3 Inversion — `SVTYPE=INV` (deferred; the biggest blind spot)
**Signature is STRAND, not distance:** a read spanning an inversion breakpoint splits into segments of
**opposite** orientation (`SA` 3rd field ≠ primary strand; FLAG `0x10`). Read pairs go **FF/RR**
(`INV3`/`INV5`; DELLY `3to3`/`5to5`; Manta `INV3`/`INV5`) — a balanced inversion has **two reciprocal
junctions**. **Copy-number-neutral** (`CVGR ≈ 1`), so depth/`AFC` give *nothing*; **`AFJ` is the only
signal**, and `AFDIFF` (AFC-vs-AFJ) will be large by construction.

**The blocker is one line:** the same-strand filter in `extract_junctions`. An INV path needs a
**separate opposite-strand branch** that (a) does **not** route into the colinear `svlen=bp3−bp5−1`
formula, (b) does **not** apply the coverage-drop gate (it would reject every true inversion), (c)
scores on `AFJ` only, (d) reports a `CIPOS/CIEND` window (breakpoints sit inside the mediating inverted
repeat with microhomology, so a point estimate is wrong).

**Posture: DETECT-AND-FLAG.** mtDNA inversions are biologically rare and usually low-VAF/artifactual
(one canonical pathogenic case, a 7-nt ND1 inversion, PMID 10775530; otherwise inverted-repeat /
replication-dependent aging accumulation, PMID 24040073). Default **non-PASS** unless high `AFJ` +
clean masks + sufficient supporting split reads on **both** segments at adequate MAPQ. Emit `SVTYPE=INV`
(with `INV3/INV5`) or paired `BND`, `SVCLAIM=J` (a junction/adjacency claim with **no** abundance
evidence — exactly true for a CN-neutral event).

**Watch:** an opposite-strand junction **with a gain** is an **inverted-DUP** (fold-back), not a
balanced inversion — classify it as `INVDUP`, not `INV`.

### 3.4 Complex — dup-del, inverted-dup, multi-event (deferred)
Compound molecules combine the above signatures and **stress clustering separation**. Rules:
- Keep the dup junction and the del junction as **separate clusters**; never let an embedded deletion
  PASS spuriously when global coverage is a **net gain**.
- `INVDUP` (Sniffles2 nested type) = opposite-strand `SA` **+** gain.
- Multi-event samples confound dosage on overlapping arcs → quantification must accept **`AFJ` *or*
  `AFC`** per event, and truth must mark shared spans.

---

## 4. Heteroplasmy quantification — make `AFC` sign-aware

One change unifies DEL/DUP dosage: `AFC = |1 − inside/outside|` with the **sign** taken from the
junction orientation —
- DEL (drop): `AFC = 1 − inside/outside`
- DUP (gain): `AFC = inside/outside − 1`
- INV (neutral): `AFC` is **undefined/NA**; report `AFJ` only and set an `INV`/CN-neutral flag.

`AFJ = JR/(JR+SR)` is event-type-agnostic and always reported. `AFDIFF = |AFC − AFJ|` stays the QC
field (a large `AFDIFF` is *expected* for INV and a red flag for a depth-only DUP without a junction).

---

## 5. The circular DEL-vs-DUP resolution (origin-preservation)

> *(This was the former `SV_DELDUP_RESOLUTION.md`; it is one case of the §1 "coverage adjudicates,
> except near the origin" principle.)*

**Problem.** A single junction is consistent with both a deletion of one arc and a duplication of the
complement (indistinguishable from short reads — MitoSAlt, Basu et al., *PLOS Genetics* 2020,
[PMC7769605](https://pmc.ncbi.nlm.nih.gov/articles/PMC7769605/)). Coverage normally adjudicates, **but
when the event spans ~half the genome the local inside/flank window can't** — and an origin-crossing
deletion of a *small* arc is linearized as its near-whole-genome complement (`svlen > MTLEN/2`).

**Today.** Such junctions are `WRAP`-flagged, non-PASS, `SVCONF='.'` — detected, never miscalled, never
resolved (`SV_METHODS.md` §8).

**Resolution (MitoSAlt rule).** Parameterize `HP_SV_ORIH` (≈110–441) and `HP_SV_ORIL` (5721–5798).
For the junction's two complementary readings, **default to DELETION**; **reclassify as DUPLICATION of
the complement iff the deletion's surviving molecule excises a replication origin** (favor the reading
where both origins are intact; failing that, where neither is deleted). Compute both arcs `mod HP_MTLEN`;
corroborate with coverage.

**Worked pair** (origins `OriH=110–441`, `OriL=5721–5798`):

| junction | deletion molecule | OriH | OriL | resolved call |
|---|---|---|---|---|
| `16400→100` | retains `100..16400` | intact | intact | **DEL** of the ~268 bp origin arc |
| `16400→200` | retains `200..16400` | **clipped** (110–199) | intact | **DUP** of the ~16 kb body |

A 100 bp breakpoint shift flips DEL↔DUP purely by sparing vs clipping OriH. These two mocks
(`sv_del_origin_spares` = the 16400→100 sparing case, and `sv_origin` = the 16400→200 OriH-clipping
case) are the origin-resolution regression fixture.

**Limits.** A heuristic, not ground truth: undecidable when *both* readings spare origins or *neither*
does — keep a withhold (`WRAP`) path for that residue. True structural resolution wants long reads.

---

## 6. Pitfalls to handle (deduplicated)

**Inversion**
- Opposite-strand `SA` is discarded at the same-strand filter → 0 records; and CN-neutral → no dosage.
  *Fix:* opposite-strand branch, `AFJ`-only scoring, no coverage-drop gate.
- Breakpoints sit *inside* the mediating inverted repeat (2–5 bp microhomology) → report `CIPOS/CIEND`,
  not a point; future INV match-tolerance ≥ IR length.
- Opposite-strand chimeras are exactly what NUMTs/palindromes mis-produce, and real INV is low-VAF
  (high FP regime) → apply `HP`/`DLOOP`/`NUMT` masks to INV breakpoints, MAPQ floors on **both**
  segments, hard min split-read count + min `AFJ`. Don't lower these to chase sensitivity.

**Duplication**
- Reverse-order junction *alone* ≠ DUP (it's the DEL-of-complement reading) → require a co-located gain.
- The DEL `AFC` formula gives a *negative* AF on a gain → make `AFC` sign-aware (§4).
- Dispersed dups decouple donor-gain from acceptor-junction → pair the two junctions; do **not** infer a
  big DEL from the inter-junction gap; today expect 0 PASS.

**Complex / circular**
- Partial-dup and its complementary deletion share an identical junction → only coverage direction
  separates them; keep clusters separate; don't let an embedded del PASS under a net gain.
- Fold-back `INVDUP` = strand-flip **+** gain → not a balanced (CN-neutral) inversion.
- `circSam.pl` wraps coords/CIGARs but **not** `SA` strand → origin events stack the wrap blind spot;
  `WRAP`-withhold them (all event types) until origin-resolution lands; keep them amber/non-PASS/out of
  plots.

**NUMT / repeats / quantification**
- NUMTs (dispersed nuclear mtDNA copies) are the single biggest false-**DUP** source → `NUMT` masks +
  competitive alignment + MAPQ floor + min-VAF; matters *more* for DUP than DEL.
- Homopolymers/STRs (D310, control-region poly-C) make spurious reverse-order splits and micro coverage
  bumps → `HP`/`HS` masks; the real `del4977` 13 bp `ACCTCCCTCACCA` repeat must **annotate** (`REPEAT`),
  not auto-reject.
- A real `>BIGDEL` non-origin deletion must PASS via **dosage** (the J path is blocked for big dels); a
  no-drop big "deletion" is the origin artifact → stays `WRAP`.
- Low VAF: dosage unreliable, `AFJ` dominates — but that's where NUMT/repeat FPs concentrate; require
  corroboration + min reads + min size (all `HP_*` thresholds); INV relies on `AFJ` alone with stricter
  gates.

---

## 7. Schema additions (append-only, when built)

- `SVTYPE=DUP` (positive `SVLEN`), `SVTYPE=INV`, or paired `BND` for inversions.
- `SVCLAIM`: `DJ` for DUP (junction + dosage gain), `J` for INV (junction/adjacency only).
- `INV3`/`INV5` orientation flags; `CIPOS`/`CIEND` homology windows.
- `AFC` becomes sign-aware (negative-of-drop convention documented in the header); `AFC=.` for INV.
- New flags: `INVDUP`, `DISPERSED`, `DUPDEL`. New `HP_SV_*`: `HP_SV_ORIH`, `HP_SV_ORIL`,
  `HP_SV_RESOLVE_ORIGIN`, INV/DUP enable + threshold knobs. All default **off / DEL-only**.

---

## 8. Phased roadmap

- **v3a — tandem DUP + origin-resolution.** Sign-aware `AFC`; `reverse-order ∧ gain → SVTYPE=DUP`;
  `HP_SV_RESOLVE_ORIGIN` (OriH/OriL) to turn today's `WRAP` set into DEL/DUP. Opt-in.
- **v3b — INV detect-and-flag.** Opposite-strand junction branch; `AFJ`-only; `SVTYPE=INV` (`INV3/INV5`);
  default non-PASS; `CIPOS/CIEND`.
- **v3c — complex.** dup-del cluster separation; `INVDUP`; dispersed-dup junction pairing; multi-event
  `AFJ`-or-`AFC` quantification.
- Throughout: NUMT/HP/DLOOP masks extended to every event type; long-read engine (`HP_SV=sniffles`) as
  the eventual ground-truth path for complex/origin events.

---

## 9. Test mapping (what each fixture pins)

The committed mock set (21 BAMs in `test/sv/bams/` — the authoritative catalog is
[`../test/sv/TEST_BAMS.md`](../test/sv/TEST_BAMS.md) §2). Each `truth.tsv` row carries a `kind` ∈
{`del`, `delwrap`, `dup`, `dupdel`, `inv`, `invdup`, `none`} and an `expect` (today's caller behavior) ∈
{`pass`, `detected`, `no_pass`, `no_record`, `wrap`, `known_fp`}. Forward-looking events assert their
**current** state; the fixture names below are the real ones:

| Fixture(s) | Class | Today's `expect` | Lands with |
|---|---|---|---|
| `sv_del4977_h30/h05`, `sv_del6000_h50`, `sv_multidel`, `sv_homoplasmy`, `sv_dloop`, `sv_lowcov`, `sv_del_500` | DEL | `pass` / `detected` | (done) |
| `sv_del_45` (45 bp) | DEL | `no_record` (< minsize; CIGAR-`D`) | (done) |
| `sv_del_13kb` (13 kb majority-arc, with drop) | DEL | `pass` (via dosage) | (done) |
| `sv_origin` (16400→200, clips OriH) / `sv_del_origin_spares` (16400→100, spares) | origin DEL | both `wrap` (indistinguishable now) | v3a origin-resolution → DUP vs DEL |
| `sv_dup` (1 kb) / `sv_dup_large` (5 kb) | tandem DUP | `no_pass` (gain-blocked DEL-shaped record) | v3a |
| `sv_dupdel` | complex dup-del | `known_fp` (embedded del spuriously PASSes) | v3a/v3c |
| `sv_invdup` | INVDUP (fold-back) | `no_record` | v3c |
| `sv_inv_small/large/origin/lowhet` | INV | `no_record` (strand-filtered) | v3b |
| `sv_wt` | control | `no_pass` | (done) |

A simulator **signature pre-assertion** (`run_test.py`'s `bam_signature_ok`) guards faithfulness: each
forward-looking mock must actually carry its intended signal (opposite-strand `SA` for INV, a coverage
gain for DUP, off-origin `SA` for wrap) **before** the test asserts the caller's behavior — so a broken
simulator fails loudly instead of a test passing for the wrong reason.

---

## 10. References

- **MitoSAlt** — Basu et al., *PLOS Genetics* 2020, [PMC7769605](https://pmc.ncbi.nlm.nih.gov/articles/PMC7769605/)
  (DEL/DUP, origin-preservation disambiguation).
- **eKLIPse** (del + soft-clip), **Sniffles2** (long-read `INVDUP`/`INVDEL` nested types), **Manta/DELLY/GRIDSS**
  (`INV3/INV5`, `BND`, `CT` orientation conventions).
- mtDNA inversions: PMID 10775530 (ND1 7-nt inversion), PMID 24040073 (IR/replication-linked accumulation).
- In-repo: `SV_METHODS.md` (as-built DEL), `SV_CALLING.md` §10 (roadmap), `callsv.py` (`extract_junctions`
  strand filter, `cvg_gain`, `ORIH`/`ORIL`, `wrapf`), `test/sv/TEST_BAMS.md` (fixtures).
