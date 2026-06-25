# Implementation note — resolving the circular DEL-vs-DUP ambiguity (origin-preservation)

**Status:** *design captured, NOT implemented. Deferred (roadmap `SV_CALLING.md` §10 v3a).* Default
caller behavior is unchanged; this note records the plan so it can be picked up cold later.

**One-line goal:** turn the origin-crossing / majority-arc junctions the caller currently `WRAP`-
*withholds* into a principled `DEL`-or-`DUP` *classification*, behind an opt-in flag — and gain the
ability to emit genuine duplication calls (today the caller is deletion-only and rejects all dups).

---

## 1. The problem (why one junction isn't enough)

On a **circular** genome a single split-read junction joining `a→b` is consistent with **both**:
- a **deletion** of one arc (retain the complementary arc), and
- a **tandem duplication** of the complementary arc.

They produce the *identical* junction and the *identical relative-coverage* pattern, so they are
**indistinguishable from short reads alone**. This is the canonical view — MitoSAlt (Basu et al.,
*PLOS Genetics* 2020, [PMC7769605](https://pmc.ncbi.nlm.nih.gov/articles/PMC7769605/)) states it
directly: *"every split read can represent either a deletion or, alternatively, a duplication of the
mtDNA arc complementary to the deletion, and these two possibilities are indistinguishable when using
short read sequencing."*

Quantitatively, a **deletion of arc D at heteroplasmy `h`** is the same observation as a
**duplication of the complement R at `h′`**, with `1/(1−h) = 1+h′` (e.g. a 40% deletion ≡ a 67%
duplication — same junction, same coverage ratio, different reported het).

What does **not** resolve it:
- **Junction orientation** alone — on a circle, both DEL and DUP can yield reverse-order ("backward
  jump") split reads; the origin-crossing case flips orientation, so a true origin deletion's junction
  looks dup-like. (Orientation *does* separate a duplication of the **same** region, which is why our
  coverage-gain check correctly rejects `sv_dup`.)
- **Coverage / copy number** alone — it's baseline-relative. A scalar genome-average mtDNA CN (mean
  chrM depth ÷ nuclear depth, what `getCN.pl` reports) can't, and even absolute per-arc coverage fits
  *both* readings (you still must declare which arc is the baseline molecule count). See the reasoning
  recorded in `SV_METHODS.md` §8.

---

## 2. Current behavior (v2, as implemented)

- The caller is **deletion-only** (`SVTYPE=DEL`); it computes `svlen` **linearly** (`bp3−bp5−1`) and so
  reports the **majority arc**.
- A **majority-arc deletion (`svlen > MTLEN/2`) with no coverage drop** is detected as the inverted
  complementary arc and **`WRAP`-flagged → non-PASS → `SVCONF='.'`** (withheld, not resolved). It never
  emits a `DUP`. (See `callsv.py` `wrapf`; `SV_METHODS.md` §8.)
- Net: origin-crossing / complementary-arc events are **detected-but-withheld**.

---

## 3. Proposed resolution — origin-preservation (MitoSAlt) + coverage corroboration

The field-standard disambiguator is **biological parsimony via replication competence**, *not*
coverage. MitoSAlt's rule:

> default to **deletion**; reclassify as a **duplication** of the complementary arc iff the deletion's
> surviving molecule is **replication-incompetent** because it lost a replication origin (**OriH** or
> **OriL**). Favor the interpretation where the origins are intact.

Then **corroborate with coverage**: the predicted-*deleted* arc should show depletion; the
predicted-*duplicated* arc should show a gain. (Belt-and-suspenders; neither signal is airtight alone.)

**Worked examples** (origins per `callsv.py`: `OriH = 110–441`, `OriL = 5721–5798`):

| junction | deletion-reading molecule | OriH | OriL | resolved call |
|---|---|---|---|---|
| `16400→100` | retains `100..16400` | intact (110–441 ⊂ body) | intact | **DELETION** of the ~268 bp origin arc |
| `16400→200` (`sv_origin`) | retains `200..16400` | **clipped** (110–199 deleted) | intact | **DUPLICATION** of the ~16 kb body |

The only difference is the breakpoint moving from 100 to 200 — crossing OriH's start (110). At 100 the
deletion spares OriH → stays a deletion; at 200 it clips OriH → flips to a duplication. This is the
whole mechanism in one example, and these two are the natural first test cases.

---

## 4. What already exists (the easy 10%)

- `ORIH=(110,441)`, `ORIL=(5721,5798)` constants in `callsv.py` (used only for `SVIMPACT` today).
- The `wrapf` majority-arc detector already *identifies the exact set of junctions to resolve*.
- Masked, wrap-aware coverage windows for the corroboration step.
- `chrMR` (rotated reference) is available for origin-region work (roadmap §3/§9).

The origin `if`-check is small. **The work is everything downstream of the classification.**

---

## 5. Scope of work (the hard 90%)

1. **DUP heteroplasmy math** — a duplication is a coverage *gain*, not a drop. Invert the `AFC`
   estimate; report the dup heteroplasmy (note the `del-h ≡ dup-h′` relation, so the number differs
   from what a DEL reading would print for the same junction).
2. **`SVTYPE=DUP` VCF representation** — which arc is duplicated, the tandem junction, **positive**
   `SVLEN`, `SVCLAIM`, `CIPOS/CIEND`; `COMMON`/`DELCLASS` are DEL-specific (N/A or redefined).
3. **Wire DUP through the cohort + viz** — `getSVSummary.sh` (`sv.tab`/merged/sites), `svReport.py`
   gallery (status/labels), `svplot.sh` (`samplot -t DUP`). All are DEL-shaped today.
4. **Gating** — `HP_SV_RESOLVE_ORIGIN` (default **OFF**), additive, frozen behavior intact. Emit
   resolved calls with an explicit **`ORIGIN_RESOLVED`** + **low-confidence** flag.
5. **Keep `WRAP`/withhold for the genuinely-undecidable residue** (see §7).

---

## 6. Validation (needed before trusting any of it)

Today the suite simulates only deletions plus **one** dup-not-called negative (`sv_dup`). To trust
classification, add:
- a **clean tandem-dup mock** (should classify as `DUP` with correct het — `sv_dup` should now be
  *called as a DUP*, not merely rejected);
- the **`16400→100` mock** (should resolve to `DEL`);
- `sv_origin` (`16400→200`) **should now resolve to `DUP`**;
- ideally a **real spiked duplication** (mirroring the del4977 spike-in).
Assertions: correct `DEL`/`DUP` label + correct (inverted) heteroplasmy on each.

---

## 7. Known limitations & open questions (read before implementing)

- **It's a heuristic, not ground truth.** Replication-competence parsimony is usually right but is
  **undecidable** when *both* readings spare origins (a junction far from OriH/OriL — deleting either
  arc is viable) or *neither* does. MitoSAlt's reviewer flagged exactly this; keep a default/withhold
  path for the residue rather than forcing a call.
- **Partial origin disruption is murky** — is "breakpoint inside OriH" = origin lost? Graded vs binary.
- **Replication-competence is simplified** — OriL-independent / alternative replication exists; CSB/TAS
  control-region elements are ignored by an OriH/OriL-only test.
- **Coverage corroboration is baseline-relative**, so it confirms rather than independently proves.
- **Payoff is concentrated in rare, often biologically-degenerate events** (near-whole-genome
  complementary arcs). The common, clinically important deletions (del4977 &c.) are already called
  correctly. `sv_origin`'s own "resolution" is a ~64% duplication of 16 kb — itself an odd event; that
  mock may be biologically implausible as *either* class (a molecule missing part of OriH wouldn't
  propagate). It remains a good **technical** test of safe origin-junction handling.
- **True structural resolution needs long/linked reads** (direct molecule length/content). Short reads
  give one junction + relative depth, which is fundamentally symmetric — this rule is a principled
  *best guess*, not a measurement.

---

## 8. Suggested phased v1 (when we return to it)

`HP_SV_RESOLVE_ORIGIN=1` (opt-in, default off → byte-identical to today). For **only** the junctions
the caller currently `WRAP`-withholds:
1. apply origin-preservation to pick `DEL` vs `DUP`;
2. corroborate with the arc coverage (deleted⇒depleted, duplicated⇒gain);
3. emit the resolved call with `ORIGIN_RESOLVED` + a low-confidence flag;
4. **withhold (keep `WRAP`) the undecidable residue** (both/neither arc origin-viable).
Plus the dup validation mocks from §6. Everything else (the frozen DEL path) is untouched.

---

## 9. References

- **MitoSAlt** — Basu et al., *"Accurate mapping of mitochondrial DNA deletions and duplications using
  deep sequencing,"* PLOS Genetics 2020, [PMC7769605](https://pmc.ncbi.nlm.nih.gov/articles/PMC7769605/).
  Origin-preservation (OriH/OriL) disambiguation; the "indistinguishable from short reads" statement.
- This repo: `docs/SV_METHODS.md` §8 (circular handling, the del-vs-dup reasoning), `docs/SV_CALLING.md`
  §10 v3a + §11 (roadmap/risks), `callsv.py` (`ORIH`/`ORIL`, `wrapf`/majority-arc detection),
  `test/sv/TEST_BAMS.md` (`sv_origin`, `sv_dup`).
