# How SV heteroplasmy (VAF) is computed

A self-contained, reproducible description of the two variant-allele-fraction (heteroplasmy) estimates
MitoHPC's SV caller reports for each deletion: a **coverage-dosage** VAF (depth-based) and a **junction**
VAF (split-read-based). All logic is in `scripts/callsv.py`; defaults below are the `--`flag defaults.
See also `docs/SV_METHODS.md` (full method) and `docs/SV_TAB_DICTIONARY.md` (column schema).

## Two estimates

| | column | based on | role |
|---|--------|----------|------|
| **Coverage-dosage VAF** | `af_coverage` (**AFC**) | read depth | **primary** heteroplasmy |
| **Junction VAF** | `af_junction` (**AFJ**) | split reads | corroborating evidence |

`afdiff = |AFJ − AFC|` is reported as an agreement/QC signal (small = the two estimates concur).
Coordinates wrap **mod `HP_MTLEN` (16569)** throughout, since chrM is circular.

---

## 1. Junction VAF (AFJ) — split-read based

```
AFJ = JR / (JR + SR)
```

- **JR** — distinct read templates whose **split-read junction** supports this `(bp5, bp3)`. Junctions are
  parsed from each read's `SA:Z:` supplementary-alignment tag and clustered to a breakpoint pair; only
  reads with **MAPQ ≥ `--minmapq` (20)** are used. JR is depth-independent in spirit: a few clean junction
  reads already give a meaningful fraction.
- **SR** — distinct **wild-type templates that span the breakpoint** (evidence *against* the deletion). A
  template spans a boundary `B` if one of its aligned blocks covers **both `B` and `B+1`** continuously; a
  read soft-clipped or SA-split *at* `B` ends its block at `B` and is correctly excluded. Counted over
  primary **and** supplementary blocks, deduplicated by read name, MAPQ ≥ 20.
  - `SR = min(spanning at bp5, spanning at bp3 − 1)` — the smaller of the two breakpoints (conservative).

*Code:* `extract_junctions` (JR + clustering), `count_spanning_boundaries` + `spanning_count` (SR),
AFJ at `callsv.py:` the per-junction loop (`afj = jr / (jr + sr)`).

---

## 2. Coverage-dosage VAF (AFC) — depth based

```
cvgr = d_inside / d_flank
AFC  = clamp( 1 − cvgr , 0, 1 )
```

**Per-base depth** `dep[p]` = number of reads with an A/C/G/T base aligned at position `p` — excludes
unmapped/secondary/QC-fail/duplicate reads, **keeps** supplementary (so origin-crossing arcs count), and
applies **no** MAPQ filter (`pysam count_coverage(quality_threshold=0)` semantics; computed via an
equivalent ~8× faster `get_blocks` difference-array — byte-identical).

**Windows** (transition pad `PADt = --trans (150)` excludes the breakpoint smear; flank width
`F = --flank (200)`):

- **inside** = `[ bp5 + PADt + 1 … bp3 − PADt − 1 ]`
- **flank** = `[ bp5 − PADt − F + 1 … bp5 − PADt ]` ∪ `[ bp3 + PADt … bp3 + PADt + F − 1 ]`

**Masking** — before aggregating, drop "fragile" positions from *both* windows: D-loop (`DLOOP.bed.gz`),
homopolymers (`HP.bed.gz`), NUMT-like sites (`NUMT.vcf.gz`), and anything within `--originpad (20)` bp of
the origin. (These regions produce coverage "bowls"/spikes with no real junction.)

**Aggregate** with a **trimmed median** — drop 15% of each tail, then take the median (robust to
NUMT/homopolymer depth spikes):

- `d_inside = trimmed_median(inside depths)`,  `d_flank = trimmed_median(flank depths)`.

**Estimability guard** — AFC is computed only if **each window has ≥ 50 unmasked bases** (`MINBASE`) and
`d_flank > 0`. Otherwise the dosage is not estimable and the call **falls back to AFJ** (`dose = False`).
For inversions (copy-number-neutral) AFC is undefined and reported as `.`.

*Code:* `per_base_depth` (depth), `masked()` (mask predicate), `masked_depths` (windowed, wrap-aware),
`trimmed_median`, and the per-junction AFC block in `callsv.py`.

---

## Default parameters

| quantity | flag / constant | default |
|----------|-----------------|---------|
| min MAPQ for JR / SR | `--minmapq` | 20 |
| transition pad `PADt` | `--trans` | 150 bp |
| flank width `F` | `--flank` | 200 bp |
| trimmed-median tail | constant | 15% per tail |
| min usable bases per window | `MINBASE` constant | 50 |
| origin mask radius | `--originpad` | 20 bp |
| min junction reads (PASS) | `--minjr` | 3 |

(Depth itself is unfiltered by MAPQ; JR and SR require MAPQ ≥ 20 — a deliberate asymmetry to reproduce.)

---

## Worked example — the common deletion (committed `del4977` call)

From the example `sv.tab` row: `JR = 154`, `SR = 401`, `cvgr = 0.732`.

- **AFJ** = 154 / (154 + 401) = 154 / 555 = **0.277**
- **AFC** = 1 − 0.732 = **0.268**  → primary heteroplasmy ≈ **27 %**
- **afdiff** = |0.277 − 0.268| = 0.009 → tight agreement (a confidence signal)

The PASS logic then uses both: a call PASSes either via the **DJ** path (a coverage drop, i.e. AFC, that a
proportional junction AFJ corroborates) or the **J** path (a clean, well-supported junction alone, since
mtDNA read depth can be unreliable). See `docs/SV_METHODS.md`.
