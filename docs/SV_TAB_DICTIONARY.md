# `sv.tab` data dictionary

Column reference for the MitoHPC structural-variant tables. The **per-sample** `$O.sv.tab` (written by
`scripts/callsv.py`) and the **cohort** `$ODIR/sv.tab` (concatenated by `scripts/getSVSummary.sh`) share
the same columns, **except** the cohort table gains one trailing column, `mitobreak` (see below).

A machine-readable copy (`column<TAB>type<TAB>description`) is emitted next to the cohort table as
`$ODIR/sv.tab.dict.tsv` and is the source for the table below (`scripts/sv.tab.dict.tsv`).

One row per **(sample, SV event)**. Tidy/long — load directly into R/pandas. Coordinates are 1-based on
the circular rCRS (`chrM`, 16569 bp).

| # | column | type | meaning |
|--:|--------|------|---------|
| 1 | `sample` | string | Sample name (from `$HP_IN`). |
| 2 | `chrom` | string | Contig (`chrM` / rCRS); circular, 16569 bp. |
| 3 | `svtype` | `DEL`\|`DUP`\|`INV` | SV type. **DEL** is default; **DUP**/**INV** are opt-in (`HP_SV_DUP`/`HP_SV_INV`). |
| 4 | `pos_bp5` | int (1-based) | 5′ breakpoint = last **retained** base before the event (VCF `POS`). |
| 5 | `end_bp3` | int (1-based) | 3′ breakpoint = last **deleted** base for DEL (VCF `END`). *(MitoBreak's 3′ bp = `end_bp3 + 1`.)* |
| 6 | `svlen` | int (bp, >0) | Event length; the VCF `SVLEN` is signed (negative for DEL). |
| 7 | `svclaim` | `DJ`\|`J` | PASS evidence path: **DJ** = coverage-drop and junction agree; **J** = junction-strong alone. |
| 8 | `jr` | int | Junction-supporting split reads. |
| 9 | `sr` | int | Wild-type reads spanning the breakpoints (denominator of AFJ). |
| 10 | `af_junction` | float [0,1] | Junction VAF (**AFJ**) = `jr/(jr+sr)`. |
| 11 | `af_coverage` | float [0,1] or `.` | Coverage-dosage VAF (**AFC**) = `1 − inside/flank` depth; **primary heteroplasmy**. `.` for INV. |
| 12 | `afdiff` | float | `|AFC − AFJ|` — agreement QC between the two heteroplasmy estimates. |
| 13 | `cvgr` | float | Masked coverage ratio inside/flank (basis of AFC). |
| 14 | `flank_dp` | int | Trimmed-median flanking read depth. |
| 15 | `homlen` | int | Breakpoint microhomology / direct-repeat length (bp). |
| 16 | `homseq` | string | Microhomology/repeat sequence; empty when `homlen=0`. |
| 17 | `delclass` | `I`\|`II`\|`III`\|`.` | Homology class: **I** = perfect repeat ≥5 bp, **II** = 1–4 bp microhomology, **III** = blunt. `.` if not DEL. |
| 18 | `common` | `0`\|`1` | 1 if the call matches the common deletion del4977 (m.8470_13447) within tolerance. |
| 19 | `ngene` | int | Number of genes/features overlapped. |
| 20 | `gene_list` | string | Comma-joined overlapped gene/feature names, or `.`. |
| 21 | `hgvs` | string | HGVS-style description, e.g. `NC_012920.1:m.8483_13446del`. |
| 22 | `filter` | string | VCF FILTER: `PASS`, or `;`-joined non-PASS reasons (`lowJR`, `no_cvg_drop`, `low_dosage`, `lowAFJ`, `WRAP`, …). |
| 23 | `flags` | string | FP-control flags (`REPEAT`, `NUMT`, `HP`, `DLOOP`, `WRAP`, `INVDUP`, …), comma-joined or `.`; mostly non-rejecting. |
| 24 | `srcons` | float [0,1] | Split-read size consistency (microhomology-invariant). |
| 25 | `srsb` | float [0,1] | Split-read strand balance. |
| 26 | `jsup` | `HIGH`\|`MOD`\|`LOW` | Split-read evidence tier (depth-independent). |
| 27 | `svconf` | int [0,100] or `.` | Per-call confidence (higher = more likely true); `.` for WRAP/origin calls. |
| 28 | `svimpact` | int | Biological-impact score (gene/origin/constraint weighted). |
| 29 | `svimpact_band` | string | Impact band label (e.g. `SEVERE`/`MODERATE`/…). |
| 30 | `mitobreak` | string | **Cohort-only.** Closest matching MitoBreak previously-reported breakpoint id (`DEL_bp5_bp3` / `DUP_bp5_bp3`) within `HP_SV_MITOBREAK_TOL` bp, or `.`. **Not present in per-sample `$O.sv.tab`.** |

## The `mitobreak` annotation (cohort-only)

Added **once at consolidation** by `getSVSummary.sh` → `scripts/svMitoBreak.py`, so individual samples never
need re-running. It also appears in the cohort VCFs as an INFO field:

```
##INFO=<ID=MITOBREAK,Number=.,Type=String,Description="Previously-reported breakpoint id(s) from the
        MitoBreak database … matched within HP_SV_MITOBREAK_TOL bp on both breakpoints …">
```

and the interactive `sv.report.html` gains a **MitoBreak** column in the recurrence table and a
**“MitoBreak-reported only”** filter.

- **Database:** `RefSeq/mitobreak.tsv.gz` — normalized from the MitoBreak database (Damas et al. 2014,
  *NAR* — [PMC3965124](https://pmc.ncbi.nlm.nih.gov/articles/PMC3965124/);
  <http://mitobreak.portugene.com>) by `scripts/mkMitoBreak.py` from the source CSVs in `resources/`
  (1369 deletions + 44 duplications). Columns: `svtype, bp5, bp3, length, location, origin_impact,
  disease, refs, mitobreak_id`.
- **Coordinate convention:** MitoBreak's 3′ breakpoint is the first **retained** base after the event
  (exclusive); our `end_bp3` is the last **deleted** base (inclusive). They differ by 1 for deletions
  (`MitoBreak_bp3 = end_bp3 + 1`); 5′ breakpoints share convention. `svMitoBreak.py` reconciles this
  before matching.
- **Tolerance:** `HP_SV_MITOBREAK_TOL` (default **20** bp, per breakpoint) — covers the ~13 bp
  direct-repeat breakpoint spread (e.g. the common deletion) plus the +1 offset while keeping spurious
  matches well under 1%. Forwarded by `mitohpc-batch-container.sh`.

See `docs/SV_METHODS.md` for the calling method and `scripts/sv.vcf` for the VCF header of the
per-sample/merged SV VCFs.
