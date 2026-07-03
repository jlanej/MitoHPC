# MitoHPC SV test BAMs — what each one is and why we test it

This catalogs every BAM used to test the structural-variant (SV) caller: the **21 simulated mock
BAMs** in [`bams/`](bams/) — spanning **deletions, duplications, inversions, and complex events** —
and the **real 1000 Genomes alignments** in [`real/`](real/). For each, it records *what the BAM
contains*, *why that scenario exists* (the specific caller behavior it pins down), and *what the
harness asserts* about it.

**The caller is deletion-only today.** Many of the duplication / inversion / complex fixtures are
**forward-looking** — they exist so that when those call paths land (design in
[`../docs/SV_EVENT_TYPES.md`](../docs/SV_EVENT_TYPES.md)) the behavior is already pinned. Until then
they assert their *current* state (mostly `0 records` / `0 PASS`), recorded per-sample in the
`expect` column of `truth.tsv` and enforced by `run_test.py`.

Sources of truth this document is derived from — keep it in sync with them if they change:
- [`make_testdata.py`](make_testdata.py) — the read simulator (`SAMPLES`, the event constructors).
- [`truth.tsv`](truth.tsv) — per (sample, event) `kind`, breakpoints, het, depth, and the `expect`
  (today's expected caller behavior: `pass` / `detected` / `no_pass` / `no_record` / `wrap` / `known_fp`).
- [`run_test.py`](run_test.py) — the expect-driven assertions (`check_sample`, `check_degenerate`,
  `check_cohort`, `check_plots`, `check_recall`).
- [`../docs/SV_EVENT_TYPES.md`](../docs/SV_EVENT_TYPES.md) — the DUP/INV/complex design these fixtures pin.
- [`real/README.md`](real/README.md) — the real-data litmus + LoD/accuracy harness.
- `../docs/SV_METHODS.md` — the as-built (deletion) caller method, fields, and thresholds.

---

## 1. How the mock BAMs are built (shared model)

[`make_testdata.py`](make_testdata.py) is a pure-stdlib, deterministic (fixed seed) paired-end read
simulator. For each sample it draws reads from a **mixture of circular genomes** at a known
heteroplasmy `h`:

- a **wild-type** chrM genome contributing depth fraction `1 − Σh`, and
- one **event** genome per simulated variant, each contributing fraction `h`. Constructors:
  `make_deletion`, `make_dup` (tandem), `make_delwrap` (origin-crossing), `make_inversion`
  (reverse-**complement** in place → opposite-strand junctions), `make_inv_dup` (fold-back),
  `make_dupdel` (partial dup-del). The event genome is **doubled** (`eg+eg`)
  before sampling so fragments wrap its (shifted) origin, exercising the circular path.

Reads are then aligned back through the pipeline's **circular path** (`gen_bams.sh`: minimap2 →
`circSam.pl` → sort) to produce a BAM that is byte-for-byte equivalent in *form* to the per-sample
circular-aware `$O.bam` the real pipeline hands the caller. Aligning a WT+DEL mixture to the WT
circular reference reproduces **both** signals a real deletion makes:

1. **split / soft-clipped reads with `SA` tags** spanning the deletion junction (the DEL molecules), and
2. a **coverage drop** across the deleted span (only WT molecules cover it).

The simulator picks read counts so that, outside the deletion, the DEL genome contributes fraction
`h` of depth. Hence the **coverage-dosage AF (`AFC`) ≈ h** and the **junction-read fraction (`AFJ`)
≈ h** — the two heteroplasmy estimates the caller reports. (Production aligns with `bwa mem`;
minimap2 is used only to mint mock data — the split-read + coverage signal the caller consumes is
aligner-agnostic, confirmed by [`real/aligntest.tsv`](real/aligntest.tsv).)

**Breakpoint convention** (same in `truth.tsv`, the caller, and the assertions): `bp5` = last
*retained* base left of the deletion, `bp3` = first *retained* base right of it; the deleted span is
`bp5+1 … bp3−1` and `SVLEN = bp3 − bp5 − 1`. VCF `POS=bp5`, `END=bp3−1`.

### 1.1 How the reads are simulated and the events created — semi-verbose

Everything is synthesized from the wild-type reference `RefSeq/chrM.fa` (rCRS, 16569 bp) — there are
**no external reads**. Per sample the simulator does two things: **(a)** build one or more *event
genomes* by editing the WT sequence in silico, then **(b)** sample paired-end reads from a WT + event
**mixture** at the target heteroplasmy and write FASTQ.

**(a) Event genomes — string surgery on the 1-based WT sequence** (`event_genome` dispatch):

| kind | constructor | in-silico operation | signal it creates |
|---|---|---|---|
| `del` | `make_deletion(bp5,bp3)` | `WT[:bp5] + WT[bp3-1:]` — fuse base `bp5` to base `bp3`, dropping `bp5+1..bp3-1` | same-strand junction **+ coverage drop** |
| `dup` | `make_dup(a,b)` | `WT[:b] + WT[a-1:b] + WT[b:]` — a 2nd tandem copy of `[a,b]` | everted junction **+ coverage gain** |
| `delwrap` | `make_delwrap(bp5,bp3)`, `bp5>bp3` | `WT[bp3-1:bp5]` — keep only arc `[bp3,bp5]`; delete the origin-crossing complement | junction whose **deleted arc crosses the origin** |
| `inv` | `make_inversion(a,b)` | `WT[:a-1] + revcomp(WT[a-1:b]) + WT[b:]` — reverse-**complement** `[a,b]` in place | two **opposite-strand** junctions, **CN-neutral** (no gain) |
| `invdup` | `make_inv_dup(a,b)` | `WT[:b] + revcomp(WT[a-1:b]) + WT[b:]` — inverted extra copy | opposite-strand junction **+ gain** (INVDUP) |
| `dupdel` | `make_dupdel(a,b,da,db)` | tandem copy of `[a,b]` whose 2nd copy omits internal `[da,db]` | dup junction **and** an embedded del junction |

**(b) Reads.** Paired-end, FR-oriented, read length `rlen=150` bp (flat Q37), fragment length uniform
in `[300, 450]`, per-base substitution error `0.1 %` (`-err 0.001`; a random base). **Deterministic** —
the RNG is seeded per sample (`42 + Σ ord(name)`), so regeneration is bit-identical. Each genome (WT
**and** each event) is **doubled** (`g+g`) and the fragment start is drawn over the first half, so a
fragment that runs off the end **wraps the origin** — that is where origin-crossing reads come from.

**Heteroplasmy → read counts.** For a sample at outside-depth `D` and event fraction `h`:
`n_WT = round(D·(1−Σh)·16569 / rlen)` and, per event, `n_event = round(D·h·len(event) / rlen)`. Because
the shorter/longer event genome redistributes its reads, **outside** the deletion the event contributes
fraction `h` of depth and WT `(1−h)`; **inside** only WT contributes → `inside/outside = 1−h` →
**`AFC = 1−ratio = h`**, and the junction-read fraction **`AFJ ≈ h`**. (Each fragment yields two mates,
so the realized per-base coverage is ≈ `2·D`; the heteroplasmy *ratio* is unaffected.)

**Alignment (`gen_bams.sh`).** `minimap2 -ax sr` → circularized reference `chrMC.fa`
(16869 bp = 16569 + the 300 bp origin overlap) → `samtools view -F 0x90C` (keep primary; drop
unmapped/mate-unmapped/secondary/supplementary) → `circSam.pl` (wraps reads aligned into the 300 bp
overlap back to `1..16569` **and** emits the `SA`-tagged split alignments for origin-crossing reads) →
`samtools sort`. The product is identical in *form* to the pipeline's per-sample `$O.bam`. (Production
aligns with `bwa mem`; minimap2 only mints mock data — the split-read + coverage signal is
aligner-agnostic.)

### 1.2 Verbose walkthrough — one sample end-to-end (teach-an-agent version)

Follow `sv_del4977_h30` (`del 8469→13447`, `h=0.30`, `depth=300`) from reference to BAM:

1. **Reference.** `WT = RefSeq/chrM.fa` = rCRS, 16569 bp, 1-based.
2. **Build the deletion genome.** `make_deletion(WT, 8469, 13447)` = `WT[:8469] + WT[13446:]` →
   `16569 − 4977 = 11592 bp`: bases `1..8469` immediately followed by `13447..16569`. The seam fuses
   base 8469 to base 13447 — the **deletion junction**. (Those breakpoints sit inside the 13 bp direct
   repeat `ACCTCCCTCACCA` at 8470–8482 / 13447–13459, so the fusion point is ambiguous *within* the
   repeat; the caller "slides" and lands near 8482/13446 — still within `BP_TOL`.)
3. **Circularize it.** `del2 = del + del` (23184 bp) so a fragment can run off the end and wrap the
   deletion genome's (shifted) origin.
4. **Choose read counts** (`h=0.30, D=300, rlen=150`):
   - WT: `round(300·0.70·16569/150) ≈ 23197` fragments.
   - DEL: `round(300·0.30·11592/150) ≈ 6955` fragments.
   The shorter DEL genome packs those 6955 fragments to the same **per-base** share as `h`: realized
   coverage ≈ 420× (WT) + 180× (DEL) = **600× outside** the deletion (≈ `2·D`; the ×2 is the two mates).
5. **Sample each fragment.** Pick `flen∈[300,450]`, `start∈[0,glen)`, `frag = template[start:start+flen]`;
   `R1 = frag[:150]`, `R2 = revcomp(frag[-150:])`; add 0.1 % base noise; write to
   `sv_del4977_h30_1.fq` / `_2.fq`. Three outcomes for DEL fragments:
   - **entirely on one side** of the seam → aligns normally to WT.
   - **straddling the seam** → contains `…(ends at ref 8469)(starts at ref 13447)…`, contiguous in the
     deletion genome but **not** in WT. Aligned to WT its two halves map to ~8469 and ~13447 → a
     **split read with an `SA` tag** = the junction evidence (`JR`).
   - **nothing covers 8470..13446** (that span doesn't exist in the deletion genome) → inside the
     deletion only WT reads land → the **coverage drop** (420× inside vs 600× outside →
     `AFC = 1 − 420/600 = 0.30`).
6. **Align → BAM** (`gen_bams.sh`): minimap2 → `chrMC` → `-F 0x90C` → `circSam.pl` → sort →
   `bams/sv_del4977_h30.bam`.
7. **What the caller now sees:** an `SA`-junction cluster at ~8482/~13446 (`JR`, with wild-type
   spanning reads as the `AFJ` denominator → `AFJ≈0.30`) and a ~30 % inside coverage drop
   (`AFC≈0.30`) → **PASS**, `COMMON=1`, `HOMLEN=13`, `HOMSEQ=ACCTCCCTCACCA`, `DELCLASS=I`.
8. **Ground truth** to `truth.tsv`: `sv_del4977_h30  del  8469  13447  4977  0.300  300  pass`.

**How the other kinds diverge from this template:**
- **`dup`** inserts a 2nd tandem copy → an **everted** junction (reference order reverses) **plus a
  coverage gain** (`CVGR ≈ 1+h`, ≈1.5 at `h=0.5`). The DEL-only caller sees the junction but the gain
  (`ratio > 1+gainpad`) blocks a deletion PASS.
- **`inv`** reverse-**complements** the segment in place: the breakpoints become **opposite-strand**
  junctions and copy number is unchanged (no gain). `extract_junctions` keeps only same-strand `SA`
  pairs, so this is invisible (`no_record`) until `--call-inv`.
- **`delwrap`** (`bp5>bp3`) keeps only the arc `[bp3,bp5]`, so the **deleted** arc crosses the origin.
  The caller computes `SVLEN` linearly, cannot represent that arc, and instead reports the
  **complementary retained arc** (a full-coverage majority-arc "deletion") with a `WRAP` flag — never a
  PASS.
- **`invdup` / `dupdel`** compose the above (opposite-strand + gain; dup junction + embedded del) for
  the forward-looking complex fixtures.

Inspect any regenerated sample's raw signal directly, e.g.
`samtools view bams/sv_del4977_h30.bam | grep -c SA:Z:` (junction reads) and
`samtools depth -a bams/sv_del4977_h30.bam` (the coverage drop).

### Thresholds the assertions use (`run_test.py`)

| Constant | Value | Meaning |
|---|---|---|
| `BP_TOL` | 30 bp | breakpoint match tolerance (≥ the del4977 13 bp repeat slide) |
| `SVLEN_TOL` | 40 bp | deletion-size tolerance (caller `svlen` vs truth; enforced on PASS/detected del events) |
| `AF_TOL` | 0.15 | heteroplasmy tolerance; **`AFC` or `AFJ`** may match (either suffices) |

Whether a deletion is *required* to PASS is decided **per fixture** by the `truth.tsv` `expect` column
(`require_pass = (expect == "pass")` in `check_sample`), not by a heteroplasmy threshold.

A genuine simulated deletion is also required to be a **clean junction**: split-read tier
`JSUP = HIGH/MOD` (not the `LOW` artifact tier) with size-consistency `SRCONS ≥ 0.7`.

---

## 2. The 21 simulated BAMs (`bams/`)

Quick reference (from `truth.tsv` / `SAMPLES`); `expect` is the **current** caller behavior the suite
asserts. New event constructors live in `make_testdata.py`: `make_inversion` (revcomp in place),
`make_inv_dup` (fold-back), `make_dupdel`. A **signature pre-assertion** confirmed
each forward-looking BAM actually carries its signal (opposite-strand `SA` for INV, a coverage gain for
DUP, off-origin `SA` for wrap) — so a broken simulator fails loudly rather than a test passing for the
wrong reason.

**Deletions (the implemented, PASS-able class)**

| BAM | event | het | depth | `expect` | what it pins |
|---|---|---|---|---|---|
| `sv_del4977_h30` | del 8469–13447 (4977 bp) | 30% | 300× | pass | common deletion; COMMON/REPEAT/HOMLEN=13/DELCLASS=I |
| `sv_del4977_h05` | del 8469–13447 | 5% | 400× | detected | low-het floor (detected, sub-PASS) |
| `sv_del6000_h50` | del 5999–10999 (4999 bp) | 50% | 300× | pass | non-repeat (not COMMON/REPEAT, Class III) |
| `sv_multidel` | del4977 **and** del6000 | 25/15% | 400× | pass | two concurrent deletions, separate records |
| `sv_homoplasmy` | del 8469–13447 | 95% | 300× | pass | near-homoplasmic; `AFJ→~1.0`, no ÷0 |
| `sv_dloop` | del 400–6000 | 40% | 300× | pass | 5′ bp in D-loop → PASS **+ `DLOOP`** flag |
| `sv_lowcov` | del 8469–13447 | 50% | 40× | pass | low depth, still detected |
| `sv_del_500` | del 8000–8501 (500 bp) | 50% | 300× | pass | small **detectable** deletion (PASS via DJ) |
| `sv_del_45` | del 9000–9046 (45 bp) | 50% | 300× | no_record | < minsize=50 **and** a CIGAR-`D` not a split → 0 records |
| `sv_del_13kb` | del 2000–15001 (13 kb) | 60% | 300× | pass | majority-arc (>BIGDEL) **with** a drop → PASS via dosage (≠ origin artifact) |

**Origin-crossing deletions (WRAP-withheld; the resolution pair)**

| BAM | event | het | depth | `expect` | what it pins |
|---|---|---|---|---|---|
| `sv_origin` | delwrap 16400→200 (368 bp, clips OriH) | 40% | 400× | wrap | reported as ~16 kb complement, `WRAP`, 0 PASS → future **DUP** |
| `sv_del_origin_spares` | delwrap 16400→100 (268 bp, spares origins) | 40% | 400× | wrap | identical WRAP today → future **DEL** (the [origin-resolution](../docs/SV_EVENT_TYPES.md#5-the-circular-del-vs-dup-resolution-origin-preservation) regression pair) |

**Duplications & complex (forward-looking — detected-but-not-PASS today)**

| BAM | event | het | depth | `expect` | what it pins |
|---|---|---|---|---|---|
| `sv_dup` | tandem dup 6000–7000 (1 kb) | 50% | 300× | no_pass | coverage *gain* blocks PASS (not a deletion) |
| `sv_dup_large` | tandem dup 4000–9000 (5 kb) | 50% | 300× | no_pass | gain guard scales to large dups |
| `sv_dupdel` | dup 5000–8000 w/ internal del 6000–6500 | 40% | 400× | **known_fp** | **documented gap**: the embedded del **spuriously PASSes** today (compound-event; fixed by the DUP-aware caller) |
| `sv_invdup` | fold-back inverted dup 7000–7400 | 40% | 400× | no_record | opposite-strand arm + gain (Sniffles2 INVDUP) → 0 records |

**Inversions (DEL-only default: invisible — strand-filtered & CN-neutral; callable via `--call-inv`)**

| BAM | event | het | depth | `expect` | what it pins |
|---|---|---|---|---|---|
| `sv_inv_small` | inv 6000–6500 (500 bp) | 50% | 300× | no_record | balanced inversion → 0 records (the blind spot) |
| `sv_inv_large` | inv 5000–9000 (4 kb) | 50% | 300× | no_record | invisibility is size-independent |
| `sv_inv_origin` | inv 16300–16560 (near origin) | 40% | 400× | no_record | strand-skip stacked near the origin |
| `sv_inv_lowhet` | inv 8000–9000 | 5% | 300× | no_record | low-het inversion (future sensitivity floor) |

**Control**

| BAM | event | het | depth | `expect` | what it pins |
|---|---|---|---|---|---|
| `sv_wt` | none (wild-type) | — | 300× | no_pass | specificity (also the degenerate-input substrate, §4) |

### Positive controls — must detect with correct biology

**`sv_del4977_h30` — the canonical "common deletion," the primary positive control.**
The 4977 bp m.8470–13447 deletion (the most clinically important mtDNA SV; Kearns-Sayre / CPEO /
Pearson) at a moderate, comfortably-detectable 30% heteroplasmy. *Why:* it is the headline call the
module exists to make, and it exercises the repeat-mediated machinery — its breakpoints sit inside a
13 bp direct repeat (`ACCTCCCTCACCA`), so the junction "slides" and the caller lands near
`8482/13447`. *Asserted:* `FILTER=PASS`; `COMMON=1`; `HOMLEN=13` / `HOMSEQ=ACCTCCCTCACCA`;
`DELCLASS=I` (direct-repeat class); `AFC`/`AFJ` ≈ 0.30; genes annotated; clean junction
(`JSUP=HIGH/MOD`). This is also the BAM the committed `example/` report and the samplot/keep-BAM
checks are built on.

**`sv_del6000_h50` — a non-repeat deletion, the negative control for "COMMON."**
A ~5 kb deletion at m.6000–11000 at 50%, with **no** direct repeat at its breakpoints. *Why:* proves
the caller detects ordinary deletions and, crucially, does **not** spuriously stamp them
`COMMON`/`REPEAT` — i.e. the common-deletion annotation is specific to the real 4977 site, not
applied to any large deletion. *Asserted:* PASS, breakpoints/`AFC` correct, and `COMMON≠1` /
`DELCLASS≠I`.

**`sv_homoplasmy` — near-homoplasmic (95%) common deletion, a numerical edge case.**
*Why:* at 95% the wild-type fraction is tiny, so the junction VAF approaches 1.0 and the
spanning-read denominator approaches 0 — the classic place a VAF formula divides by zero or
overflows. *Asserted:* PASS with `AFJ → ~1.0` and no crash/NaN. Guards the heteroplasmy arithmetic
at the top of its range.

### Sensitivity / robustness — detection under stress

**`sv_del4977_h05` — the low-heteroplasmy floor (5%).**
The same common deletion at just 5%. *Why:* the hardest *true positive* — it sets the lower
sensitivity bound. At 5% the coverage drop is marginal, so the call is expected to be **detected with
correct breakpoints/annotation but to fall short of PASS** (typically `no_cvg_drop`), consistent with
the empirical PASS limit-of-detection of ≈ 8% established in [`real/`](real/README.md). *Asserted:*
the deletion is matched (junction present) with `AFC`/`AFJ` ≈ 0.05 and a clean junction; PASS is
**not** required (its `expect` is `detected`, not `pass`). Tests that low-level events surface for
curation rather than being silently dropped.

**`sv_lowcov` — low sequencing depth (40×).**
The common deletion at 50% het but only 40× depth (vs the ~300–400× of the others, and the
thousands-fold real mtDNA depth). *Why:* cohorts have wildly variable chrM depth; this confirms the
caller still recovers a clear high-heteroplasmy deletion when reads are scarce, so depth variability
across a cohort does not drop real calls. *Asserted:* detected and PASS.

**`sv_multidel` — two independent concurrent deletions in one sample.**
del4977 at 25% **and** del6000 at 15%, simulated independently. *Why:* a single mtDNA can carry
multiple deletions; this checks the caller emits **separate** records per junction (clustering does
not merge distinct events) and that the overlapping events don't corrupt each other's heteroplasmy —
where the coverage dosage is confounded by the other deletion's overlap, the **junction** estimate
must stay specific. *Asserted:* both truth deletions matched as distinct records with correct
breakpoints; for each, `AFC` *or* `AFJ` is within tolerance (the harness accepts either precisely
because dosage is confounded here).

### Negative / withheld — must NOT produce a PASS deletion

(Note the three differ in *why*: `sv_wt` and `sv_dup` are **true negatives** — no real deletion is
present — so 0 PASS is *specificity*. `sv_origin` **does contain a real deletion**, but an
origin-crossing one the caller cannot yet resolve, so its 0 PASS is a **deliberate withholding** of a
real event — a documented known-limitation false-negative, not specificity.)

**`sv_wt` — pure wild-type (specificity).**
No event at all. *Why:* the bedrock false-positive guard — a healthy genome must yield no deletion.
*Asserted:* **0 PASS**. (This BAM is also the substrate the degenerate-input subtests mutate; see §4.)

**`sv_dup` — a tandem duplication (must not masquerade as a deletion).**
A tandem duplication of m.6000–7000 at 50%: a coverage *gain* over that span plus a junction where
reference order reverses. *Why:* duplications produce split reads too, and a naïve caller could
report the reverse-order junction as a deletion. The caller must reject it — a coverage *increase*
(`CVGR > 1`) blocks the deletion (`no_cvg_drop`), and the gain tolerance blocks any junction-only
PASS. *Asserted:* **0 PASS**, while still emitting non-PASS record(s) (it is detected as *something*,
just not called a deletion). Pins down deletion-vs-duplication discrimination.

**`sv_origin` — an origin-crossing deletion (circularization safety).**
A deletion whose *deleted arc crosses the artificial linear origin* (retains m.200–16400, deletes the
368 bp arc spanning position 1/16569). *Why:* the circular-genome corner case (§2 of `CLAUDE.md` /
`docs/SV_METHODS.md` §8). The caller computes junction `SVLEN` **linearly**, so it cannot represent the
true 368 bp origin-crossing arc; instead it reports the **complementary (retained) arc** — a 16,044 bp
"deletion" `m.301_16344del`. The giveaway that this is the *retained* arc, not a real deletion, is that
it has **full coverage** (`CVGR≈1.0`, no dosage drop). The requirement is **never a wrong (PASS) call
and never an invalid VCF**. Two guards enforce it: (1) the junction-only PASS path is blocked for
`svlen ≥ BIGDEL`, and the absence of a coverage drop fails the dosage path → `no_cvg_drop`; and (2) a
*majority-arc deletion (`svlen > MTLEN/2`) with no coverage drop* is detected as the inverted
origin-crossing arc and **`WRAP`-flagged**, which also blanks its confidence (`SVCONF='.'`). So the
call is emitted as `FILTER=no_cvg_drop;WRAP`, `flags=WRAP,HP,DLOOP`, `SVCONF=.`. *Asserted:* **0 PASS**,
and every emitted coordinate ≤ contig length (valid VCF). Under `HP_SV_PLOT_ALL` this call *is* plotted,
where samplot's linear view renders it misleadingly as a genome-spanning "duplication" — which is why
the gallery's status column (here showing the `WRAP` non-PASS status) matters.

### Annotation flag — passes, but flagged

**`sv_dloop` — a real deletion with a breakpoint in the control region (D-loop).**
A genuine 5599 bp deletion whose 5′ breakpoint (m.400) lies in the D-loop / control region. *Why:*
the D-loop is a known artifact hotspot, so breakpoints there get a `DLOOP` flag — but a *strong,
clean* junction there is still a real deletion and must PASS. This BAM proves the flag is set **without
suppressing a well-supported call** (flagging ≠ filtering). *Asserted:* `FILTER=PASS` **and** `DLOOP`
present in `flags`. It does double duty in the visualization tests: in **default** samplot mode it is
correctly **skipped** (the `HP_SV_PLOT_SKIP=HP,DLOOP,NUMT` artifact filter — `check_plots`'
`samplot_filter`), and in **all-plots** mode (`HP_SV_PLOT_ALL`) it **is** shown (`samplot_all_mode`),
carrying its PASS status and `DLOOP` flag in the gallery.

### Size-range deletions — the detection bounds

**`sv_del_500` / `sv_del_45` / `sv_del_13kb`** bracket the deletion size range.
- **`sv_del_500`** (500 bp) is the **small detectable** positive control: large enough that the aligner
  splits it (190 `SA` reads) so it PASSes via the DJ path with `AFC≈0.50` — the lower bound of clean
  detection. *Asserted:* `pass`.
- **`sv_del_45`** (45 bp) is a deliberate **hard-negative**: below `minsize=50`, *and* small enough that
  the aligner emits it as an internal CIGAR `D` (no `SA` split), so the caller — which seeds only from
  split-read junctions — sees nothing. *Asserted:* `no_record`. (Together with `sv_del_500` these pin
  the split-vs-gap boundary the research flagged: an event smaller than a read is jumped, not split.)
- **`sv_del_13kb`** (13 kb) is a genuine **majority-arc deletion** (`svlen > MTLEN/2`) that — unlike the
  origin artifact — has a real coverage drop, so it PASSes via **dosage**. *Asserted:* `pass`. This is
  the regression boundary that keeps the `WRAP` rule honest: a big "deletion" *without* a drop stays
  `WRAP` (the origin complement), one *with* a drop is a real call.

### Forward-looking — duplications, inversions & complex events

By **default** (deletion-only) these assert their DEL-only outcome (below). The **opt-in** tandem-DUP
(`HP_SV_DUP`) and inversion (`HP_SV_INV`) paths are now implemented
([`../docs/SV_EVENT_TYPES.md`](../docs/SV_EVENT_TYPES.md) §3.2/§3.3); `run_test.py`'s **`check_dup_inv`**
re-runs these fixtures with the flags on and asserts the real DUP/INV calls (`sv_dup`→DUP PASS with
`AFC`≈het, `sv_inv_small/large`→INV PASS, `sv_inv_lowhet`/`sv_inv_origin`→withheld, `sv_invdup`→INVDUP
flagged) plus that a real deletion is unaffected by the flags. Each BAM was verified to carry its real
signal first.

**Duplications — `sv_dup` (1 kb tandem) / `sv_dup_large` (5 kb tandem).** A tandem dup is the *mirror*
of a deletion: a reverse-order junction **plus a coverage gain** (verified `CVGR≈1.5` over the 5 kb arc,
i.e. `1+h`, not `0.5`). Today the reverse-order junction surfaces as a DEL-shaped record but the gain
guard (`ratio > 1+gainpad`) blocks PASS. *Asserted:* `no_pass`. They pin that the gain guard rejects
dups at both scales — and become the DUP positive controls once `SVTYPE=DUP` lands.

**Compound — `sv_dupdel` (partial duplication / dup-del).** One molecule carries **both** a dup junction
and an internal deletion (the clinically dominant mtDNA dup class — KSS/Pearson). It exposes a real
**KNOWN GAP**: the embedded deletion (m.6000–6500) shows a *relative* coverage drop against the
duplicated flanks and so **spuriously PASSes today** as a standalone deletion. *Asserted (as a
characterization test):* `known_fp` — `≥1 PASS`, labelled a documented gap that the DUP-aware caller
will close (it must not let an embedded del PASS under a net gain). When that lands, this assertion
flips and signals the fix.

**Inverted dup — `sv_invdup` (fold-back).** An inverted extra copy: opposite-strand `SA` (the inverted
arm) **co-located with a gain** — a Sniffles2 `INVDUP`, *not* a balanced inversion. The inverted arm is
strand-filtered and there is no same-strand junction, so today → `no_record`.

**Inversions — `sv_inv_small/large/origin/lowhet`.** Balanced inversions are the single biggest blind
spot: their junctions are **opposite-strand** (verified: 324/350/390/37 opposite-strand `SA` reads
respectively) which `extract_junctions` discards, and they are **copy-number-neutral** so dosage sees
nothing either. *Asserted:* `no_record` for all four — demonstrating the invisibility is independent of
size (`small`/`large`), heteroplasmy (`lowhet` @5%), and position (`origin`, near the artificial origin).
They become the INV detect-and-flag fixtures once the opposite-strand branch lands.

Small committed **real** chrM alignments from 1000 Genomes 30× high-coverage (GRCh38 chrM == rCRS ==
`RefSeq/chrM.fa`), extracted and realigned through the same circular path, subsampled to the
pipeline's working depth (~2000–2400×). They provide what simulation cannot: real NUMT-derived reads,
real D-loop complexity, and real error profiles. See [`real/README.md`](real/README.md).

| BAM | sample | role |
|---|---|---|
| `NA12718.chrM.bam` | NA12718 (CEU), ENA `ERR3239480`, ~2400× | **specificity** — healthy, no large deletion → **0 PASS** |
| `NA12748.chrM.bam` | NA12748 (CEU), ENA `ERR3239481`, ~2400× | specificity → **0 PASS** |
| `NA12775.chrM.bam` | NA12775 (CEU), ENA `ERR3239482`, ~2200× | specificity → **0 PASS** |
| `spike_del4977_h20.chrM.bam` | del4977 spiked into a real WT background @20% | **positive control on real data** → recovered **PASS + COMMON** |

*Why these matter:* the three healthy backgrounds are the real-world false-positive guard — a clean
caller returns zero PASS calls on real, error-laden, NUMT-bearing data, not just on tidy simulated
reads. The spike-in is the matching positive control: a known del4977 dropped into a *real* WT
background must still be recovered with the correct COMMON annotation. These same backgrounds form the
**REAL arm** of the quantitative LoD/accuracy sweep ([`real/lod_sweep.tsv`](real/lod_sweep.tsv) +
[`real/lod_report/`](real/lod_report/)): a heteroplasmy × depth grid plus an `HP_ARTIFACT` hard-negative
and an `ORIGIN` suppression check, establishing the ≈ 8% PASS LoD, 0/62 negatives PASS, and AUPRC ≈ 0.97
(see `docs/SV_METHODS.md` §10).

---

## 4. Derived negative / robustness inputs (not committed)

`run_test.py`'s `check_degenerate` mints malformed inputs **at test time** from `sv_wt.bam`, asserting
the caller fails *cleanly* (clear error, non-zero exit) and **never** throws a Python traceback:

- **empty BAM** (valid header, 0 reads) → header-only VCF, 0 records, exit 0.
- **unindexed BAM** → clean ERROR, non-zero exit.
- **wrong contig name** → clean ERROR.
- **`mtlen` mismatch** → clean ERROR.

These guard the input-validation paths so a bad cohort sample produces a legible failure, not a crash.

---

## 5. Regenerating

The simulator is deterministic, so regenerated BAMs are reproducible. Requires `minimap2` +
`samtools` + `perl` (and `RefSeq/chrMC.*`):

```bash
python3 test/sv/make_testdata.py -ref RefSeq/chrM.fa -out test/sv/fastq
bash    test/sv/gen_bams.sh        test/sv/fastq test/sv/bams
cp      test/sv/fastq/truth.tsv    test/sv/truth.tsv
```

Real-data BAMs: [`real/gen_real.sh`](real/gen_real.sh) (healthy backgrounds) and
[`real/gen_spike.sh`](real/gen_spike.sh) (the del4977 spike-in) — need `samtools` + `minimap2` +
network access. Run the whole evaluation with `bash test/sv/run_test.sh`.
