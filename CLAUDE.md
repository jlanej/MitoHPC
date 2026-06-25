# CLAUDE.md — MitoHPC working notes & guardrails

This file is durable context for Claude (and humans) working in this repo. Read it before
touching anything related to **structural variant (SV) calling**.

---

## 0. THE NON-NEGOTIABLE CONSTRAINT (read this first)

We are adding **mitochondrial structural-variant (SV) calling** to MitoHPC as a **brand-new,
standalone, purely-additive module**.

> **It is absolutely critical that the SV module NEVER changes the output of any existing
> deliverable.** Existing SNV/heteroplasmy/copy-number/haplogroup outputs are byte-for-byte
> frozen. The SV module is built *on top of* the existing machinery (especially the circular-
> genome handling) and writes *new* files only.

Concretely, this means:

1. **Additive only.** New behavior must be gated behind a flag and must default to **off**, so
   that a default run produces exactly the same files it does today. Follow the existing
   precedent: the `HP_V` ("SV caller") environment variable already gates the optional GRIDSS
   path and is empty by default (see `scripts/init.sh`, `scripts/run.sh`, `scripts/filter.sh`).
2. **Never edit an existing output file or its schema.** No new columns, no reordered rows, no
   new INFO/FORMAT fields in the existing `{mutect2,mutserve,freebayes}.*` VCFs, `count.tab`,
   `cvg.tab`, `*.summary`, `*.tab`, `*.haplogroup*`, `*.haplocheck*`, `*.fa`, `*.merge.bed`, etc.
3. **Read, don't mutate.** The SV module may *consume* existing intermediate artifacts
   (the realigned circular-aware BAM, per-base coverage, the split-read BED) but must not alter,
   move, or delete them in a way that changes existing behavior.
4. **New outputs live under a clearly separate namespace** (e.g. `*.sv.*` files and/or an
   `sv/` subdirectory), aggregated by a *separate* summary step — NOT folded into
   `getSummary.sh`'s existing tables. If aggregation is wanted, add it conditionally (only when
   the SV flag is set), mirroring how `getSummary.sh` only emits `$V.00.concat.vcf` when
   `HP_V` is set.
5. **Verification of "no regression" is part of done.** Before/after a change, a default-config
   run must produce identical existing outputs (diff the output dir). See `examples1/`,
   `examples2/` for fixtures.

If a requested change cannot be done additively, **stop and surface the conflict** rather than
modifying a frozen deliverable.

---

## 1. What MitoHPC is

A bash + perl pipeline that estimates mitochondrial **copy number** and **SNV heteroplasmy** from
short-read WGS (BAM/CRAM pre-aligned to a whole-genome reference). Reference: Battle et al., NAR
2022 (PMC9112767). Driven entirely by `HP_*` environment variables.

**Control flow:**
- `scripts/init.sh` (and `init.hs38DH.sh`, `init.mm39.sh`) set every `HP_*` variable and build
  `$HP_IN` (a TSV: `sampleName <tab> inputFile <tab> outputPrefix`).
- `scripts/run.sh` reads `$HP_IN` and **emits a shell script** (`run.all.sh`) that calls
  `scripts/filter.sh` once per sample, then `scripts/getSummary.sh` once at the end.
- `scripts/filter.sh` is the **per-sample workhorse** (subsample → realign → call SNVs →
  consensus → 2nd-iteration realign/recall). This is where an SV step would hook in per sample.
- `scripts/getSummary.sh` **aggregates** per-sample outputs into the cohort-level tables.

Key vars (defaults in `scripts/init.sh`): `HP_M` (SNV caller: mutect2|mutserve|freebayes),
`HP_I` (0/1/2 iterations), `HP_L` (subsample read cap, 222000 ≈ 2000× chrM), `HP_DP` (min depth),
`HP_E` (circularization extension, **300**), `HP_MTLEN` (**16569**), `HP_V` (SV caller, empty).

---

## 2. Circular-genome machinery (the part the SV module must reuse)

chrM is circular but stored linearly (1..16569). Reads and variants spanning the artificial
origin/breakpoint need special handling. MitoHPC already solves this; **reuse it, do not
reinvent it.**

- **`HP_E=300`** — extension/overlap length used everywhere for circularization.
- **Circularized reference `chrMC` (`$HP_MTC`)** — chrM with its first `HP_E` bp appended to the
  end (built by `scripts/circFasta.sh`). **Alignment target**: `filter.sh` runs `bwa mem` against
  `$HP_RDIR/$HP_MTC` so reads crossing the origin map contiguously.
- **`scripts/circSam.pl`** — post-processes alignments to the circular reference: wraps
  coordinates `> MTLEN` back into `1..MTLEN` and splits CIGARs at the boundary, emitting
  supplementary alignments (`flag |= 2048`). It is run twice with `-offset 0` (→ `$O.bam`) and
  `-offset $HP_E` (→ `$OR.bam`).
- **Rotated reference `chrMR` (`$HP_MTR`)** — reference rotated by `HP_E` (built by
  `scripts/rotateFasta.sh`), used to **call variants in the origin region** that a linear
  reference would split. Coordinates are then mapped back with
  `POS_final = (POS - HP_E) % HP_MTLEN` (see `filter.sh` and the `fixsnpPos.pl` logic).
- Net effect available to an SV module: a coordinate-correct, circular-aware, deduplicated,
  subsampled **`$O.bam`** on a 16569 bp linear chrM, where origin-spanning reads are represented
  as supplementary/split alignments.

**Implication for SV calling:** an SV / junction that crosses the origin must be detectable. The
existing split-read representation already encodes origin-crossing reads, so a split-read-based
caller inherits circular-correct **detection input** "for free" — the wrapped split reads plus
coverage windows that wrap modulo `HP_MTLEN` (both true in `callsv.py`). This is **not** the same as
resolving every origin event: `callsv.py` computes junction `svlen` linearly (`bp3-bp5-1`), so a
deletion whose *deleted arc* spans the artificial origin is detected-but-suppressed (negative svlen →
dropped, or `WRAP`-flagged out of PASS), never a wrong call but also not yet resolved (del-vs-
complementary-arc-dup; see `docs/SV_METHODS.md` §8). Any coverage/segmentation method must treat
position 1 and 16569 as adjacent.

---

## 3. SV-relevant primitives that ALREADY EXIST (available to build on — reuse is OPTIONAL)

Reusing these is a convenience, **not a requirement**. There is nothing wrong with building our
own, better SV signal from scratch if it improves on what is partially baked in today — the only
hard rules are §0 (additive, freeze existing outputs) and §2 (reuse the *circular* coordinate
machinery so origin-crossing events are handled correctly).

- **`$O.sa.bed` — latent split-read breakpoints, already computed but UNUSED downstream.**
  `filter.sh:120` runs `samtools view $O.bam | sam2bedSA.pl | uniq.pl -i 3 | sort > $O.sa.bed`
  (and `$OS.sa.bed` in iteration 2). `scripts/sam2bedSA.pl` parses each read's `SA:Z:`
  supplementary-alignment tag into a **pair of BED intervals** (the two segments of a split
  read) — i.e. candidate breakpoint junctions. It is one possible starting signal, but it is
  coarse (it only uses `SA` tags, no soft-clip clustering, no breakpoint refinement, no
  coverage). We are free to ignore it and compute richer breakpoints ourselves.
- **`$O.cvg` — per-base coverage** (`bedtools genomecov -d`) and `$O.cvg.stat` summary. The basis
  for read-depth / coverage-drop deletion evidence and heteroplasmy estimation. **Persists** to
  the output dir.
- **Common-deletion assets already in `RefSeq/`**: `chrM.8455-8482_13460-13474.fa` (and `_33mer`,
  `_23mer` siblings) are **synthetic junction sequences** spanning the canonical 4977 bp common
  deletion (`del4977`), embedding the 13 bp direct repeat `ACCTCCCTCACCA`. These are ready-made
  positive controls / targeted-detection references.
- **False-positive masking resources in `RefSeq/`**: `HP.bed.gz` (homopolymers),
  `HS.bed.gz` (hotspots), `MCC.bed.gz` (per-region mappability/complexity scores),
  `NUMT.vcf.gz` (NUMT-like sites), `$HP_RNUMT` (nuclear NUMT loci).
- **GRIDSS stub** — `HP_V=gridss` runs GRIDSS on `$O.bam` and post-processes to `$O.gridss.00.vcf`
  via the `scripts/gridss.vcf` header template (`filter.sh:162-171`). Off by default,
  experimental, partially baked. Useful mainly as the **precedent for how an SV path was wired in
  additively and placed BEFORE the BAM cleanup** — see §3a. GRIDSS on a tiny high-copy circular
  genome is heavyweight; see the research doc before relying on it.

## 3a. ⚠️ The realigned BAM is TRANSIENT — it is deleted at the end of every per-sample run

This is the single most important practical constraint for the SV module.

`filter.sh` builds the circular-aware, subsampled, deduplicated chrM alignment `$O.bam` (and the
rotated `$OR.bam`, and iter-2 `$OS.bam`/`$OSR.bam`) but **deletes all of them before it exits**:
- `filter.sh:232` → `rm -f $O.bam* $OR.bam*`
- `filter.sh:237` and `filter.sh:282` → `rm -f $OS.bam* $OSR.bam*`

The subsampled reads `$O.fq` are also deleted (`filter.sh:229`). The user's invocation path
(`mitohpc-batch-container.sh $base 25 docker://…` → in-container `run_parallel_mitohpc.sh` →
`init.sh` → `run.sh` → `filter.sh` → `getSummary.sh`) therefore leaves **no BAM and no FASTQ** in
`out/`. What persists per sample: `$O.cvg`, `$O.cvg.stat`, `$O.sa.bed`, the `*.vcf`s, `$OS.fa`,
`$OS.merge.bed`, `$OS.idxstats`, `$O.count`, haplogroup/haplocheck.

Consequence — the SV computation cannot run as a naive *post-hoc* pass over `out/`, because a
real SV caller needs the alignment (or the reads). Two viable patterns (pick per §4):
- **(A) Hook the SV step inside `filter.sh` while `$O.bam` is still alive** — i.e. in a gated
  block placed *before* the `rm -f $O.bam*` at line 232, exactly where the GRIDSS stub already
  sits. The output is still standalone/additive; only the *computation* is co-located. Preferred.
- **(B) Conditionally preserve the alignment** — when the SV flag is set, copy/keep `$O.bam` (e.g.
  as `$O.sv.bam`) past cleanup so a separate `callSV.sh` can consume it. More I/O, but keeps the
  SV logic fully out of `filter.sh`.

"Standalone" here means the **output** is separate and additive (per the user) — it does **not**
require the SV code to run as a detached pass with no access to the pipeline's alignment.

---

## 4. Integration contract for the new SV module

When implementing, honor these rules (they operationalize §0):

1. **Gate behind a flag**, default off (extend/parallel the `HP_V` pattern; validate it in
   `run.sh` exactly like the existing `HP_V`/`HP_M` checks).
2. **Get access to the alignment before it is deleted (§3a).** The SV logic should live in a
   **separate script** (e.g. `callSV.sh`) for isolation, but it must be **invoked from within
   `filter.sh` before the `rm -f $O.bam*` at line 232** (i.e. where the GRIDSS stub sits), so it
   receives the live `$O.bam` — OR the flag-gated path must preserve `$O.bam` for a later pass.
   Do **not** wire the SV step as a post-`filter.sh` step in `run.sh`: by then the BAM is gone.
   Keep the touch to `filter.sh` to a single small gated block (one `if [ $HP_SV ]` calling
   `callSV.sh`) to minimize blast radius on the frozen path. Reusing `$O.sa.bed` is optional —
   `callSV.sh` may recompute its own (better) breakpoint signal directly from `$O.bam`.
3. **Write only new files** under an `*.sv.*` namespace (and/or `sv/` subdir). Suggested:
   `$O.sv.vcf` (per-sample SV calls) + `$O.sv.tab` (flat table). Aggregate with a **separate**
   `getSVSummary.sh` invoked only when the flag is set; do not edit `getSummary.sh`'s existing
   tables.
4. **Reuse the circular machinery** (`HP_E`, `circSam.pl`'s split-read output, the `% HP_MTLEN`
   wrap). Origin-crossing junctions must be handled; do not introduce a new linearization scheme.
5. **Be defensible**: prefer split-read + coverage corroboration; apply explicit, documented
   false-positive controls (NUMT regions via `$HP_RNUMT`, homopolymers via `HP.bed.gz`, minimum
   supporting reads, minimum deletion size, mapping quality). Record thresholds as `HP_*` vars.
6. **Quantify heteroplasmy explicitly** and document the formula in the VCF header
   (e.g. `AF = junction_reads / (junction_reads + spanning_reads)`, and/or a coverage-ratio
   estimate) so the number is interpretable and reproducible.
7. **Keep dependencies light and best-fit.** The SV caller uses **Python 3 + `pysam`** (the one
   added dependency — a pip manylinux wheel that bundles htslib; no compiler, installed in Docker).
   Language is chosen for fit, not repo tradition: `pysam` gives robust in-process BAM/CIGAR/SA/depth
   handling. `samtools`/`bedtools`/`bcftools` remain for the rest of the pipeline.

---

## 5. Pointers

- Design & literature: **`docs/SV_CALLING.md`** (state-of-the-art review, tool comparison,
  algorithm primitives, and the phased implementation plan for MitoHPC).
- Method as implemented (keep in sync with the code): **`docs/SV_METHODS.md`** — intuitive but
  precise description of `callsv.py`/`callSV.sh`, formulas, schema, parameters.
- Deferred-feature design note: **`docs/SV_DELDUP_RESOLUTION.md`** — plan to resolve the circular
  DEL-vs-complementary-arc-DUP ambiguity via MitoSAlt-style origin (OriH/OriL) preservation, turning
  today's `WRAP`-withheld origin-crossing junctions into opt-in DEL/DUP calls. NOT implemented.
- Existing pipeline outputs/legend: `README.md` (the `## OUTPUT ##` section is the list of frozen
  deliverables).
- Test fixtures: `examples1/`, `examples2/`.

### SV module (v1 — implemented, default off via `HP_SV`)

- `scripts/callsv.py` — **the caller** (Python 3 + `pysam`): split-read junction extraction +
  clustering from `SA:Z:` tags, in-process per-base depth (`count_coverage`), coverage
  corroboration, two heteroplasmy estimates (AFJ/AFC) + AFDIFF QC, FP flags
  (REPEAT/NUMT/HP/DLOOP/WRAP), PASS/FILTER logic, VCF+tab. Replaces the former perl
  `sa2del.pl`/`svCall.pl` (field-for-field parity verified).
- `scripts/callSV.sh` — thin per-sample driver; resolves `HP_SV_*` thresholds + masks and runs
  `$HP_PYTHON(=python3) callsv.py` on the live `$O.bam`, writing only `$O.sv.vcf` + `$O.sv.tab`.
  Invoked by the gated block in `filter.sh` (after the GRIDSS block, before the BAM `rm`).
- `scripts/sv.vcf` — VCF header template (mirrors `gridss.vcf`).
- `scripts/svplot.sh` — optional **samplot** visualization (default-off via `HP_SV_PLOT`). Invoked by
  `callSV.sh` while `$O.bam` is still alive; runs `samplot plot` on the *visualizable* subset (PASS,
  `AFC>=HP_SV_PLOT_MINAF`[0.03], breakpoints NOT in `HP`/`DLOOP`/`NUMT` artifact regions — all
  configurable; `REPEAT` is NOT skipped so the real common deletion shows), with a **gene annotation
  track** (`samplot -A genes.bed.gz` by default; `HP_SV_PLOT_ANNOT` is a configurable comma-separated list
  of tabixed `$HP_RDIR` BEDs), writing `${O}.sv.<bp5>_<end>.png` + a manifest `${O}.sv.plots.tsv`. samplot
  is installed in the image (Dockerfile, pip from a pinned GitHub commit `2929e4a` — PyPI's `0.0.1` is broken;
  `2929e4a` = the v1.3.0 release plus Python-3.11/numpy forward-compat fixes — with jinja2; no conda).
  **`matplotlib` is pinned `==3.6.3` (and `numpy<2`):** matplotlib `>=3.7` regressed samplot's axes — a
  spurious `0..1` normalized axis is overprinted on the genomic x-axis and coverage/insert-size y-axes,
  making breakpoints unreadable so calls appear not to line up with the reads (samplot issues #189/#201;
  glaring on narrow/small-deletion windows, subtle on wide ones). This — not the samplot commit — was the
  real cause of bad plots; an earlier unpinned `matplotlib` + `--no-deps` install let `>=3.7` in. A build-time
  assert fails the image if matplotlib ever resolves `>=3.7`. `svplot.sh` degrades gracefully if samplot is
  absent. CI smoke-tests it in the image, asserting `matplotlib<3.7` (`.github/workflows/docker-publish.yml`).
  Visualization reuses the SAME circular-aware `$O.bam` the caller reads, so for the plotted (PASS, non-`WRAP`)
  subset samplot's coordinates match the caller's `1..MTLEN` frame; only genuinely origin-crossing reads carry
  un-wrapped SA tags (circSam.pl does not rewrite them), and those calls are `WRAP`-flagged out of the plotted
  set, so the circularization blind spot never reaches a rendered plot.
- `scripts/getSVSummary.sh` — cohort aggregator (tidy `$ODIR/sv.tab`, `bcftools merge` matrix
  `$ODIR/sv.merged.vcf.gz` with `NS` recurrence, sites union `$ODIR/sv.sites.vcf.gz`, and the
  interactive `$ODIR/sv.report.html`); SEPARATE from `getSummary.sh`, gated on `HP_SV`.
- `scripts/svReport.py` — builds a **self-contained, offline, interactive HTML report** (vanilla
  SVG/JS, no deps) from `sv.tab` + `genes.bed.gz`: circular mtDNA map + linear genome browser with
  gene/OXPHOS-complex annotation, a per-position deletion-frequency track, VAF-coloured calls, live
  filtering, summary stats, a recurrence table, and (when `--plots` is given) an **interactive samplot
  gallery** — a table whose rows reveal the base64-embedded samplot PNG on click (still single-file/offline).
  The gallery is **subsampled to one representative call per breakpoint cluster** (`--plot-dedup` bp,
  default 25 = the recurrence-table rounding; highest-heteroplasmy call kept, `samples` column shows the
  cluster size), and the report spells this out; only representative PNGs are embedded. Each gallery row
  carries a **status** column (PASS green / non-PASS amber, via theme-aware `--ok`/`--warn` CSS vars) read
  from the manifest's appended `filter`/`flags` (png stays column 7 for back-compat); the caption shows the
  flags. The gallery prose is JS-built from `meta` so it stays accurate in both modes.
- **`HP_SV_PLOT_ALL=1` — "show everything" mode** (rarely used; ON for the committed example so it
  visualizes ALL samples). It flips `svplot.sh`'s filter DEFAULTS to permissive (plot EVERY call: no
  PASS / heteroplasmy / `HP`/`DLOOP`/`NUMT` filtering — an explicit `HP_SV_PLOT_PASS`/`MINAF`/`SKIP`
  still overrides) AND makes `getSVSummary.sh` default `--plot-dedup` to 0 (no representative collapsing)
  and pass `--plot-all` (report wording). Non-PASS / `WRAP` calls are then shown — the status column and
  the §2 origin caveat matter here: a `WRAP`/origin-crossing call renders MISLEADINGLY under samplot's
  linear view (e.g. a genome-spanning "duplication"), so its amber non-PASS status is the disambiguator.
  Forwarded by `mitohpc-batch-container.sh`.
- Output is general VCF/SV best practice (NOT this repo's other VCFs): `##contig`/`##reference`/
  provenance headers, sample-named genotype column, `HOMLEN`/`HOMSEQ`/`DELCLASS`/`CIPOS`/`CIEND`,
  `SVCLAIM`, `COMMON`, `GENE`/`NGENE`, `HGVS`, `FORMAT GT:DP:AD:AF:SR`; VCF 4.2 (negative `SVLEN`).
- `scripts/recallSV.sh` — **fast SV re-run**: re-runs ONLY the SV caller (+ samplot + cohort report)
  on persisted circular-aware BAMs, skipping the expensive subsample→realign→SNV front end. Enabled by
  `HP_SV_KEEPBAM=1` on a full run (gated copy in `callSV.sh`: `$O.bam`→`$O.sv.bam`, which `filter.sh`'s
  `rm -f $O.bam*` does NOT match, so it survives — no frozen-file edit; ~10 MB/sample). Then
  `HP_SV_RECALL=1 mitohpc-batch-container.sh <dir> …` (or `recallSV.sh <out> <jobs>`) refreshes `sv.*`
  in minutes (~3 s/sample) instead of a multi-hour full re-run — byte-identical to a fresh SV call; never
  touches the frozen SNV/CN deliverables. The self-copy guard in `callSV.sh` lets it re-run on the
  persisted BAM safely.
- Wiring: `HP_SV` + `HP_SV_*` tunables in `init.sh`; validation + exports + gated summary in
  `run.sh`; one gated block in `filter.sh`. With `HP_SV` empty the pipeline is unchanged. SV
  visualization defaults ON and is configured/forwarded by `mitohpc-batch-container.sh`
  (`HP_SV_PLOT` enable + `HP_SV_PLOT_*` filter), so a default batch run also produces the samplot
  gallery; `HP_SV_KEEPBAM`/`HP_SV_RECALL` (re-run speedup) are forwarded too; the bare pipeline is unchanged.
- Tests/mock data: **`test/sv/`** — `run_test.py` (via `run_test.sh`) runs 24 checks (+2 samplot
  visualization checks when `samplot` is on PATH, e.g. inside the image): 10 mock
  scenarios (multi-deletion, near-homoplasmy, tandem-dup-not-called, origin-crossing, D-loop,
  low-coverage, …) against committed mock BAMs (`test/sv/bams/`, ~13 MB), degenerate inputs, cohort
  recurrence, a `bcftools` VCF-spec gate, a schema check on the committed example outputs, plus
  real-data vetting under `test/sv/real/`: 3 healthy 1000G high-cov chrM → 0 PASS (specificity) and a
  del4977 spiked into a real WT background → recovered PASS+COMMON (positive control). The v2 caller
  reports coverage-dosage AF (AFC) as primary heteroplasmy with a corrected junction VAF (AFJ) and a
  junction-strong PASS path; see `docs/SV_METHODS.md`.
- Committed example outputs: **`test/sv/example/`** — a representative per-sample VCF/tab, the cohort
  `sv.tab`/`sv.merged.vcf.gz`/`sv.sites.vcf.gz`, and the interactive `sv.report.html`; regenerate with
  `test/sv/make_example.sh` (paths sanitized to repo-relative). See `test/sv/README.md`.

## 6. Conventions

- Scripts live in `scripts/`; shell stages are `*.sh`, legacy helpers are `*.pl` (perl), newer
  logic is Python 3 (`*.py`, e.g. the SV caller `callsv.py`). Use the best-fit language. Match existing
  naming (`fix*Vcf.pl`, `*2*.pl`, `filter*.sh`) and the env-var-driven style.
- Keep edits surgical and reversible. When in doubt about whether something is a frozen
  deliverable, treat it as frozen and ask.
