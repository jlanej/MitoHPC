#!/usr/bin/env python3
"""
LoD + accuracy SWEEP for the MitoHPC SV caller — the data-generation half of the in-silico
evaluation (the report/figures are built by lod_report.py from the TSV this writes).

Design (CLSI EP17-A2 LoD surface + spike-in benchmark; see docs/SV_METHODS.md §10.1):
  TWO ARMS at a shared heteroplasmy(VAF) × depth lattice, each replicate independently seeded:
    SIM  — simulate a WT+EVENT circular-genome mixture, realign through the pipeline's circular
           path (minimap2 -ax sr chrMC → -F 0x90C → circSam.pl), call. Cheap → carries the grid + CIs.
    REAL — spike EVENT reads into a *real* 1000G WT chrM background (extract WT FASTQ from the real
           BAM once, add D event read-pairs, realign, call). Real error/NUMT/coverage structure →
           validates the simulated grid (sim↔real concordance).
  VARIANTS: del4977 (repeat-mediated, COMMON), NONREP (~5 kb non-repeat → recovery is not repeat-
    specific), HP_ARTIFACT (control-region del(305,965): both ends in the D-loop poly-C homopolymer
    so nfragile≥2 — the dominant real-cohort FALSE POSITIVE, the hard-negative that defends the
    SVCONF fragile penalty; clears HP_SV_MINSIZE so it is actually emitted), ORIGIN (origin-crossing
    delwrap → must be WRAP-suppressed, SVCONF='.').
  NEGATIVES: the VAF=0 column of each arm (LoB) + the 3 real healthy BAMs called directly.
  TRUTH labels: is_true=1 only for a genuinely recoverable deletion (del4977/NONREP); is_true=0 for
    WT, HP_ARTIFACT, ORIGIN — so sensitivity, specificity and the SVCONF artifact-penalty are
    measured on the right populations.

Reproducibility: every replicate's exact integer seed is a deterministic crc32 of
(variant, vaf, depth, rep) — injective over the grid (asserted), independent of PYTHONHASHSEED, and
written into the TSV. Re-running with the same grid reproduces every row.

Output: a single tidy TSV (one row per replicate-call). Needs samtools + minimap2 + python3/pysam.
  HP_PYTHON=/venv/bin/python python3 lod_sweep.py --quick            # fast smoke grid (committed data)
  HP_PYTHON=/venv/bin/python python3 lod_sweep.py --full             # publication grid (run on a cluster)
"""
import argparse
import importlib.util
import os
import subprocess
import sys
import tempfile
import time
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SDIR = os.environ.get("HP_SDIR", os.path.join(ROOT, "scripts"))
RDIR = os.environ.get("HP_RDIR", os.path.join(ROOT, "RefSeq"))
PY = os.environ.get("HP_PYTHON", sys.executable)
MTLEN = 16569
RLEN = 150
BP_TOL = 30            # summed |bp5_err|+|bp3_err| tolerance to call a deletion "detected"

spec = importlib.util.spec_from_file_location("mt", os.path.join(HERE, "make_testdata.py"))
mt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mt)

# variant -> retained breakpoints + kind + truth label. is_true=1 ONLY for a recoverable deletion.
VARIANTS = {
    "del4977":     dict(bp5=8469, bp3=13447, kind="del", is_true=1),
    "NONREP":      dict(bp5=5999, bp3=10999, kind="del", is_true=1),
    "HP_ARTIFACT": dict(bp5=305,  bp3=965,   kind="del", is_true=0),   # D-loop homopolymer hard-negative
    "ORIGIN":      dict(bp5=16400, bp3=200,  kind="delwrap", is_true=0),  # origin-crossing -> suppress
}

TSV_COLUMNS = [
    "arm", "del_variant", "is_true", "target_vaf", "target_depth", "background", "rep_index", "seed",
    "truth_bp5", "truth_bp3", "truth_svlen", "detected", "passed", "matched_to_truth",
    "n_calls_total", "n_calls_pass", "filter", "svclaim", "called_bp5", "called_bp3",
    "bp5_err", "bp3_err", "svlen_err", "bp_err_abs", "afc", "afj", "afc_err", "afj_err",
    "jr", "sr", "srcons", "srsb", "jsup", "cvgr", "nfragile", "flags",
    "svconf", "svimpact", "svimpact_band", "runtime_s",
]

NA = "NA"


def seed_for(variant, vaf, depth, rep):
    """Deterministic, PYTHONHASHSEED-independent, injective over the grid (asserted in main)."""
    return zlib.crc32(("%s|%.4f|%d|%d" % (variant, vaf, depth, rep)).encode()) & 0x7FFFFFFF


def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def realign(work, r1, r2, sample):
    """minimap2 -ax sr chrMC -> -F 0x90C -> circSam.pl -> sorted $O.bam (the caller's input)."""
    bam = os.path.join(work, sample + ".bam")
    sh("minimap2 -ax sr '%s/chrMC.fa' '%s' '%s' 2>/dev/null "
       "| samtools view -h -F 0x90C - "
       "| '%s/circSam.pl' -ref_len '%s/chrM.fa.fai' -offset 0 "
       "| samtools sort -o '%s' --write-index - 2>/dev/null" % (RDIR, r1, r2, SDIR, RDIR, bam))
    return bam


def call_rows(sample, bam, work):
    """Run callSV.sh and return ALL emitted call rows as dicts (empty list if none)."""
    pref = os.path.join(work, sample)
    subprocess.run(["bash", os.path.join(SDIR, "callSV.sh"), sample, bam, pref],
                   env=dict(os.environ, HP_SDIR=SDIR, HP_RDIR=RDIR, HP_PYTHON=PY, HP_MT="chrM",
                            HP_MTLEN=str(MTLEN)), capture_output=True, text=True)
    tab = pref + ".sv.tab"
    if not os.path.exists(tab):
        return []
    with open(tab) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if len(lines) < 2:
        return []
    hdr = lines[0].lstrip("#").split("\t")
    return [dict(zip(hdr, ln.split("\t"))) for ln in lines[1:]]


def best_match(rows, bp5, bp3):
    """Closest call to truth by summed breakpoint distance, or None."""
    best, bestd = None, None
    for r in rows:
        try:
            d = abs(int(r["pos_bp5"]) - bp5) + abs(int(r["end_bp3"]) - (bp3 - 1))
        except (ValueError, KeyError):
            continue
        if bestd is None or d < bestd:
            best, bestd = r, d
    return best, bestd


def n_fragile(flags):
    s = set((flags or "").split(","))
    return sum(1 for f in ("DLOOP", "HP", "NUMT", "WRAP") if f in s)


def record(arm, variant, vaf, depth, background, rep, seed, rows, truth_bp5, truth_bp3, kind, t0):
    """Build one TSV row from the caller output for one replicate."""
    truth_svlen = (truth_bp3 - truth_bp5 - 1) if kind == "del" else mt.event_svlen(kind, truth_bp5, truth_bp3, MTLEN)
    n_total = len(rows)
    n_pass = sum(1 for r in rows if r.get("filter") == "PASS")
    m, md = best_match(rows, truth_bp5, truth_bp3)
    matched = 1 if (m is not None and md is not None and md <= BP_TOL) else 0
    row = {c: NA for c in TSV_COLUMNS}
    row.update(arm=arm, del_variant=variant, is_true=VARIANTS.get(variant, {}).get("is_true", 0),
               target_vaf="%.4f" % vaf, target_depth=depth, background=background, rep_index=rep,
               seed=seed, truth_bp5=truth_bp5, truth_bp3=truth_bp3, truth_svlen=truth_svlen,
               detected=matched, n_calls_total=n_total, n_calls_pass=n_pass,
               runtime_s="%.2f" % (time.time() - t0))
    if matched:
        cb5, cb3 = int(m["pos_bp5"]), int(m["end_bp3"])
        b5e, b3e = cb5 - truth_bp5, cb3 - (truth_bp3 - 1)
        afc = _f(m.get("af_coverage")); afj = _f(m.get("af_junction"))
        row.update(passed=1 if m["filter"] == "PASS" else 0, matched_to_truth=1,
                   filter=m["filter"], svclaim=m.get("svclaim"), called_bp5=cb5, called_bp3=cb3,
                   bp5_err=b5e, bp3_err=b3e, svlen_err=int(m["svlen"]) - truth_svlen,
                   bp_err_abs=abs(b5e) + abs(b3e),
                   afc=m.get("af_coverage"), afj=m.get("af_junction"),
                   afc_err=("%.4f" % (afc - vaf)) if afc is not None else NA,
                   afj_err=("%.4f" % (afj - vaf)) if afj is not None else NA,
                   jr=m.get("jr"), sr=m.get("sr"), srcons=m.get("srcons"), srsb=m.get("srsb"),
                   jsup=m.get("jsup"), cvgr=m.get("cvgr"), nfragile=n_fragile(m.get("flags")),
                   flags=m.get("flags"), svconf=m.get("svconf"), svimpact=m.get("svimpact"),
                   svimpact_band=m.get("svimpact_band"))
    else:
        row.update(passed=0, matched_to_truth=0, filter="UNDETECTED")
    return row


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


# --------------------------------------------------------------------------- #
def sim_run(ref, variant, vaf, depth, rep, work):
    v = VARIANTS[variant]
    seed = seed_for(variant, vaf, depth, rep)
    import random
    rng = random.Random(seed)
    glen = len(ref)
    r1, r2 = os.path.join(work, "r1.fq"), os.path.join(work, "r2.fq")
    with open(r1, "w") as f1, open(r2, "w") as f2:
        n_wt = round(depth * (1 - vaf) * glen / RLEN)
        mt.emit_reads(f1, f2, ref + ref, n_wt, RLEN, 300, 450, 0.001, rng, "wt")
        if vaf > 0:
            eg = mt.event_genome(ref, v["kind"], v["bp5"], v["bp3"])
            n_e = round(depth * vaf * len(eg) / RLEN)
            mt.emit_reads(f1, f2, eg + eg, n_e, RLEN, 300, 450, 0.001, rng, "ev")
    bam = realign(work, r1, r2, "s")
    t0 = time.time()
    rows = call_rows("s", bam, work)
    return record("SIM", variant, vaf, depth, "sim", rep, seed, rows, v["bp5"], v["bp3"], v["kind"], t0)


def extract_wt_fastq(real_bam, work):
    """Collate a real chrM BAM to FR FASTQ (WT background). Returns (r1, r2, n_pairs)."""
    w1, w2 = os.path.join(work, "wt1.fq"), os.path.join(work, "wt2.fq")
    sh("samtools collate -u -O '%s' 2>/dev/null | samtools fastq -1 '%s' -2 '%s' -0 /dev/null "
       "-s /dev/null -n - 2>/dev/null" % (real_bam, w1, w2))
    with open(w1) as fh:
        n = sum(1 for _ in fh) // 4
    return w1, w2, n


def real_run(ref, variant, vaf, bg_name, wt1, wt2, w_pairs, rep, work):
    """Spike `variant` reads into pre-extracted real WT FASTQ at heteroplasmy `vaf`."""
    v = VARIANTS[variant]
    seed = seed_for(variant + ":" + bg_name, vaf, 0, rep)
    r1, r2 = os.path.join(work, "r1.fq"), os.path.join(work, "r2.fq")
    sh("cp '%s' '%s' && cp '%s' '%s'" % (wt1, r1, wt2, r2))
    if vaf > 0:
        import random
        rng = random.Random(seed)
        eg = mt.event_genome(ref, v["kind"], v["bp5"], v["bp3"])
        d = int((vaf / (1 - vaf)) * w_pairs * len(eg) / MTLEN)
        d1, d2 = os.path.join(work, "d1.fq"), os.path.join(work, "d2.fq")
        with open(d1, "w") as f1, open(d2, "w") as f2:
            mt.emit_reads(f1, f2, eg + eg, d, RLEN, 300, 450, 0.001, rng, "spike")
        sh("cat '%s' >> '%s' && cat '%s' >> '%s'" % (d1, r1, d2, r2))
    bam = realign(work, r1, r2, "r")
    # achieved depth bucket (mean coverage) for the real arm (depth is whatever the BAM yields)
    cov = sh("samtools depth -a '%s' 2>/dev/null | awk '{s+=$3;n++}END{if(n)printf \"%%.0f\",s/n}'" % bam).stdout
    depth = int(cov) if cov.strip().isdigit() else 0
    t0 = time.time()
    rows = call_rows("r", bam, work)
    return record("REAL", variant, vaf, depth, bg_name, rep, seed, rows, v["bp5"], v["bp3"], v["kind"], t0)


def real_wt_run(real_bam, bg_name, work):
    """Call a real healthy chrM BAM directly (a real negative / specificity control)."""
    bam = os.path.join(work, "rw.bam")
    sh("samtools sort -o '%s' --write-index '%s' 2>/dev/null" % (bam, real_bam))
    t0 = time.time()
    rows = call_rows("rw", bam, work)
    return record("REAL", "WT", 0.0, 0, bg_name, 0, 0, rows, VARIANTS["del4977"]["bp5"],
                  VARIANTS["del4977"]["bp3"], "del", t0)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=os.path.join(HERE, "real", "lod_sweep.tsv"))
    ap.add_argument("--quick", action="store_true", help="small grid for a committed-data smoke run")
    ap.add_argument("--full", action="store_true", help="publication grid (slow; run on a cluster)")
    ap.add_argument("--vafs", default="")
    ap.add_argument("--depths", default="")
    ap.add_argument("--reps", type=int, default=0)
    ap.add_argument("--real-bams", default="NA12718,NA12748,NA12775")
    ap.add_argument("--arm", choices=["both", "sim", "real"], default="both",
                    help="run only one arm (e.g. --arm real to append the real arm to an existing sim TSV)")
    ap.add_argument("--append", action="store_true", help="append to --out instead of overwriting (no header)")
    args = ap.parse_args()

    if args.full:
        vafs = [0, .005, .01, .02, .03, .04, .05, .06, .08, .10, .15, .20, .30, .50]
        depths = [250, 500, 1000, 2000, 4000]
        reps = args.reps or 30
        real_vafs = [0, .02, .03, .05, .08, .10, .20, .30, .50]
        real_reps = 5
    else:   # --quick (default): tractable, still multi-replicate for real CIs
        vafs = [0, .02, .03, .05, .08, .10, .20, .30]
        depths = [500, 1000, 2000]
        reps = args.reps or 8
        real_vafs = [0, .03, .05, .08, .10, .20]
        real_reps = 3
    if args.vafs:
        vafs = [float(x) for x in args.vafs.split(",")]
    if args.depths:
        depths = [int(x) for x in args.depths.split(",")]

    ref = mt.read_fasta_single(os.path.join(RDIR, "chrM.fa"))

    # --- enumerate the run plan, assert injective seeds, then execute ---
    plan = []   # (kind, variant, vaf, depth, bg, rep)
    for variant in ("del4977", "NONREP"):
        for depth in depths:
            for vaf in vafs:
                for rep in range(reps):
                    plan.append(("SIM", variant, vaf, depth, "sim", rep))
    for vaf in [v for v in vafs if v > 0]:             # HP_ARTIFACT hard-negative (one depth, spiked levels)
        for rep in range(reps):
            plan.append(("SIM", "HP_ARTIFACT", vaf, 2000, "sim", rep))
    for vaf in (0.10, 0.30):                            # ORIGIN suppression (two VAF)
        for rep in range(max(2, reps // 2)):
            plan.append(("SIM", "ORIGIN", vaf, 2000, "sim", rep))

    seeds = [seed_for(v, vf, d, r) for (k, v, vf, d, b, r) in plan]
    assert len(seeds) == len(set(seeds)), "seed collision — grid not injective"

    real_bams = [b for b in args.real_bams.split(",") if b]
    out = open(args.out, "a" if args.append else "w")
    if not args.append:
        out.write("\t".join(TSV_COLUMNS) + "\n")
    n_done = 0
    t_start = time.time()

    def emit(rec):
        out.write("\t".join(str(rec[c]) for c in TSV_COLUMNS) + "\n")
        out.flush()

    # SIMULATED arm
    if args.arm in ("both", "sim"):
        for (k, variant, vaf, depth, bg, rep) in plan:
            with tempfile.TemporaryDirectory() as work:
                rec = sim_run(ref, variant, vaf, depth, rep, work)
            emit(rec)
            n_done += 1
            sys.stderr.write("[lod] SIM %-11s vaf=%.3f depth=%4d rep=%d -> det=%s pass=%s svconf=%s (%d/%d)\n"
                             % (variant, vaf, depth, rep, rec["detected"], rec["passed"], rec["svconf"],
                                n_done, len(plan)))

    # REAL arm: call EVERY healthy background directly as a WT negative (cheap); SPIKE del4977 into a
    # subset (all backgrounds in --full as biological replicates; the first in --quick for tractability).
    spike_bgs = real_bams if args.full else real_bams[:1]
    for bg in (real_bams if args.arm in ("both", "real") else []):
        rb = os.path.join(HERE, "real", bg + ".chrM.bam")
        if not os.path.exists(rb):
            sys.stderr.write("[lod] REAL background %s absent — skipping\n" % bg)
            continue
        with tempfile.TemporaryDirectory() as work:
            emit(real_wt_run(rb, bg, work))                       # direct real-WT negative
        if bg not in spike_bgs:
            continue
        with tempfile.TemporaryDirectory() as work:
            wt1, wt2, wp = extract_wt_fastq(rb, work)
            for vaf in real_vafs:
                for rep in range(real_reps):
                    with tempfile.TemporaryDirectory() as w2:
                        sh("cp '%s' '%s/wt1.fq' && cp '%s' '%s/wt2.fq'" % (wt1, w2, wt2, w2))
                        rec = real_run(ref, "del4977", vaf, bg,
                                       os.path.join(w2, "wt1.fq"), os.path.join(w2, "wt2.fq"),
                                       wp, rep, w2)
                    emit(rec)
                    sys.stderr.write("[lod] REAL del4977 bg=%s vaf=%.3f rep=%d -> det=%s pass=%s svconf=%s\n"
                                     % (bg, vaf, rep, rec["detected"], rec["passed"], rec["svconf"]))

    out.close()
    sys.stderr.write("[lod] wrote %s (%d sim + real runs, %.0fs)\n"
                     % (args.out, n_done, time.time() - t_start))


if __name__ == "__main__":
    main()
