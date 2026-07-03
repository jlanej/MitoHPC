#!/usr/bin/env python3
"""
Evaluate the MitoHPC SV caller (scripts/callsv.py via scripts/callSV.sh) against the
committed mock BAMs in bams/ using the ground truth in truth.tsv, plus degenerate-input
robustness subtests and (when bcftools is available) a VCF-spec-compliance gate.

Self-contained: needs only python3 + pysam (the caller does everything in-process).
Run:  bash test/sv/run_test.sh   (wrapper)   or   python3 test/sv/run_test.py
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SDIR = os.path.join(ROOT, "scripts")
RDIR = os.path.join(ROOT, "RefSeq")
BAMS = os.path.join(HERE, "bams")
TRUTH = os.path.join(HERE, "truth.tsv")
PYEXE = os.environ.get("HP_PYTHON", sys.executable)   # interpreter that has pysam

MTLEN = 16569     # chrM length (the mocks are all chrM/16569)
BP_TOL = 30       # breakpoint tolerance (>= the del4977 13bp repeat ambiguity)
SVLEN_TOL = 40    # deletion-size tolerance (>= the repeat slide); enforced on PASS/detected del events
AF_TOL = 0.15     # heteroplasmy tolerance (AFC or AFJ may match)
# (PASS is required per-fixture via the truth.tsv `expect` column, not a heteroplasmy threshold.)

results = []      # (name, detail, ok)


def record(name, ok, detail=""):
    results.append((name, ok, detail))
    print("  %-26s %s  %s" % (name, "PASS" if ok else "FAIL", detail))


def run_caller(sample, bam, prefix, mtlen=16569, chrom="chrM"):
    # Pin contig + length so the suite is HERMETIC w.r.t. any ambient HP_MT/HP_MTLEN. The committed
    # BAMs are all chrM / 16569, but callSV.sh reads HP_MT to choose the contig AND its reference, so
    # a stray HP_MT inherited from the environment (e.g. a Docker build RUN that exported HP_MT=RSRS
    # for the reference-install loop) would make it fetch a non-existent contig and every call would
    # error — exactly the all-checks-fail signature that masquerades as a caller regression.
    env = dict(os.environ, HP_SDIR=SDIR, HP_RDIR=RDIR, HP_PYTHON=PYEXE,
               HP_MT=chrom, HP_MTLEN=str(mtlen))
    p = subprocess.run(["bash", os.path.join(SDIR, "callSV.sh"), sample, bam, prefix],
                       env=env, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def run_callsv_direct(bam, chrom, mtlen, prefix):
    """Invoke callsv.py directly (bypassing callSV.sh) to exercise its own input validation."""
    cmd = [PYEXE, os.path.join(SDIR, "callsv.py"), "--bam", bam,
           "--ref", os.path.join(RDIR, "chrM.fa"), "--header", os.path.join(SDIR, "sv.vcf"),
           "--sample", "x", "--out", prefix + ".vcf", "--tab", prefix + ".tab",
           "--chrom", chrom, "--mtlen", str(mtlen)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def read_tab(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if not lines:
        return []
    hdr = lines[0].lstrip("#").split("\t")
    return [dict(zip(hdr, ln.split("\t"))) for ln in lines[1:]]


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def valid_report(path):
    """A self-contained svReport HTML: data substituted (no placeholder), both views + calls present."""
    if not os.path.exists(path):
        return False
    h = open(path).read()
    return ("/*__DATA__*/" not in h and 'id="circ"' in h and 'id="lin"' in h
            and '"calls":' in h and len(h) > 5000)


def load_truth():
    """truth.tsv -> {name: {"expect": <str>, "events": [<event dict>...]}}. `expect` (col 8) is the
    behavior the CURRENT deletion-only caller should show today; forward-looking DUP/INV/complex
    fixtures carry no_pass/no_record/known_fp. See make_testdata.py / docs/SV_EVENT_TYPES.md."""
    samples = {}
    with open(TRUTH) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            name, kind = f[0], f[1]
            expect = f[7] if len(f) > 7 else "pass"
            s = samples.setdefault(name, {"expect": expect, "events": []})
            s["expect"] = expect
            if kind == "none":                 # wild-type marker, not an event
                continue
            s["events"].append({"kind": kind, "bp5": f[2], "bp3": f[3], "svlen": f[4],
                                "het": fnum(f[5]) or 0.0, "depth": int(f[6])})
    return samples


def match_del(rows, bp5, bp3):
    """Best called row whose breakpoints are within tolerance of (bp5, bp3)."""
    best, bestscore = None, None
    for r in rows:
        d5 = abs(int(r["pos_bp5"]) - bp5)
        d3 = abs(int(r["end_bp3"]) - (bp3 - 1))
        if d5 <= BP_TOL and d3 <= BP_TOL and (bestscore is None or d5 + d3 < bestscore):
            best, bestscore = r, d5 + d3
    return best


def is_del4977(bp5, bp3):
    return 8460 <= bp5 <= 8490 and 13440 <= bp3 <= 13460


SIG_BASELINE = (11000, 12000)   # coverage region untouched by any event (DUP-gain baseline)


def _meancov(bam, a, b):
    with pysam.AlignmentFile(bam, "rb") as f:
        c = f.count_coverage("chrM", a - 1, b, quality_threshold=0)
    return sum(c[0][i] + c[1][i] + c[2][i] + c[3][i] for i in range(b - a + 1)) / max(1, b - a + 1)


def _sa_first(r):
    return r.get_tag("SA").split(";")[0].split(",") if r.has_tag("SA") else None


def _opp_strand_sa(bam):                # the INVERSION signature: SA segment on the opposite strand
    n = 0
    with pysam.AlignmentFile(bam, "rb") as f:
        for r in f.fetch():
            if r.is_supplementary or r.is_secondary or r.is_unmapped:
                continue
            sa = _sa_first(r)
            if sa and len(sa) >= 3 and sa[2] != ("-" if r.is_reverse else "+"):
                n += 1
    return n


def _offorigin_sa(bam):                 # the origin-WRAP signature: a split read linking the two ends
    n = 0
    with pysam.AlignmentFile(bam, "rb") as f:
        for r in f.fetch():
            if r.is_supplementary:
                continue
            sa = _sa_first(r)
            if not sa or not sa[1].isdigit():
                continue
            p = int(sa[1])
            if (r.reference_start < 300 and p > 16000) or (r.reference_start > 16000 and p < 300):
                n += 1
    return n


def bam_signature_ok(bam, kind, bp5, bp3):
    """Independently confirm the mock BAM carries the SIGNAL its kind implies, BEFORE asserting the
    caller's behavior — so a `no_record`/`no_pass` assertion can't pass merely because the simulator
    silently produced nothing (e.g. an inversion that lost its opposite-strand reads). (ok, detail)."""
    if kind in ("inv", "invdup"):
        n = _opp_strand_sa(bam)
        return n > 0, "opposite-strand SA=%d" % n
    if kind in ("dup", "dupdel"):
        g = _meancov(bam, bp5 + 50, bp3 - 50) / max(1e-9, _meancov(bam, *SIG_BASELINE))
        return g > 1.15, "arc/baseline coverage=%.2f" % g
    if kind == "delwrap":
        n = _offorigin_sa(bam)
        return n > 0, "off-origin SA=%d" % n
    return True, ""   # del (incl. the sub-minsize CIGAR-D negative): no independent signal to assert


# --------------------------------------------------------------------------- #
def check_sample(name, info, outdir):
    expect, events = info["expect"], info["events"]
    prefix = os.path.join(outdir, name)
    bam = os.path.join(BAMS, name + ".bam")
    if not os.path.exists(bam):
        record(name, False, "missing BAM")
        return
    rc, log = run_caller(name, bam, prefix)
    if rc != 0 or "Traceback" in log:
        record(name, False, "caller failed rc=%d" % rc)
        return
    rows = read_tab(prefix + ".sv.tab")
    npass = sum(1 for r in rows if r["filter"] == "PASS")

    # forward-looking fixtures: first confirm the BAM actually carries its intended SV signal, so a
    # 0-record / 0-PASS result can't pass for the wrong reason (a silently broken simulator).
    if expect in ("no_record", "no_pass", "wrap", "known_fp") and events:
        e0 = events[0]
        sok, sdet = bam_signature_ok(bam, e0["kind"], int(e0["bp5"]), int(e0["bp3"]))
        if not sok:
            record(name, False, "simulator did not produce the %s signal (%s)" % (e0["kind"], sdet))
            return

    # --- expectation-driven outcomes for the non-deletion / forward-looking classes ---
    if expect == "no_record":            # INV (strand-filtered, CN-neutral) / sub-minsize del / fold-back
        record(name, len(rows) == 0, "want 0 records, got %d" % len(rows))
        return
    if expect == "no_pass":              # tandem DUP / controls: detected-but-not-PASS or nothing
        record(name, npass == 0, "want 0 PASS, got %d (%d records)" % (npass, len(rows)))
        return
    if expect == "wrap":                 # origin-crossing: 0 PASS + a WRAP-flagged MAJORITY-ARC record
        wraprec = [r for r in rows if "WRAP" in r.get("flags", "") and int(r["svlen"]) > MTLEN // 2]
        record(name, npass == 0 and len(wraprec) >= 1,
               "want 0 PASS + WRAP majority-arc record; %d PASS, wrap-recs=%d" % (npass, len(wraprec)))
        return
    if expect == "known_fp":             # KNOWN GAP: dup-del's EMBEDDED deletion spuriously PASSes today
        fp = match_del(rows, 6000, 6501)   # the internal del [6000..6500] (see make_testdata sv_dupdel)
        ok = fp is not None and fp["filter"] == "PASS"
        record(name, ok, "KNOWN-GAP: embedded del m.6000_6500 PASSes (to be fixed by the DUP-aware "
               "caller); matched_PASS=%s" % ok)
        return

    # --- expect in {pass, detected}: each truth deletion must be detected with correct fields ---
    require_pass = (expect == "pass")
    all_ok, details = True, []
    for e in events:
        bp5, bp3, het = int(e["bp5"]), int(e["bp3"]), e["het"]
        m = match_del(rows, bp5, bp3)
        if m is None:
            all_ok = False
            details.append("del@%d/%d not detected" % (bp5, bp3))
            continue
        sub = []
        # heteroplasmy: the caller reports a coverage-dosage AF (af_coverage, PRIMARY) and a
        # corrected junction VAF (af_junction). For an isolated deletion the dosage tracks truth;
        # for OVERLAPPING deletions (sv_multidel) the dosage is confounded by the other event's
        # overlap, but the junction estimate stays specific. Accept if EITHER is within tolerance.
        afc = fnum(m["af_coverage"])
        afj = fnum(m["af_junction"])
        errs = [abs(a - het) for a in (afc, afj) if a is not None]
        if not errs or min(errs) > AF_TOL:
            all_ok = False
            sub.append("AFC=%s/AFJ=%s vs het=%.2f" % (m["af_coverage"], m["af_junction"], het))
        if abs(int(m["svlen"]) - int(e["svlen"])) > SVLEN_TOL:
            all_ok = False
            sub.append("svlen=%s vs truth %s" % (m["svlen"], e["svlen"]))
        if require_pass and m["filter"] != "PASS":
            all_ok = False
            sub.append("want PASS got %s" % m["filter"])
        # split-read evidence lens: a true simulated deletion is a clean, consistent junction, so it
        # must be JSUP HIGH/MOD (not the LOW artifact tier) with high size-consistency.
        if m.get("jsup") == "LOW" or (fnum(m.get("srcons")) is not None and fnum(m["srcons"]) < 0.7):
            all_ok = False
            sub.append("JSUP=%s SRCONS=%s (want HIGH/MOD)" % (m.get("jsup"), m.get("srcons")))
        if is_del4977(bp5, bp3):
            if m["common"] != "1":
                all_ok = False; sub.append("COMMON!=1")
            if m["homlen"] != "13" or m["homseq"] != "ACCTCCCTCACCA":
                all_ok = False; sub.append("homlen/seq=%s/%s" % (m["homlen"], m["homseq"]))
            if m["delclass"] != "I":
                all_ok = False; sub.append("delclass=%s" % m["delclass"])
        else:
            if m["common"] == "1":
                all_ok = False; sub.append("unexpected COMMON")
        if int(m["ngene"] or 0) < 1:
            all_ok = False; sub.append("no genes annotated")
        if name == "sv_dloop" and "DLOOP" not in m["flags"]:
            all_ok = False; sub.append("missing DLOOP flag")
        details.append("del@%d/%d[%s]" % (bp5, bp3, ",".join(sub) if sub else "ok"))
    record(name, all_ok, "; ".join(details))


# --------------------------------------------------------------------------- #
def check_degenerate(outdir):
    print("\n[degenerate inputs — must fail cleanly, never a traceback]")
    src = os.path.join(BAMS, "sv_wt.bam")

    # (a) empty BAM (valid header, 0 reads) -> header-only VCF, 0 records, clean exit
    empty = os.path.join(outdir, "empty.bam")
    with pysam.AlignmentFile(src, "rb") as a:
        with pysam.AlignmentFile(empty, "wb", header=a.header):
            pass
    pysam.index(empty)
    rc, log = run_caller("empty", empty, os.path.join(outdir, "empty"))
    rows = read_tab(os.path.join(outdir, "empty.sv.tab"))
    record("empty_bam", rc == 0 and len(rows) == 0 and "Traceback" not in log,
           "rc=%d records=%d" % (rc, len(rows)))

    # (b) unindexed BAM -> clean ERROR, non-zero exit, no traceback
    noidx = os.path.join(outdir, "noindex.bam")
    shutil.copy(src, noidx)                       # copy BAM only, no .csi/.bai
    rc, log = run_caller("noindex", noidx, os.path.join(outdir, "noindex"))
    record("unindexed_bam", rc != 0 and "Traceback" not in log and "ERROR" in log,
           "rc=%d" % rc)

    # (c) wrong contig name -> clean ERROR (test callsv.py directly: callSV.sh would fail
    #     earlier on the missing per-contig reference, also cleanly)
    rc, log = run_callsv_direct(src, "NOPE", 16569, os.path.join(outdir, "badchrom"))
    record("wrong_contig", rc != 0 and "Traceback" not in log and "ERROR" in log, "rc=%d" % rc)

    # (d) mtlen mismatch -> clean ERROR
    rc, log = run_callsv_direct(src, "chrM", 99999, os.path.join(outdir, "badlen"))
    record("wrong_mtlen", rc != 0 and "Traceback" not in log and "ERROR" in log, "rc=%d" % rc)


def check_cohort(outdir, names):
    for tool in ("bcftools", "bgzip", "tabix"):
        if not shutil.which(tool):
            print("\n[cohort gate — SKIPPED (%s not on PATH)]" % tool)
            return
    print("\n[cohort aggregation — getSVSummary.sh: merge matrix + sites + recurrence]")
    intxt = os.path.join(outdir, "in.txt")
    with open(intxt, "w") as fh:
        for n in names:
            if os.path.exists(os.path.join(outdir, n + ".sv.vcf")):
                fh.write("%s\t%s\t%s\n" % (n, os.path.join(BAMS, n + ".bam"),
                                           os.path.join(outdir, n)))
    env = dict(os.environ, HP_SDIR=SDIR, HP_RDIR=RDIR, HP_IN=intxt, HP_ODIR=outdir)
    p = subprocess.run(["bash", os.path.join(SDIR, "getSVSummary.sh"), outdir],
                       env=env, capture_output=True, text=True)
    merged = os.path.join(outdir, "sv.merged.vcf.gz")
    sites = os.path.join(outdir, "sv.sites.vcf.gz")
    ok = p.returncode == 0 and os.path.exists(merged) and os.path.exists(sites)
    recur = False
    if ok:
        q = subprocess.run(["bcftools", "query", "-f", "%INFO/NS\n", merged],
                           capture_output=True, text=True)
        recur = any(int(x) >= 2 for x in q.stdout.split() if x.isdigit())  # del6000 shared
        v1 = subprocess.run(["bcftools", "view", merged], capture_output=True).returncode == 0
        v2 = subprocess.run(["bcftools", "view", sites], capture_output=True).returncode == 0
        ok = ok and v1 and v2
    record("cohort_getSVSummary", ok and recur,
           "merge+sites valid, recurrence(NS>=2)=%s" % recur)

    # the interactive HTML report (svReport.py; emitted by getSVSummary.sh)
    record("html_report", valid_report(os.path.join(outdir, "sv.report.html")),
           "self-contained interactive report")


def check_examples(outdir):
    """The committed example outputs (test/sv/example/) must reflect the CURRENT schema —
    if the VCF/tab/report format changes, regenerate them: bash test/sv/make_example.sh."""
    EX = os.path.join(HERE, "example")
    if not os.path.isdir(EX):
        print("\n[committed examples — SKIPPED (test/sv/example/ not present)]")
        return
    print("\n[committed examples — schema in sync with current code]")
    record("example_report", valid_report(os.path.join(EX, "sv.report.html")),
           "interactive report present + valid")

    # The committed example sv.tab is the COHORT table, whose header = the per-sample header plus the
    # cohort-only `mitobreak` column appended by getSVSummary.sh (see SV_TAB_DICTIONARY.md).
    fresh_tab = os.path.join(outdir, "sv_del4977_h30.sv.tab")
    ex_tab = os.path.join(EX, "sv.tab")
    tok = False
    if os.path.exists(fresh_tab) and os.path.exists(ex_tab):
        expected = open(fresh_tab).readline().strip() + "\tmitobreak"
        tok = open(ex_tab).readline().strip() == expected
    record("example_tab_schema", tok,
           "columns match (+ cohort mitobreak)" if tok else "DRIFT — run: bash test/sv/make_example.sh")

    fresh_vcf = os.path.join(outdir, "sv_del4977_h30.sv.vcf")
    ex_vcf = os.path.join(EX, "sv_del4977_h30.sv.vcf")
    vok = False
    if os.path.exists(fresh_vcf) and os.path.exists(ex_vcf):
        ids = lambda p: set(re.findall(r"##INFO=<ID=([^,]+)", open(p).read()))
        h = open(ex_vcf).read()
        chrom_ok = any(l.startswith("#CHROM") and l.rstrip().endswith("sv_del4977_h30")
                       for l in h.splitlines())
        vok = (ids(fresh_vcf) == ids(ex_vcf)
               and {"HOMLEN", "SVCLAIM", "GENE", "CIPOS", "DELCLASS"} <= ids(ex_vcf)
               and chrom_ok and "##contig=" in h)
    record("example_vcf_schema", vok,
           "INFO/FORMAT + sample column match" if vok else "DRIFT — run: bash test/sv/make_example.sh")


def check_vcf_spec(outdir):
    bcftools = shutil.which("bcftools")
    if not bcftools:
        print("\n[vcf spec gate — SKIPPED (bcftools not on PATH)]")
        return
    print("\n[vcf spec gate — bcftools view/sort, END<=contig]")
    ok_all = True
    for fn in sorted(os.listdir(outdir)):
        if not fn.endswith(".sv.vcf"):
            continue
        path = os.path.join(outdir, fn)
        p = subprocess.run([bcftools, "view", path], capture_output=True, text=True)
        bad = p.returncode != 0 or "not defined" in p.stderr.lower() or "undefined" in p.stderr.lower()
        if bad:
            ok_all = False
            print("    %s: %s" % (fn, p.stderr.strip().splitlines()[:1]))
    record("bcftools_spec", ok_all, "all per-sample VCFs valid")


def check_real(outdir):
    """Real-data vetting on committed real/*.chrM.bam (skipped if absent):
      * healthy 1000G high-coverage samples (NA*.chrM.bam) -> ZERO PASS (specificity guard the
        simulated mocks cannot give: real NUMT / D-loop / error structure);
      * a del4977 SPIKED INTO a real WT background (spike_del4977_h*.chrM.bam) -> the deletion is
        recovered (PASS + COMMON) with the right breakpoints and no off-target PASS (a realistic
        positive control: real error/coverage + known truth). See real/README.md, real/gen_spike.sh."""
    realdir = os.path.join(HERE, "real")
    if not os.path.isdir(realdir):
        return
    bams = sorted(f for f in os.listdir(realdir) if f.endswith(".chrM.bam"))
    if not bams:
        return
    print("[real-data vetting — healthy (0 PASS) + del4977-into-real-background (positive control)]")
    for fn in bams:
        s = fn[:-len(".chrM.bam")]
        rc, _ = run_caller(s, os.path.join(realdir, fn), os.path.join(outdir, "real_" + s))
        rows = read_tab(os.path.join(outdir, "real_" + s + ".sv.tab"))
        npass = sum(1 for r in rows if r.get("filter") == "PASS")
        if s.startswith("spike_del4977"):
            m = match_del(rows, 8469, 13447)               # the spiked common deletion
            ok = (rc == 0 and m is not None and m["filter"] == "PASS"
                  and m["common"] == "1" and npass == 1)   # recovered, COMMON, no off-target PASS
            record("real_%s" % s, ok, "spiked del4977 -> %s common=%s AFC=%s; %d PASS total (want 1)"
                   % (m["filter"] if m else "MISSING", m["common"] if m else "-",
                      m["af_coverage"] if m else "-", npass))
        else:
            record("real_%s" % s, rc == 0 and npass == 0,
                   "healthy real chrM: %d PASS (want 0), %d total calls" % (npass, len(rows)))


# --------------------------------------------------------------------------- #
def check_plots(outdir):
    """Optional: exercise the HP_SV_PLOT samplot visualization end-to-end (skipped if samplot absent).
    A visualizable call (del4977: PASS, high AFC, no HP/DLOOP/NUMT flag) -> a PNG + manifest that
    svReport embeds; a DLOOP-flagged call -> NOT visualized (the artifact filter)."""
    if not shutil.which("samplot"):
        print("\n[samplot visualization — SKIPPED (samplot not on PATH)]")
        return
    print("\n[samplot visualization — HP_SV_PLOT end-to-end]")
    env = dict(os.environ, HP_SDIR=SDIR, HP_RDIR=RDIR, HP_PYTHON=PYEXE, HP_MT="chrM",
               HP_MTLEN="16569", HP_SV_PLOT="1")

    def run_plot(name):
        pref = os.path.join(outdir, "plot_" + name)
        subprocess.run(["bash", os.path.join(SDIR, "callSV.sh"), name, os.path.join(BAMS, name + ".bam"), pref],
                       env=env, capture_output=True, text=True)
        man = pref + ".sv.plots.tsv"
        return pref, man, (os.path.exists(man) and os.path.getsize(man) > 0)

    pref, man, has = run_plot("sv_del4977_h30")
    png_ok = has and os.path.exists(open(man).readline().rstrip("\n").split("\t")[6])  # png is column 7
    rep = os.path.join(outdir, "plot_report.html")
    subprocess.run([PYEXE, os.path.join(SDIR, "svReport.py"), "--tab", pref + ".sv.tab",
                    "--plots", man, "--nsamples", "1", "--out", rep], capture_output=True, text=True)
    h = open(rep).read() if os.path.exists(rep) else ""
    record("samplot_plot", png_ok and '"png":' in h and 'id="plotsection"' in h,
           "del4977 -> PNG + manifest embedded in report")
    _, dman, dhas = run_plot("sv_dloop")
    record("samplot_filter", not dhas, "DLOOP-flagged call correctly NOT visualized")

    # HP_SV_PLOT_ALL: the DLOOP call (filtered out by default, above) IS visualized in all-plots mode,
    # the manifest carries the filter/flags columns, and the report marks plotAll (unfiltered gallery).
    apref = os.path.join(outdir, "plotall_sv_dloop")
    subprocess.run(["bash", os.path.join(SDIR, "callSV.sh"), "sv_dloop", os.path.join(BAMS, "sv_dloop.bam"), apref],
                   env=dict(env, HP_SV_PLOT_ALL="1"), capture_output=True, text=True)
    aman = apref + ".sv.plots.tsv"
    cols = open(aman).readline().rstrip("\n").split("\t") if os.path.exists(aman) else []
    arep = os.path.join(outdir, "plotall_report.html")
    subprocess.run([PYEXE, os.path.join(SDIR, "svReport.py"), "--tab", apref + ".sv.tab", "--plots", aman,
                    "--plot-dedup", "0", "--plot-all", "--nsamples", "1", "--out", arep],
                   capture_output=True, text=True)
    ah = open(arep).read() if os.path.exists(arep) else ""
    record("samplot_all_mode", len(cols) >= 9 and '"plotAll":true' in ah and 'status' in ah,
           "HP_SV_PLOT_ALL visualizes the DLOOP call; manifest carries filter/flags; report marks plotAll")

    # .filterpass naming (ALL mode only): a call that would ALSO pass the default filter gets
    # `.filterpass` embedded in its PNG name; one that would not (the DLOOP call above — skipped by the
    # default filter) stays plain. The common deletion clears the default filter, so it IS marked.
    dloop_png = open(aman).readline().split("\t")[6] if os.path.exists(aman) else ""
    dloop_unmarked = dloop_png.endswith(".png") and ".filterpass." not in dloop_png
    fpref = os.path.join(outdir, "plotall_sv_del4977_h30")
    subprocess.run(["bash", os.path.join(SDIR, "callSV.sh"), "sv_del4977_h30",
                    os.path.join(BAMS, "sv_del4977_h30.bam"), fpref],
                   env=dict(env, HP_SV_PLOT_ALL="1"), capture_output=True, text=True)
    fman = fpref + ".sv.plots.tsv"
    fp_png = open(fman).readline().split("\t")[6] if os.path.exists(fman) else ""
    record("samplot_filterpass_name", fp_png.endswith(".filterpass.png") and dloop_unmarked,
           "ALL mode: default-filter-passing call -> .filterpass.png; DLOOP call -> plain .png")


def check_dup_inv(outdir):
    """The OPT-IN tandem-DUP (HP_SV_DUP) and INV (HP_SV_INV) call paths. Re-run the forward-looking
    fixtures with the flags ON and assert the correct SVTYPE=DUP/INV calls; and prove the deletion
    path is unaffected by the flags (no spurious DUP/INV on a real deletion)."""
    print("\n[opt-in DUP/INV call paths — HP_SV_DUP / HP_SV_INV]")
    base = dict(os.environ, HP_SDIR=SDIR, HP_RDIR=RDIR, HP_PYTHON=PYEXE, HP_MT="chrM", HP_MTLEN="16569")

    def run(name, **extra):
        pref = os.path.join(outdir, "di_" + name)
        subprocess.run(["bash", os.path.join(SDIR, "callSV.sh"), name, os.path.join(BAMS, name + ".bam"), pref],
                       env=dict(base, **extra), capture_output=True, text=True)
        return read_tab(pref + ".sv.tab")

    def find(rows, svtype, bp5=None, bp3=None, passed=None, flag=None):
        for r in rows:
            if r["svtype"] != svtype:
                continue
            if bp5 is not None and abs(int(r["pos_bp5"]) - bp5) > BP_TOL:
                continue
            if bp3 is not None and abs(int(r["end_bp3"]) - bp3) > BP_TOL + 5:
                continue
            if passed is not None and (r["filter"] == "PASS") != passed:
                continue
            if flag is not None and flag not in r["flags"]:
                continue
            return r
        return None

    # tandem DUP -> SVTYPE=DUP PASS, AFC = the GAIN fraction (~ the dup heteroplasmy). Breakpoints are
    # approximate (the reverse-order junction comes from the shared deletion extractor), so we bound them
    # loosely to the duplicated arc rather than matching exactly.
    def near_arc(r, lo, hi, pad=400):
        return lo - pad <= int(r["pos_bp5"]) and int(r["end_bp3"]) <= hi + pad
    dup = [r for r in run("sv_dup", HP_SV_DUP="1")
           if r["svtype"] == "DUP" and r["filter"] == "PASS" and near_arc(r, 6000, 7000)]
    record("dup_tandem", len(dup) >= 1 and abs(fnum(dup[0]["af_coverage"]) - 0.5) <= AF_TOL,
           "sv_dup -> DUP PASS in arc, AFC=%s (want ~0.50)" % (dup[0]["af_coverage"] if dup else "none"))
    dupL = [r for r in run("sv_dup_large", HP_SV_DUP="1")
            if r["svtype"] == "DUP" and r["filter"] == "PASS" and near_arc(r, 4000, 9000)]
    record("dup_large", len(dupL) >= 1, "sv_dup_large -> DUP PASS in arc (%d)" % len(dupL))

    # balanced INV -> SVTYPE=INV PASS; low-het INV detected but withheld for the RIGHT reason; origin
    # INV WRAP-withheld (assert the specific FILTER token, not just "not PASS", so a wrong-reason fail
    # can't false-green).
    record("inv_balanced", find(run("sv_inv_small", HP_SV_INV="1"), "INV", 6000, 6500, passed=True) is not None,
           "sv_inv_small -> INV PASS")
    record("inv_large", find(run("sv_inv_large", HP_SV_INV="1"), "INV", 5000, 9000, passed=True) is not None,
           "sv_inv_large -> INV PASS")
    rl = find(run("sv_inv_lowhet", HP_SV_INV="1"), "INV", 8000, 9000)
    record("inv_lowhet", rl is not None and "lowAFJ" in rl["filter"],
           "sv_inv_lowhet -> INV withheld with lowAFJ (%s)" % (rl["filter"] if rl else "none"))
    ro = [r for r in run("sv_inv_origin", HP_SV_INV="1") if r["svtype"] == "INV"]
    record("inv_origin", len(ro) >= 1 and all("WRAP" in r["filter"] for r in ro),
           "sv_inv_origin -> INV record(s), all WRAP-withheld (%d)" % len(ro))

    # fold-back INVDUP -> flagged INVDUP, withheld with the not_balanced filter (not a balanced inversion)
    invdup = find(run("sv_invdup", HP_SV_INV="1"), "INV", flag="INVDUP")
    record("inv_dup_flag", invdup is not None and "not_balanced" in invdup["filter"],
           "sv_invdup -> INVDUP flag, not_balanced (%s)" % (invdup["filter"] if invdup else "none"))

    # FREEZE: the deletion VCF record is BYTE-IDENTICAL with the flags off — vs the committed baseline
    # (whose record predates the DUP/INV refactor), and unchanged when the flags are ON.
    pref = os.path.join(outdir, "di_frozen")
    subprocess.run(["bash", os.path.join(SDIR, "callSV.sh"), "sv_del4977_h30",
                    os.path.join(BAMS, "sv_del4977_h30.bam"), pref], env=base, capture_output=True, text=True)
    recs = lambda p: [l for l in open(p).read().splitlines() if l and not l.startswith("#")] if os.path.exists(p) else []
    fresh = recs(pref + ".sv.vcf")
    baseline = recs(os.path.join(HERE, "example", "sv_del4977_h30.sv.vcf"))
    record("del_frozen_vcf", len(fresh) >= 1 and fresh == baseline,
           "DEL VCF record byte-identical to the committed baseline (flags off)")
    rdel = run("sv_del4977_h30", HP_SV_DUP="1", HP_SV_INV="1")
    deld = find(rdel, "DEL", 8469, 13447, passed=True)
    record("del_unaffected_by_flags", deld is not None and not any(r["svtype"] in ("DUP", "INV") for r in rdel),
           "del4977 with DUP+INV on -> DEL PASS, no spurious DUP/INV (%d records)" % len(rdel))


def check_recall(outdir):
    """HP_SV_KEEPBAM persists the circular-aware BAM; re-calling on it (the fast re-run path, no
    realign) must reproduce the original calls byte-for-byte."""
    print("\n[SV re-run: keep-BAM persistence + recall fidelity]")
    bam = os.path.join(BAMS, "sv_del4977_h30.bam")
    env = dict(os.environ, HP_SDIR=SDIR, HP_RDIR=RDIR, HP_PYTHON=PYEXE, HP_MT="chrM", HP_MTLEN="16569")
    p1 = os.path.join(outdir, "keep")
    subprocess.run(["bash", os.path.join(SDIR, "callSV.sh"), "sv_del4977_h30", bam, p1],
                   env=dict(env, HP_SV_KEEPBAM="1"), capture_output=True, text=True)
    kept = p1 + ".sv.bam"
    persisted = os.path.exists(kept) and os.path.getsize(kept) > 0
    record("keepbam_persist", persisted, "HP_SV_KEEPBAM wrote $O.sv.bam")
    p2 = os.path.join(outdir, "recall")
    if persisted:
        subprocess.run(["bash", os.path.join(SDIR, "callSV.sh"), "sv_del4977_h30", kept, p2],
                       env=dict(env, HP_SV_KEEPBAM=""), capture_output=True, text=True)
        a, b = read_tab(p1 + ".sv.tab"), read_tab(p2 + ".sv.tab")
        keys = ("pos_bp5", "end_bp3", "svlen", "af_coverage", "svconf", "svimpact", "filter")
        same = (len(a) == len(b) and len(a) >= 1
                and all([r1[k] for k in keys] == [r2[k] for k in keys] for r1, r2 in zip(a, b)))
        record("recall_fidelity", same, "re-call on persisted BAM reproduces the call(s)")
    else:
        record("recall_fidelity", False, "no persisted BAM to re-call")


def check_mitobreak(outdir):
    """Cohort-only MitoBreak known-breakpoint annotation: the matcher (convention + tolerance), the
    sv.tab column, the merged-VCF MITOBREAK INFO, the data dictionary emitted alongside, the report
    column+filter, and that PER-SAMPLE outputs stay frozen (no MitoBreak fields). DB =
    RefSeq/mitobreak.tsv.gz; depends on check_cohort having run getSVSummary.sh into outdir."""
    db = os.path.join(RDIR, "mitobreak.tsv.gz")
    if not os.path.exists(db):
        print("\n[MitoBreak annotation — SKIPPED (RefSeq/mitobreak.tsv.gz absent)]")
        return
    print("\n[MitoBreak annotation — cohort-only previously-reported-breakpoint flagging]")
    import importlib.util
    spec = importlib.util.spec_from_file_location("svMitoBreak", os.path.join(SDIR, "svMitoBreak.py"))
    mb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mb)
    D = mb.load_db(db)

    # (a) matcher on the real DB: the common deletion matches; a clearly novel breakpoint does not
    common = mb.best_match(D, "DEL", 8482, 13446, 20)
    novel = mb.best_match(D, "DEL", 3000, 9000, 20)
    record("mitobreak_match", bool(common) and not novel,
           "common-del matches (%s); novel breakpoint -> no match" % (common or "-"))
    # (b) hermetic mini-DB: the +1 (DEL) convention and the exact tolerance boundary
    mini = {"DEL": [(1000, 2000, "DEL_1000_2000")], "DUP": []}
    m_in = mb.best_match(mini, "DEL", 1000, 1999, 20)   # bp3 1999+1=2000, bp5 exact -> match
    m_tol = mb.best_match(mini, "DEL", 1020, 1999, 20)  # bp5 off by 20 (== tol) -> match
    m_out = mb.best_match(mini, "DEL", 1021, 1999, 20)  # bp5 off by 21 (tol+1) -> NO match
    record("mitobreak_convention_tol",
           m_in == "DEL_1000_2000" and m_tol == "DEL_1000_2000" and m_out == "",
           "+1 bp3 convention + exact per-breakpoint tolerance boundary")

    # cohort-dependent checks below need getSVSummary.sh to have built the cohort sv.tab (bcftools/
    # bgzip/tabix present — see check_cohort). Skip cleanly otherwise; the matcher tests above stand.
    if not os.path.exists(os.path.join(outdir, "sv.tab")):
        print("  (cohort sv.tab not built — skipping cohort annotation checks)")
        return

    # (c) cohort sv.tab gained the column; del4977 rows annotated, a non-DB deletion is '.'
    tab = os.path.join(outdir, "sv.tab")
    hdr, rows = [], []
    if os.path.exists(tab):
        L = open(tab).read().splitlines()
        if L:
            hdr = L[0].split("\t")
            rows = [r.split("\t") for r in L[1:] if r]
    ix = {h: i for i, h in enumerate(hdr)}
    has_col = bool(hdr) and hdr[-1] == "mitobreak"
    del4977_ok = has_col and any(r[ix["mitobreak"]] not in ("", ".")
                                 for r in rows if len(r) > ix["mitobreak"] and r[ix["pos_bp5"]] == "8482")
    record("mitobreak_tab", has_col and del4977_ok, "sv.tab has mitobreak col; del4977 rows annotated")

    # (d) merged VCF INFO header + a tagged record; sites VCF must INHERIT the field (annotation runs
    #     before the sites derivation — guard against a future reorder that would silently drop it).
    merged = os.path.join(outdir, "sv.merged.vcf.gz")
    sites = os.path.join(outdir, "sv.sites.vcf.gz")
    hdrok = recok = sitesok = False
    if shutil.which("bcftools"):
        if os.path.exists(merged):
            h = subprocess.run(["bcftools", "view", "-h", merged], capture_output=True, text=True).stdout
            v = subprocess.run(["bcftools", "view", merged], capture_output=True, text=True).stdout
            hdrok = "ID=MITOBREAK" in h
            recok = any("MITOBREAK=" in ln for ln in v.splitlines() if not ln.startswith("#"))
        if os.path.exists(sites):
            sh = subprocess.run(["bcftools", "view", "-h", sites], capture_output=True, text=True).stdout
            sitesok = "ID=MITOBREAK" in sh
    record("mitobreak_vcf", hdrok and recok and sitesok,
           "merged VCF has MITOBREAK header+record; sites VCF inherits the field")

    # (e) data dictionary emitted alongside sv.tab, incl. the mitobreak row
    dictf = os.path.join(outdir, "sv.tab.dict.tsv")
    record("mitobreak_dict", os.path.exists(dictf) and "mitobreak" in open(dictf).read(),
           "sv.tab.dict.tsv emitted alongside sv.tab")

    # (f) per-sample outputs FROZEN — no MitoBreak field leaked into a sample's tab/vcf
    ps_tab = os.path.join(outdir, "sv_del4977_h30.sv.tab")
    ps_vcf = os.path.join(outdir, "sv_del4977_h30.sv.vcf")
    frozen = ((not os.path.exists(ps_tab) or "mitobreak" not in open(ps_tab).readline())
              and (not os.path.exists(ps_vcf) or "MITOBREAK" not in open(ps_vcf).read()))
    record("mitobreak_persample_frozen", frozen, "per-sample tab/vcf carry no MitoBreak annotation")

    # (g) report has the filter + the column
    rep = os.path.join(outdir, "sv.report.html")
    h = open(rep).read() if os.path.exists(rep) else ""
    record("mitobreak_report", 'id="fmb"' in h and ">MitoBreak<" in h,
           "report has MitoBreak-reported filter + column")


def check_depth_identity(outdir):
    """The per_base_depth fast path (get_blocks + N-base correction) must stay BYTE-IDENTICAL to the
    pysam count_coverage(quality_threshold=0) it replaced — on simulated mocks (no N) AND on a real
    BAM (N bases exercise the correction). Imports callsv directly and compares the two depth arrays."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("callsv", os.path.join(SDIR, "callsv.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    import pysam

    def ref_count_coverage(bam, chrom, m):
        dep = [0] * (m + 1)
        with pysam.AlignmentFile(bam, "rb") as af:
            a, c, g, t = af.count_coverage(chrom, 0, m, quality_threshold=0)
            for i in range(m):
                dep[i + 1] = a[i] + c[i] + g[i] + t[i]
        return dep

    cands = [os.path.join(BAMS, "sv_del4977_h30.bam")]
    rb = os.path.join(HERE, "real", "NA12718.chrM.bam")
    if os.path.isfile(rb):
        cands.append(rb)                      # real reads exercise the N-base correction
    for bam in cands:
        fast = mod.per_base_depth(bam, "chrM", 16569)
        ref = ref_count_coverage(bam, "chrM", 16569)
        ndiff = sum(1 for p in range(1, 16570) if fast[p] != ref[p])
        record("depth_identity:" + os.path.basename(bam), ndiff == 0,
               "per_base_depth == count_coverage (%d/16569 differ)" % ndiff)


def main():
    samples = load_truth()
    outdir = tempfile.mkdtemp()
    print("[scenario checks]")
    try:
        for name in sorted(samples):
            check_sample(name, samples[name], outdir)
        check_cohort(outdir, sorted(samples))
        check_mitobreak(outdir)
        check_examples(outdir)
        check_degenerate(outdir)
        check_vcf_spec(outdir)
        check_real(outdir)
        check_plots(outdir)
        check_dup_inv(outdir)
        check_depth_identity(outdir)
        check_recall(outdir)
    finally:
        shutil.rmtree(outdir, ignore_errors=True)

    nfail = sum(1 for _, ok, _ in results if not ok)
    print("\n%s  (%d checks, %d failed)" % (
        "ALL TESTS PASSED" if nfail == 0 else "SOME TESTS FAILED", len(results), nfail))
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
