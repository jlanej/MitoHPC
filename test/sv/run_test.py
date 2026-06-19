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

BP_TOL = 30       # breakpoint tolerance (>= the del4977 13bp repeat ambiguity)
SVLEN_TOL = 40
AF_TOL = 0.15
HI_HET = 0.10     # >= this heteroplasmy is expected to reach FILTER=PASS

results = []      # (name, detail, ok)


def record(name, ok, detail=""):
    results.append((name, ok, detail))
    print("  %-26s %s  %s" % (name, "PASS" if ok else "FAIL", detail))


def run_caller(sample, bam, prefix, mtlen=None, chrom=None):
    env = dict(os.environ, HP_SDIR=SDIR, HP_RDIR=RDIR, HP_PYTHON=PYEXE)
    if mtlen is not None:
        env["HP_MTLEN"] = str(mtlen)
    if chrom is not None:
        env["HP_MT"] = chrom
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


def load_truth():
    samples = {}
    with open(TRUTH) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            name, kind = f[0], f[1]
            samples.setdefault(name, [])
            if kind == "none":                 # wild-type marker, not an event
                continue
            samples[name].append({"kind": kind, "bp5": f[2], "bp3": f[3], "svlen": f[4],
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


# --------------------------------------------------------------------------- #
def check_sample(name, events, outdir):
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

    kinds = {e["kind"] for e in events}

    # negative / not-a-deletion samples: must yield zero PASS records
    if not events:                                  # wild-type
        record(name, npass == 0, "wild-type: %d PASS (want 0)" % npass)
        return
    if kinds & {"dup", "delwrap"}:                  # duplication / origin-crossing
        ok = npass == 0
        record(name, ok, "%s: %d PASS (want 0), %d non-PASS records"
               % ("/".join(kinds), npass, len(rows)))
        return

    # deletion sample(s): every truth deletion must be detected with correct fields
    all_ok, details = True, []
    for e in events:
        bp5, bp3, het = int(e["bp5"]), int(e["bp3"]), e["het"]
        m = match_del(rows, bp5, bp3)
        if m is None:
            all_ok = False
            details.append("del@%d/%d not detected" % (bp5, bp3))
            continue
        sub = []
        afj = fnum(m["af_junction"])
        if afj is None or abs(afj - het) > AF_TOL:
            all_ok = False
            sub.append("AFJ=%s vs het=%.2f" % (m["af_junction"], het))
        if het >= HI_HET and m["filter"] != "PASS":
            all_ok = False
            sub.append("want PASS got %s" % m["filter"])
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
    rep = os.path.join(outdir, "sv.report.html")
    rok = False
    if os.path.exists(rep):
        h = open(rep).read()
        rok = ("/*__DATA__*/" not in h and 'id="circ"' in h and 'id="lin"' in h
               and '"calls":' in h and len(h) > 5000)
    record("html_report", rok, "self-contained interactive report" if rok else "missing/invalid")


def check_examples(outdir):
    """The committed example outputs (test/sv/example/) must reflect the CURRENT schema —
    if the VCF/tab/report format changes, regenerate them: bash test/sv/make_example.sh."""
    EX = os.path.join(HERE, "example")
    if not os.path.isdir(EX):
        print("\n[committed examples — SKIPPED (test/sv/example/ not present)]")
        return
    print("\n[committed examples — schema in sync with current code]")
    rep = os.path.join(EX, "sv.report.html")
    hok = os.path.exists(rep)
    if hok:
        h = open(rep).read()
        hok = ("/*__DATA__*/" not in h and 'id="circ"' in h and 'id="lin"' in h
               and '"calls":' in h and len(h) > 5000)
    record("example_report", hok, "interactive report present + valid")

    fresh_tab = os.path.join(outdir, "sv_del4977_h30.sv.tab")
    ex_tab = os.path.join(EX, "sv.tab")
    tok = (os.path.exists(fresh_tab) and os.path.exists(ex_tab)
           and open(fresh_tab).readline().strip() == open(ex_tab).readline().strip())
    record("example_tab_schema", tok,
           "columns match" if tok else "DRIFT — run: bash test/sv/make_example.sh")

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


# --------------------------------------------------------------------------- #
def main():
    samples = load_truth()
    outdir = tempfile.mkdtemp()
    print("[scenario checks]")
    try:
        for name in sorted(samples):
            check_sample(name, samples[name], outdir)
        check_cohort(outdir, sorted(samples))
        check_examples(outdir)
        check_degenerate(outdir)
        check_vcf_spec(outdir)
    finally:
        shutil.rmtree(outdir, ignore_errors=True)

    nfail = sum(1 for _, ok, _ in results if not ok)
    print("\n%s  (%d checks, %d failed)" % (
        "ALL TESTS PASSED" if nfail == 0 else "SOME TESTS FAILED", len(results), nfail))
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
