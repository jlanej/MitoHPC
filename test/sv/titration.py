#!/usr/bin/env python3
"""
Heteroplasmy x depth titration for the MitoHPC SV caller — measures SENSITIVITY (detection &
PASS rate vs heteroplasmy and depth), the limit of detection (LoD), and AF accuracy, and exercises
the two PASS paths (DJ = depth+junction, J = junction-strong). Not part of run_test.py (slow:
simulate -> realign -> call per grid point); run manually to calibrate the HP_SV_* thresholds.

  python3 titration.py [--out titration.tsv] [--hets 0.01,...] [--depths 500,...] [--del DEL|NONREP]

For each (het, depth) it simulates reads from a WT+DEL circular-genome mixture, realigns through the
pipeline's circular path (minimap2 -> circSam.pl), runs callSV.sh, and records whether the truth
deletion is detected / PASS, by which evidence path, and AFC/AFJ vs the spiked het.
Needs: samtools, minimap2, python3+pysam (HP_PYTHON). Reuses make_testdata.py's simulator.
"""
import argparse
import os
import random
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SDIR = os.environ.get("HP_SDIR", os.path.join(ROOT, "scripts"))
RDIR = os.environ.get("HP_RDIR", os.path.join(ROOT, "RefSeq"))
PY = os.environ.get("HP_PYTHON", "python3")
sys.path.insert(0, HERE)
from make_testdata import read_fasta_single, make_deletion, emit_reads   # noqa: E402

DELS = {                       # name -> (bp5, bp3): retained breakpoints (SVLEN = bp3-bp5-1)
    "DEL": (8469, 13447),      # del4977 common deletion (13 bp direct repeat -> Class I)
    "NONREP": (5999, 10999),   # non-repeat ~5 kb deletion (Class III)
}


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)


def simulate(ref_seq, bp5, bp3, het, depth, rlen, work, seed):
    """Write WT+DEL mixture FASTQ at the given heteroplasmy/outside-depth (make_testdata model)."""
    rng = random.Random(seed)
    glen = len(ref_seq)
    eg = make_deletion(ref_seq, bp5, bp3)
    f1 = open(os.path.join(work, "r1.fq"), "w")
    f2 = open(os.path.join(work, "r2.fq"), "w")
    n_wt = round(depth * (1 - het) * glen / rlen)
    n_e = round(depth * het * len(eg) / rlen)
    emit_reads(f1, f2, ref_seq + ref_seq, n_wt, rlen, 300, 450, 0.001, rng, "wt")
    emit_reads(f1, f2, eg + eg, n_e, rlen, 300, 450, 0.001, rng, "del")
    f1.close()
    f2.close()


def realign(work, sample):
    """minimap2 -ax sr chrMC -> -F 0x90C -> circSam.pl -> sorted $O.bam (the caller's input)."""
    bam = os.path.join(work, sample + ".bam")
    cmd = ("minimap2 -ax sr '%s/chrMC.fa' '%s/r1.fq' '%s/r2.fq' 2>/dev/null "
           "| samtools view -h -F 0x90C - "
           "| '%s/circSam.pl' -ref_len '%s/chrM.fa.fai' -offset 0 "
           "| samtools sort -o '%s' --write-index - 2>/dev/null") % (RDIR, work, work, SDIR, RDIR, bam)
    sh(cmd)
    return bam


def call_and_parse(sample, bam, work, bp5, bp3, bp_tol=30):
    env = dict(os.environ, HP_SDIR=SDIR, HP_RDIR=RDIR, HP_PYTHON=PY)
    pref = os.path.join(work, sample)
    subprocess.run(["bash", os.path.join(SDIR, "callSV.sh"), sample, bam, pref],
                   env=env, capture_output=True, text=True)
    tab = pref + ".sv.tab"
    if not os.path.exists(tab):
        return None
    with open(tab) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if len(lines) < 2:
        return None
    hdr = lines[0].lstrip("#").split("\t")
    best = None
    for ln in lines[1:]:
        r = dict(zip(hdr, ln.split("\t")))
        d = abs(int(r["pos_bp5"]) - bp5) + abs(int(r["end_bp3"]) - (bp3 - 1))
        if d <= bp_tol and (best is None or d < best[0]):
            best = (d, r)
    return best[1] if best else None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=os.path.join(HERE, "titration.tsv"))
    ap.add_argument("--hets", default="0.01,0.02,0.03,0.05,0.10,0.20,0.30,0.50")
    ap.add_argument("--depths", default="500,1000,2000,4000")
    ap.add_argument("--del", dest="delname", default="DEL", choices=list(DELS))
    ap.add_argument("--rlen", type=int, default=150)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    bp5, bp3 = DELS[args.delname]
    ref = read_fasta_single(os.path.join(RDIR, "chrM.fa"))
    hets = [float(x) for x in args.hets.split(",")]
    depths = [int(x) for x in args.depths.split(",")]

    rows = []
    out = open(args.out, "w")
    out.write("del\thet\tdepth\tdetected\tfilter\tsvclaim\tafc\tafj\tjr\tsr\n")
    for depth in depths:
        for het in hets:
            with tempfile.TemporaryDirectory() as work:
                simulate(ref, bp5, bp3, het, depth, args.rlen, work,
                         args.seed + int(het * 1000) + depth)
                bam = realign(work, "t")
                r = call_and_parse("t", bam, work, bp5, bp3)
            if r is None:
                det, flt, svc, afc, afj, jr, sr = 0, "-", "-", "-", "-", "-", "-"
            else:
                det, flt, svc = 1, r["filter"], r["svclaim"]
                afc, afj, jr, sr = r["af_coverage"], r["af_junction"], r["jr"], r["sr"]
            out.write("%s\t%.3f\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n"
                      % (args.delname, het, depth, det, flt, svc, afc, afj, jr, sr))
            rows.append((het, depth, det, flt, svc, afc, afj))
            sys.stderr.write("[titration] %s het=%.2f depth=%d -> det=%d %s %s afc=%s afj=%s\n"
                             % (args.delname, het, depth, det, flt, svc, afc, afj))
    out.close()

    # summary grid: PASS path per (het, depth) + per-depth LoD (lowest het reaching PASS)
    print("\n=== %s: PASS grid (DJ / J / . =non-PASS / x=undetected) ===" % args.delname)
    print("het\\depth  " + "  ".join("%6d" % d for d in depths))
    for het in hets:
        cells = []
        for depth in depths:
            r = next(x for x in rows if x[0] == het and x[1] == depth)
            cells.append("x" if not r[2] else (r[4] if r[3] == "PASS" else "."))
        print("%-9.2f  " % het + "  ".join("%6s" % c for c in cells))
    print("\nPASS LoD (lowest het reaching PASS) per depth:")
    for depth in depths:
        passing = [h for (h, d, det, flt, svc, afc, afj) in rows if d == depth and flt == "PASS"]
        print("  depth %5d : %s" % (depth, ("%.2f" % min(passing)) if passing else "none"))
    print("\nwrote %s" % args.out)


if __name__ == "__main__":
    main()
