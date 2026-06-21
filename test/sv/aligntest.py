#!/usr/bin/env python3
"""
Alignment-sensitivity harness — does the ALIGNER (or its settings) stop the SV caller from SEEING
deletions via split reads?

A junction-spanning read only becomes usable evidence if the aligner emits it as a split read
(an `SA:Z:` supplementary alignment) or at least soft-clips it at the breakpoint (which the caller
then harvests, §4.4). The two knobs that govern this in `bwa mem` are `-T` (min score to report an
alignment, incl. the supplementary arm; default 30 → the shorter arm needs ~30 bp) and `-L` (clipping
penalty; default 5 → higher penalties force reads to align THROUGH a junction instead of clipping).
Production MitoHPC uses `bwa mem -Y` with both at default (scripts/filter.sh).

This harness simulates reads from a WT+DEL mixture at several heteroplasmies, aligns each set under
several aligner configurations, runs callSV, and reports per config: the deletion's JR / AFJ / AFC /
breakpoint accuracy / PASS, AND the raw junction-read accounting (SA-tagged vs soft-clip-only at the
breakpoint) so you can see how many junction reads each setting EXPOSES. Not in run_test.py (slow;
needs bwa + minimap2). Reuses make_testdata.py.

Usage: HP_SDIR=../scripts HP_RDIR=../RefSeq python3 aligntest.py [--del DEL|NONREP] [--hets ...] [--depth N]
"""
import argparse
import importlib.util
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
spec = importlib.util.spec_from_file_location("mt", os.path.join(HERE, "make_testdata.py"))
mt = importlib.util.module_from_spec(spec); spec.loader.exec_module(mt)

DELS = {"DEL": (8469, 13447), "NONREP": (5999, 10999)}

# (label, argv prefix). The reference + r1 + r2 are appended. bwa needs an indexed ref (handled below).
ALIGNERS = [
    ("bwa_prod",    ["bwa", "mem", "-Y", "-v", "1"]),                 # == filter.sh (defaults: -T 30 -L 5)
    ("bwa_T20",     ["bwa", "mem", "-Y", "-v", "1", "-T", "20"]),     # lower supplementary-score threshold
    ("bwa_T15",     ["bwa", "mem", "-Y", "-v", "1", "-T", "15"]),     # even lower (shorter arms reported)
    ("bwa_L2",      ["bwa", "mem", "-Y", "-v", "1", "-L", "2"]),      # cheaper clipping (clip vs align-through)
    ("bwa_T15_L2",  ["bwa", "mem", "-Y", "-v", "1", "-T", "15", "-L", "2"]),
    ("minimap2_sr", ["minimap2", "-ax", "sr"]),                       # == gen_bams.sh / the mock path
]


def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def simulate(ref_seq, bp5, bp3, het, depth, work, seed):
    rng = random.Random(seed)
    glen = len(ref_seq)
    eg = mt.make_deletion(ref_seq, bp5, bp3)
    with open(f"{work}/r1.fq", "w") as f1, open(f"{work}/r2.fq", "w") as f2:
        mt.emit_reads(f1, f2, ref_seq + ref_seq, round(depth * (1 - het) * glen / 150), 150, 300, 450, 0.001, rng, "wt")
        mt.emit_reads(f1, f2, eg + eg, round(depth * het * len(eg) / 150), 150, 300, 450, 0.001, rng, "del")


def align_and_call(label, argv, ref_for_aligner, work, fai, bp5, bp3):
    """Align r1/r2 with `argv`, wrap to circular, call, and return metrics + raw junction accounting."""
    a = " ".join(argv) + f" '{ref_for_aligner}' '{work}/r1.fq' '{work}/r2.fq'"
    bam = f"{work}/{label}.bam"
    sh(f"{a} 2>/dev/null | samtools view -h -F 0x90C - 2>/dev/null "
       f"| '{SDIR}/circSam.pl' -ref_len '{fai}' -offset 0 2>/dev/null "
       f"| samtools sort -o '{bam}' --write-index - 2>/dev/null")
    # run the caller
    pref = f"{work}/{label}"
    subprocess.run(["bash", os.path.join(SDIR, "callSV.sh"), label, bam, pref],
                   env=dict(os.environ, HP_SDIR=SDIR, HP_RDIR=RDIR, HP_PYTHON=PY),
                   capture_output=True, text=True)
    row = _best(pref + ".sv.tab", bp5, bp3)
    # raw junction accounting at the CALLED bp5 (the 13bp repeat shifts it off the truth): how many
    # junction-spanning reads did the aligner expose as SA-tagged split reads vs soft-clip-only?
    cbp = int(row["pos_bp5"]) if row else bp5
    sa, clip = _junction_reads(bam, cbp)
    return sa, clip, row


def _junction_reads(bam, bp, pad=8, minclip=10):
    import pysam
    sa = clip = 0
    try:
        af = pysam.AlignmentFile(bam, "rb")
    except Exception:
        return 0, 0
    for r in af.fetch("chrM", max(0, bp - 50), bp + 50):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        cg = r.cigartuples or []
        clipped = ((cg and cg[-1][0] in (4, 5) and cg[-1][1] >= minclip and abs((r.reference_end or 0) - bp) <= pad)
                   or (cg and cg[0][0] in (4, 5) and cg[0][1] >= minclip and abs((r.reference_start + 1) - bp) <= pad))
        if not clipped:
            continue
        if r.has_tag("SA"):
            sa += 1
        else:
            clip += 1
    return sa, clip


def _best(tab, bp5, bp3, bptol=30):
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
        if d <= bptol and (best is None or d < best[0]):
            best = (d, r)
    return best[1] if best else None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--del", dest="delname", default="DEL", choices=list(DELS))
    ap.add_argument("--hets", default="0.05,0.10,0.30")
    ap.add_argument("--depth", type=int, default=3000)
    ap.add_argument("--out", default=os.path.join(HERE, "real", "aligntest.tsv"))
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    bp5, bp3 = DELS[args.delname]
    ref = mt.read_fasta_single(os.path.join(RDIR, "chrM.fa"))
    fai = os.path.join(RDIR, "chrM.fa.fai")
    hets = [float(x) for x in args.hets.split(",")]

    # bwa needs an indexed reference; index a private copy of chrMC so RefSeq stays clean
    idxdir = tempfile.mkdtemp()
    bwa_ref = os.path.join(idxdir, "chrMC.fa")
    sh(f"cp '{RDIR}/chrMC.fa' '{bwa_ref}' && bwa index '{bwa_ref}' 2>/dev/null")
    mm_ref = os.path.join(RDIR, "chrMC.fa")   # minimap2 indexes on the fly

    out = open(args.out, "w")
    out.write("del\thet\taligner\tdetected\tfilter\tsvclaim\tjsup\tjr\tafj\tafc\tbp_err\tsa_reads\tclip_only\n")
    rows = []
    for het in hets:
        with tempfile.TemporaryDirectory() as work:
            simulate(ref, bp5, bp3, het, args.depth, work, args.seed + int(het * 1000))
            for (label, argv) in ALIGNERS:
                ref_for = bwa_ref if argv[0] == "bwa" else mm_ref
                sa, clip, r = align_and_call(label, argv, ref_for, work, fai, bp5, bp3)
                if r is None:
                    det, flt, svc, jsup, jr, afj, afc, bperr = 0, "-", "-", "-", 0, "-", "-", "-"
                else:
                    det, flt, svc, jsup = 1, r["filter"], r["svclaim"], r.get("jsup", "-")
                    jr, afj, afc = r["jr"], r["af_junction"], r["af_coverage"]
                    bperr = abs(int(r["pos_bp5"]) - bp5)
                out.write(f"{args.delname}\t{het:.2f}\t{label}\t{det}\t{flt}\t{svc}\t{jsup}\t{jr}\t{afj}\t{afc}\t{bperr}\t{sa}\t{clip}\n")
                rows.append((het, label, det, flt, jr, afj, sa, clip))
                sys.stderr.write(f"[aligntest] {args.delname} het={het:.2f} {label:12s} det={det} jr={jr} "
                                 f"sa={sa} clip={clip} {flt}\n")
    out.close()
    sh(f"rm -rf '{idxdir}'")

    print(f"\n=== {args.delname}: JR (junction reads recovered) per aligner × het ===")
    labels = [a[0] for a in ALIGNERS]
    print("het\\aligner   " + "  ".join(f"{l:>12s}" for l in labels))
    for het in hets:
        cells = []
        for l in labels:
            r = next((x for x in rows if x[0] == het and x[1] == l), None)
            cells.append(f"{r[4]}" if r and r[2] else "x")
        print(f"{het:<12.2f}  " + "  ".join(f"{c:>12s}" for c in cells))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
