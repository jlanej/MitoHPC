#!/usr/bin/env python3
"""
MitoHPC mtDNA structural-variant (large deletion) caller.

Single self-contained module (pysam, stdlib only) that replaces the v1 perl cores
(sa2del.pl + svCall.pl). It reads the circular-aware chrM alignment ($O.bam) in
process and writes ONLY new files:
    --out  $O.sv.vcf   (VCFv4.2, SVTYPE=DEL; header from scripts/sv.vcf)
    --tab  $O.sv.tab   (flat table)

Method (v1): cluster split-read deletion junctions (from SA:Z tags), require a
corroborating coverage drop for PASS, and report two heteroplasmy estimates
(junction fraction AFJ + coverage ratio AFC) with an AFDIFF QC field. See
docs/SV_METHODS.md (kept in sync with this code).

Driven by scripts/callSV.sh; all thresholds come from HP_SV_* env vars there.
"""
import argparse
import gzip
import sys
from collections import Counter

import pysam

REF_CONSUMING = frozenset("MDN=X")  # CIGAR ops that advance the reference


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def ref_len_from_cigar(cigar):
    """Reference bases consumed by a CIGAR string (M/D/N/=/X)."""
    n = 0
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            if ch in REF_CONSUMING and num:
                n += int(num)
            num = ""
    return n


def mode(values):
    """Most frequent value; ties broken by the smallest value (matches sa2del.pl)."""
    c = Counter(values)
    top = max(c.values())
    return min(v for v, k in c.items() if k == top)


def median(vals):
    if not vals:
        return 0
    s = sorted(vals)
    n = len(s)
    if n % 2:
        return s[n // 2]
    return (s[n // 2 - 1] + s[n // 2]) / 2


def rnd(x):
    """Round half up for non-negative values (matches perl int(x+0.5))."""
    return int(x + 0.5)


def load_bed_gz(path):
    """BED(.gz) -> list of [start1, end1] 1-based inclusive intervals."""
    iv = []
    if not path:
        return iv
    try:
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split()
                if len(f) >= 3 and f[1].isdigit():
                    iv.append((int(f[1]) + 1, int(f[2])))
    except OSError:
        pass
    return iv


def load_vcf_pos(path):
    """VCF(.gz) -> set of POS (1-based), any contig."""
    s = set()
    if not path:
        return s
    try:
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split()
                if len(f) >= 2 and f[1].isdigit():
                    s.add(int(f[1]))
    except OSError:
        pass
    return s


def in_iv(iv, p):
    return any(a <= p <= b for a, b in iv)


def near(p, a, b, pad):
    return (a - pad) <= p <= (b + pad)


# --------------------------------------------------------------------------- #
# Stage A: split-read deletion junctions (replaces sa2del.pl)
# --------------------------------------------------------------------------- #
def extract_junctions(bam, chrom, minmapq, minsize, maxsize, pad, minsupport):
    """Return clustered junctions: list of (bp5, bp3, svlen, JR, strand)."""
    pts = []  # (up, dn, rid, strand)
    with pysam.AlignmentFile(bam, "rb") as af:
        for r in af.fetch(chrom):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            if r.mapping_quality < minmapq:
                continue
            if not r.has_tag("SA"):
                continue
            strand = "-" if r.is_reverse else "+"
            abeg = r.reference_start + 1          # 1-based start
            aend = r.reference_end                # 1-based inclusive end
            if aend is None:
                continue

            # first SA segment only (by design, matching the original sa2del.pl): a read
            # with multiple supplementary alignments contributes just its first junction.
            sa = r.get_tag("SA").split(";")[0]
            s = sa.split(",")
            if len(s) < 4:
                continue
            sref, spos, sstrand, scig = s[0], s[1], s[2], s[3]
            if sref != chrom or sstrand != strand:   # same contig + same strand => DEL
                continue
            if not spos.isdigit() or scig == "*" or not scig:   # skip degenerate SA (cf. sa2del.pl)
                continue
            sbeg = int(spos)
            send = sbeg + ref_len_from_cigar(scig) - 1

            if abeg <= sbeg:
                up, dn = aend, sbeg
            else:
                up, dn = send, abeg
            svlen = dn - up - 1
            if svlen < minsize or svlen > maxsize:
                continue

            rid = r.query_name + ("/1" if r.is_read1 else "") + ("/2" if r.is_read2 else "")
            pts.append((up, dn, rid, strand))

    # greedy single-linkage clustering to the cluster seed (matches sa2del.pl)
    pts.sort(key=lambda x: (x[0], x[1]))
    clusters = []
    for p in pts:
        if clusters and abs(p[0] - clusters[-1]["su"]) <= pad and abs(p[1] - clusters[-1]["sd"]) <= pad:
            clusters[-1]["pts"].append(p)
        else:
            clusters.append({"su": p[0], "sd": p[1], "pts": [p]})

    out = []
    for c in clusters:
        ups = [p[0] for p in c["pts"]]
        dns = [p[1] for p in c["pts"]]
        strands = [p[3] for p in c["pts"]]
        jr = len({p[2] for p in c["pts"]})
        if jr < minsupport:
            continue
        bp5, bp3, strand = mode(ups), mode(dns), mode(strands)
        svlen = bp3 - bp5 - 1
        if svlen < minsize or svlen > maxsize:
            continue
        out.append((bp5, bp3, svlen, jr, strand))
    return out


# --------------------------------------------------------------------------- #
# per-base read depth (pysam count_coverage)
# --------------------------------------------------------------------------- #
def per_base_depth(bam, chrom, mtlen):
    """Per-base read depth over `chrom`, as a 1-based array of length mtlen+1.

    count_coverage sums A/C/G/T base counts, so it excludes deletions/ref-skips and (with
    read_callback="all") skips unmapped/secondary/qcfail/dup reads; quality_threshold=0
    counts all bases regardless of base quality. On the deduplicated, primary chrM $O.bam
    this equals `samtools depth -a` (verified field-for-field on the mock BAMs).
    """
    dep = [0] * (mtlen + 1)  # 1-based
    with pysam.AlignmentFile(bam, "rb") as af:
        if chrom not in af.references:
            raise ValueError("contig %r not found in %s" % (chrom, bam))
        reflen = af.get_reference_length(chrom)
        if reflen != mtlen:
            raise ValueError("--mtlen %d != %s length %d in %s" % (mtlen, chrom, reflen, bam))
        a, c, g, t = af.count_coverage(chrom, 0, mtlen, quality_threshold=0)
        for i in range(mtlen):
            dep[i + 1] = a[i] + c[i] + g[i] + t[i]
    return dep


# --------------------------------------------------------------------------- #
# Stage B: coverage corroboration, heteroplasmy, flags (replaces svCall.pl)
# --------------------------------------------------------------------------- #
def call(args):
    maxsize = args.maxsize if args.maxsize else args.mtlen - 1
    m = args.mtlen

    def wrap(p):
        return ((p - 1) % m) + 1

    def med_range(dep, a, b):
        if b < a:
            return 0
        return median([dep[wrap(p)] for p in range(a, b + 1)])

    def med_flank(dep, bp5, bp3, flank):
        vals = [dep[wrap(p)] for p in range(bp5 - flank + 1, bp5 + 1)]
        vals += [dep[wrap(p)] for p in range(bp3, bp3 + flank)]
        return median(vals)

    junctions = extract_junctions(args.bam, args.chrom, args.minmapq,
                                  args.minsize, maxsize, args.pad, args.minsupport)
    dep = per_base_depth(args.bam, args.chrom, m)

    fa = pysam.FastaFile(args.ref)
    seq = fa.fetch(args.chrom)
    fa.close()

    hp = load_bed_gz(args.hp)
    dloop = load_bed_gz(args.dloop)
    numt = load_vcf_pos(args.numt)

    rep = (args.rep5a, args.rep5b, args.rep3a, args.rep3b)

    vcf_records = []
    tab_rows = []
    for (bp5, bp3, svlen, jr, strand) in junctions:
        med_in = med_range(dep, bp5 + 1, bp3 - 1)
        med_fl = med_flank(dep, bp5, bp3, args.flank)
        ratio = med_in / med_fl if med_fl > 0 else 1.0
        span = rnd((dep[wrap(bp5)] + dep[wrap(bp3)]) / 2.0)
        sr = max(span - jr, 0)
        afj = jr / (jr + sr) if (jr + sr) > 0 else 0.0
        afc = max(0.0, min(1.0, 1.0 - ratio))
        afdiff = abs(afj - afc)

        flags = []
        repeat = (near(bp5, rep[0], rep[1], args.pad) or near(bp3, rep[2], rep[3], args.pad)
                  or near(bp5, rep[2], rep[3], args.pad) or near(bp3, rep[0], rep[1], args.pad))
        wrapf = (bp5 <= args.originpad or bp5 >= m - args.originpad
                 or bp3 <= args.originpad or bp3 >= m - args.originpad)
        inhp = in_iv(hp, bp5) or in_iv(hp, bp3)
        indl = in_iv(dloop, bp5) or in_iv(dloop, bp3)
        innumt = (bp5 in numt) or (bp3 in numt)
        if repeat:
            flags.append("REPEAT")
        if wrapf:
            flags.append("WRAP")
        if inhp:
            flags.append("HP")
        if indl:
            flags.append("DLOOP")
        if innumt:
            flags.append("NUMT")

        fil = []
        if jr < args.minjr:
            fil.append("lowJR")
        if ratio > args.drop:
            fil.append("no_cvg_drop")
        if wrapf:
            fil.append("WRAP")
        if args.mindepth and med_fl < args.mindepth:
            fil.append("lowDP")
        flt = ";".join(fil) if fil else "PASS"

        refbase = seq[bp5 - 1].upper() if 1 <= bp5 <= len(seq) else "N"
        end = bp3 - 1
        info = ("SM=%s;SVTYPE=DEL;END=%d;SVLEN=%d;JR=%d;SR=%d;AFJ=%.3f;AFC=%.3f;AFDIFF=%.3f;CVGR=%.3f"
                % (args.sample, end, -svlen, jr, sr, afj, afc, afdiff, ratio))
        for fl in flags:
            info += ";" + fl
        vcf_records.append("%s\t%d\t.\t%s\t<DEL>\t.\t%s\t%s\tGT:DP:AF\t0/1:%d:%.3f"
                           % (args.chrom, bp5, refbase, flt, info, rnd(med_fl), afj))
        tab_rows.append("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%s\t%s"
                        % (args.sample, args.chrom, bp5, bp3, svlen, jr, sr, afj, afc,
                           afdiff, ratio, rnd(med_fl), flt, ",".join(flags) if flags else "."))

    write_vcf(args, vcf_records)
    if args.tab:
        with open(args.tab, "w") as t:
            t.write("#sample\tchrom\tbp5\tbp3\tsvlen\tJR\tSR\tAFJ\tAFC\tAFDIFF\tCVGR\tFLANKDP\tFILTER\tflags\n")
            for row in tab_rows:
                t.write(row + "\n")
    sys.stderr.write("[callsv] %s -> %s (%d records)\n" % (args.sample, args.out, len(vcf_records)))


def write_vcf(args, records):
    out = open(args.out, "w") if args.out else sys.stdout
    with open(args.header) as h:
        for line in h:
            if line.startswith("#CHROM"):
                out.write("##sample=%s\n" % args.sample)
            out.write(line)
    for rec in records:
        out.write(rec + "\n")
    if args.out:
        out.close()


# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True)
    p.add_argument("--ref", required=True, help="chrM FASTA (indexed)")
    p.add_argument("--header", required=True, help="VCF header template (scripts/sv.vcf)")
    p.add_argument("--sample", default="SAMPLE")
    p.add_argument("--out", help="output VCF path (default stdout)")
    p.add_argument("--tab", help="output flat-table path")
    p.add_argument("--chrom", default="chrM")
    p.add_argument("--mtlen", type=int, default=16569)
    p.add_argument("--minmapq", type=int, default=20)
    p.add_argument("--minjr", type=int, default=3)
    p.add_argument("--minsize", type=int, default=50)
    p.add_argument("--maxsize", type=int, default=0)
    p.add_argument("--pad", type=int, default=25)
    p.add_argument("--drop", type=float, default=0.9)
    p.add_argument("--flank", type=int, default=200)
    p.add_argument("--mindepth", type=int, default=0)
    p.add_argument("--hp")
    p.add_argument("--numt")
    p.add_argument("--dloop")
    # Fixed v1 constants (callSV.sh does not expose these as HP_SV_*; the canonical defaults
    # live here, used for standalone/test invocation): cluster floor, origin guard, and the
    # del4977 13bp direct-repeat windows (m.8470-8482 / 13447-13459).
    p.add_argument("--minsupport", type=int, default=2)
    p.add_argument("--originpad", type=int, default=20)
    p.add_argument("--rep5a", type=int, default=8470)
    p.add_argument("--rep5b", type=int, default=8482)
    p.add_argument("--rep3a", type=int, default=13447)
    p.add_argument("--rep3b", type=int, default=13459)
    args = p.parse_args()
    try:
        call(args)
    except (ValueError, OSError) as e:
        sys.exit("[callsv] ERROR: %s (bam=%s, ref=%s, chrom=%s)"
                 % (e, args.bam, args.ref, args.chrom))


if __name__ == "__main__":
    main()
