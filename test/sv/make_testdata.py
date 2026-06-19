#!/usr/bin/env python3
"""
Realistic mock-read simulator for the MitoHPC structural-variant (SV) module.

It draws paired-end short reads from a mixture of WILD-TYPE and DELETED circular
mitochondrial genomes at a known heteroplasmy, so that aligning the reads back to
the wild-type circular reference (chrMC) reproduces BOTH signals a real deletion
produces:
  * split / soft-clipped reads with SA tags spanning the deletion junction
  * a coverage drop across the deleted span (only WT molecules cover it)

Ground truth (breakpoints + heteroplasmy) is written to truth.tsv so the caller
can be evaluated directly.

Heteroplasmy model (h = fraction of mtDNA molecules carrying the deletion):
  Outside the deletion both genomes contribute; inside, only WT contributes.
  We pick read counts so that, OUTSIDE the deletion, the DEL genome contributes a
  fraction h of depth and WT contributes (1-h). Then inside/outside depth ratio
  = (1-h), i.e. coverage-based AF = 1 - ratio = h, and the junction-read fraction
  also ~ h. See docs/SV_CALLING.md sec 3 & 9.2.

Pure stdlib; deterministic (fixed seed). Writes <out>/<sample>_1.fq, _2.fq.
"""
import argparse
import os
import random
import sys

RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(RC)[::-1]


def read_fasta_single(path: str) -> str:
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq).upper()


def make_deletion(seq: str, bp5: int, bp3: int) -> str:
    """Remove 1-based positions (bp5, bp3) exclusive of the breakpoints, i.e.
    delete bases bp5+1 .. bp3-1, joining base bp5 directly to base bp3.
    Returns the deleted linear genome. (bp5 = last retained left base,
    bp3 = first retained right base; SVLEN = bp3 - bp5 - 1.)"""
    # 1-based -> 0-based slice: keep [0, bp5) + [bp3-1, end)
    return seq[:bp5] + seq[bp3 - 1:]


def emit_reads(out1, out2, template_circular, n_frags, rlen, fmin, fmax,
               err, rng, name_prefix):
    """Sample n_frags paired-end fragments from a circular template (passed as a
    doubled linear string) and write FR-oriented reads to out1/out2."""
    glen = len(template_circular) // 2
    qual = "F" * rlen  # Q37
    bases = "ACGT"
    k = 0
    for _ in range(n_frags):
        flen = rng.randint(fmin, fmax)
        if flen < rlen:
            flen = rlen
        start = rng.randrange(glen)
        frag = template_circular[start:start + flen]
        if len(frag) < rlen:
            continue
        r1 = frag[:rlen]
        r2 = revcomp(frag[-rlen:])
        if err > 0:
            r1 = mutate(r1, err, rng, bases)
            r2 = mutate(r2, err, rng, bases)
        name = "%s:%d" % (name_prefix, k)
        out1.write("@%s/1\n%s\n+\n%s\n" % (name, r1, qual))
        out2.write("@%s/2\n%s\n+\n%s\n" % (name, r2, qual))
        k += 1


def mutate(s, err, rng, bases):
    out = []
    for c in s:
        if rng.random() < err:
            out.append(rng.choice(bases))
        else:
            out.append(c)
    return "".join(out)


# sample definitions: (name, bp5, bp3, heteroplasmy, outside_depth)
# bp5/bp3 are 1-based retained breakpoint bases; svlen = bp3-bp5-1.
# A heteroplasmy of 0 / bp5==0 means wild-type only (negative control).
SAMPLES = [
    ("sv_del4977_h30", 8469, 13447, 0.30, 300),  # canonical common deletion, 30%
    ("sv_del4977_h05", 8469, 13447, 0.05, 400),  # common deletion, 5% (sensitivity)
    ("sv_del6000_h50", 5999, 10999, 0.50, 300),  # non-repeat deletion, 50% (generality)
    ("sv_wt",          0,    0,     0.00, 300),  # wild-type only (specificity)
]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-ref", required=True, help="wild-type chrM FASTA (e.g. RefSeq/chrM.fa)")
    ap.add_argument("-out", required=True, help="output directory for FASTQ + truth.tsv")
    ap.add_argument("-rlen", type=int, default=150)
    ap.add_argument("-fmin", type=int, default=300)
    ap.add_argument("-fmax", type=int, default=450)
    ap.add_argument("-err", type=float, default=0.001, help="per-base substitution error")
    ap.add_argument("-seed", type=int, default=42)
    args = ap.parse_args()

    seq = read_fasta_single(args.ref)
    glen = len(seq)
    os.makedirs(args.out, exist_ok=True)
    truth = open(os.path.join(args.out, "truth.tsv"), "w")
    truth.write("#sample\tbp5\tbp3\tsvlen\theteroplasmy\toutside_depth\n")

    for (name, bp5, bp3, het, depth) in SAMPLES:
        rng = random.Random(args.seed + sum(ord(c) for c in name))
        wt2 = seq + seq  # circularized WT for wrap-around fragments
        f1 = open(os.path.join(args.out, name + "_1.fq"), "w")
        f2 = open(os.path.join(args.out, name + "_2.fq"), "w")

        if bp5 == 0:  # wild-type only
            n_wt = round(depth * glen / args.rlen)
            emit_reads(f1, f2, wt2, n_wt, args.rlen, args.fmin, args.fmax,
                       args.err, rng, name + "_wt")
            truth.write("%s\t.\t.\t.\t0\t%d\n" % (name, depth))
        else:
            del_seq = make_deletion(seq, bp5, bp3)
            dlen = len(del_seq)
            del2 = del_seq + del_seq
            svlen = bp3 - bp5 - 1
            # depth contributions: WT -> depth*(1-h), DEL -> depth*h (outside)
            n_wt = round(depth * (1 - het) * glen / args.rlen)
            n_del = round(depth * het * dlen / args.rlen)
            emit_reads(f1, f2, wt2, n_wt, args.rlen, args.fmin, args.fmax,
                       args.err, rng, name + "_wt")
            emit_reads(f1, f2, del2, n_del, args.rlen, args.fmin, args.fmax,
                       args.err, rng, name + "_del")
            truth.write("%s\t%d\t%d\t%d\t%.3f\t%d\n"
                        % (name, bp5, bp3, svlen, het, depth))
        f1.close()
        f2.close()
        sys.stderr.write("[make_testdata] %s: wrote FASTQ (het=%.2f)\n" % (name, het))

    truth.close()
    sys.stderr.write("[make_testdata] truth.tsv written to %s\n" % args.out)


if __name__ == "__main__":
    main()
