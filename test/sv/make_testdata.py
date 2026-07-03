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


def make_dup(seq, a, b):
    """Tandem duplication of 1-based [a, b]: ...[a..b][a..b]... (coverage GAIN, junction
    where reference order reverses). Used to confirm the caller does NOT PASS a DEL on it."""
    return seq[:b] + seq[a - 1:b] + seq[b:]


def make_delwrap(seq, bp5, bp3):
    """Origin-crossing deletion: retain the arc [bp3, bp5] (bp5 > bp3), delete the
    complementary arc that crosses the origin. The junction joins bp5 -> bp3."""
    return seq[bp3 - 1:bp5]


def make_inversion(seq, a, b):
    """Balanced inversion of 1-based [a, b]: the segment is REVERSE-COMPLEMENTED in place
    (revcomp, NOT just reversed — a plain reverse would still map to the forward strand and
    produce no opposite-strand SA, silently no-op'ing the test). Copy-number-NEUTRAL: the two
    breakpoints (a-1|a and b|b+1) become OPPOSITE-STRAND junctions — the signal the caller's
    same-strand SA filter currently discards (docs/SV_EVENT_TYPES.md §3.3)."""
    return seq[:a - 1] + revcomp(seq[a - 1:b]) + seq[b:]


def make_inv_dup(seq, a, b):
    """Fold-back inverted duplication of 1-based [a, b]: ...[a..b][revcomp(a..b)]... — an extra,
    INVERTED copy inserted right after the segment. Presents as an opposite-strand junction (the
    inverted arm) CO-LOCATED with a coverage GAIN over [a,b]; must NOT be read as a balanced
    inversion (which is CN-neutral). Sniffles2 INVDUP."""
    return seq[:b] + revcomp(seq[a - 1:b]) + seq[b:]


def make_dupdel(seq, a, b, da, db):
    """Compound duplication-with-internal-deletion (partial duplication / 'dup-del'): a tandem extra
    copy of [a, b] whose second copy is missing its internal [da, db]. One molecule yields BOTH a
    dup junction (b -> a) AND a del junction (da-1 -> db+1), the clinically dominant mtDNA dup class
    (KSS/Pearson). da,db are 1-based, a < da <= db < b."""
    second = seq[a - 1:da - 1] + seq[db:b]          # the duplicated copy, internal [da,db] deleted
    return seq[:b] + second + seq[b:]


# Sample definitions: (name, outside_depth, events, expect). An event is a tuple starting
# (kind, p1, p2, het); compound kinds carry extra params after het. Kinds:
#   del      p1=bp5, p2=bp3 retained breakpoints (deleted span p1+1..p2-1)
#   dup      p1=a, p2=b tandem-duplicated segment (coverage GAIN)
#   delwrap  p1=bp5, p2=bp3 (p1>p2): origin-crossing deletion
#   inv      p1=a, p2=b balanced inversion (revcomp in place; CN-neutral; OPPOSITE-strand junctions)
#   invdup   p1=a, p2=b fold-back inverted duplication (opposite-strand junction + gain)
#   dupdel   (kind,a,b,het,da,db): tandem dup of [a,b] whose 2nd copy lacks internal [da,db]
# Wild-type fraction per sample = 1 - sum(event hets). Empty events => pure wild-type.
# `expect` is what the CURRENT deletion-only caller should do TODAY (run_test asserts it):
#   pass        a PASS record matching the truth deletion (with the DEL field checks)
#   detected    detected/matched but PASS not required (low-het edge)
#   no_pass     >=0 records but ZERO PASS (dup/complex/control)
#   no_record   ZERO records emitted (sub-minsize del; INV — strand-filtered & CN-neutral)
#   wrap        zero PASS and a WRAP-flagged record (origin/majority-arc)
#   known_fp    a documented compound-event gap: a record PASSes today that the DUP-aware
#               caller will reject (used by sv_dupdel's spuriously-PASSing embedded deletion)
# Forward-looking DUP/INV/complex are mostly no_pass/no_record until those paths land
# (docs/SV_EVENT_TYPES.md); the kinds beyond del/dup/delwrap are NOT yet callable.
SAMPLES = [
    # --- deletions: the implemented, PASS-able class ---
    ("sv_del4977_h30", 300, [("del", 8469, 13447, 0.30)], "pass"),       # common deletion 30% (positive control)
    ("sv_del4977_h05", 400, [("del", 8469, 13447, 0.05)], "detected"),   # common deletion 5% (low-het floor)
    ("sv_del6000_h50", 300, [("del", 5999, 10999, 0.50)], "pass"),       # non-repeat deletion 50% (Class III)
    ("sv_multidel",    400, [("del", 8469, 13447, 0.25),                 # TWO concurrent deletions
                             ("del", 5999, 10999, 0.15)], "pass"),
    ("sv_homoplasmy",  300, [("del", 8469, 13447, 0.95)], "pass"),       # near-homoplasmic common deletion
    ("sv_dloop",       300, [("del", 400, 6000, 0.40)], "pass"),         # 5' breakpoint in the D-loop (DLOOP flag)
    ("sv_lowcov",      40,  [("del", 8469, 13447, 0.50)], "pass"),       # low coverage (cohort depth variability)
    ("sv_del_500",     300, [("del", 8000, 8501, 0.50)], "pass"),        # 500bp del -> small detectable PASS
    ("sv_del_45",      300, [("del", 9000, 9046, 0.50)], "no_record"),   # 45bp < minsize=50; also a CIGAR-D not a split
    ("sv_del_13kb",    300, [("del", 2000, 15001, 0.60)], "pass"),       # 13kb (>BIGDEL, majority-arc) PASS via dosage
    # --- duplications + complex: detected-but-not-PASS today (forward-looking) ---
    ("sv_dup",         300, [("dup", 6000, 7000, 0.50)], "no_pass"),     # tandem dup 1kb (must NOT PASS as DEL)
    ("sv_dup_large",   300, [("dup", 4000, 9000, 0.50)], "no_pass"),     # tandem dup 5kb (gain guard scales)
    ("sv_dupdel",      400, [("dupdel", 5000, 8000, 0.40, 6000, 6500)], "known_fp"),  # partial dup-del: embedded del spuriously PASSes today (known compound-event gap)
    ("sv_invdup",      400, [("invdup", 7000, 7400, 0.40)], "no_record"),# fold-back inverted dup
    # --- inversions: architecturally invisible today (strand-filtered, CN-neutral) ---
    ("sv_inv_small",   300, [("inv", 6000, 6500, 0.50)], "no_record"),   # 500bp balanced inversion
    ("sv_inv_large",   300, [("inv", 5000, 9000, 0.50)], "no_record"),   # 4kb inversion (invisibility is size-indep)
    ("sv_inv_origin",  400, [("inv", 16300, 16560, 0.40)], "no_record"), # near-origin inversion (strand + origin)
    ("sv_inv_lowhet",  300, [("inv", 8000, 9000, 0.05)], "no_record"),   # low-het inversion (future sensitivity floor)
    # --- origin-crossing deletions: WRAP-withheld; the resolution regression pair ---
    ("sv_origin",            400, [("delwrap", 16400, 200, 0.40)], "wrap"),  # clips OriH -> future DUP
    ("sv_del_origin_spares", 400, [("delwrap", 16400, 100, 0.40)], "wrap"),  # spares origins -> future DEL
    # --- controls ---
    ("sv_wt",          300, [], "no_pass"),                              # wild-type only (specificity + degenerate substrate)
]


def event_genome(seq, ev):
    kind, p1, p2 = ev[0], ev[1], ev[2]
    if kind == "del":           return make_deletion(seq, p1, p2)
    if kind == "dup":           return make_dup(seq, p1, p2)
    if kind == "delwrap":       return make_delwrap(seq, p1, p2)
    if kind == "inv":           return make_inversion(seq, p1, p2)
    if kind == "invdup":        return make_inv_dup(seq, p1, p2)
    if kind == "dupdel":        return make_dupdel(seq, p1, p2, ev[4], ev[5])
    raise ValueError("unknown event kind %r" % kind)


def event_svlen(ev, mtlen):
    """Truth svlen: a real deleted length for del/delwrap; for dup/inv/complex the AFFECTED
    (duplicated/inverted) segment length, recorded for reference (these are not yet called)."""
    kind, p1, p2 = ev[0], ev[1], ev[2]
    if kind == "del":
        return p2 - p1 - 1
    if kind == "delwrap":            # origin-crossing deleted arc length
        return (mtlen - p1) + (p2 - 1)
    return p2 - p1 + 1               # dup/inv/invdup/dupdel: affected segment span


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
    truth.write("#sample\tkind\tbp5\tbp3\tsvlen\thet\tdepth\texpect\n")

    wt2 = seq + seq  # circularized WT for wrap-around fragments
    for (name, depth, events, expect) in SAMPLES:
        rng = random.Random(args.seed + sum(ord(c) for c in name))
        f1 = open(os.path.join(args.out, name + "_1.fq"), "w")
        f2 = open(os.path.join(args.out, name + "_2.fq"), "w")

        hetsum = sum(e[3] for e in events)
        n_wt = round(depth * (1 - hetsum) * glen / args.rlen)
        emit_reads(f1, f2, wt2, n_wt, args.rlen, args.fmin, args.fmax, args.err, rng, name + "_wt")

        if not events:
            truth.write("%s\tnone\t.\t.\t.\t0\t%d\t%s\n" % (name, depth, expect))
        for i, ev in enumerate(events):
            kind, het = ev[0], ev[3]
            eg = event_genome(seq, ev)
            eg2 = eg + eg                # double the EVENT genome so fragments span its (shifted) origin
            n_e = round(depth * het * len(eg) / args.rlen)
            emit_reads(f1, f2, eg2, n_e, args.rlen, args.fmin, args.fmax,
                       args.err, rng, "%s_%s%d" % (name, kind, i))
            truth.write("%s\t%s\t%d\t%d\t%d\t%.3f\t%d\t%s\n"
                        % (name, kind, ev[1], ev[2], event_svlen(ev, glen), het, depth, expect))
        f1.close()
        f2.close()
        sys.stderr.write("[make_testdata] %s: %d event(s), depth=%d, expect=%s\n"
                         % (name, len(events), depth, expect))

    truth.close()
    sys.stderr.write("[make_testdata] truth.tsv written to %s\n" % args.out)


if __name__ == "__main__":
    main()
