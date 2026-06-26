#!/usr/bin/env python3
"""Annotate cohort SV summaries with previously-reported MitoBreak breakpoints — COHORT-ONLY.

Run once at consolidation (scripts/getSVSummary.sh). It NEVER touches per-sample outputs; it only
enriches the cohort sv.tab and the merged/sites VCFs:
  - tab mode:  appends a `mitobreak` column = the closest matching MitoBreak id within tolerance,
               or "." when none (so the column never has empty cells).
  - vcf mode:  inserts ##INFO=<ID=MITOBREAK,...> and appends MITOBREAK=<id> to matching records.

Coordinate convention (the subtle part): MitoBreak's 3' breakpoint is the FIRST RETAINED base after
the event (exclusive); MitoHPC's end_bp3 / VCF END is the LAST DELETED base (inclusive). For
deletions they differ by 1 (MitoBreak_bp3 == our_bp3 + 1); the 5' breakpoints share convention (last
retained base before the event), so they compare directly. The matcher converts our bp3 to MitoBreak's
convention before comparing. Duplications are matched directly (best-effort; DUP is opt-in and rarely
emitted — origin-wrapping DUP breakpoints are not yet wrap-corrected here).

Matching: same svtype AND |bp5 - mb_bp5| <= tol AND |bp3' - mb_bp3| <= tol. When several DB rows
match, the closest (minimum summed breakpoint distance) id is reported. Default tol = 20 bp, which
covers the ~13 bp direct-repeat breakpoint spread (e.g. the common deletion) plus the +1 offset while
keeping spurious matches well under 1%.
"""
import argparse
import gzip
import sys

INFO_HDR = ('##INFO=<ID=MITOBREAK,Number=.,Type=String,Description="Previously-reported breakpoint '
            'id(s) from the MitoBreak database (Damas 2014, PMC3965124) matched within '
            'HP_SV_MITOBREAK_TOL bp on both breakpoints, coordinate-convention reconciled; '
            'see RefSeq/mitobreak.tsv.gz">')


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_db(path):
    """RefSeq/mitobreak.tsv(.gz) -> {svtype: [(bp5, bp3, id), ...]}."""
    db = {"DEL": [], "DUP": []}
    with _open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                t = f[ix["svtype"]]
                b5 = int(f[ix["bp5"]])
                b3 = int(f[ix["bp3"]])
                mid = f[ix["mitobreak_id"]]
            except (KeyError, ValueError, IndexError):
                continue
            db.setdefault(t, []).append((b5, b3, mid))
    return db


def best_match(db, svtype, bp5, bp3, tol):
    """Closest MitoBreak id whose 5' and 3' breakpoints are both within `tol`, else "" (no match)."""
    rows = db.get(svtype)
    if not rows:
        return ""
    q3 = bp3 + 1 if svtype == "DEL" else bp3      # our last-deleted -> MitoBreak first-retained
    best = None
    for (m5, m3, mid) in rows:
        d5, d3 = abs(bp5 - m5), abs(q3 - m3)
        if d5 <= tol and d3 <= tol:
            d = d5 + d3
            if best is None or d < best[0]:
                best = (d, mid)
    return best[1] if best else ""


def _info_get(info, key):
    for kv in info.split(";"):
        if kv.startswith(key + "="):
            return kv[len(key) + 1:]
    return None


def annotate_tab(db, tol, src, dst):
    """Append a `mitobreak` column to a tidy sv.tab (matched-id or '.')."""
    with open(src) as fi, open(dst, "w") as fo:
        head = fi.readline().rstrip("\n")
        hdr = head.split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        fo.write(head + "\tmitobreak\n")
        have = all(k in ix for k in ("svtype", "pos_bp5", "end_bp3"))
        for line in fi:
            line = line.rstrip("\n")
            mid = ""
            if have and line:
                f = line.split("\t")
                try:
                    mid = best_match(db, f[ix["svtype"]], int(f[ix["pos_bp5"]]),
                                     int(f[ix["end_bp3"]]), tol)
                except (ValueError, IndexError):
                    mid = ""
            fo.write(line + "\t" + (mid or ".") + "\n")


def annotate_vcf(db, tol, src, dst):
    """Insert the MITOBREAK INFO header and tag matching records (POS=bp5, INFO END=bp3, SVTYPE)."""
    with _open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    fo.write(INFO_HDR + "\n")     # ##INFO is valid anywhere before #CHROM
                fo.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                fo.write(line)
                continue
            info = f[7]
            svtype = _info_get(info, "SVTYPE") or "DEL"
            end = _info_get(info, "END")
            mid = ""
            try:
                mid = best_match(db, svtype, int(f[1]), int(end), tol)
            except (TypeError, ValueError):
                mid = ""
            if mid and _info_get(info, "MITOBREAK") is None:   # idempotent: never double-tag
                f[7] = ("" if info in ("", ".") else info + ";") + "MITOBREAK=" + mid
                fo.write("\t".join(f) + "\n")
            else:
                fo.write(line)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True, help="RefSeq/mitobreak.tsv(.gz)")
    ap.add_argument("--tol", type=int, default=20, help="per-breakpoint match tolerance (bp)")
    ap.add_argument("--tab")
    ap.add_argument("--tab-out", dest="tab_out")
    ap.add_argument("--vcf")
    ap.add_argument("--vcf-out", dest="vcf_out")
    a = ap.parse_args(argv)
    db = load_db(a.db)
    if a.tab and a.tab_out:
        annotate_tab(db, a.tol, a.tab, a.tab_out)
    if a.vcf and a.vcf_out:
        annotate_vcf(db, a.tol, a.vcf, a.vcf_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
