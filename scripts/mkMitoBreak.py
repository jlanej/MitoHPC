#!/usr/bin/env python3
"""Normalize the committed MitoBreak database CSVs (resources/) into RefSeq/mitobreak.tsv.gz.

MitoBreak (http://mitobreak.portugene.com, Damas et al. 2014, NAR — PMC3965124) is a curated
database of previously reported mitochondrial DNA breakpoints. The portal export ships as two
quoted CSVs whose names contain spaces and an em-dash:
  - "MitoBreak — The mitochondrial DNA breakpoints database.csv"      (deletions)
  - "MitoBreak — The mitochondrial DNA breakpoints database 1.csv"    (duplications)

This one-off generator parses both (robustly, via csv) and emits a single clean, version-controlled,
tab-delimited table that `svMitoBreak.py` consumes at COHORT CONSOLIDATION to flag whether a call's
breakpoints were previously reported. Coordinates are 1-based rCRS (NC_012920.1, 16569 bp), copied
VERBATIM from MitoBreak (5'/3' breakpoint columns) — the caller-vs-MitoBreak coordinate convention
is reconciled in svMitoBreak.py, not here.

Usage:
  python3 scripts/mkMitoBreak.py resources/ | bgzip > RefSeq/mitobreak.tsv.gz
  # (provenance — source files, row counts — is printed to stderr)

Output columns (header included): svtype, bp5, bp3, length, location, origin_impact, disease, refs,
mitobreak_id   (see docs/SV_TAB_DICTIONARY.md).
"""
import csv
import glob
import os
import sys

# disease/tissue columns shared by both CSVs -> compact tokens (non-empty cell => token present)
DISEASE_COLS = [("Healthy tissues", "Healthy"), ("Parkinson Disease", "PD"),
                ("Inclusion Body Myositis", "IBM"), ("Tumour", "Tumour"),
                ("Other clinical features", "Other")]
HEADER = ["svtype", "bp5", "bp3", "length", "location", "origin_impact",
          "disease", "refs", "mitobreak_id"]


def _clean(s):
    return (s or "").strip().strip('"').strip()


def _disease_tokens(row):
    toks = [tok for col, tok in DISEASE_COLS if _clean(row.get(col))]
    return ";".join(toks) if toks else "."


def _emit(rows, svtype, len_col, loc_col, origin_col, out):
    n = 0
    for r in rows:
        b5, b3 = _clean(r.get("5' breakpoint")), _clean(r.get("3' breakpoint"))
        if not (b5.isdigit() and b3.isdigit()):
            continue
        length = _clean(r.get(len_col)) or "."
        location = _clean(r.get(loc_col)) or "."
        origin = _clean(r.get(origin_col)) or "."
        disease = _disease_tokens(r)
        refs = _clean(r.get("References")) or "."
        mbid = "%s_%s_%s" % (svtype, b5, b3)
        out.append([svtype, b5, b3, length, location, origin, disease, refs, mbid])
        n += 1
    return n


def main(argv):
    resdir = argv[1] if len(argv) > 1 else "resources"
    csvs = sorted(glob.glob(os.path.join(resdir, "*.csv")))
    # the deletions file is the big one; the "… 1.csv" sibling is duplications
    dels = [p for p in csvs if "1.csv" not in os.path.basename(p)]
    dups = [p for p in csvs if "1.csv" in os.path.basename(p)]

    out = []
    nd = nu = 0
    for p in dels:
        with open(p, newline="") as fh:
            nd += _emit(list(csv.DictReader(fh)), "DEL",
                        "Deletion length - bp", "Location of the deleted region",
                        "Deletion of replication origins", out)
        sys.stderr.write("[mkMitoBreak] deletions  <- %s\n" % os.path.basename(p))
    for p in dups:
        with open(p, newline="") as fh:
            nu += _emit(list(csv.DictReader(fh)), "DUP",
                        "Duplication length - bp", "Location of the duplicated region",
                        "Duplicates replication origin", out)
        sys.stderr.write("[mkMitoBreak] duplications <- %s\n" % os.path.basename(p))

    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(HEADER)
    for row in out:
        w.writerow(row)
    sys.stderr.write("[mkMitoBreak] wrote %d rows (%d DEL, %d DUP)\n" % (len(out), nd, nu))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
