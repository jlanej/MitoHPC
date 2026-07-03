#!/usr/bin/env python3
"""
Build a self-contained, interactive HTML report of mitochondrial structural-variant (large
deletion) calls across a cohort: a circular mtDNA overview + a linear genome browser with
gene/complex annotations, a per-position deletion-frequency map, VAF-encoded calls, live
filtering, and summary statistics. No external dependencies — the output opens offline in any
browser.

Inputs:
  --tab       cohort long table ($ODIR/sv.tab from getSVSummary.sh; or a single $O.sv.tab)
  --genes     6-col BED(.gz) of mtDNA features (RefSeq/genes.bed.gz)
  --nsamples  total cohort size (samples processed; default = distinct samples in --tab)
  --mtlen     mitochondrial genome length (default 16569)
  --out       output HTML path ($ODIR/sv.report.html)
"""
import argparse
import base64
import datetime
import gzip
import json
import os
import sys


# --- mtDNA feature -> functional category (OXPHOS complex / RNA / control) ---
CAT_LABEL = {
    "ci": "Complex I (NADH dehydrogenase)", "ciii": "Complex III (cytb)",
    "civ": "Complex IV (cyt c oxidase)", "cv": "Complex V (ATP synthase)",
    "rrna": "rRNA", "trna": "tRNA", "dloop": "D-loop / control", "other": "other",
}
CAT_COLOR = {
    "ci": "#378ADD", "ciii": "#7F77DD", "civ": "#1D9E75", "cv": "#EF9F27",
    "rrna": "#97C459", "trna": "#B4B2A9", "dloop": "#ED93B1", "other": "#C9C7BD",
}
CAT_ORDER = ["ci", "ciii", "civ", "cv", "rrna", "trna", "dloop"]


def category(name):
    n = name.upper()
    if n.startswith("ND"):
        return "ci"
    if n.startswith("COX") or (n.startswith("CO") and n[2:3].isdigit()):
        return "civ"
    if n.startswith("ATP"):
        return "cv"
    if n.startswith("CYB") or n.startswith("CYTB") or n.startswith("MT-CYB"):
        return "ciii"
    if n.startswith("RNR") or n.endswith("RRNA"):
        return "rrna"
    if n.startswith("TRN") or n.startswith("MT-T"):
        return "trna"
    if n.startswith("DLOOP") or n.startswith("HV") or n.startswith("OL") or n == "CR":
        return "dloop"
    return "other"


def load_features(path, chrom_hint="chrM"):
    feats = []
    if not path or not os.path.exists(path):
        return feats
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 4 or not f[1].isdigit():
                continue
            if f[0] != chrom_hint:
                continue
            name = f[3]
            feats.append({"name": name, "s": int(f[1]) + 1, "e": int(f[2]),
                          "cat": category(name)})
    return feats


def load_calls(path):
    calls = []
    if not path or not os.path.exists(path):
        return calls
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if not lines:
        return calls
    hdr = lines[0].lstrip("#").split("\t")
    for ln in lines[1:]:
        r = dict(zip(hdr, ln.split("\t")))
        try:
            calls.append({
                "smp": r.get("sample", ""),
                "bp5": int(r["pos_bp5"]), "end": int(r["end_bp3"]),
                "len": int(r.get("svlen", 0)),
                "vaf": float(r.get("af_coverage", 0) or 0),   # PRIMARY heteroplasmy = coverage-dosage AFC
                "afj": float(r.get("af_junction", 0) or 0),   # junction-fraction evidence
                "jr": int(r.get("jr", 0) or 0),               # split (junction) read count
                "jsup": r.get("jsup", ""),                    # split-read evidence tier HIGH/MOD/LOW
                "srsb": float(r.get("srsb", 0) or 0),         # strand balance
                "claim": r.get("svclaim", ""),                # DJ (depth+junction) vs J (junction-only)
                "pass": r.get("filter", "") == "PASS",
                "filter": r.get("filter", ""),
                "cls": r.get("delclass", ""),
                "common": r.get("common", "0") == "1",
                "ngene": int(r.get("ngene", 0) or 0),
                "genes": r.get("gene_list", ".") if r.get("gene_list", ".") != "." else "",
                "flags": r.get("flags", ".") if r.get("flags", ".") != "." else "",
                # cohort-only MitoBreak annotation (added by getSVSummary.sh); "" when absent/no match
                "mb": r.get("mitobreak", ".") if r.get("mitobreak", ".") != "." else "",
            })
        except (KeyError, ValueError):
            continue
    return calls


def load_plots(path, dedup_bp=25):
    """samplot manifest TSV (sample, bp5, end, afc, svlen, svconf, png_path) -> representative plot
    dicts with the PNG base64-EMBEDDED (so the report stays single-file/offline).

    SUBSAMPLING: to keep the gallery legible, calls whose breakpoints fall in the same cluster
    (rounded to `dedup_bp`, the same rounding the recurrence table uses; 0 disables) are collapsed to
    ONE representative — the highest-heteroplasmy call (tie-broken by SVCONF). Each representative
    records how many distinct samples shared that site (`nsmp`). Returns (representatives, n_total)."""
    raw = []
    if path and os.path.exists(path):
        seen = set()
        with open(path) as fh:
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if len(f) < 7:
                    continue
                smp, bp5, end, afc, svlen, svconf, png = f[:7]
                filt = f[7] if len(f) > 7 else ""      # optional (appended after png for back-compat)
                flags = f[8] if len(f) > 8 else ""
                key = (smp, bp5, end)
                if key in seen or not os.path.exists(png):
                    continue
                seen.add(key)
                try:
                    raw.append({"smp": smp, "bp5": int(bp5), "end": int(end), "afc": float(afc or 0),
                                "svlen": int(svlen or 0), "svconf": svconf, "png": png,
                                "filter": (filt if filt not in ("", ".") else ""),
                                "flags": (flags if flags not in ("", ".") else "")})
                except ValueError:
                    continue
    n_total = len(raw)

    def conf(p):
        try:
            return float(p["svconf"])
        except (TypeError, ValueError):
            return -1.0
    if dedup_bp and dedup_bp > 0:
        groups = {}
        for p in raw:
            k = (round(p["bp5"] / dedup_bp) * dedup_bp, round(p["end"] / dedup_bp) * dedup_bp)
            groups.setdefault(k, []).append(p)
        reps = []
        for members in groups.values():
            rep = dict(max(members, key=lambda m: (m["afc"], conf(m))))
            rep["nsmp"] = len({m["smp"] for m in members})
            reps.append(rep)
    else:
        reps = [dict(p, nsmp=1) for p in raw]

    plots = []                                   # base64 ONLY the representatives
    for rep in sorted(reps, key=lambda p: (-p["afc"], p["smp"])):
        try:
            with open(rep["png"], "rb") as ph:
                rep["png"] = base64.b64encode(ph.read()).decode()
            plots.append(rep)
        except OSError:
            continue
    return plots, n_total


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tab", required=True)
    ap.add_argument("--genes")
    ap.add_argument("--nsamples", type=int, default=0)
    ap.add_argument("--mtlen", type=int, default=16569)
    ap.add_argument("--chrom", default="chrM")
    ap.add_argument("--plots", help="samplot manifest TSV (sample,bp5,end,afc,svlen,svconf,png) to embed")
    ap.add_argument("--plot-dedup", type=int, default=25,
                    help="collapse plots whose breakpoints cluster within N bp to one representative (0=off)")
    ap.add_argument("--plot-all", action="store_true",
                    help="gallery shows EVERY call unfiltered/unsubsampled (sets the report's wording; "
                         "pair with the permissive svplot.sh filter + --plot-dedup 0)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    feats = load_features(args.genes, args.chrom)
    calls = load_calls(args.tab)
    plots, plots_total = load_plots(args.plots, args.plot_dedup)
    nsamp = args.nsamples or len({c["smp"] for c in calls}) or 1

    data = {
        "meta": {
            "mtlen": args.mtlen, "nsamples": nsamp, "ncalls": len(calls),
            "generated": datetime.date.today().strftime("%Y-%m-%d"),
            "plotsTotal": plots_total, "plotDedup": args.plot_dedup,
            "plotAll": bool(args.plot_all),
        },
        "features": feats,
        "calls": calls,
        "plots": plots,
        "catColor": CAT_COLOR, "catLabel": CAT_LABEL, "catOrder": CAT_ORDER,
    }
    html = TEMPLATE.replace("/*__DATA__*/", json.dumps(data, separators=(",", ":")))
    with open(args.out, "w") as fh:
        fh.write(html)
    sys.stderr.write("[svReport] %d calls, %d samples, %d plots -> %s\n"
                     % (len(calls), nsamp, len(plots), args.out))


TEMPLATE = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>MitoHPC — mtDNA structural-variant cohort report</title>
<style>
:root{--bg:#fbfbf9;--surf:#ffffff;--ink:#1d1d1b;--mut:#56554f;--hint:#84837a;--line:#e3e1d8;--accent:#534ab7;--ok:#137333;--warn:#9a5b00}
:root[data-svtheme=dark]{--bg:#13120c;--surf:#211f16;--ink:#f3f2e9;--mut:#bdbcb0;--hint:#92917f;--line:#403d31;--accent:#bcb6ef;--ok:#5fbf85;--warn:#e6a85c}
@media(prefers-color-scheme:dark){:root:not([data-svtheme=light]){--bg:#13120c;--surf:#211f16;--ink:#f3f2e9;--mut:#bdbcb0;--hint:#92917f;--line:#403d31;--accent:#bcb6ef;--ok:#5fbf85;--warn:#e6a85c}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
.wrap{max-width:1180px;margin:0 auto;padding:28px 22px 60px}
h1{font-size:22px;font-weight:500;margin:0 0 2px}
h2{font-size:17px;font-weight:500;margin:30px 0 12px}
.sub{color:var(--mut);font-size:14px;margin:0 0 18px}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:18px 0}
.card{background:var(--surf);border:1px solid var(--line);border-radius:12px;padding:14px 16px}
.card .lab{font-size:12px;color:var(--mut);margin:0 0 4px;text-transform:uppercase;letter-spacing:.03em}
.card .val{font-size:26px;font-weight:500}
.card .val small{font-size:14px;color:var(--mut);font-weight:400}
.panel{background:var(--surf);border:1px solid var(--line);border-radius:14px;padding:16px 18px;margin:14px 0}
.controls{display:flex;flex-wrap:wrap;gap:18px 26px;align-items:center}
.ctl{display:flex;flex-direction:column;gap:4px;font-size:13px;color:var(--mut)}
.ctl .row{display:flex;align-items:center;gap:8px}
input[type=range]{width:150px}
input[type=text]{height:32px;padding:0 10px;border:1px solid var(--line);border-radius:8px;background:var(--surf);color:var(--ink);font-size:13px}
.chk{display:inline-flex;align-items:center;gap:5px;font-size:13px;color:var(--ink);cursor:pointer}
.legend{display:flex;flex-wrap:wrap;gap:6px 14px;font-size:12px;color:var(--mut);margin:6px 0 0}
.legend span{display:inline-flex;align-items:center;gap:5px}
.sw{width:11px;height:11px;border-radius:3px;display:inline-block}
.viz{width:100%;overflow:visible}
.two{display:grid;grid-template-columns:1fr 1fr;gap:16px}
@media(max-width:780px){.two{grid-template-columns:1fr}}
table{width:100%;border-collapse:collapse;font-size:13px}
th,td{text-align:left;padding:6px 8px;border-bottom:1px solid var(--line)}
th{color:var(--mut);font-weight:500}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
.tip{position:absolute;pointer-events:none;background:var(--surf);border:1px solid var(--line);border-radius:8px;padding:8px 10px;font-size:12px;max-width:280px;opacity:0;transition:opacity .08s;z-index:9;box-shadow:0 2px 8px rgba(0,0,0,.12)}
.tip b{font-weight:500}
details{margin:8px 0}summary{cursor:pointer;color:var(--accent);font-size:14px}
.gloss{margin:8px 0 0;font-size:13px}.gloss dt{font-weight:600;color:var(--ink);margin-top:9px}.gloss dd{margin:2px 0 0 0;color:var(--mut);line-height:1.5}
.muted{color:var(--mut)}.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
.sr-only{position:absolute;width:1px;height:1px;overflow:hidden;clip:rect(0 0 0 0)}
.vafbar{height:10px;border-radius:5px;background:linear-gradient(90deg,#FDE3AE,#F0997B,#D85A30,#A32D2D)}
.pill{font-size:11px;padding:1px 7px;border-radius:10px;border:1px solid var(--line);color:var(--mut)}
</style></head>
<body><div class="wrap">
<h2 class="sr-only">Interactive cohort report of mitochondrial DNA deletions: a circular and linear genome map of where deletions occur, how frequent they are, and their heteroplasmy.</h2>
<div style="display:flex;justify-content:space-between;align-items:flex-start;gap:12px">
<h1>mtDNA structural-variant cohort report</h1>
<button id="theme" aria-label="toggle colour theme" style="height:32px;padding:0 12px;border:1px solid var(--line);border-radius:8px;background:var(--surf);color:var(--ink);cursor:pointer;font-size:13px;white-space:nowrap">theme: auto</button>
</div>
<p class="sub" id="subtitle"></p>

<div class="cards" id="cards"></div>

<div class="panel">
  <div class="controls">
    <div class="ctl"><label for="fpass" class="chk"><input type="checkbox" id="fpass" checked> PASS calls only</label>
      <label for="fcommon" class="chk"><input type="checkbox" id="fcommon"> common deletion (del4977) only</label>
      <label for="fmb" class="chk" title="Show only calls whose breakpoints were previously reported in the MitoBreak database (within HP_SV_MITOBREAK_TOL bp)."><input type="checkbox" id="fmb"> MitoBreak-reported only</label></div>
    <div class="ctl"><div class="row"><span>min heteroplasmy (VAF)</span></div>
      <div class="row"><input type="range" id="fvaf" min="0" max="1" step="0.01" value="0"><span class="mono" id="fvafv" style="min-width:34px">0%</span></div></div>
    <div class="ctl"><span>class</span><div class="row" id="fcls"></div></div>
    <div class="ctl"><span title="Split-read evidence — independent of read depth. Filter low-level events by junction support tier (JSUP) and read count.">split-read evidence (JSUP)</span>
      <div class="row" id="fjsup"></div></div>
    <div class="ctl"><div class="row"><span>min split reads (JR)</span></div>
      <div class="row"><input type="range" id="fjr" min="0" max="50" step="1" value="0"><span class="mono" id="fjrv" style="min-width:34px">0</span></div></div>
    <div class="ctl"><label for="fsmp">sample contains</label><input type="text" id="fsmp" placeholder="all samples" size="14"></div>
    <div class="ctl"><span>&nbsp;</span><button id="reset" style="height:32px;padding:0 12px;border:1px solid var(--line);border-radius:8px;background:var(--surf);color:var(--ink);cursor:pointer">reset</button></div>
  </div>
  <p class="sub" style="margin:8px 0 0">Tip: uncheck <b>PASS calls only</b> and set <b>split-read evidence</b> to <span class="mono">HIGH/MOD</span> to curate credible <i>low-level</i> junctions (clean split reads, no depth signal) apart from depth-supported calls.</p>
</div>

<h2>cohort overview — circular mtDNA</h2>
<p class="sub">Outer ring: genes coloured by OXPHOS (oxidative phosphorylation) complex. Inner arcs: each deletion, coloured by heteroplasmy (VAF, variant allele fraction). The shaded inner band is the per-position deletion frequency across the cohort.</p>
<div id="circ"></div>
<div class="legend" id="legend"></div>

<h2>genome browser — linear track</h2>
<p class="sub">Top: deletions overlapping each position. Middle: gene / feature annotation. Bottom: individual deletions (width = span, colour = VAF). Hover for details.</p>
<div style="position:relative"><div id="lin"></div><div class="tip" id="tip"></div></div>

<div class="two">
  <div><h2>heteroplasmy distribution</h2><div id="vafhist"></div></div>
  <div><h2>deletion size distribution</h2><div id="szhist"></div></div>
</div>

<h2>recurrent deletions</h2>
<p class="sub">Distinct deletion sites (breakpoints rounded to 25 bp), ranked by the number of samples carrying them.</p>
<table id="rec"><thead><tr><th>breakpoints (m.)</th><th class="n">size (bp)</th><th class="n">samples</th><th class="n">cohort %</th><th class="n">median VAF</th><th>genes</th><th title="Closest previously-reported MitoBreak breakpoint id within tolerance">MitoBreak</th><th>tags</th></tr></thead><tbody></tbody></table>

<style>
#plotsection .plotwrap{display:flex;gap:18px;align-items:flex-start;flex-wrap:wrap;margin-top:6px}
#plotsection #plottbl{flex:0 0 auto;display:block;max-height:520px;overflow:auto}
#plotsection #plottbl tbody tr{cursor:pointer}
#plotsection #plottbl tbody tr.sel{background:rgba(10,140,90,.16);font-weight:600}
#plotsection .plotfig{flex:1 1 460px;min-width:340px;margin:0}
#plotsection .plotfig img{max-width:100%;border:1px solid #ddd;border-radius:6px;background:#fff}
#plotsection figcaption{margin-top:6px}
</style>
<section id="plotsection" hidden>
<h2>structural-variant plots <span class="muted" style="font-weight:400;font-size:13px">(samplot)</span></h2>
<p class="sub" id="plotdesc"></p>
<div class="plotwrap">
  <table id="plottbl"><thead><tr><th>sample</th><th>breakpoints (m.)</th><th class="n">size</th><th class="n">VAF</th><th class="n">SVCONF</th><th class="n">samples</th><th>status</th></tr></thead><tbody></tbody></table>
  <figure class="plotfig"><img id="samplotimg" alt="samplot image"><figcaption id="plotcap" class="muted">select a row to view its samplot image</figcaption></figure>
</div>
</section>

<details><summary>how these calls are made</summary>
<p class="muted" style="font-size:14px">Each sample's circular-aware chrM alignment is scanned for <b>split reads</b> (reads whose two halves map across a deletion junction). Junctions are clustered to base-pair breakpoints. A deletion is reported <b>PASS</b> either when a <b>coverage drop</b> and the junction agree (<span class="mono">SVCLAIM=DJ</span>) <i>or</i> when the split-read junction alone is strong and clean (<span class="mono">SVCLAIM=J</span> — valuable because mtDNA read depth is finicky). Heteroplasmy (VAF) is the <b>coverage-dosage</b> fraction (<span class="mono">AFC</span>, shown here), with the junction-read fraction (<span class="mono">AFJ</span>) as corroborating evidence; their agreement is a QC signal (<span class="mono">AFDIFF</span>). The breakpoint microhomology / direct repeat (e.g. the 13&nbsp;bp repeat of the common deletion) is reported as <span class="mono">HOMLEN</span>/<span class="mono">DELCLASS</span>. See <span class="mono">docs/SV_METHODS.md</span>.</p></details>
<details><summary>glossary &mdash; definitions of terms used in this report</summary>
<dl class="gloss">
  <dt>VAF (variant allele fraction) / heteroplasmy</dt><dd>Fraction of mtDNA molecules carrying the deletion (0&ndash;100%). The colour scale and the "VAF" shown here are the <b>coverage-dosage</b> estimate <span class="mono">AFC</span> (below) &mdash; the standard heteroplasmy measure for large mtDNA deletions; the junction estimate <span class="mono">AFJ</span> is reported as corroborating evidence.</dd>
  <dt>Class&nbsp;I / II / III <span class="muted">(<span class="mono">DELCLASS</span>)</span></dt><dd>Breakpoint-homology class, a clue to the deletion mechanism: <b>Class&nbsp;I</b> = the breakpoints sit in a <b>perfect direct repeat &ge;5&nbsp;bp</b> (e.g. the 13&nbsp;bp repeat of the common deletion &mdash; slipped-strand mispairing); <b>Class&nbsp;II</b> = <b>1&ndash;4&nbsp;bp microhomology</b>; <b>Class&nbsp;III</b> = <b>no breakpoint homology</b> (blunt). Larger homology = the breakpoint can sit anywhere within the repeat, so its position is reported as imprecise.</dd>
  <dt><span class="mono">HOMLEN</span> / <span class="mono">HOMSEQ</span></dt><dd>Length and sequence of that breakpoint microhomology / direct repeat (drives the Class above).</dd>
  <dt><span class="mono">AFC</span> &mdash; coverage VAF <span class="muted">(primary)</span></dt><dd>Heteroplasmy from the dosage loss: <span class="mono">AFC = 1 &minus; trimmed-median(depth inside) / trimmed-median(depth in flanks)</span>, computed over D-loop/origin/homopolymer/NUMT-masked, transition-excluded windows. This is the reported VAF.</dd>
  <dt><span class="mono">AFJ</span> &mdash; junction VAF <span class="muted">(evidence)</span></dt><dd><span class="mono">AFJ = JR / (JR + SR)</span>, where SR counts wild-type reads aligned contiguously across the breakpoint. Used to confirm the deletion and gate calls (the coverage drop must be backed by a proportional junction), not as the primary load.</dd>
  <dt><span class="mono">AFDIFF</span></dt><dd>|AFJ &minus; AFC| &mdash; how far the two estimates disagree. Large values flag amplification bias or a duplication/artifact masquerading as a deletion (a QC signal).</dd>
  <dt><span class="mono">JR</span> / <span class="mono">SR</span></dt><dd><b>JR</b> = number of distinct <b>split (junction) reads</b> spanning the deletion breakpoint; <b>SR</b> = <b>wild-type spanning reads</b> (the AFJ denominator at the breakpoints).</dd>
  <dt><span class="mono">JSUP</span> &mdash; split-read evidence tier</dt><dd>A depth-independent curation lens for the junction itself: <b>HIGH</b> = consistent breakpoint, &ge; the read threshold, two-strand support; <b>MOD</b> = a clean, consistent junction but low read count or one-strand (a <i>credible low-level event</i> &mdash; a handful of reads all clipping at the same base); <b>LOW</b> = scattered breakpoint sizes (likely a mapping/NUMT artifact). Filter to <span class="mono">HIGH</span>/<span class="mono">MOD</span> (with <b>PASS only</b> off) to review low-heteroplasmy junctions that have strong split-read evidence but no read-depth signal. Backed by <span class="mono">SRCONS</span> (size consistency, microhomology-invariant) and <span class="mono">SRSB</span> (strand balance).</dd>
  <dt><span class="mono">CVGR</span></dt><dd>Coverage ratio = median depth inside the deletion / median depth in the flanks. <span class="mono">&le;0.9</span> means a &ge;10% coverage drop (the corroboration gate).</dd>
  <dt><span class="mono">SVCLAIM</span></dt><dd>Which evidence supports the call: <b>DJ</b> = the split-read junction <i>and</i> the coverage drop agree; <b>J</b> = split-read junction only (no confirming coverage drop).</dd>
  <dt><b>PASS</b> &amp; filters</dt><dd><b>PASS</b> via either the <b>DJ</b> path (coverage drop corroborated by a proportional junction) or the <b>J</b> path (junction-strong: clean, well-supported split reads, no coverage drop required). Non-PASS reasons: <span class="mono">lowJR</span> (too few junction reads), <span class="mono">no_cvg_drop</span> (no &ge;10% coverage drop), <span class="mono">low_dosage</span> (coverage-dosage AF below threshold), <span class="mono">lowAFJ</span> (corrected junction VAF below the floor), <span class="mono">unexplained_drop</span> (coverage drop not backed by a proportional junction &mdash; a likely artifact), <span class="mono">fragile_weakJ</span> (D-loop/NUMT/origin breakpoint without strong junction support), <span class="mono">bigdel_weakJ</span> (very large deletion with weak junction support), <span class="mono">WRAP</span> (breakpoint at the artificial origin / origin-crossing; deletion-vs-duplication unresolved), <span class="mono">lowDP</span> (flanking depth below threshold).</dd>
  <dt>tags / flags</dt><dd><span class="mono">COMMON</span> = matches the common deletion del4977; <span class="mono">REPEAT</span> = breakpoint in the del4977 13&nbsp;bp direct repeat; <span class="mono">NUMT</span> / <span class="mono">HP</span> / <span class="mono">DLOOP</span> = breakpoint overlaps a known NUMT-like site / homopolymer run / control region (D-loop); <span class="mono">WRAP</span> = breakpoint near the origin. Flags annotate; they do not by themselves reject a call (except <span class="mono">WRAP</span>).</dd>
  <dt>common deletion (<span class="mono">del4977</span>)</dt><dd>The canonical ~4977&nbsp;bp "common" mtDNA deletion (m.8470_13447), flanked by a 13&nbsp;bp direct repeat; accumulates with age in post-mitotic tissue.</dd>
  <dt><span class="mono">HGVS</span></dt><dd>Approximate deletion span in Human Genome Variation Society (HGVS) notation on the revised Cambridge Reference Sequence (rCRS, NC_012920.1).</dd>
  <dt><span class="mono">NUMT</span></dt><dd>Nuclear mitochondrial DNA segment &mdash; a mtDNA-derived paralogue in the nuclear genome; reads from NUMTs are the classic source of spurious mtDNA structural calls (suppressed upstream by competitive alignment, and flagged here).</dd>
</dl></details>
<p class="sub" id="foot"></p>

<script>
const DATA=/*__DATA__*/;
const M=DATA.meta, MT=M.mtlen, NS=M.nsamples, F=DATA.features, ALL=DATA.calls;
const PLOTS=DATA.plots||[];
const CC=DATA.catColor, CL=DATA.catLabel, CO=DATA.catOrder;
const $=id=>document.getElementById(id);
const rnd=(x,d=0)=>{const p=Math.pow(10,d);return Math.round(x*p)/p};
const pct=x=>rnd(x*100)+'%';
const SVGNS='http://www.w3.org/2000/svg';
function el(t,a){const e=document.createElementNS(SVGNS,t);for(const k in(a||{}))e.setAttribute(k,a[k]);return e}
const accent=()=>getComputedStyle(document.documentElement).getPropertyValue('--accent').trim()||'#534ab7';
// heteroplasmy colour ramp (light->deep red)
const RAMP=[[0,'#FDE3AE'],[.33,'#F0997B'],[.66,'#D85A30'],[1,'#A32D2D']];
function vcol(v){v=Math.max(0,Math.min(1,v));for(let i=1;i<RAMP.length;i++){if(v<=RAMP[i][0]){const[a,ca]=RAMP[i-1],[b,cb]=RAMP[i];return lerp(ca,cb,(v-a)/(b-a||1))}}return RAMP[RAMP.length-1][1]}
function hex(c){return[parseInt(c.slice(1,3),16),parseInt(c.slice(3,5),16),parseInt(c.slice(5,7),16)]}
function lerp(a,b,t){const x=hex(a),y=hex(b);return'#'+[0,1,2].map(i=>Math.round(x[i]+(y[i]-x[i])*t).toString(16).padStart(2,'0')).join('')}

// ---- filter state ----
const st={pass:true,common:false,mb:false,vaf:0,cls:{I:true,II:true,III:true,'':true},smp:'',jsup:{HIGH:true,MOD:true,LOW:true,'':true},minjr:0};
function filtered(){return ALL.filter(c=>(!st.pass||c.pass)&&(!st.common||c.common)&&(!st.mb||c.mb)&&c.vaf>=st.vaf&&(st.cls[c.cls]!==false)&&(st.jsup[c.jsup]!==false)&&(c.jr>=st.minjr)&&(!st.smp||c.smp.toLowerCase().includes(st.smp)))}

// ---- summary cards ----
function cards(cs){
  const samp=new Set(cs.map(c=>c.smp)), pas=cs.filter(c=>c.pass);
  const psamp=new Set(pas.map(c=>c.smp)), com=new Set(cs.filter(c=>c.common).map(c=>c.smp));
  const vafs=pas.map(c=>c.vaf).sort((a,b)=>a-b), med=vafs.length?vafs[vafs.length>>1]:0;
  const card=(lab,val,sub)=>`<div class="card"><p class="lab">${lab}</p><div class="val">${val}${sub?' <small>'+sub+'</small>':''}</div></div>`;
  $('cards').innerHTML=
    card('cohort samples',NS)+
    card('deletions shown',cs.length,pas.length+' PASS')+
    card('samples affected',psamp.size,pct(psamp.size/NS))+
    card('common deletion',com.size,'samples · '+pct(com.size/NS))+
    card('median VAF (PASS)',pct(med))+
    card('distinct sites',recurrence(cs).length);
}

// ---- circular overview ----
function circle(cs){
  const W=560,R=250,cx=W/2,cy=W/2+6;const ang=p=>(p/MT)*2*Math.PI-Math.PI/2;
  const pc=(r,a)=>[cx+r*Math.cos(a),cy+r*Math.sin(a)];
  function arc(r,p1,p2,w){const a1=ang(p1),a2=ang(p2);const large=(a2-a1)>Math.PI?1:0;
    const o1=pc(r+w/2,a1),o2=pc(r+w/2,a2),i2=pc(r-w/2,a2),i1=pc(r-w/2,a1);
    return`M${o1[0]} ${o1[1]} A${r+w/2} ${r+w/2} 0 ${large} 1 ${o2[0]} ${o2[1]} L${i2[0]} ${i2[1]} A${r-w/2} ${r-w/2} 0 ${large} 0 ${i1[0]} ${i1[1]} Z`}
  const svg=el('svg',{viewBox:`0 0 ${W} ${W+12}`,class:'viz',width:W,height:W+12,role:'img'});
  svg.appendChild(el('title',{})).textContent='Circular map of mtDNA genes and cohort deletions';
  // gene ring
  F.forEach(f=>{const p=el('path',{d:arc(R,f.s,f.e,16),fill:CC[f.cat]||CC.other,opacity:.92});p.appendChild(el('title',{})).textContent=f.name+' ('+CL[f.cat]+')';svg.appendChild(p)});
  // position ticks every 2kb
  for(let p=0;p<MT;p+=2000){const a=ang(p),o=pc(R+10,a),i=pc(R+16,a);svg.appendChild(el('line',{x1:i[0],y1:i[1],x2:o[0],y2:o[1],stroke:'var(--line)','stroke-width':1}));
    const t=pc(R+26,a),tx=el('text',{x:t[0],y:t[1],fill:'var(--hint)','font-size':10,'text-anchor':'middle','dominant-baseline':'middle'});tx.textContent=(p/1000)+'k';svg.appendChild(tx)}
  // frequency band (calls covering each base) just inside the genes
  const cov=coverage(cs),mx=Math.max(1,...cov);const fr=R-26;
  const AC=accent();for(let p=0;p<MT;p+=40){const v=cov[p]/mx;if(v<=0)continue;svg.appendChild(el('path',{d:arc(fr,p,Math.min(p+40,MT),3+18*v),fill:AC,opacity:.16+.32*v}))}
  // deletion arcs (inner), colour by VAF, ordered by size so big ones sit inside
  const dl=cs.slice().sort((a,b)=>b.len-a.len);let r=fr-30;const step=Math.max(1.6,Math.min(7,(r-70)/Math.max(1,dl.length)));
  dl.forEach(c=>{svg.appendChild(el('path',{d:arc(r,c.bp5,c.end,Math.max(1.4,step*.8)),fill:vcol(c.vaf),opacity:c.pass?.95:.4}));r-=step;if(r<70)r=fr-30});
  const cap=el('text',{x:cx,y:cy,fill:'var(--mut)','font-size':12,'text-anchor':'middle'});cap.textContent=MT.toLocaleString()+' bp';svg.appendChild(cap);
  const cap2=el('text',{x:cx,y:cy+16,fill:'var(--hint)','font-size':11,'text-anchor':'middle'});cap2.textContent='chrM';svg.appendChild(cap2);
  $('circ').innerHTML='';$('circ').appendChild(svg);
}

// ---- linear genome browser ----
function coverage(cs){const a=new Float64Array(MT+2);cs.forEach(c=>{const s=Math.max(1,c.bp5+1),e=Math.min(MT,c.end);if(e>=s){a[s]++;a[e+1]--}});const cov=new Float64Array(MT+1);let run=0;for(let p=1;p<=MT;p++){run+=a[p];cov[p]=run}return cov}
function lanes(cs){const so=cs.slice().sort((a,b)=>a.bp5-b.bp5||a.end-b.end);const ends=[];so.forEach(c=>{let L=ends.findIndex(e=>e<c.bp5);if(L<0){L=ends.length;ends.push(0)}ends[L]=c.end;c._lane=L});return Math.max(1,ends.length)}
function linear(cs){
  const W=1140,ML=46,MR=14,iw=W-ML-MR,xs=p=>ML+(p-1)/(MT-1)*iw;
  const hFreq=110,hFeat=54,laneH=7,nl=lanes(cs),hCalls=Math.max(40,nl*laneH+8);
  const H=20+hFreq+16+hFeat+12+hCalls+26;
  const svg=el('svg',{viewBox:`0 0 ${W} ${H}`,class:'viz',width:'100%',height:H,role:'img'});
  svg.appendChild(el('title',{})).textContent='Linear mtDNA browser of cohort deletions';
  let y=18;
  // frequency area
  const cov=coverage(cs),mx=Math.max(1,...cov);
  let d='M'+xs(1)+' '+(y+hFreq);for(let p=1;p<=MT;p+=20)d+=' L'+rnd(xs(p),1)+' '+rnd(y+hFreq-cov[p]/mx*hFreq,1);d+=' L'+xs(MT)+' '+(y+hFreq)+' Z';
  const AC=accent();svg.appendChild(el('path',{d,fill:AC,opacity:.22}));
  svg.appendChild(el('path',{d:d.replace(' Z','').replace('M'+xs(1)+' '+(y+hFreq),'M'+xs(1)+' '+rnd(y+hFreq-cov[1]/mx*hFreq,1)),fill:'none',stroke:AC,'stroke-width':1.4,opacity:.85}));
  const yl=el('text',{x:ML,y:y-4,fill:'var(--hint)','font-size':11});yl.textContent='deletions overlapping (max '+mx+')';svg.appendChild(yl);
  y+=hFreq+16;
  // feature track
  F.forEach(f=>{const x1=xs(f.s),x2=xs(f.e),w=Math.max(1.2,x2-x1);const tr=f.cat==='trna';
    const r=el('rect',{x:x1,y:tr?y+14:y,width:w,height:tr?20:34,rx:tr?1:3,fill:CC[f.cat]||CC.other,opacity:.92});
    r.appendChild(el('title',{})).textContent=f.name+' · '+CL[f.cat]+' · m.'+f.s+'-'+f.e;svg.appendChild(r);
    if(!tr&&w>26){const t=el('text',{x:(x1+x2)/2,y:y+22,fill:'#241f12','font-size':10,'text-anchor':'middle'});t.textContent=f.name.length>7?f.name.slice(0,7):f.name;svg.appendChild(t)}});
  y+=hFeat+12;
  // call lanes
  cs.slice().sort((a,b)=>a.bp5-b.bp5).forEach(c=>{const x1=xs(c.bp5),x2=xs(c.end);
    const r=el('rect',{x:x1,y:y+c._lane*laneH,width:Math.max(1.5,x2-x1),height:laneH-2,rx:2,fill:vcol(c.vaf),opacity:c.pass?.95:.45,stroke:c.common?'#241f12':'none','stroke-width':c.common?.8:0,'data-c':1});
    r.__c=c;svg.appendChild(r)});
  y+=hCalls+6;
  // axis
  for(let p=0;p<=MT;p+=2000){const x=xs(Math.max(1,p));svg.appendChild(el('line',{x1:x,y1:y,x2:x,y2:y+4,stroke:'var(--line)'}));
    const t=el('text',{x:x,y:y+16,fill:'var(--hint)','font-size':10,'text-anchor':'middle'});t.textContent=(p/1000)+'k';svg.appendChild(t)}
  const c=$('lin');c.innerHTML='';c.appendChild(svg);
  // tooltip
  const tip=$('tip');svg.addEventListener('mousemove',e=>{const t=e.target;if(t.__c){const c=t.__c;
    tip.innerHTML=`<b>${c.smp}</b><br>m.${c.bp5+1}_${c.end}del · ${c.len.toLocaleString()} bp<br>VAF ${pct(c.vaf)} (junction ${pct(c.afj)}) · class ${c.cls} · ${c.filter}<br><span class="muted">split reads JR=${c.jr} · evidence ${c.jsup}${c.claim?'/'+c.claim:''} · strand-bal ${c.srsb.toFixed(2)}</span>`+(c.genes?'<br><span class="muted">'+c.genes+'</span>':'')+(c.flags?'<br><span class="muted">'+c.flags+'</span>':'');
    const b=c.getBoundingClientRect?null:null;tip.style.opacity=1;tip.style.left=(e.offsetX+14)+'px';tip.style.top=(e.offsetY-10)+'px'}else tip.style.opacity=0});
  svg.addEventListener('mouseleave',()=>tip.style.opacity=0);
}

// ---- histograms ----
function hist(id,vals,bins,fmt,acc){const W=540,H=150,ML=38,MB=26,iw=W-ML-12,ih=H-MB-10;
  const cnt=new Array(bins).fill(0);vals.forEach(v=>{let b=Math.min(bins-1,Math.floor(acc(v)*bins));if(b<0)b=0;cnt[b]++});
  const mx=Math.max(1,...cnt),bw=iw/bins;
  const svg=el('svg',{viewBox:`0 0 ${W} ${H}`,class:'viz',width:'100%',height:H,role:'img'});
  svg.appendChild(el('title',{})).textContent='histogram';
  const AC=accent();for(let i=0;i<bins;i++){const h=cnt[i]/mx*ih,x=ML+i*bw;svg.appendChild(el('rect',{x:x+1,y:10+ih-h,width:bw-2,height:h,rx:2,fill:AC,opacity:.85}));}
  svg.appendChild(el('line',{x1:ML,y1:10+ih,x2:W-12,y2:10+ih,stroke:'var(--line)'}));
  for(let i=0;i<=4;i++){const x=ML+i/4*iw;const t=el('text',{x:x,y:H-8,fill:'var(--hint)','font-size':10,'text-anchor':'middle'});t.textContent=fmt(i/4);svg.appendChild(t)}
  const ym=el('text',{x:ML,y:8,fill:'var(--hint)','font-size':10});ym.textContent='n (max '+mx+')';svg.appendChild(ym);
  $(id).innerHTML='';$(id).appendChild(svg);
}

// ---- recurrence ----
function recurrence(cs){const m=new Map();cs.forEach(c=>{const k=Math.round(c.bp5/25)*25+'_'+Math.round(c.end/25)*25;
  if(!m.has(k))m.set(k,{bp5:c.bp5,end:c.end,len:c.len,smp:new Set(),vaf:[],genes:c.genes,common:c.common,mb:c.mb||''});
  const g=m.get(k);g.smp.add(c.smp);g.vaf.push(c.vaf);if(c.common)g.common=true;if(c.mb&&!g.mb)g.mb=c.mb});
  return[...m.values()].map(g=>{g.vaf.sort((a,b)=>a-b);g.n=g.smp.size;g.med=g.vaf[g.vaf.length>>1];return g}).sort((a,b)=>b.n-a.n||b.med-a.med)}
function recTable(cs){const rows=recurrence(cs).slice(0,12).map(g=>{
  const tags=(g.common?'<span class="pill" style="color:#A32D2D;border-color:#F0997B">del4977</span> ':'');
  return`<tr><td class="mono">${g.bp5+1}_${g.end}</td><td class="n">${g.len.toLocaleString()}</td><td class="n">${g.n}</td><td class="n">${pct(g.n/NS)}</td><td class="n">${pct(g.med)}</td><td class="muted" style="font-size:12px">${(g.genes||'').split(',').slice(0,4).join(', ')}${(g.genes||'').split(',').length>4?'…':''}</td><td class="mono" style="font-size:12px">${g.mb||''}</td><td>${tags}</td></tr>`}).join('');
  document.querySelector('#rec tbody').innerHTML=rows||'<tr><td colspan="8" class="muted">no calls match the current filters</td></tr>'}

function legend(){$('legend').innerHTML=CO.map(c=>`<span><span class="sw" style="background:${CC[c]}"></span>${CL[c]}</span>`).join('')+
  ' <span style="margin-left:10px">VAF <span class="vafbar" style="width:90px;display:inline-block;vertical-align:-1px"></span> 0→100%</span>'}

// ---- samplot gallery: an interactive table; selecting a row shows its embedded PNG ----
function showPlot(i,row){document.querySelectorAll('#plottbl tbody tr').forEach(r=>r.classList.remove('sel'));
  if(row)row.classList.add('sel');const p=PLOTS[i];
  $('samplotimg').src='data:image/png;base64,'+p.png;
  $('plotcap').textContent=`${p.smp} · m.${p.bp5+1}_${p.end}del · ${p.svlen.toLocaleString()} bp · VAF ${pct(p.afc)} · SVCONF ${p.svconf}`+(p.filter&&p.filter!=='PASS'?` · ${p.filter}`:'')+(p.flags?` · flags: ${p.flags}`:'')+((p.nsmp||1)>1?` · representative of ${p.nsmp} samples`:'')}
function statusCell(p){if(!p.filter)return '';
  if(p.filter==='PASS')return `<span style="color:var(--ok);font-weight:600">PASS</span>`;
  return `<span style="color:var(--warn);font-weight:600" title="flags: ${p.flags||p.filter}">${p.filter}</span>`}
function buildPlots(){if(!PLOTS.length)return;$('plotsection').hidden=false;
  const dd=M.plotDedup||0, tot=M.plotsTotal||PLOTS.length, n=PLOTS.length;
  const filt=M.plotAll
    ?`<b>All</b> calls &mdash; <b>no</b> PASS / heteroplasmy / artifact-region filtering (non-PASS and HP / D-loop / NUMT / WRAP calls included)`
    :`<b>Curated</b> calls &mdash; <b>PASS</b>, heteroplasmy above the configured floor, and breakpoints <b>outside</b> the homopolymer / D-loop / NUMT artifact regions`;
  const sub=dd>0
    ?`<b>Subsampling:</b> calls sharing a breakpoint site (rounded to ${dd} bp) are collapsed to <b>one representative</b> (the highest-heteroplasmy call); the <b>samples</b> column is how many samples carried that site. Showing ${n} representative site${n==1?'':'s'} from ${tot} visualizable call${tot==1?'':'s'}.`
    :`<b>No subsampling</b> &mdash; every visualizable call is shown (the <b>samples</b> column is the per-site sample count). Showing all ${n} call${n==1?'':'s'}.`;
  $('plotdesc').innerHTML=`${filt}, rendered with <span class="mono">samplot</span> (read depth, split reads, the deletion span). ${sub} Click a row to view its image. (Configurable via <span class="mono">HP_SV_PLOT_*</span> / <span class="mono">HP_SV_PLOT_ALL</span> / <span class="mono">--plot-dedup</span>.)`;
  const tb=document.querySelector('#plottbl tbody');
  tb.innerHTML=PLOTS.map((p,i)=>`<tr data-i="${i}"><td>${p.smp}</td><td class="mono">m.${p.bp5+1}_${p.end}del</td><td class="n">${p.svlen.toLocaleString()}</td><td class="n">${pct(p.afc)}</td><td class="n">${p.svconf}</td><td class="n">${p.nsmp||1}</td><td>${statusCell(p)}</td></tr>`).join('');
  tb.querySelectorAll('tr').forEach(r=>r.onclick=()=>showPlot(+r.dataset.i,r));
  showPlot(0,tb.querySelector('tr'))}

function render(){const cs=filtered();cards(cs);circle(cs);linear(cs);
  hist('vafhist',cs.filter(c=>c.pass).map(c=>c.vaf),20,v=>pct(v),v=>v);
  const mlen=Math.max(1,...ALL.map(c=>c.len));hist('szhist',cs.map(c=>c.len),20,v=>rnd(v*mlen/1000,1)+'k',v=>v/mlen);
  recTable(cs);
  $('subtitle').textContent=`${NS} samples · ${ALL.length} deletion calls · generated ${M.generated}`;
  $('foot').textContent=`Generated by MitoHPC svReport. VAF = coverage-dosage heteroplasmy (AFC); AFJ = junction fraction (evidence). Frequency track counts deletions (PASS by default) overlapping each base.`}

// ---- wire controls ----
function ui(){
  $('fcls').innerHTML=['I','II','III'].map(k=>`<label class="chk"><input type="checkbox" data-cls="${k}" checked> ${k}</label>`).join('');
  $('fjsup').innerHTML=['HIGH','MOD','LOW'].map(k=>`<label class="chk" title="${k=='HIGH'?'consistent, >=MINJR reads, two-strand':k=='MOD'?'clean junction but low-count / one-strand (credible low-level)':'scattered breakpoint sizes (likely artifact)'}"><input type="checkbox" data-jsup="${k}" checked> ${k}</label>`).join('');
  $('fpass').onchange=e=>{st.pass=e.target.checked;render()};
  $('fcommon').onchange=e=>{st.common=e.target.checked;render()};
  $('fmb').onchange=e=>{st.mb=e.target.checked;render()};
  $('fvaf').oninput=e=>{st.vaf=+e.target.value;$('fvafv').textContent=pct(st.vaf);render()};
  $('fjr').oninput=e=>{st.minjr=+e.target.value;$('fjrv').textContent=e.target.value;render()};
  $('fsmp').oninput=e=>{st.smp=e.target.value.trim().toLowerCase();render()};
  document.querySelectorAll('[data-cls]').forEach(b=>b.onchange=e=>{st.cls[e.target.dataset.cls]=e.target.checked;render()});
  document.querySelectorAll('[data-jsup]').forEach(b=>b.onchange=e=>{st.jsup[e.target.dataset.jsup]=e.target.checked;render()});
  $('reset').onclick=()=>{st.pass=true;st.common=false;st.mb=false;st.vaf=0;st.smp='';st.cls={I:true,II:true,III:true,'':true};st.jsup={HIGH:true,MOD:true,LOW:true,'':true};st.minjr=0;
    $('fpass').checked=true;$('fcommon').checked=false;$('fmb').checked=false;$('fvaf').value=0;$('fvafv').textContent='0%';$('fsmp').value='';$('fjr').value=0;$('fjrv').textContent='0';
    document.querySelectorAll('[data-cls]').forEach(b=>b.checked=true);document.querySelectorAll('[data-jsup]').forEach(b=>b.checked=true);render()};
  legend();
  buildPlots();
  const TH=['auto','light','dark'];let ti=0;
  try{const s=localStorage.getItem('svtheme');if(s){const i=TH.indexOf(s);if(i>=0)ti=i}}catch(e){}
  function applyTheme(){const t=TH[ti];if(t==='auto')document.documentElement.removeAttribute('data-svtheme');else document.documentElement.setAttribute('data-svtheme',t);$('theme').textContent='theme: '+t;try{localStorage.setItem('svtheme',t)}catch(e){}render()}
  $('theme').onclick=()=>{ti=(ti+1)%3;applyTheme()};
  applyTheme();
}
ui();
</script>
</div></body></html>
"""


if __name__ == "__main__":
    main()
