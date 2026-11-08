#!/usr/bin/env python3
"""Build report/dashboard.html from the Task 3 and Task 4 result files.

Inputs (read-only):
    results/mutation_analysis.json
    results/statistics_summary.json
    results/mutation_types.csv
    results/summary.csv
    results/figures/*.png

Output:
    report/dashboard.html  (self-contained: figures are embedded as base64)

Usage:
    python3 report/make_dashboard.py
"""
import base64
import csv
import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUTPUT = ROOT / "report" / "dashboard.html"

TYPE_COLORS = {
    "substitution": "#5eead4",
    "deletion": "#fbbf24",
    "insertion": "#a78bfa",
    "ambiguous": "#64748b",
}


def number(value):
    return f"{int(value):,}"


def embed(path):
    return base64.b64encode(path.read_bytes()).decode()


def gather():
    analysis = json.loads((RESULTS / "mutation_analysis.json").read_text())
    stats = json.loads((RESULTS / "statistics_summary.json").read_text())
    with (RESULTS / "mutation_types.csv").open(newline="") as handle:
        types = list(csv.DictReader(handle))
    with (RESULTS / "summary.csv").open(newline="") as handle:
        top = list(csv.DictReader(handle))[:20]
    return analysis, stats, types, top


def render(analysis, stats, types, top):
    counts = analysis["counts"]
    max_type = max(int(row["n_mutations"]) for row in types)
    type_rows = "".join(
        '<div class="trow"><span class="tname">{name}</span>'
        '<span class="tbar"><i style="width:{width:.1f}%;background:{color}"></i></span>'
        '<span class="tval">{count} <em>{share:.1f}%</em></span></div>'.format(
            name=html.escape(row["type"]),
            width=int(row["n_mutations"]) / max_type * 100,
            color=TYPE_COLORS.get(row["type"], "#94a3b8"),
            count=number(row["n_mutations"]),
            share=int(row["n_mutations"]) / counts["mutation_rows"] * 100,
        )
        for row in types
    )
    max_freq = max(float(row["frequency"]) for row in top)
    top_rows = "".join(
        '<tr><td class="mono">{mutation}</td><td>{type}</td><td class="num">{seqs}</td>'
        '<td class="fcell"><span class="fbar"><i style="width:{width:.1f}%"></i></span>'
        '<span class="fval">{freq:.1f}%</span></td></tr>'.format(
            mutation=html.escape(row["mutation"]),
            type=html.escape(row["type"]),
            seqs=number(row["n_sequences"]),
            width=float(row["frequency"]) / max_freq * 100,
            freq=float(row["frequency"]) * 100,
        )
        for row in top
    )
    funnel = (
        ("NCBI download", 142254, "families"),
        ("HBV retained", counts["input_sequences"], "Task 2"),
        ("Eligible ≥100 nt", counts["eligible_sequences"], "Task 3"),
        ("Analysed", counts["analyzed_sequences"], "aligned"),
    )
    funnel_html = "".join(
        '<div class="step"><div class="snum">{v}</div>'
        '<div class="slabel">{label}</div><div class="stag">{tag}</div></div>'.format(
            v=number(value), label=html.escape(label), tag=html.escape(tag)
        )
        for label, value, tag in funnel
    )
    top_fig = embed(RESULTS / "figures" / "top_mutations.png")
    genome_fig = embed(RESULTS / "figures" / "mutations_along_genome.png")

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>HBV Genome Pipeline — Dashboard</title>
<style>
:root{{--bg:#0b1220;--card:#121c2e;--line:#1f2f47;--text:#e6edf7;--dim:#8ba1bd;--teal:#5eead4;--amber:#fbbf24}}
*{{box-sizing:border-box}} body{{margin:0;background:radial-gradient(1200px 600px at 15% -10%,#16294a 0%,var(--bg) 55%);color:var(--text);
font:15px/1.5 ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;padding:34px 38px 60px}}
.wrap{{max-width:1180px;margin:0 auto}}
h1{{font-size:27px;margin:0 0 4px;letter-spacing:-.02em}} h1 span{{color:var(--teal)}}
.sub{{color:var(--dim);font-size:13px;margin-bottom:26px}} .sub b{{color:#cfe0f5;font-weight:600}}
.kpis{{display:grid;grid-template-columns:repeat(5,1fr);gap:14px;margin-bottom:26px}}
.kpi{{background:linear-gradient(180deg,#16233a,#111b2c);border:1px solid var(--line);border-radius:14px;padding:16px 17px}}
.kpi .v{{font-size:25px;font-weight:700;letter-spacing:-.02em}} .kpi .l{{color:var(--dim);font-size:11.5px;text-transform:uppercase;letter-spacing:.09em;margin-top:3px}}
.kpi.hl .v{{color:var(--teal)}}
.grid{{display:grid;grid-template-columns:1.15fr .85fr;gap:16px;margin-bottom:16px}}
.card{{background:var(--card);border:1px solid var(--line);border-radius:14px;padding:18px 20px}}
.card h2{{font-size:13px;text-transform:uppercase;letter-spacing:.1em;color:var(--dim);margin:0 0 16px;font-weight:600}}
.funnel{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px}}
.step{{background:#0e1830;border:1px solid #1c2b44;border-radius:11px;padding:12px}}
.snum{{font-size:19px;font-weight:700}} .slabel{{font-size:12px;color:#cfe0f5}} .stag{{font-size:10.5px;color:var(--dim);text-transform:uppercase;letter-spacing:.08em}}
.trow{{display:grid;grid-template-columns:96px 1fr auto;align-items:center;gap:12px;margin:9px 0}}
.tname{{font-size:13px;text-transform:capitalize}} .tbar,.fbar{{display:block;background:#0d1729;border-radius:6px;height:12px;overflow:hidden}}
.tbar i,.fbar i{{display:block;height:100%;border-radius:6px;background:var(--teal)}}
.tval{{font-size:12.5px;font-variant-numeric:tabular-nums;color:#cfe0f5}} .tval em{{color:var(--dim);font-style:normal;font-size:11px}}
table{{width:100%;border-collapse:collapse;font-size:13px}} th{{text-align:left;color:var(--dim);font-size:10.5px;text-transform:uppercase;letter-spacing:.09em;padding:6px 8px;border-bottom:1px solid var(--line)}}
td{{padding:7px 8px;border-bottom:1px solid #16233a}} tr:last-child td{{border-bottom:0}}
.mono{{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;color:#d7e6fb}} .num{{text-align:right;font-variant-numeric:tabular-nums}}
.fcell{{display:flex;align-items:center;gap:9px;min-width:150px}} .fbar{{flex:1;height:8px}} .fval{{width:46px;text-align:right;font-variant-numeric:tabular-nums;color:#cfe0f5;font-size:12px}}
img{{width:100%;border-radius:10px;border:1px solid var(--line);background:#fff}}
.note{{font-size:12px;color:var(--dim);margin-top:10px}}
.caveats{{display:grid;grid-template-columns:repeat(3,1fr);gap:12px;margin-top:16px}}
.cav{{background:#0e1830;border:1px solid #1c2b44;border-left:3px solid var(--amber);border-radius:9px;padding:11px 13px;font-size:12.3px;color:#c9d7e9}}
.cav b{{color:#f3d38a;display:block;margin-bottom:2px;font-size:12px}}
footer{{margin-top:28px;color:#5f7592;font-size:11.5px;border-top:1px solid #16233a;padding-top:14px;display:flex;justify-content:space-between;gap:20px;flex-wrap:wrap}}
@media(max-width:900px){{.kpis{{grid-template-columns:repeat(2,1fr)}}.grid{{grid-template-columns:1fr}}.caveats{{grid-template-columns:1fr}}body{{padding:22px}}}}
</style></head><body><div class="wrap">
<h1>HBV Genome <span>Pipeline</span> — Results Dashboard</h1>
<div class="sub">Reference <b>NC_003977.2</b> · MAFFT <b>{html.escape(analysis["alignment"]["version"])}</b> · run <b>{html.escape(analysis["created_at"][:10])}</b></div>

<div class="kpis">
 <div class="kpi"><div class="v">{number(counts["input_sequences"])}</div><div class="l">Input records</div></div>
 <div class="kpi"><div class="v">{number(counts["analyzed_sequences"])}</div><div class="l">Analysed</div></div>
 <div class="kpi hl"><div class="v">{number(counts["mutation_rows"])}</div><div class="l">Variant events</div></div>
 <div class="kpi"><div class="v">{number(stats["n_distinct_mutations"])}</div><div class="l">Distinct mutations</div></div>
 <div class="kpi"><div class="v">{number(stats["n_uncertain_calls_excluded"])}</div><div class="l">Uncertain excluded</div></div>
</div>

<div class="grid">
 <div class="card"><h2>Pipeline funnel</h2><div class="funnel">{funnel_html}</div>
  <div class="note">1,399 non-HBV hepadnaviruses removed in preprocessing · 12,482 records skipped in Task 3 (10,017 below 100 nt, 2,464 below 70% identity, 1 insufficient aligned bases).</div></div>
 <div class="card"><h2>Variant types</h2>{type_rows}
  <div class="note">Share of {number(counts["mutation_rows"])} events. Ambiguous IUPAC calls are reported but excluded from frequencies.</div></div>
</div>

<div class="card" style="margin-bottom:16px"><h2>Top 20 mutations · frequency of analysed sequences</h2>
<table><thead><tr><th>Mutation</th><th>Type</th><th class="num">Sequences</th><th>Frequency</th></tr></thead><tbody>{top_rows}</tbody></table>
<div class="note">Denominator: {number(stats["denominator"])} sequences ({html.escape(stats["denominator_basis"])}).</div></div>

<div class="grid">
 <div class="card"><h2>Most frequent mutations</h2><img src="data:image/png;base64,{top_fig}" alt="Top mutations"></div>
 <div class="card"><h2>Mutations along the genome</h2><img src="data:image/png;base64,{genome_fig}" alt="Mutations along genome">
  <div class="note">Calls per 50-nt bin over the 3,182 bp reference.</div></div>
</div>

<div class="caveats">
 <div class="cav"><b>Partial genomes</b>127,395 of 140,855 records are PARTIAL; terminal unobserved sequence is not called as deletion.</div>
 <div class="cav"><b>Duplicate sequences</b>31,288 records duplicate another accession's sequence; frequencies count records, not unique genomes.</div>
 <div class="cav"><b>No genotype or gene layer</b>Metadata has no genotype column and the table has no codon annotation, so neither comparison is made.</div>
</div>

<footer><span>Generated from results/statistics_summary.json, results/summary.csv and results/mutation_analysis.json.</span>
<span>{number(counts["mutation_rows"])} events · {number(stats["n_distinct_mutations"])} distinct · {number(counts["analyzed_sequences"])} sequences</span></footer>
</div></body></html>"""


def main():
    for name in ("mutation_analysis.json", "statistics_summary.json", "mutation_types.csv", "summary.csv"):
        if not (RESULTS / name).is_file():
            raise SystemExit(f"Missing {RESULTS / name}; run the pipeline first.")
    for name in ("top_mutations.png", "mutations_along_genome.png"):
        if not (RESULTS / "figures" / name).is_file():
            raise SystemExit(f"Missing figure {RESULTS / 'figures' / name}; run Task 4 first.")
    analysis, stats, types, top = gather()
    OUTPUT.write_text(render(analysis, stats, types, top), encoding="utf-8")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
