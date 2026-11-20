#!/usr/bin/env python3
"""5.4 Report generation: read the report_inputs.json produced by 5.3 and render report.html.

Position: 5.3 (risk scoring → results.tsv / species_risk.tsv / report_inputs.json)
          → this script (render the scores into a self-contained HTML report).

Input:
  - report_inputs.json (step-5.3 output; rows / species / relations / qc_dir)
Output:
  - report.html (self-contained HTML: QC (fastp links) → scoring → rationale (KG + recommendations))

Usage:
  python 5.4-report.py --bundle report_inputs.json
  python 5.4-report.py --config config.yaml          # reads results_dir/00_results/report_inputs.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path


def _esc(s) -> str:
    """HTML-escape to prevent special characters from breaking the page."""
    return (str(s)
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
            .replace('"', "&quot;"))


REL_CN = {
    "virus_detected_in_host": "Virus detected in host",
    "virus_has_gene_or_segment": "Virus has gene/segment",
    "virus_coinfects_with": "Coinfection",
    "virus_related_to": "Related virus",
    "virus_reassortment": "Reassortment",
    "virus_host_switching": "Host switching",
    "virus_associated_with_host": "Virus associated with host",
    "variant_reported": "Variant reported",
}


# ---------- Circular bar chart (one per L-segment sequence: inner = 3 detection items, outer = 8 species-risk items) ----------
CBAR_INNER = [
    ("DIAMOND", "#2166ac"),
    ("Tree", "#67a9cf"),
    ("LucaProt", "#d1e5f0"),
]
# Outer ring 8 items: Host (orange) + Phylogeny (purple) + 4 literature items (green) + Human infection (red)
CBAR_OUTER = [
    ("Host type", "host_type_score", 10.0, "#e08214"),
    ("Host breadth", "host_breadth_score", 5.0, "#f4a582"),
    ("Phylogeny", "phylogeny_score", 15.0, "#7570b3"),
    ("Reassortment", "reassortment", 10.0, "#1b7837"),
    ("Host switching", "host_switching", 10.0, "#5aae61"),
    ("Coinfection", "coinfection", 5.0, "#a6dba0"),
    ("Related", "related_viruses", 5.0, "#d9f0d3"),
    ("Human infection", "human_infection_score", 40.0, "#d73027"),
]

# 4 literature traits -> (relation type, whether direct_target_match is required), kept in sync with 5.3 literature_traits
LIT_TRAIT = {
    "reassortment": ("virus_reassortment", False),
    "host_switching": ("virus_host_switching", False),
    "coinfection": ("virus_coinfects_with", True),
    "related_viruses": ("virus_related_to", True),
}
LIT_TRAIT_CN = {
    "reassortment": "Reassortment",
    "host_switching": "Host switching",
    "coinfection": "Coinfection",
    "related_viruses": "Related",
}
LIT_PLACEHOLDER = {"", "related virus", "co-infecting virus", "unspecified variant"}
MISS_COLOR = "#d9d9d9"

# Host annotation source (host_source field from 5.3)
HOST_SOURCE_CN = {"VMR": "VMR annotation", "literature": "Literature search", "none": "None"}

# The upstream host_desc values are a mix of a Chinese family name and a Latin
# name. The report shows the Latin name only. These keys are match patterns, not
# display strings, so the Chinese fragments must stay verbatim.
# Example input: "\u9f20\u79d1 Muridae" -> "Muridae".
HOST_DESC_EN = {
    "鼠科 Muridae": "Muridae",
    "鼩鼱科 Soricidae": "Soricidae",
    "仓鼠科 Cricetidae": "Cricetidae",
    "蝙蝠科 Vespertilionidae": "Vespertilionidae",
    "菊头蝠科 Rhinolophidae": "Rhinolophidae",
    "刺猬科 Erinaceidae": "Erinaceidae",
    "鸟类": "Bird",
    "昆虫": "Insect",
}


def _en_host(desc) -> str:
    d = str(desc or "").strip()
    return HOST_DESC_EN.get(d, d)


_cbar_uid = 0


def _num(v, default=0.0):
    try:
        if v is None or v == "":
            return default
        return float(v)
    except (TypeError, ValueError):
        return default


def _lit_evidence(sp: str, rels_by_sp: dict) -> dict:
    """Supporting relations (PMID + evidence sentence) for the 4 literature traits of one species."""
    out = {t: [] for t in LIT_TRAIT}
    for rel in rels_by_sp.get(sp, []):
        relt = rel.get("relation")
        for trait, (want, need_direct) in LIT_TRAIT.items():
            if relt != want:
                continue
            if rel.get("polarity") != "positive":
                continue
            if need_direct and rel.get("direct_target_match") is not True:
                continue
            if trait == "related_viruses" and str(rel.get("object") or "") in LIT_PLACEHOLDER:
                continue
            out[trait].append(rel)
    return out


def _lit_subscore(sr: dict, trait: str, maxv: float) -> float:
    """Literature trait score: full mark if hit, 0 otherwise."""
    detail = sr.get("literature_detail") or {}
    return maxv if detail.get(trait) is True else 0.0


def _hover_builder(o: dict, sr: dict, lit_ev: dict, full: dict) -> dict:
    """Generate hover HTML for each bar (raw alignment / tree / host / literature evidence)."""
    def kv(label, v, suffix=""):
        s = str(v or "").strip()
        return f"{label}: {_esc(s)}{suffix}" if s else None

    def pct(score, maxv):
        return round(min(1.0, max(0.0, score / maxv)) * 100, 1) if maxv else 0.0

    h = {}

    # ---- Inner ring (detection confidence) ----
    d = _num(o.get("detection_diamond"))
    dlines = [f"Score: {d} / 1 ({pct(d, 1)}%)"]
    for ln in [kv("Reference", full.get("sseqid")),
               kv("Identity (pident)", full.get("pident"), "%"),
               kv("Alignment length", full.get("length"), " bp"),
               kv("evalue", full.get("evalue")),
               kv("bitscore", full.get("bitscore"))]:
        if ln:
            dlines.append(ln)
    if full.get("qstart") and full.get("qend"):
        dlines.append(f"Query position: {full.get('qstart')}-{full.get('qend')}")
    if full.get("sstart") and full.get("send"):
        dlines.append(f"Reference position: {full.get('sstart')}-{full.get('send')}")
    h["DIAMOND"] = "<br/>".join(dlines)

    t = _num(o.get("detection_tree"))
    tlines = [f"Score: {t} / 1 ({pct(t, 1)}%)"]
    for ln in [kv("Tree position", full.get("tree_position")),
               kv("bootstrap", full.get("tree_bootstrap")),
               kv("Nearest neighbour", full.get("tree_nn_sseqid")),
               kv("Distance", full.get("tree_nn_distance"))]:
        if ln:
            tlines.append(ln)
    h["Tree"] = "<br/>".join(tlines)

    l = _num(o.get("detection_lucaprot"))
    h["LucaProt"] = f"Score: {l} / 1 (RdRp probability, {pct(l, 1)}%)"

    # ---- Outer ring (species risk) ----
    ht = _num(sr.get("host_type_score"))
    hlines = [f"Score: {ht} / 10 ({pct(ht, 10)}%)"]
    for ln in [kv("Host", full.get("host")), kv("Class", sr.get("host_class"))]:
        if ln:
            hlines.append(ln)
    h["Host type"] = "<br/>".join(hlines)

    br = sr.get("host_breadth") or 0
    h["Host breadth"] = f"Score: {br} / 5 ({pct(br, 5)}%)<br/>Host species count: {br}"

    phy = _num(sr.get("phylogeny_score"))
    plines = [f"Score: {phy} / 15 ({pct(phy, 15)}%)"]
    for ln in [kv("Genus", full.get("genus")), kv("Nearest-neighbour distance", full.get("tree_nn_distance"))]:
        if ln:
            plines.append(ln)
    h["Phylogeny"] = "<br/>".join(plines)

    for trait, cn, maxv in [("reassortment", "Reassortment", 10), ("host_switching", "Host switching", 10),
                            ("coinfection", "Coinfection", 5), ("related_viruses", "Related", 5)]:
        hit = bool((sr.get("literature_detail") or {}).get(trait))
        if not hit:
            h[cn] = f"Not detected: 0 / {maxv}"
            continue
        rels = lit_ev.get(trait) or []
        pmids = list(dict.fromkeys(str(r.get("pmid")) for r in rels if r.get("pmid")))[:3]
        sents = [str(r.get("evidence_sentence") or "").strip()
                 for r in rels if r.get("evidence_sentence")]
        lines = [f"Score: {maxv} / {maxv} (100%)"]
        if pmids:
            lines.append("Source PMID: " + ", ".join(_esc(p) for p in pmids))
        if sents:
            lines.append("Evidence sentence: " + _esc(sents[0][:120]))
        h[cn] = "<br/>".join(lines)

    if sr.get("human_infection_vmr"):
        h["Human infection"] = "Score: 40 / 40 (100%)<br/>Source: VMR official human isolate"
    elif sr.get("human_infection_literature"):
        hpm = [str(p) for p in (sr.get("human_pmids") or [])]
        lines = ["Score: 40 / 40 (100%)"]
        if hpm:
            lines.append("Source PMID: " + ", ".join(_esc(p) for p in hpm))
        h["Human infection"] = "<br/>".join(lines)
    else:
        h["Human infection"] = "Not detected: 0 / 40"

    return h


def _cbar_option(o: dict, sr: dict, lit_ev: dict, full: dict) -> dict:
    """ECharts option for the two-layer circular bar chart of one sequence (inner = 3 detection items, outer = 8 species-risk items)."""
    hover = _hover_builder(o, sr, lit_ev, full)
    inner = [
        ("DIAMOND", _num(o.get("detection_diamond")), 1.0, "#2166ac", None),
        ("Tree", _num(o.get("detection_tree")), 1.0, "#67a9cf", None),
        ("LucaProt", _num(o.get("detection_lucaprot")), 1.0, "#d1e5f0", None),
    ]
    outer_raw = [
        ("Host type", _num(sr.get("host_type_score")), 10.0, "#e08214", None),
        ("Host breadth", _num(sr.get("host_breadth_score")), 5.0, "#f4a582", None),
        ("Phylogeny", _num(sr.get("phylogeny_score")), 15.0, "#7570b3", None),
        ("Reassortment", _lit_subscore(sr, "reassortment", 10.0), 10.0, "#1b7837", "reassortment"),
        ("Host switching", _lit_subscore(sr, "host_switching", 10.0), 10.0, "#5aae61", "host_switching"),
        ("Coinfection", _lit_subscore(sr, "coinfection", 5.0), 5.0, "#a6dba0", "coinfection"),
        ("Related", _lit_subscore(sr, "related_viruses", 5.0), 5.0, "#d9f0d3", "related_viruses"),
        ("Human infection", _num(sr.get("human_infection_score")), 40.0, "#d73027", "human"),
    ]

    def ring(items):
        labels = [it[0] for it in items]
        data = []
        for name, val, maxv, color, trait in items:
            prop = min(1.0, max(0.0, val / maxv)) if maxv else 0.0
            item = {"value": round(prop, 4), "itemStyle": {"color": color},
                    "hover": hover.get(name, name)}
            if trait is not None:
                hit = (bool(sr.get("human_infection_vmr") or sr.get("human_infection_literature"))
                       if trait == "human"
                       else bool((sr.get("literature_detail") or {}).get(trait)))
                if not hit:
                    item["value"] = 0.04
                    item["itemStyle"] = {"color": MISS_COLOR}
            data.append(item)
        return labels, data

    inner_labels, inner_data = ring(inner)
    outer_labels, outer_data = ring(outer_raw)
    angle_off = {"axisLabel": {"show": False}, "axisLine": {"show": False},
                 "axisTick": {"show": False}, "splitLine": {"show": False}}
    radius_off = {"axisLabel": {"show": False}, "axisLine": {"show": False},
                  "axisTick": {"show": False},
                  "splitLine": {"show": True, "lineStyle": {"color": "#eeeeee"}}}
    return {
        "animation": False,
        "tooltip": {"trigger": "item", "formatter": "@@CBARTIP@@", "confine": True, "appendToBody": True,
                    "extraCssText": "max-width:340px;white-space:normal;box-shadow:0 6px 20px rgba(0,0,0,.14);border-radius:8px;padding:8px 10px;"},
        "polar": [
            {"center": ["50%", "50%"], "radius": ["26%", "52%"]},
            {"center": ["50%", "50%"], "radius": ["60%", "92%"]},
        ],
        "angleAxis": [
            {"type": "category", "data": inner_labels, "polarIndex": 0, "startAngle": 90, **angle_off},
            {"type": "category", "data": outer_labels, "polarIndex": 1, "startAngle": 90, **angle_off},
        ],
        "radiusAxis": [
            {"type": "value", "min": 0, "max": 1, "polarIndex": 0, **radius_off},
            {"type": "value", "min": 0, "max": 1, "polarIndex": 1, **radius_off},
        ],
        "series": [
            {"type": "bar", "coordinateSystem": "polar", "polarIndex": 0,
             "data": inner_data, "barWidth": 75},
            {"type": "bar", "coordinateSystem": "polar", "polarIndex": 1,
             "data": outer_data, "barWidth": 30,
             "label": {"show": True, "position": "outside", "color": "#444", "fontSize": 10}},
        ],
    }


def _cbar_html(o: dict, sr: dict, esc, lit_ev: dict, full: dict) -> str:
    """Render the circular-bar-chart card for one sequence (div + centre final score + click jump + ECharts script)."""
    global _cbar_uid
    _cbar_uid += 1
    uid = _cbar_uid
    tier = str(o.get("tier") or "Low")
    key = f"{o.get('sample', '')}_{o.get('id', '')}"
    opt = json.dumps(_cbar_option(o, sr, lit_ev, full), ensure_ascii=False).replace('"@@CBARTIP@@"', "cbarTip")
    head = (
        f'<div class="cbar-head"><b>{esc(o.get("id", ""))}</b> '
        f'<span class="tag {esc(tier)}">{esc(tier)}</span>'
        f'<div class="cbar-sub">{esc(o.get("sample", ""))} · L segment</div></div>'
    )
    center = (
        f'<div class="cbar-center"><div class="cbar-score">{esc(o.get("final_score", ""))}</div>'
        f'<div class="cbar-score-label">Final score</div></div>'
    )
    links = (
        f'<div class="cbar-links">'
        f'<a href="#seq-{esc(key)}">Section 5 detail →</a> '
        f'<a href="#kg-{esc(key)}">Section 6 KG →</a></div>'
    )
    script = (
        f'<script>var c=echarts.init(document.getElementById("cbar_{uid}"));'
        f'c.setOption({opt});'
        f'c.on("click",function(){{location.hash="seq-{key}";}});</script>'
    )
    return (
        f'<div class="cbar-card" id="chart-{esc(key)}">{head}'
        f'<div class="cbar-wrap"><div id="cbar_{uid}" style="width:250px;height:250px;"></div>{center}</div>'
        f'{links}'
        f'{script}'
        f'</div>'
    )


def write_report_html(path: Path, rows_out: list[dict], species_risk: dict,
                      relations: list[dict] | None, qc_dir: Path | None,
                      rows_full: list[dict] | None = None,
                      species_lit: dict | None = None) -> None:
    """Generate a self-contained HTML report: QC (fastp links) → scoring → rationale (KG + recommendations)."""
    esc = _esc

    # Index full rows (raw alignment / tree columns) by sample_id for hover detail
    full_by_key: dict[str, dict] = {}
    for f in (rows_full or []):
        key = f"{f.get('sample', '')}_{f.get('id', '')}"
        if key:
            full_by_key[key] = f

    # ---- Group KG relations by species (record_id is sample_contig, matching the row key) ----
    key2sp: dict[str, str] = {}
    for o in rows_out:
        key = f"{o.get('sample', '')}_{o.get('id', '')}"
        if key and key not in key2sp:
            key2sp[key] = o.get("species", "Unclassified")
    rels_by_sp: dict[str, list[dict]] = defaultdict(list)
    for rel in (relations or []):
        sp = key2sp.get(str(rel.get("record_id") or ""))
        if sp:
            rels_by_sp[sp].append(rel)
    # Literature evidence (per species): supporting relations for the 4 traits (PMID + evidence sentence)
    lit_ev_by_sp = {sp: _lit_evidence(sp, rels_by_sp) for sp in species_risk}
    # Contigs needing a Section-6 anchor: contigs with relations + all L-segment contigs
    rel_keys = {str(rel.get("record_id")) for rel in (relations or [])}
    l_keys = {f"{o.get('sample', '')}_{o.get('id', '')}" for o in rows_out
              if str(o.get("segment", "")).strip().lower() == "l"}
    kg_keys = rel_keys | l_keys

    # ---- QC samples (deduplicated in result-table order; only samples with detections) ----
    sample_order: list[str] = []
    seen: set[str] = set()
    for o in rows_out:
        s = (o.get("sample") or "").strip()
        if s and s not in seen:
            seen.add(s)
            sample_order.append(s)
    qc_rows: list[dict] = []
    for s in sample_order:
        exists = False
        href = ""
        if qc_dir:
            html = qc_dir / f"{s}_fastp.html"
            if html.is_file():
                exists = True
                try:
                    href = os.path.relpath(html, path.parent)
                except ValueError:
                    href = str(html)
        cnt = sum(1 for o in rows_out if (o.get("sample") or "").strip() == s)
        qc_rows.append({"sample": s, "exists": exists, "href": href, "count": cnt})

    # ---- Basic statistics ----
    n_rows = len(rows_out)
    n_species = len(species_risk)
    tier_counts: dict[str, int] = defaultdict(int)
    for sr in species_risk.values():
        tier_counts[sr.get("tier", "Low")] += 1

    by_sp: dict[str, list[dict]] = defaultdict(list)
    for o in rows_out:
        by_sp[o.get("species", "Unclassified")].append(o)
    ordered = sorted(species_risk.values(), key=lambda x: x.get("species_risk_score", 0.0), reverse=True)

    def tag(tier: str) -> str:
        return f'<span class="tag {esc(tier)}">{esc(tier)}</span>'

    def _cap(x: str) -> str:
        return (x[0].upper() + x[1:]) if x else x

    def format_host_cell(sr: dict) -> str:
        """Host-type cell of table 03: class name + specific species (italic binomials) + human-isolate marker.

        - show full binomials only; fall back to genus/generic names when needed
        - list up to 2 in full, collapse 3+ to "first et al. (N species)"
        """
        out = [esc(_en_host(sr.get("host_desc", "")) or "—")]
        species = [s for s in str(sr.get("host_species", "") or "").split("; ") if s]
        full = [_cap(s) for s in species if " " in s]
        short = [_cap(s) for s in species if " " not in s]
        shown = full if full else short
        if shown:
            if len(shown) == 1:
                txt = f"<i>{esc(shown[0])}</i>"
            elif len(shown) == 2:
                txt = ", ".join(f"<i>{esc(s)}</i>" for s in shown)
            else:
                txt = f"<i>{esc(shown[0])}</i> et al. ({len(species)} species)"
            out.append(txt)
        if str(sr.get("human_isolation", "")).lower() in ("true", "1"):
            out.append('<span class="note">(human isolate)</span>')
        return "<br>".join(out)

    L: list[str] = []
    L.append("<!DOCTYPE html>")
    L.append('<html lang="en">')
    L.append("<head>")
    L.append('<meta charset="utf-8">')
    L.append('<meta name="viewport" content="width=device-width, initial-scale=1">')
    L.append("<title>Hantavirus Detection and Risk Scoring Report</title>")
    L.append("<style>")
    L.append(":root{--accent:#2166ac;--text:#1d1d1f;--muted:#6e6e73;--border:#e5e5e7;}")
    L.append("html{scroll-behavior:smooth;}")
    L.append("*{box-sizing:border-box;}")
    L.append("body{font-family:Arial,'Helvetica Neue','PingFang SC','Microsoft YaHei',sans-serif;background:#fff;color:var(--text);margin:0;line-height:1.6;-webkit-font-smoothing:antialiased;}")
    L.append("a{color:var(--accent);text-decoration:none;}")
    L.append("a:hover{text-decoration:underline;}")
    L.append(".container{max-width:1200px;margin:0 auto;padding:0 28px;}")
    L.append(".topnav{position:fixed;top:0;left:0;right:0;z-index:100;background:rgba(255,255,255,.72);backdrop-filter:saturate(180%) blur(20px);-webkit-backdrop-filter:saturate(180%) blur(20px);border-bottom:1px solid rgba(0,0,0,.08);}")
    L.append(".topnav-inner{max-width:1200px;margin:0 auto;display:flex;align-items:center;justify-content:space-between;height:52px;padding:0 28px;}")
    L.append(".topnav-title{font-size:14px;font-weight:700;color:var(--text);}")
    L.append(".topnav-links{display:flex;gap:20px;}")
    L.append(".topnav-links a{font-size:12px;color:#1d1d1f;font-weight:500;transition:color .2s;}")
    L.append(".topnav-links a:hover{color:var(--accent);text-decoration:none;}")
    L.append(".hero{text-align:center;padding:150px 24px 96px;}")
    L.append(".hero-badge{display:inline-block;font-size:12px;color:var(--accent);font-weight:700;letter-spacing:.08em;text-transform:uppercase;margin-bottom:16px;}")
    L.append(".hero h1{font-size:54px;font-weight:800;letter-spacing:-.02em;margin:0 0 16px;line-height:1.08;}")
    L.append(".hero-sub{font-size:20px;color:var(--muted);font-weight:400;margin:0 0 10px;}")
    L.append(".hero-meta{font-size:13px;color:#aaa;margin:0;}")
    L.append(".section{padding:72px 0 36px;}")
    L.append(".section .sec-no{color:#c7c7cc;font-weight:800;margin-right:12px;}")
    L.append("h2{font-size:32px;font-weight:800;letter-spacing:-.01em;margin:0 0 14px;}")
    L.append("h3{font-size:19px;margin:34px 0 12px;font-weight:700;}")
    L.append("h4{font-size:14px;margin:22px 0 6px;font-weight:700;color:#333;}")
    L.append(".sec-desc{color:var(--muted);font-size:15px;max-width:780px;margin:0 0 26px;}")
    L.append("table{border-collapse:collapse;width:100%;margin:16px 0;font-size:13px;}")
    L.append("th,td{padding:10px 12px;text-align:left;vertical-align:top;border-bottom:1px solid var(--border);}")
    L.append("th{font-weight:600;color:#6e6e73;font-size:12px;background:#fafafa;}")
    L.append("tr:hover td{background:#fbfbfd;}")
    L.append("ul{margin:8px 0 8px 22px;padding:0;}")
    L.append(".tag{display:inline-block;padding:2px 10px;border-radius:20px;font-size:11px;font-weight:600;color:#fff;line-height:1.5;}")
    L.append(".High{background:#d73027;}.Medium{background:#f0a000;}.Low{background:#4daf4a;}")
    L.append(".note{color:var(--muted);font-size:12px;}")
    L.append(".evidence{margin:6px 0;padding:10px 14px;background:#fafafa;border-left:3px solid #d0d0d5;border-radius:8px;font-size:12px;}")
    L.append(".cbar-grid{display:flex;flex-wrap:wrap;gap:18px;margin:18px 0;}")
    L.append(".cbar-card{border:1px solid var(--border);border-radius:18px;padding:14px;width:286px;background:#fff;box-shadow:0 2px 14px rgba(0,0,0,.05);transition:transform .25s ease,box-shadow .25s ease;}")
    L.append(".cbar-card:hover{transform:translateY(-4px);box-shadow:0 12px 32px rgba(0,0,0,.10);}")
    L.append(".cbar-head{font-size:13px;word-break:break-all;margin-bottom:4px;}")
    L.append(".cbar-sub{color:var(--muted);font-size:11px;}")
    L.append(".cbar-wrap{position:relative;width:250px;height:250px;margin:6px auto 0;}")
    L.append(".cbar-center{position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);text-align:center;pointer-events:none;line-height:1.15;}")
    L.append(".cbar-score{font-size:26px;font-weight:800;color:var(--text);}")
    L.append(".cbar-score-label{font-size:10px;color:var(--muted);}")
    L.append(".cbar-links{font-size:11px;text-align:center;margin-top:4px;}")
    L.append(".cbar-links a{margin:0 4px;}")
    L.append("details{margin:12px 0;border:1px solid var(--border);border-radius:14px;padding:4px 20px 12px;}")
    L.append("summary{cursor:pointer;font-size:20px;font-weight:700;padding:12px 0;outline:none;}")
    L.append("tr:target{background:#fff3cd;}")
    L.append(".cbar-card:target{border:1px solid #e7298a;box-shadow:0 0 0 3px rgba(231,41,138,.22);}")
    L.append(".fade{opacity:0;transform:translateY(24px);transition:opacity .7s ease,transform .7s ease;}")
    L.append(".fade.fade-in{opacity:1;transform:none;}")
    L.append(".totop{position:fixed;right:28px;bottom:28px;width:48px;height:48px;border-radius:50%;border:none;background:rgba(29,29,31,.85);color:#fff;font-size:20px;cursor:pointer;opacity:0;pointer-events:none;transition:opacity .3s,transform .3s;z-index:99;backdrop-filter:blur(10px);}")
    L.append(".totop.show{opacity:1;pointer-events:auto;}")
    L.append(".totop:hover{transform:translateY(-3px);background:#1d1d1f;}")
    L.append("@media (prefers-reduced-motion: reduce){.fade{opacity:1;transform:none;transition:none;}html{scroll-behavior:auto;}}")
    L.append("@media (max-width:760px){.hero h1{font-size:36px;}.topnav-links{display:none;}}")
    L.append("</style>")
    L.append('<script src="https://cdn.jsdelivr.net/npm/echarts@5/dist/echarts.min.js"></script>')
    L.append("<script>")
    L.append("function cbarTip(params){var d=params.data||{};return d.hover||params.name;}")
    L.append("</script>")
    L.append("</head>")
    L.append("<body>")

    # Top navigation
    L.append('<nav class="topnav"><div class="topnav-inner">'
             '<div class="topnav-title">Hantavirus Risk Scoring</div>'
             '<div class="topnav-links">'
             '<a href="#s0">Contents</a><a href="#s1">QC</a><a href="#s2">Summary</a><a href="#s3">Overview</a>'
             '<a href="#s4">Details</a><a href="#s5">Evidence</a><a href="#s6">Assessment</a>'
             '</div></div></nav>')

    # Hero
    L.append('<header class="hero">')
    L.append('<div class="hero-badge">Hantavirus Risk Assessment</div>')
    L.append('<h1>Hantavirus Detection and Risk<br>Scoring Report</h1>')
    L.append(f'<p class="hero-sub">{n_rows} detected sequences · {n_species} species · detection confidence × species risk</p>')
    L.append(f'<p class="hero-meta">High {tier_counts.get("High", 0)} | Medium {tier_counts.get("Medium", 0)} | Low {tier_counts.get("Low", 0)} | Generated on {date.today()}</p>')
    L.append('</header>')
    L.append('<main class="container">')

    # 0 Contents
    L.append('<div class="section fade" id="s0">')
    L.append('<h2><span class="sec-no">00</span>Contents</h2>')
    L.append('<p class="note">Description of the final result files under 00_results; start with report.html (this report).</p>')
    L.append("<table><tr><th>File</th><th>Description</th></tr>")
    result_desc = {
        "all_vmr_final_classified.tsv": "Combined summary table",
        "l.fasta": "Extracted hantavirus L-segment protein sequences",
        "l.fasta.csv": "LucaProt predictions (id + probability)",
        "RdRPs_only_using_threshold0.500000.csv": "LucaProt results: RdRps above the 0.5 threshold",
        "results.tsv": "Final results",
        "species_risk.tsv": "Species risk scores (High/Medium/Low)",
        "report_inputs.json": "Intermediate file for step 5.3 → 5.4 (not needed)",
        "report.html": "Result report (read first)",
    }
    existing_files = {p.name for p in path.parent.iterdir() if p.is_file()}
    ordered_names = [n for n in result_desc if n in existing_files]
    ordered_names += sorted(existing_files - set(result_desc))
    for name in ordered_names:
        desc = result_desc.get(name, "Other output files")
        L.append(f"<tr><td>{esc(name)}</td><td>{esc(desc)}</td></tr>")
    L.append("</table>")
    L.append("</div>")

    # 1 Summary: tier counts + high-risk species table
    L.append('<div class="section fade" id="s1">')
    L.append('<h2><span class="sec-no">01</span>Summary</h2>')
    L.append(f"<p>Species tiers: High {tier_counts.get('High', 0)}, Medium {tier_counts.get('Medium', 0)}, Low {tier_counts.get('Low', 0)}.</p>")
    L.append("<table><tr><th>High-risk species</th><th>Contigs</th><th>Positive samples</th><th>Representative samples</th></tr>")
    for sr in [x for x in ordered if x.get("tier") == "High"]:
        sp = sr.get("species", "")
        rows_sp = by_sp.get(sp, [])
        samp_cnt = defaultdict(int)
        for o in rows_sp:
            s_ = (o.get("sample") or "").strip()
            if s_:
                samp_cnt[s_] += 1
        top2 = ", ".join(k for k, _ in sorted(samp_cnt.items(), key=lambda kv: -kv[1])[:2])
        L.append(f"<tr><td>{esc(sp)}</td><td>{len(rows_sp)}</td><td>{len(samp_cnt)}</td><td>{esc(top2)}</td></tr>")
    L.append("</table>")
    L.append("</div>")

    # 2 Quality control
    L.append('<div class="section fade" id="s2">')
    L.append('<h2><span class="sec-no">02</span>Quality control</h2>')
    L.append('<p class="note">Only samples with hantavirus detections are listed; click a sample name to open its fastp QC report.</p>')
    L.append("<table><tr><th>Sample</th><th>QC report</th><th>Detected contigs</th></tr>")
    for q in qc_rows:
        name = esc(q["sample"])
        if q["exists"]:
            cell = f'<a href="{esc(q["href"])}" target="_blank">{name}</a>'
        else:
            cell = f'{name}<span class="note">(no QC report)</span>'
        L.append(f"<tr><td>{name}</td><td>{cell}</td><td>{q['count']}</td></tr>")
    L.append("</table>")
    L.append("</div>")

    # 3 Species-level verdict
    L.append('<div class="section fade" id="s3">')
    L.append('<h2><span class="sec-no">03</span>Species-level verdict</h2>')
    L.append("<table><tr><th>Species</th><th>Virus</th><th>Sequences</th><th>Host type</th><th>Species risk</th><th>Verdict</th></tr>")
    for sr in ordered:
        sp = sr.get("species", "")
        rows_sp = by_sp.get(sp, [])
        virus_name = next((o.get("virus_name", "") for o in rows_sp if o.get("virus_name")), "")
        L.append(
            f"<tr><td>{esc(sp)}</td><td>{esc(virus_name)}</td><td>{len(rows_sp)}</td>"
            f"<td>{format_host_cell(sr)}</td><td>{sr.get('species_risk_score', '')}</td>"
            f"<td>{tag(sr.get('tier', 'Low'))}</td></tr>"
        )
    L.append("</table>")
    L.append("</div>")

    # 4 Species-level assessment details
    L.append('<div class="section fade" id="s4">')
    L.append('<h2><span class="sec-no">04</span>Species-level assessment details</h2>')
    L.append('<div class="note" style="margin:8px 0;">Circular bar charts: inner ring = detection confidence ('
             '<span style="color:#2166ac;">●DIAMOND</span> '
             '<span style="color:#67a9cf;">●Tree</span> '
             '<span style="color:#d1e5f0;">●LucaProt</span>); outer ring = species risk: host ('
             '<span style="color:#e08214;">●Host type</span> '
             '<span style="color:#f4a582;">●Host breadth</span>), phylogeny ('
             '<span style="color:#7570b3;">●Phylogeny</span>), literature ('
             '<span style="color:#1b7837;">●Reassortment</span> '
             '<span style="color:#5aae61;">●Host switching</span> '
             '<span style="color:#a6dba0;">●Coinfection</span> '
             '<span style="color:#d9f0d3;">●Related</span>), human infection ('
             '<span style="color:#d73027;">●Human infection</span>). Bar length = item score ÷ full mark; missed items are drawn grey; hover for values and PMIDs.</div>')
    for idx, sr in enumerate(ordered, 1):
        sp = sr.get("species", "")
        rows_sp = by_sp.get(sp, [])
        virus_name = next((o.get("virus_name", "") for o in rows_sp if o.get("virus_name")), "")
        pmids = "; ".join(sr.get("human_pmids", [])) or "None"
        if sr.get("human_infection_vmr"):
            human_src = "Yes (VMR human isolate)"
        elif sr.get("human_infection_literature"):
            human_src = f"Yes (literature report, PMID: {pmids})"
        else:
            human_src = "None"
        L.append(f"<h3>4.{idx} {esc(sp)}</h3>")
        L.append("<ul>")
        if virus_name:
            L.append(f"<li>Virus name: {esc(virus_name)}</li>")
        L.append(f"<li>Host: {esc(_en_host(sr.get('host_desc', '')))}(class {esc(sr.get('host_class', ''))}, "
                 f"source {esc(HOST_SOURCE_CN.get(str(sr.get('host_source', '')), str(sr.get('host_source', '')) or '-'))})</li>")
        L.append(f"<li>Natural host: {esc(sr.get('host_species', '') or '—')}</li>")
        human_iso_txt = "Yes" if str(sr.get("human_isolation", "")).lower() in ("true", "1") else "No"
        L.append(f"<li>Human isolate: {human_iso_txt}</li>")
        L.append(f"<li>Detected sequences: {len(rows_sp)}</li>")
        L.append(
            f"<li>Components: Host type {sr.get('host_type_score', 0)} + Host breadth {sr.get('host_breadth_score', 0)} + "
            f"Phylogeny {sr.get('phylogeny_score', 0)} + Literature {sr.get('literature_score', 0)} + "
            f"Human infection {sr.get('human_infection_score', 0)}</li>"
        )
        L.append(f"<li>Species risk score: <b>{sr.get('species_risk_score', '')}</b> → Verdict: {tag(sr.get('tier', 'Low'))}</li>")
        L.append(f"<li>Human-infection evidence: {esc(human_src)}</li>")
        L.append("</ul>")

        # Literature / human-infection evidence sources (species-level, PMIDs clickable)
        lit_ev = lit_ev_by_sp.get(sp, {})
        src_parts = []
        for trait, cn in LIT_TRAIT_CN.items():
            if (sr.get("literature_detail") or {}).get(trait) is True:
                pms = sorted({str(r.get("pmid")) for r in lit_ev.get(trait, []) if r.get("pmid")})
                if pms:
                    links = ", ".join(
                        f'<a href="https://pubmed.ncbi.nlm.nih.gov/{esc(p)}/" target="_blank">PMID:{esc(p)}</a>'
                        for p in pms)
                    src_parts.append(f"{cn}({links})")
                else:
                    src_parts.append(cn)
        if sr.get("human_infection_vmr"):
            src_parts.append("Human infection (VMR official human isolate)")
        elif sr.get("human_infection_literature") and sr.get("human_pmids"):
            links = ", ".join(
                f'<a href="https://pubmed.ncbi.nlm.nih.gov/{esc(p)}/" target="_blank">PMID:{esc(p)}</a>'
                for p in sr.get("human_pmids", []))
            src_parts.append(f"Human infection ({links})")
        if src_parts:
            L.append(f'<p class="note">Literature / human-infection evidence sources: {"; ".join(src_parts)}</p>')
        else:
            L.append('<p class="note">Literature / human-infection evidence sources: None.</p>')

        l_rows = [o for o in rows_sp if str(o.get("segment", "")).strip().lower() == "l"]
        if l_rows:
            L.append('<p class="note">Scoring sources for the L-segment sequences of this species (hover for evidence and PMIDs; click to jump to Section 5/6):</p>')
            L.append('<div class="cbar-grid">')
            for o in l_rows:
                key = f"{o.get('sample', '')}_{o.get('id', '')}"
                L.append(_cbar_html(o, sr, esc, lit_ev, full_by_key.get(key, {})))
            L.append('</div>')
    L.append("</div>")

    # 5 Per-contig evidence (expanded by default, collapsible; chart click jumps here and highlights the row)
    L.append('<div class="section fade" id="s5">')
    L.append('<h2><span class="sec-no">05</span>Per-contig evidence</h2>')
    L.append("<details open>")
    L.append("<summary>Click to collapse / expand the table</summary>")
    L.append("<table><tr><th>Sample</th><th>Contig ID</th><th>Segment</th><th>DIAMOND</th><th>Tree</th><th>LucaProt</th><th>Detection confidence</th><th>Species risk</th><th>Final score</th><th>Verdict</th><th>Jump</th></tr>")
    for o in rows_out:
        seg = (o.get("segment") or "").strip().upper() or "-"
        d_tree = o.get("detection_tree") or "-"
        d_luca = o.get("detection_lucaprot") or "-"
        key = f"{o.get('sample', '')}_{o.get('id', '')}"
        kg_link = f'<a href="#kg-{esc(key)}">KG →</a>' if key in kg_keys else "—"
        chart_link = f'<a href="#chart-{esc(key)}">Back to chart →</a>' if key in l_keys else "—"
        L.append(
            f'<tr id="seq-{esc(key)}">'
            f"<td>{esc(o.get('sample', ''))}</td><td>{esc(o.get('id', ''))}</td><td>{esc(seg)}</td>"
            f"<td>{o.get('detection_diamond', '')}</td><td>{d_tree}</td><td>{d_luca}</td>"
            f"<td>{o.get('detection_confidence', '')}</td><td>{o.get('species_risk_score', '')}</td>"
            f"<td>{o.get('final_score', '')}</td><td>{tag(o.get('tier', 'Low'))}</td>"
            f"<td>{chart_link} {kg_link}</td></tr>"
        )
    L.append("</table>")
    L.append("</details>")
    L.append("</div>")

    # 6 Evidence and recommendations (KG + score components + recommendation; relations split by contig, clickable anchors)
    L.append('<div class="section fade" id="s6">')
    L.append('<h2><span class="sec-no">06</span>Evidence and recommendations</h2>')
    L.append('<p class="note">For each species, the score components, recommendation, and contig-level knowledge-graph literature evidence (with PMIDs and evidence sentences) are listed below.</p>')

    def render_rel(rel):
        relt = rel.get("relation") or "relation"
        rel_cn = REL_CN.get(relt, relt)
        subject = rel.get("subject") or ""
        obj = rel.get("object") or ""
        pmid = rel.get("pmid") or ""
        sent = (rel.get("evidence_sentence") or "").strip()
        direct = rel.get("direct_target_match")
        direct_note = "Direct evidence" if direct is True else "Context/related"
        pmid_html = (f'<a href="https://pubmed.ncbi.nlm.nih.gov/{esc(pmid)}/" target="_blank">PMID:{esc(pmid)}</a>'
                     if pmid else "—")
        obj_part = f" → {esc(obj)}" if obj else ""
        out = [f'<div class="evidence"><b>{esc(rel_cn)}</b>({esc(relt)}) | {esc(subject)}{obj_part} | '
               f'{pmid_html} | <span class="note">{direct_note}</span>']
        if sent:
            out.append(f'<br><span class="note">Evidence sentence: {esc(sent[:300])}</span>')
        out.append("</div>")
        return "".join(out)

    for idx, sr in enumerate(ordered, 1):
        sp = sr.get("species", "")
        rels = rels_by_sp.get(sp, [])
        tier = sr.get("tier", "Low")
        suggest = {
            "High": "Recommend adding more sequences of the same species, sampling metadata, and re-checking the literature evidence.",
            "Medium": "Recommend continuous tracking of new literature and host-association evidence.",
            "Low": "Routine monitoring is recommended, with attention to host distribution.",
        }.get(tier, "Recommend continuous tracking of new literature and host-association evidence.")
        L.append(f"<h3>6.{idx} {esc(sp)}</h3>")
        L.append("<ul>")
        L.append(
            f"<li>Score components: Host type {sr.get('host_type_score', 0)} + Host breadth {sr.get('host_breadth_score', 0)} + "
            f"Phylogeny {sr.get('phylogeny_score', 0)} + Literature {sr.get('literature_score', 0)} + "
            f"Human infection {sr.get('human_infection_score', 0)} = Species risk score {sr.get('species_risk_score', '')} ({tag(tier)})</li>"
        )
        # KG (contig-level) relations: state honestly when none were produced, so "0 relations" is not misread as "no literature"
        if rels:
            L.append(f"<li>Knowledge-graph evidence: {len(rels)} relations</li>")
        else:
            L.append('<li>Knowledge-graph evidence: no contig-level relations were produced in this run (that step was not executed); species-level literature results follow.</li>')
        # Species-level literature evidence (PMIDs of human-infection reports from PubTator3; contig-level KG is a separate step)
        lit = (species_lit or {}).get(sp) or {}
        lit_pmids = [str(p) for p in (lit.get("pmids") or []) if str(p).strip()]
        if lit_pmids:
            links = ", ".join(
                f'<a href="https://pubmed.ncbi.nlm.nih.gov/{esc(p)}/" target="_blank">PMID:{esc(p)}</a>'
                for p in lit_pmids)
            L.append(f"<li>Literature evidence (species-level): {links} ({len(lit_pmids)} reports of human infection)</li>")
        elif lit:
            L.append("<li>Literature evidence (species-level): no human-infection reports found</li>")
        L.append(f"<li>Recommendation: {suggest}</li>")
        L.append("</ul>")

        # Split relations by contig (anchors only for "contigs with relations + L-segment contigs")
        rels_by_contig = defaultdict(list)
        for rel in rels:
            rels_by_contig[str(rel.get("record_id") or "")].append(rel)
        contig_keys = []
        for o in rows_out:
            if o.get("species", "") == sp:
                k = f"{o.get('sample', '')}_{o.get('id', '')}"
                if k in kg_keys and k not in contig_keys:
                    contig_keys.append(k)
        for k in rels_by_contig:
            if k not in contig_keys:
                contig_keys.append(k)
        for ck in contig_keys:
            crels = rels_by_contig.get(ck, [])
            back = f' <a class="note" href="#chart-{esc(ck)}">[Back to chart]</a>' if ck in l_keys else ""
            L.append(f'<h4 id="kg-{esc(ck)}">{esc(ck)}({len(crels)} relations){back}</h4>')
            if crels:
                L.append("<div>")
                for rel in crels[:20]:
                    L.append(render_rel(rel))
                if len(crels) > 20:
                    L.append(f'<p class="note">... showing the first 20 of {len(crels)}.</p>')
                L.append("</div>")
            else:
                L.append('<p class="note">No literature relation evidence.</p>')

    L.append("</div>")
    L.append("</main>")
    L.append('<button class="totop" id="totop" aria-label="Back to top">↑</button>')
    L.append("<script>")
    L.append("(function(){var b=document.getElementById('totop');window.addEventListener('scroll',function(){b.classList.toggle('show',window.scrollY>600);});b.addEventListener('click',function(){window.scrollTo({top:0,behavior:'smooth'});});})();")
    L.append("(function(){var reveal=function(){document.querySelectorAll('.fade').forEach(function(el){el.classList.add('fade-in');});};if('IntersectionObserver' in window){var io=new IntersectionObserver(function(es){es.forEach(function(e){if(e.isIntersecting){e.target.classList.add('fade-in');io.unobserve(e.target);}});},{threshold:0.05});document.querySelectorAll('.fade').forEach(function(el){io.observe(el);});}else{reveal();}window.addEventListener('hashchange',reveal);setTimeout(reveal,3000);})();")
    L.append("</script>")
    L.append("</body>")
    L.append("</html>")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="5.4 Report generation: read report_inputs.json from step 5.3 → report.html")
    p.add_argument("--config", type=str, default="",
                   help="path to config.yaml; reads results_dir/00_results/report_inputs.json")
    p.add_argument("--bundle", type=Path, default=None,
                   help="path to report_inputs.json (used without --config)")
    p.add_argument("--output", type=Path, default=None,
                   help="output HTML path; defaults to report.html in the bundle directory")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if args.config:
        import yaml
        cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
        results_dir = Path(cfg["paths"]["results_dir"])
        bundle_path = results_dir / "00_results" / "report_inputs.json"
        outdir = results_dir / "00_results"
    elif args.bundle:
        bundle_path = args.bundle
        outdir = bundle_path.parent
    else:
        print("Requires --bundle (or --config)", file=sys.stderr)
        return 2

    if not bundle_path.is_file():
        print(f"Step-5.3 intermediate file not found: {bundle_path}", file=sys.stderr)
        return 2

    data = json.loads(bundle_path.read_text(encoding="utf-8"))
    rows = data.get("rows") or []
    rows_full = data.get("rows_full") or []
    species = data.get("species") or {}
    relations = data.get("relations") or None
    qc_dir = Path(data["qc_dir"]) if data.get("qc_dir") else None

    # Species-level literature evidence (pre-fetched PubTator3 results, optional)
    species_lit: dict = {}
    sp_lit_path = bundle_path.parent / "species_pmids.json"
    if sp_lit_path.is_file():
        try:
            species_lit = (json.loads(sp_lit_path.read_text(encoding="utf-8")) or {}).get("species", {})
            print(f"[info] Loaded species-level literature evidence for {len(species_lit)} species", file=sys.stderr)
        except Exception as e:
            print(f"[warning] Failed to read species_pmids.json: {e}", file=sys.stderr)

    out_html = args.output or (outdir / "report.html")
    write_report_html(out_html, rows, species, relations, qc_dir, rows_full, species_lit)
    print(json.dumps({"html": str(out_html),
                      "row_count": len(rows),
                      "species_count": len(species)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
