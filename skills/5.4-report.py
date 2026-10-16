#!/usr/bin/env python3
"""5.4 报告生成：读取 5.3 输出的中间文件 report_inputs.json，生成 report.html。

流程位置：5.3（风险打分，产出 results.tsv / species_risk.tsv / report_inputs.json）
          → 本脚本（把打分结果渲染成自包含 HTML 报告）。

输入：
  - report_inputs.json（5.3 输出；含 rows / species / relations / qc_dir）
输出：
  - report.html（自包含 HTML：质控（fastp 超链接）→ 打分 → 判断依据（知识图谱 + 建议））

用法：
  python 5.4-report.py --bundle report_inputs.json
  python 5.4-report.py --config config.yaml          # 从 results_dir/00_results/report_inputs.json 读
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
    """HTML 转义，防止数据中的特殊字符破坏页面。"""
    return (str(s)
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
            .replace('"', "&quot;"))


REL_CN = {
    "virus_detected_in_host": "病毒在宿主中检出",
    "virus_has_gene_or_segment": "病毒含基因/片段",
    "virus_coinfects_with": "共感染",
    "virus_related_to": "近缘病毒",
    "virus_reassortment": "重配",
    "virus_host_switching": "宿主转换",
    "virus_associated_with_host": "病毒与宿主关联",
    "variant_reported": "变异报道",
}


# ---------- 圆形柱状图（每条 L 段序列一张：内环=检出可信度3项，外环=物种风险5项） ----------
CBAR_INNER = [
    ("DIAMOND", "#2166ac"),
    ("Tree", "#67a9cf"),
    ("LucaProt", "#d1e5f0"),
]
# 外环 8 项：宿主(橙) + 系统发育(紫) + 文献4项(绿) + 人类感染(红)
CBAR_OUTER = [
    ("宿主类型", "host_type_score", 10.0, "#e08214"),
    ("宿主广度", "host_breadth_score", 5.0, "#f4a582"),
    ("系统发育", "phylogeny_score", 15.0, "#7570b3"),
    ("重配", "reassortment", 10.0, "#1b7837"),
    ("宿主转换", "host_switching", 10.0, "#5aae61"),
    ("共感染", "coinfection", 5.0, "#a6dba0"),
    ("近缘", "related_viruses", 5.0, "#d9f0d3"),
    ("人类感染", "human_infection_score", 40.0, "#d73027"),
]

# 文献 4 性状 -> (关系类型, 是否要求 direct_target_match)，与 5.3 literature_traits 保持一致
LIT_TRAIT = {
    "reassortment": ("virus_reassortment", False),
    "host_switching": ("virus_host_switching", False),
    "coinfection": ("virus_coinfects_with", True),
    "related_viruses": ("virus_related_to", True),
}
LIT_TRAIT_CN = {
    "reassortment": "重配",
    "host_switching": "宿主转换",
    "coinfection": "共感染",
    "related_viruses": "近缘",
}
LIT_PLACEHOLDER = {"", "related virus", "co-infecting virus", "unspecified variant"}
MISS_COLOR = "#d9d9d9"

_cbar_uid = 0


def _num(v, default=0.0):
    try:
        if v is None or v == "":
            return default
        return float(v)
    except (TypeError, ValueError):
        return default


def _lit_evidence(sp: str, rels_by_sp: dict) -> dict:
    """某物种文献 4 性状各自的支撑关系（PMID + 证据句）。"""
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
    """文献性状得分：命中=满分，未命中=0。"""
    detail = sr.get("literature_detail") or {}
    return maxv if detail.get(trait) is True else 0.0


def _hover_builder(o: dict, sr: dict, lit_ev: dict, full: dict) -> dict:
    """为每根柱生成悬停 HTML 字符串（含原始比对/树/宿主/文献依据）。"""
    def kv(label, v, suffix=""):
        s = str(v or "").strip()
        return f"{label}：{_esc(s)}{suffix}" if s else None

    def pct(score, maxv):
        return round(min(1.0, max(0.0, score / maxv)) * 100, 1) if maxv else 0.0

    h = {}

    # ---- 内环（检出可信度）----
    d = _num(o.get("detection_diamond"))
    dlines = [f"得分：{d} / 1（占比 {pct(d, 1)}%）"]
    for ln in [kv("参考序列", full.get("sseqid")),
               kv("一致率 pident", full.get("pident"), "%"),
               kv("比对长度", full.get("length"), " bp"),
               kv("evalue", full.get("evalue")),
               kv("bitscore", full.get("bitscore"))]:
        if ln:
            dlines.append(ln)
    if full.get("qstart") and full.get("qend"):
        dlines.append(f"query 位置：{full.get('qstart')}-{full.get('qend')}")
    if full.get("sstart") and full.get("send"):
        dlines.append(f"参考位置：{full.get('sstart')}-{full.get('send')}")
    h["DIAMOND"] = "<br/>".join(dlines)

    t = _num(o.get("detection_tree"))
    tlines = [f"得分：{t} / 1（占比 {pct(t, 1)}%）"]
    for ln in [kv("树位置", full.get("tree_position")),
               kv("bootstrap", full.get("tree_bootstrap")),
               kv("最近邻", full.get("tree_nn_sseqid")),
               kv("距离", full.get("tree_nn_distance"))]:
        if ln:
            tlines.append(ln)
    h["Tree"] = "<br/>".join(tlines)

    l = _num(o.get("detection_lucaprot"))
    h["LucaProt"] = f"得分：{l} / 1（RdRp 概率，占比 {pct(l, 1)}%）"

    # ---- 外环（物种风险分）----
    ht = _num(sr.get("host_type_score"))
    hlines = [f"得分：{ht} / 10（占比 {pct(ht, 10)}%）"]
    for ln in [kv("宿主", full.get("host")), kv("分类", sr.get("host_class"))]:
        if ln:
            hlines.append(ln)
    h["宿主类型"] = "<br/>".join(hlines)

    br = sr.get("host_breadth") or 0
    h["宿主广度"] = f"得分：{br} / 5（占比 {pct(br, 5)}%）<br/>宿主种数：{br}"

    phy = _num(sr.get("phylogeny_score"))
    plines = [f"得分：{phy} / 15（占比 {pct(phy, 15)}%）"]
    for ln in [kv("属", full.get("genus")), kv("最近邻距离", full.get("tree_nn_distance"))]:
        if ln:
            plines.append(ln)
    h["系统发育"] = "<br/>".join(plines)

    for trait, cn, maxv in [("reassortment", "重配", 10), ("host_switching", "宿主转换", 10),
                            ("coinfection", "共感染", 5), ("related_viruses", "近缘", 5)]:
        hit = bool((sr.get("literature_detail") or {}).get(trait))
        if not hit:
            h[cn] = f"未命中：0 / {maxv}"
            continue
        rels = lit_ev.get(trait) or []
        pmids = list(dict.fromkeys(str(r.get("pmid")) for r in rels if r.get("pmid")))[:3]
        sents = [str(r.get("evidence_sentence") or "").strip()
                 for r in rels if r.get("evidence_sentence")]
        lines = [f"得分：{maxv} / {maxv}（占比 100%）"]
        if pmids:
            lines.append("来源 PMID：" + "、".join(_esc(p) for p in pmids))
        if sents:
            lines.append("证据句：" + _esc(sents[0][:120]))
        h[cn] = "<br/>".join(lines)

    if sr.get("human_infection_vmr"):
        h["人类感染"] = "得分：40 / 40（占比 100%）<br/>来源：VMR 官方人类分离株"
    elif sr.get("human_infection_literature"):
        hpm = [str(p) for p in (sr.get("human_pmids") or [])]
        lines = ["得分：40 / 40（占比 100%）"]
        if hpm:
            lines.append("来源 PMID：" + "、".join(_esc(p) for p in hpm))
        h["人类感染"] = "<br/>".join(lines)
    else:
        h["人类感染"] = "未命中：0 / 40"

    return h


def _cbar_option(o: dict, sr: dict, lit_ev: dict, full: dict) -> dict:
    """单条序列的双层圆形柱状图 ECharts option（内环=检出3项，外环=物种风险8项）。"""
    hover = _hover_builder(o, sr, lit_ev, full)
    inner = [
        ("DIAMOND", _num(o.get("detection_diamond")), 1.0, "#2166ac", None),
        ("Tree", _num(o.get("detection_tree")), 1.0, "#67a9cf", None),
        ("LucaProt", _num(o.get("detection_lucaprot")), 1.0, "#d1e5f0", None),
    ]
    outer_raw = [
        ("宿主类型", _num(sr.get("host_type_score")), 10.0, "#e08214", None),
        ("宿主广度", _num(sr.get("host_breadth_score")), 5.0, "#f4a582", None),
        ("系统发育", _num(sr.get("phylogeny_score")), 15.0, "#7570b3", None),
        ("重配", _lit_subscore(sr, "reassortment", 10.0), 10.0, "#1b7837", "reassortment"),
        ("宿主转换", _lit_subscore(sr, "host_switching", 10.0), 10.0, "#5aae61", "host_switching"),
        ("共感染", _lit_subscore(sr, "coinfection", 5.0), 5.0, "#a6dba0", "coinfection"),
        ("近缘", _lit_subscore(sr, "related_viruses", 5.0), 5.0, "#d9f0d3", "related_viruses"),
        ("人类感染", _num(sr.get("human_infection_score")), 40.0, "#d73027", "human"),
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
    """生成一条序列的圆形柱状图卡片（div + 中心最终分 + 点击跳转 + ECharts 脚本）。"""
    global _cbar_uid
    _cbar_uid += 1
    uid = _cbar_uid
    tier = str(o.get("tier") or "Low")
    key = f"{o.get('sample', '')}_{o.get('id', '')}"
    opt = json.dumps(_cbar_option(o, sr, lit_ev, full), ensure_ascii=False).replace('"@@CBARTIP@@"', "cbarTip")
    head = (
        f'<div class="cbar-head"><b>{esc(o.get("id", ""))}</b> '
        f'<span class="tag {esc(tier)}">{esc(tier)}</span>'
        f'<div class="cbar-sub">{esc(o.get("sample", ""))} · L段</div></div>'
    )
    center = (
        f'<div class="cbar-center"><div class="cbar-score">{esc(o.get("final_score", ""))}</div>'
        f'<div class="cbar-score-label">最终分</div></div>'
    )
    links = (
        f'<div class="cbar-links">'
        f'<a href="#seq-{esc(key)}">第5节明细→</a> '
        f'<a href="#kg-{esc(key)}">第6节知识图谱→</a></div>'
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
                      rows_full: list[dict] | None = None) -> None:
    """生成自包含 HTML 报告：质控（fastp 超链接）→ 打分 → 判断依据（知识图谱 + 建议）。"""
    esc = _esc

    # 完整行（含原始比对/树列）按 sample_id 建索引，供悬停细化
    full_by_key: dict[str, dict] = {}
    for f in (rows_full or []):
        key = f"{f.get('sample', '')}_{f.get('id', '')}"
        if key:
            full_by_key[key] = f

    # ---- 知识图谱关系按物种分组（record_id 为 sample_contig，与行级 key 对应）----
    key2sp: dict[str, str] = {}
    for o in rows_out:
        key = f"{o.get('sample', '')}_{o.get('id', '')}"
        if key and key not in key2sp:
            key2sp[key] = o.get("species", "未分类")
    rels_by_sp: dict[str, list[dict]] = defaultdict(list)
    for rel in (relations or []):
        sp = key2sp.get(str(rel.get("record_id") or ""))
        if sp:
            rels_by_sp[sp].append(rel)
    # 文献证据（按物种）：4 性状各自支撑关系（PMID + 证据句）
    lit_ev_by_sp = {sp: _lit_evidence(sp, rels_by_sp) for sp in species_risk}
    # 需要建第6节锚点的 contig：有关系的 contig + 全部 L 段 contig
    rel_keys = {str(rel.get("record_id")) for rel in (relations or [])}
    l_keys = {f"{o.get('sample', '')}_{o.get('id', '')}" for o in rows_out
              if str(o.get("segment", "")).strip().lower() == "l"}
    kg_keys = rel_keys | l_keys

    # ---- 质控样本（按结果表出现顺序去重；只列有检出的样本）----
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

    # ---- 基本统计 ----
    n_rows = len(rows_out)
    n_species = len(species_risk)
    tier_counts: dict[str, int] = defaultdict(int)
    for sr in species_risk.values():
        tier_counts[sr.get("tier", "Low")] += 1

    by_sp: dict[str, list[dict]] = defaultdict(list)
    for o in rows_out:
        by_sp[o.get("species", "未分类")].append(o)
    ordered = sorted(species_risk.values(), key=lambda x: x.get("species_risk_score", 0.0), reverse=True)

    def tag(tier: str) -> str:
        return f'<span class="tag {esc(tier)}">{esc(tier)}</span>'

    L: list[str] = []
    L.append("<!DOCTYPE html>")
    L.append('<html lang="zh-CN">')
    L.append("<head>")
    L.append('<meta charset="utf-8">')
    L.append('<meta name="viewport" content="width=device-width, initial-scale=1">')
    L.append("<title>汉坦病毒检出与风险评分评估报告</title>")
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

    # 顶部导航
    L.append('<nav class="topnav"><div class="topnav-inner">'
             '<div class="topnav-title">汉坦病毒风险评分</div>'
             '<div class="topnav-links">'
             '<a href="#s1">质控</a><a href="#s2">摘要</a><a href="#s3">总览</a>'
             '<a href="#s4">详情</a><a href="#s5">明细</a><a href="#s6">依据</a>'
             '</div></div></nav>')

    # Hero
    L.append('<header class="hero">')
    L.append('<div class="hero-badge">Hantavirus Risk Assessment</div>')
    L.append('<h1>汉坦病毒检出与风险评分<br>评估报告</h1>')
    L.append(f'<p class="hero-sub">{n_rows} 条检出序列 · {n_species} 个物种 · 检出可信度 × 物种风险分</p>')
    L.append(f'<p class="hero-meta">高 {tier_counts.get("High", 0)} ｜ 中 {tier_counts.get("Medium", 0)} ｜ 低 {tier_counts.get("Low", 0)} ｜ 生成日期 {date.today()}</p>')
    L.append('</header>')
    L.append('<main class="container">')

    # 1 质控情况
    L.append('<div class="section fade" id="s1">')
    L.append('<h2><span class="sec-no">01</span>质控情况</h2>')
    L.append('<p class="note">仅列出有汉坦病毒检出的样本；点击样本名打开 fastp 质控报告。</p>')
    L.append("<table><tr><th>样本</th><th>质控报告</th><th>检出序列数</th></tr>")
    for q in qc_rows:
        name = esc(q["sample"])
        if q["exists"]:
            cell = f'<a href="{esc(q["href"])}" target="_blank">{name}</a>'
        else:
            cell = f'{name}<span class="note">（无质控报告）</span>'
        L.append(f"<tr><td>{name}</td><td>{cell}</td><td>{q['count']}</td></tr>")
    L.append("</table>")
    L.append("</div>")

    # 2 摘要
    L.append('<div class="section fade" id="s2">')
    L.append('<h2><span class="sec-no">02</span>摘要</h2>')
    L.append(f"<p>本报告对 {n_rows} 条检出序列（涉及 {n_species} 个物种/病毒）计算「检出可信度 × 物种风险分」，并划分高/中/低三档。</p>")
    L.append(f"<p>物种三档分布：高 {tier_counts.get('High', 0)} 个、中 {tier_counts.get('Medium', 0)} 个、低 {tier_counts.get('Low', 0)} 个。</p>")
    L.append("</div>")

    # 3 判定结果总览
    L.append('<div class="section fade" id="s3">')
    L.append('<h2><span class="sec-no">03</span>判定结果总览</h2>')
    L.append("<table><tr><th>物种</th><th>病毒名</th><th>序列数</th><th>宿主类型</th><th>物种风险分</th><th>判定</th></tr>")
    for sr in ordered:
        sp = sr.get("species", "")
        rows_sp = by_sp.get(sp, [])
        virus_name = next((o.get("virus_name", "") for o in rows_sp if o.get("virus_name")), "")
        L.append(
            f"<tr><td>{esc(sp)}</td><td>{esc(virus_name)}</td><td>{len(rows_sp)}</td>"
            f"<td>{esc(sr.get('host_desc', ''))}</td><td>{sr.get('species_risk_score', '')}</td>"
            f"<td>{tag(sr.get('tier', 'Low'))}</td></tr>"
        )
    L.append("</table>")
    L.append("</div>")

    # 4 各物种评估详情
    L.append('<div class="section fade" id="s4">')
    L.append('<h2><span class="sec-no">04</span>各物种评估详情</h2>')
    L.append('<div class="note" style="margin:8px 0;">圆形柱状图：内环=检出可信度（'
             '<span style="color:#2166ac;">●DIAMOND</span> '
             '<span style="color:#67a9cf;">●树</span> '
             '<span style="color:#d1e5f0;">●LucaProt</span>）；外环=物种风险分：宿主（'
             '<span style="color:#e08214;">●宿主类型</span> '
             '<span style="color:#f4a582;">●宿主广度</span>）、系统发育（'
             '<span style="color:#7570b3;">●系统发育</span>）、文献（'
             '<span style="color:#1b7837;">●重配</span> '
             '<span style="color:#5aae61;">●宿主转换</span> '
             '<span style="color:#a6dba0;">●共感染</span> '
             '<span style="color:#d9f0d3;">●近缘</span>）、人类感染（'
             '<span style="color:#d73027;">●人类感染</span>）。柱长=该项÷满分，未命中画灰柱；悬停看数值与文献 PMID。</div>')
    for idx, sr in enumerate(ordered, 1):
        sp = sr.get("species", "")
        rows_sp = by_sp.get(sp, [])
        virus_name = next((o.get("virus_name", "") for o in rows_sp if o.get("virus_name")), "")
        pmids = "; ".join(sr.get("human_pmids", [])) or "无"
        human_src = "有" if (sr.get("human_infection_vmr") or sr.get("human_infection_literature")) else "无"
        L.append(f"<h3>4.{idx} {esc(sp)}</h3>")
        L.append("<ul>")
        if virus_name:
            L.append(f"<li>病毒名：{esc(virus_name)}</li>")
        L.append(f"<li>宿主：{esc(sr.get('host_desc', ''))}（类型 {esc(sr.get('host_class', ''))}）</li>")
        L.append(f"<li>检出序列数：{len(rows_sp)}</li>")
        L.append(
            f"<li>分项：宿主类型 {sr.get('host_type_score', 0)} + 宿主广度 {sr.get('host_breadth_score', 0)} + "
            f"系统发育 {sr.get('phylogeny_score', 0)} + 文献 {sr.get('literature_score', 0)} + "
            f"人类感染 {sr.get('human_infection_score', 0)}</li>"
        )
        L.append(f"<li>物种风险分：<b>{sr.get('species_risk_score', '')}</b> → 判定：{tag(sr.get('tier', 'Low'))}</li>")
        L.append(f"<li>人类感染证据：{human_src}（PMID：{esc(pmids)}）</li>")
        L.append("</ul>")

        # 文献/人类感染证据来源（物种级，PMID 可点击）
        lit_ev = lit_ev_by_sp.get(sp, {})
        src_parts = []
        for trait, cn in LIT_TRAIT_CN.items():
            if (sr.get("literature_detail") or {}).get(trait) is True:
                pms = sorted({str(r.get("pmid")) for r in lit_ev.get(trait, []) if r.get("pmid")})
                if pms:
                    links = "、".join(
                        f'<a href="https://pubmed.ncbi.nlm.nih.gov/{esc(p)}/" target="_blank">PMID:{esc(p)}</a>'
                        for p in pms)
                    src_parts.append(f"{cn}（{links}）")
                else:
                    src_parts.append(cn)
        if sr.get("human_infection_vmr"):
            src_parts.append("人类感染（VMR 官方人类分离株）")
        elif sr.get("human_infection_literature") and sr.get("human_pmids"):
            links = "、".join(
                f'<a href="https://pubmed.ncbi.nlm.nih.gov/{esc(p)}/" target="_blank">PMID:{esc(p)}</a>'
                for p in sr.get("human_pmids", []))
            src_parts.append(f"人类感染（{links}）")
        if src_parts:
            L.append(f'<p class="note">文献/人类感染证据来源：{"；".join(src_parts)}</p>')
        else:
            L.append('<p class="note">文献/人类感染证据来源：无。</p>')

        l_rows = [o for o in rows_sp if str(o.get("segment", "")).strip().lower() == "l"]
        if l_rows:
            L.append('<p class="note">该物种 L 段序列打分来源（悬停看原始依据与文献 PMID，点图跳转第5/6节）：</p>')
            L.append('<div class="cbar-grid">')
            for o in l_rows:
                key = f"{o.get('sample', '')}_{o.get('id', '')}"
                L.append(_cbar_html(o, sr, esc, lit_ev, full_by_key.get(key, {})))
            L.append('</div>')
    L.append("</div>")

    # 5 检出可信度明细（默认展开，可折叠；点圆图跳到这里并高亮该行）
    L.append('<div class="section fade" id="s5">')
    L.append('<h2><span class="sec-no">05</span>检出可信度明细</h2>')
    L.append("<details open>")
    L.append("<summary>点击折叠 / 展开明细表</summary>")
    L.append("<table><tr><th>样本</th><th>序列ID</th><th>片段</th><th>DIAMOND</th><th>树</th><th>LucaProt</th><th>检出可信度</th><th>物种风险分</th><th>最终分</th><th>判定</th><th>跳转</th></tr>")
    for o in rows_out:
        seg = (o.get("segment") or "").strip().upper() or "-"
        d_tree = o.get("detection_tree") or "-"
        d_luca = o.get("detection_lucaprot") or "-"
        key = f"{o.get('sample', '')}_{o.get('id', '')}"
        kg_link = f'<a href="#kg-{esc(key)}">知识图谱→</a>' if key in kg_keys else "—"
        chart_link = f'<a href="#chart-{esc(key)}">返回图→</a>' if key in l_keys else "—"
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

    # 6 判断依据与建议（知识图谱 + 打分分项 + 建议；关系按 contig 拆分，可点图/链接跳转定位）
    L.append('<div class="section fade" id="s6">')
    L.append('<h2><span class="sec-no">06</span>判断依据与建议</h2>')
    L.append('<p class="note">以下列出每个物种的打分分项、建议，以及按 contig 拆分的知识图谱文献证据（附 PMID 与证据句）。</p>')

    def render_rel(rel):
        relt = rel.get("relation") or "relation"
        rel_cn = REL_CN.get(relt, relt)
        subject = rel.get("subject") or ""
        obj = rel.get("object") or ""
        pmid = rel.get("pmid") or ""
        sent = (rel.get("evidence_sentence") or "").strip()
        direct = rel.get("direct_target_match")
        direct_note = "直接证据" if direct is True else "背景/相关"
        pmid_html = (f'<a href="https://pubmed.ncbi.nlm.nih.gov/{esc(pmid)}/" target="_blank">PMID:{esc(pmid)}</a>'
                     if pmid else "—")
        obj_part = f" → {esc(obj)}" if obj else ""
        out = [f'<div class="evidence"><b>{esc(rel_cn)}</b>（{esc(relt)}）｜{esc(subject)}{obj_part}｜'
               f'{pmid_html}｜<span class="note">{direct_note}</span>']
        if sent:
            out.append(f'<br><span class="note">证据句：{esc(sent[:300])}</span>')
        out.append("</div>")
        return "".join(out)

    for idx, sr in enumerate(ordered, 1):
        sp = sr.get("species", "")
        rels = rels_by_sp.get(sp, [])
        tier = sr.get("tier", "Low")
        suggest = {
            "High": "建议补充同种病毒多样本序列、采样元数据并复核文献证据。",
            "Medium": "建议持续跟踪新增文献与宿主关联证据。",
            "Low": "建议常规监测，结合宿主分布持续关注。",
        }.get(tier, "建议持续跟踪新增文献与宿主关联证据。")
        L.append(f"<h3>6.{idx} {esc(sp)}</h3>")
        L.append("<ul>")
        L.append(
            f"<li>打分分项：宿主类型 {sr.get('host_type_score', 0)} + 宿主广度 {sr.get('host_breadth_score', 0)} + "
            f"系统发育 {sr.get('phylogeny_score', 0)} + 文献 {sr.get('literature_score', 0)} + "
            f"人类感染 {sr.get('human_infection_score', 0)} = 物种风险分 {sr.get('species_risk_score', '')}（{tag(tier)}）</li>"
        )
        L.append(f"<li>知识图谱证据：共 {len(rels)} 条关系</li>")
        L.append(f"<li>建议：{suggest}</li>")
        L.append("</ul>")

        # 按 contig 拆关系（只给「有关系的 contig + L 段 contig」建锚点）
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
            back = f' <a class="note" href="#chart-{esc(ck)}">[返回图]</a>' if ck in l_keys else ""
            L.append(f'<h4 id="kg-{esc(ck)}">{esc(ck)}（{len(crels)} 条关系）{back}</h4>')
            if crels:
                L.append("<div>")
                for rel in crels[:20]:
                    L.append(render_rel(rel))
                if len(crels) > 20:
                    L.append(f'<p class="note">……仅显示前 20 条，共 {len(crels)} 条。</p>')
                L.append("</div>")
            else:
                L.append('<p class="note">无文献关系证据。</p>')

    L.append("</div>")
    L.append("</main>")
    L.append('<button class="totop" id="totop" aria-label="回到顶部">↑</button>')
    L.append("<script>")
    L.append("(function(){var b=document.getElementById('totop');window.addEventListener('scroll',function(){b.classList.toggle('show',window.scrollY>600);});b.addEventListener('click',function(){window.scrollTo({top:0,behavior:'smooth'});});})();")
    L.append("(function(){var reveal=function(){document.querySelectorAll('.fade').forEach(function(el){el.classList.add('fade-in');});};if('IntersectionObserver' in window){var io=new IntersectionObserver(function(es){es.forEach(function(e){if(e.isIntersecting){e.target.classList.add('fade-in');io.unobserve(e.target);}});},{threshold:0.05});document.querySelectorAll('.fade').forEach(function(el){io.observe(el);});}else{reveal();}window.addEventListener('hashchange',reveal);setTimeout(reveal,3000);})();")
    L.append("</script>")
    L.append("</body>")
    L.append("</html>")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="5.4 报告生成：读取 5.3 的 report_inputs.json → report.html")
    p.add_argument("--config", type=str, default="",
                   help="config.yaml 路径；从 results_dir/00_results/report_inputs.json 读取")
    p.add_argument("--bundle", type=Path, default=None,
                   help="report_inputs.json 路径（不配合 --config 时使用）")
    p.add_argument("--output", type=Path, default=None,
                   help="输出 HTML 路径；默认 bundle 同目录 report.html")
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
        print("需要 --bundle（或 --config）", file=sys.stderr)
        return 2

    if not bundle_path.is_file():
        print(f"未找到 5.3 中间文件：{bundle_path}", file=sys.stderr)
        return 2

    data = json.loads(bundle_path.read_text(encoding="utf-8"))
    rows = data.get("rows") or []
    rows_full = data.get("rows_full") or []
    species = data.get("species") or {}
    relations = data.get("relations") or None
    qc_dir = Path(data["qc_dir"]) if data.get("qc_dir") else None

    out_html = args.output or (outdir / "report.html")
    write_report_html(out_html, rows, species, relations, qc_dir, rows_full)
    print(json.dumps({"html": str(out_html),
                      "row_count": len(rows),
                      "species_count": len(species)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
