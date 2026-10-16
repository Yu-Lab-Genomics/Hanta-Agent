#!/usr/bin/env python3
"""5.3 风险打分：检出可信度 × 物种风险分，高中低三档，并用金标准标签验证。

流程位置：4 号（建树）→ 5 号（PubTator3 文献证据）→ 本脚本（打分 + 验证）。
本脚本取代 5.2-early_warning_score.py（旧规则分只作为敏感性对照，未实现时跳过）。

输入：
  - 带树列的 *.vmr_final_classified.tsv（4 号脚本输出；含 final_classification /
    classification_source / tree_position / tree_nn_sseqid / tree_nn_distance）
  - extracted_relations.json（5 号输出；缺失时文献性状记 0 并警告）
  - retrieval_queries.tsv（5 号输出；用于检出偏倚协变量）
  - risk_labels.tsv（金标准标签；用于验证）
  - LucaProt 输出（可选 TSV：序列id + RdRp 概率；缺失时按 DIAMOND+树归一）

打分规则（与 results/figure3/scoring_rules.png 一致）：
  检出可信度（0–1）= 0.35×DIAMOND + 0.30×树 + 0.35×LucaProt（缺失项按剩余权重归一）
  物种风险分（0–100）= 宿主 15 + 系统发育 15 + 文献性状 30 + 传染人 40
    宿主：类型（啮齿10/蝙蝠5/鼩鼱·鼹鼠2/鱼·爬行0）+ 广度（1种0/2–3种2/≥4种5）
    系统发育：已知致病种 15；啮齿 Orthohantavirus 10；其余 0；按树距离微调
    文献性状：重配+10、宿主转换+10、共感染+5、近缘病毒+5（不含人类感染证据）
    传染人：文献明确记载人类感染 +40（PubTator3 实时检索；离线时跳过）
  最终分 = 检出可信度 × 物种风险分；三档：高 ≥60 / 中 30–59 / 低 <30
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

PUBTATOR_SEARCH = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/?text={text}&page=1"
HDR = {"User-Agent": "Yu-hanta-risk-score/0.1 (contact: local)"}

# ---------- 宿主清洗与分类 ----------
RODENT_KW = (
    "rattus", "ruttus", "niviventer", "peromyscus", "oligoryzomys", "microtus",
    "myodes", "clethrionomys", "apodemus", "mus ", "bandicota", "hylomyscus",
    "reithrodontomys", "oryzomys", "sigmodon", "neotom", "calomys", "eothenomys",
    "stenocephalemys", "eliurus", "akodon", "lemmus", "ochrotomys", "cricet",
    "murid", "arvicol", "rodent", "vole", "gerbil", "hamster", "rat", "mice",
)
BAT_KW = (
    "chiroptera", "rhinolophus", "pteropus", "nyctalus", "aselliscus", "taphozous",
    "rousettus", "hipposideros", "pipistrellus", "eptesicus", "miniopterus",
    "molossid", "vespertilionid", "pteropodid", "emballonurid", "bat",
)
SHREW_MOLE_KW = (
    "sorex", "crocidura", "chodsigoa", "anourosorex", "suncus", "talpa",
    "urotrichus", "scalopus", "shrew", "mole", "eulipotyphla", "insectivor",
)
FISH_REPTILE_KW = (
    "perca", "halieutaea", "triacanthodes", "eptatretus", "luposicya", "lophius",
    "hemidactylus", "gecko", "reptil", "actinopteryg", "agnathan", "perciform",
    "gobiid", "percid", "fish", "hagfish", "goosefish", "batfish", "spikefish",
)
HUMAN_KW = ("homo sapiens", "human", "patient")

HOST_TYPE_POINTS = {"rodent": 10, "bat": 5, "shrew_mole": 2, "fish_reptile": 0, "other": 0}


def clean_host(raw: str) -> str:
    h = (raw or "").casefold().strip()
    h = h.replace("ruttus", "rattus")
    h = re.sub(r"\(.*?\)", "", h)                 # 去括号俗名
    h = re.sub(r"strain.*$", "", h).strip()        # 去 strain 后缀
    h = re.sub(r"msb:mamm.*$", "", h).strip()      # 去标本馆编号
    h = re.sub(r"\s+", " ", h)
    return h


def host_class(hosts: list[str]) -> tuple[str, str]:
    cleaned = sorted({clean_host(h) for h in hosts if clean_host(h)})
    if not cleaned:
        return "other", "无宿主注释"
    if any(any(k in h for k in HUMAN_KW) for h in cleaned):
        return "human", "人类分离株"
    if any(any(k in h for k in RODENT_KW) for h in cleaned):
        return "rodent", "啮齿"
    if any(any(k in h for k in BAT_KW) for h in cleaned):
        return "bat", "蝙蝠"
    if any(any(k in h for k in SHREW_MOLE_KW) for h in cleaned):
        return "shrew_mole", "鼩鼱/鼹鼠"
    if any(any(k in h for k in FISH_REPTILE_KW) for h in cleaned):
        return "fish_reptile", "鱼/爬行"
    return "other", "其他"


# ---------- 检出可信度 ----------
def diammond_score(row: dict[str, str]) -> float:
    try:
        ev = float(row.get("evalue") or "1")
    except ValueError:
        ev = 1.0
    ev_s = min(1.0, max(0.0, -math.log10(max(ev, 1e-300)) / 100.0))
    try:
        length = float(row.get("length") or 0)
        sstart, send = float(row.get("sstart") or 0), float(row.get("send") or 0)
        span = abs(send - sstart) + 1
        cov = min(1.0, length / span) if span > 1 else 0.5
    except ValueError:
        cov = 0.5
    return 0.6 * ev_s + 0.4 * cov


def tree_score(row: dict[str, str]) -> float | None:
    pos = row.get("tree_position") or row.get("tree_classification") or ""
    if not pos:
        return None
    pos_s = {"nested": 1.0, "sister": 0.8}.get(pos.casefold(), 0.5)
    boot = row.get("tree_bootstrap") or ""
    try:
        sup = min(1.0, float(boot) / 100.0) if boot else None
    except ValueError:
        sup = None
    if sup is None:
        src = (row.get("classification_source") or "").casefold()
        sup = {"tree": 0.9, "blast": 0.6, "blast_weak": 0.4, "low_confidence": 0.3}.get(src, 0.5)
    d = row.get("tree_nn_distance") or ""
    try:
        dd = float(d) if d else None
    except ValueError:
        dd = None
    df = min(1.0, 0.05 / max(dd, 0.005)) if dd is not None else 0.5
    return 0.5 * pos_s + 0.3 * sup + 0.2 * df


def load_lucaprot(path: Path | None) -> dict[str, float]:
    if not path or not path.is_file():
        return {}
    out: dict[str, float] = {}
    with path.open(encoding="utf-8-sig", newline="", errors="ignore") as fh:
        first = fh.readline()
        fh.seek(0)
        delim = "," if first.count(",") >= first.count("\t") else "\t"
        reader = csv.DictReader(fh, delimiter=delim)
        fields = [f.strip() for f in (reader.fieldnames or [])]
        if not fields:
            return {}
        fieldmap = {f.casefold(): f for f in fields}
        id_col = fieldmap.get("protein_id") or fieldmap.get("proteinid") or fields[0]
        prob_col = fieldmap.get("prob") or fieldmap.get("probability")
        if prob_col is None:
            return {}
        for row in reader:
            # LucaProt 输出的 protein_id 可能带 FASTA 标记符 ">"，先剥掉再拼 key
            sid = (row.get(id_col) or "").strip().lstrip(">")
            raw = (row.get(prob_col) or "").strip()
            if not sid or not raw:
                continue
            try:
                prob = float(raw)
            except ValueError:
                continue
            key = re.sub(r"_\d+$", "", sid)  # 只剥 prodigal 后缀 _\d+，保留 样本_contig
            if key:
                out[key] = max(0.0, min(1.0, prob))
    return out


def row_key(r: dict[str, str]) -> str:
    """合并总表里用 sample+id 作为唯一序列 key（配 LucaProt 与文献 record_id）。"""
    sample = (r.get("sample") or "").strip()
    rid = (r.get("id") or "").strip()
    return f"{sample}_{rid}" if sample else rid


# ---------- 金标准标签 ----------
def load_labels(path: Path) -> dict[str, str]:
    if not path or not path.is_file():
        return {}
    out: dict[str, str] = {}
    with path.open(encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sp = (r.get("species") or "").strip()
            lab = (r.get("label") or "").strip()
            if sp and lab in {"0", "1", "borderline"}:
                out[sp] = lab
    return out


# ---------- 文献性状（不含人类感染证据） ----------
PLACEHOLDER = {"", "related virus", "co-infecting virus", "unspecified variant"}


def literature_traits(relations: list[dict], contig_ids: set[str]) -> tuple[float, dict]:
    rels = [r for r in relations if str(r.get("record_id")) in contig_ids]
    reassort = any(r.get("relation") == "virus_reassortment" and r.get("polarity") == "positive" for r in rels)
    switching = any(r.get("relation") == "virus_host_switching" and r.get("polarity") == "positive" for r in rels)
    coinf = any(
        r.get("relation") == "virus_coinfects_with" and r.get("polarity") == "positive"
        and r.get("direct_target_match") is True for r in rels
    )
    related = any(
        r.get("relation") == "virus_related_to" and r.get("direct_target_match") is True
        and str(r.get("object")) not in PLACEHOLDER for r in rels
    )
    score = (10 if reassort else 0) + (10 if switching else 0) + (5 if coinf else 0) + (5 if related else 0)
    detail = {"reassortment": reassort, "host_switching": switching, "coinfection": coinf, "related_viruses": related}
    return float(score), detail


# ---------- 传染人（人类感染证据） ----------
HUMAN_RE = re.compile(
    r"\bhuman[- ]to[- ]human\b|\bhuman cases?\b|\bhuman infections?\b|\binfected humans?\b|\bpatients?\b", re.I
)


def pubtator_human_evidence(virus: str) -> tuple[bool, list[str]]:
    query = f'"{virus}" AND human'
    url = PUBTATOR_SEARCH.format(text=urllib.parse.quote(query, safe=""))
    try:
        req = urllib.request.Request(url, headers=HDR)
        with urllib.request.urlopen(req, timeout=45) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        results = (data or {}).get("results", [])[:6]
        pmids: list[str] = []
        found = False
        for r in results:
            pmid = str(r.get("pmid") or "")
            hl = re.sub(r"@\w+_[\w:]+|@@@|</?m>", " ", r.get("text_hl") or "")
            if HUMAN_RE.search(hl):
                found = True
                if pmid:
                    pmids.append(pmid)
        return found, pmids
    except Exception:
        return False, []


# ---------- 排名工具（纯标准库） ----------
def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    rk = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            rk[order[k]] = avg
        i = j + 1
    return rk


def spearman(a: list[float], b: list[float]) -> float:
    if len(a) < 3:
        return float("nan")
    ra, rb = ranks(a), ranks(b)
    ma = sum(ra) / len(ra)
    mb = sum(rb) / len(rb)
    cov = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    va = sum((x - ma) ** 2 for x in ra)
    vb = sum((y - mb) ** 2 for y in rb)
    if va == 0 or vb == 0:
        return float("nan")
    return cov / math.sqrt(va * vb)


def auc(pos: list[float], neg: list[float]) -> float:
    if not pos or not neg:
        return float("nan")
    allv = pos + neg
    rk = ranks(allv)
    rp = sum(rk[: len(pos)])
    u = rp - len(pos) * (len(pos) + 1) / 2.0
    return u / (len(pos) * len(neg))


# ---------- 主流程 ----------
def score_table(input_tsv: Path, outdir: Path, relations: list[dict] | None,
                labels: dict[str, str], lucaprot: dict[str, float], offline: bool,
                query_pmids: dict[str, int]) -> dict:
    rows = list(csv.DictReader(input_tsv.open(encoding="utf-8-sig", newline=""), delimiter="\t"))
    header = list(rows[0].keys()) if rows else []
    if not rows:
        raise ValueError(f"空表：{input_tsv}")

    # 唯一序列 key（sample_id）-> 病毒名/物种
    id_to_virus: dict[str, str] = {}
    id_to_species: dict[str, str] = {}
    for r in rows:
        key = row_key(r)
        if key:
            id_to_virus[key] = (r.get("virus_name") or r.get("species") or "").strip()
            id_to_species[key] = (r.get("species") or r.get("virus_name") or "").strip()

    # 按物种聚合
    by_species: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        sp = (r.get("species") or r.get("virus_name") or "").strip() or "未分类"
        by_species[sp].append(r)

    # 物种级特征
    species_risk: dict[str, dict] = {}
    human_cache: dict[str, tuple[bool, list[str]]] = {}
    for sp, srows in by_species.items():
        hosts = [(r.get("host") or "").strip() for r in srows if (r.get("host") or "").strip()]
        hcls, hdesc = host_class(hosts)
        hscore = HOST_TYPE_POINTS.get(hcls, 0)
        breadth = len({clean_host(h) for h in hosts})
        bscore = 0 if breadth <= 1 else 2 if breadth <= 3 else 5

        # 系统发育：致病种 15；啮齿 Orthohantavirus 10；其余 0；树距离微调
        genus = (srows[0].get("genus") or "").strip()
        phylo = 0.0
        if labels.get(sp) == "1":
            phylo = 15.0
        elif genus == "Orthohantavirus" and hcls == "rodent":
            phylo = 10.0
        dists = []
        for r in srows:
            d = r.get("tree_nn_distance") or ""
            try:
                if d:
                    dists.append(float(d))
            except ValueError:
                pass
        if dists:
            dmin = min(dists)
            phylo *= min(1.0, 0.05 / max(dmin, 0.005))

        # 文献性状
        contigs = {row_key(r) for r in srows if row_key(r)}
        lit, lit_detail = literature_traits(relations or [], contigs)

        # 传染人：VMR 官方人类分离株 或 文献人类感染证据
        human_from_vmr = hcls == "human"
        human_pmids: list[str] = []
        if not offline and not human_from_vmr:
            viruses = sorted({id_to_virus.get(row_key(r), "") for r in srows}) 
            viruses = [v for v in viruses if v]
            viruses = viruses[:2]
            found = False
            for v in viruses:
                if v not in human_cache:
                    human_cache[v] = pubtator_human_evidence(v)
                    time.sleep(0.35)
                f, pms = human_cache[v]
                human_pmids.extend(pms)
                found = found or f
            human_ev = human_from_vmr or found
        else:
            human_ev = human_from_vmr

        hscore_total = hscore + bscore
        human_score = 40.0 if human_ev else 0.0
        risk = hscore_total + phylo + lit + human_score
        risk = min(100.0, max(0.0, risk))
        tier = "High" if risk >= 60 else "Medium" if risk >= 30 else "Low"

        species_risk[sp] = {
            "species": sp, "host_class": hcls, "host_desc": hdesc,
            "host_type_score": hscore, "host_breadth": breadth, "host_breadth_score": bscore,
            "phylogeny_score": round(phylo, 2), "literature_score": lit,
            "literature_detail": lit_detail,
            "human_infection_vmr": human_from_vmr,
            "human_infection_literature": bool(human_ev and not human_from_vmr),
            "human_pmids": sorted(set(human_pmids)),
            "human_infection_score": human_score,
            "species_risk_score": round(risk, 2), "tier": tier,
            "pmid_count": query_pmids.get(sp, 0),
        }

    # 行级：检出可信度 × 物种风险
    out_rows = []
    for r in rows:
        rid = (r.get("id") or "").strip()
        sample = (r.get("sample") or "").strip()
        key = row_key(r)
        sp = (r.get("species") or r.get("virus_name") or "").strip() or "未分类"
        sr = species_risk.get(sp, {})
        risk = sr.get("species_risk_score", 0.0)
        is_L = bool((r.get("vmr_l") or "").strip())
        d = diammond_score(r)
        t = tree_score(r) if is_L else None
        l = lucaprot.get(key) if is_L and lucaprot else None
        w_d, w_t, w_l = 0.35, 0.30, 0.35
        parts = [(w_d, d)]
        flags = []
        if t is not None:
            parts.append((w_t, t))
        else:
            flags.append("no_tree")
        if l is not None:
            parts.append((w_l, l))
        else:
            flags.append("no_lucaprot")
        wsum = sum(w for w, _ in parts)
        conf = sum(w * v for w, v in parts) / wsum if wsum else 0.0
        final = round(conf * risk, 2)
        tier = "High" if final >= 60 else "Medium" if final >= 30 else "Low"
        out_rows.append({
            "sample": sample, "id": rid, "species": sp,
            "virus_name": (r.get("virus_name") or "").strip(),
            "segment": (r.get("source_vmr") or r.get("segment") or "").strip(),
            "detection_diamond": round(d, 3),
            "detection_tree": round(t, 3) if t is not None else "",
            "detection_lucaprot": round(l, 3) if l is not None else "",
            "detection_confidence": round(conf, 3),
            "species_risk_score": risk,
            "final_score": final, "tier": tier,
            "flags": ";".join(flags),
        })

    # 写出
    outdir.mkdir(parents=True, exist_ok=True)
    out_header = header + ["detection_diamond", "detection_tree", "detection_lucaprot",
                           "detection_confidence", "species_risk_score", "final_score",
                           "tier", "flags"]
    rows_full: list[dict] = []
    with (outdir / "results.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_header, delimiter="\t")
        w.writeheader()
        for r, o in zip(rows, out_rows):
            merged = dict(r)
            for k, v in o.items():
                merged[k] = v
            rows_full.append(merged)
            w.writerow({k: merged.get(k, "") for k in out_header})

    with (outdir / "species_risk.tsv").open("w", encoding="utf-8", newline="") as fh:
        fields = ["species", "host_class", "host_desc", "host_type_score", "host_breadth",
                  "host_breadth_score", "phylogeny_score", "literature_score",
                  "human_infection_vmr", "human_infection_literature", "human_pmids",
                  "human_infection_score", "species_risk_score", "tier", "pmid_count"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for sp in sorted(species_risk):
            rec = dict(species_risk[sp])
            rec["human_pmids"] = ";".join(rec["human_pmids"])
            w.writerow({k: rec.get(k, "") for k in fields})

    return {"rows": out_rows, "rows_full": rows_full, "species": species_risk, "input": str(input_tsv)}


def validate(species_risk: dict[str, dict], labels: dict[str, str]) -> dict:
    """用金标准标签验证：致病种召回、AUC、敏感性分析（纯标准库）。"""
    if not labels:
        return {"note": "risk_labels.tsv 缺失或为空，跳过验证"}
    pos, neg = [], []
    pos_sp, neg_sp = [], []
    for sp, lab in labels.items():
        sr = species_risk.get(sp)
        if sr is None:
            continue
        if lab == "1":
            pos.append(sr["species_risk_score"]); pos_sp.append(sp)
        elif lab == "0":
            neg.append(sr["species_risk_score"]); neg_sp.append(sp)
    base_scores = {sp: species_risk[sp]["species_risk_score"] for sp in labels if sp in species_risk}

    # 预测版：去掉传染人项，归一化到 0–100
    pred_scores = {}
    for sp, lab in labels.items():
        sr = species_risk.get(sp)
        if sr is None:
            continue
        base = sr["host_type_score"] + sr["host_breadth_score"] + sr["phylogeny_score"] + sr["literature_score"]
        pred_scores[sp] = round(base * 100.0 / 60.0, 2)

    pred_pos = [pred_scores[sp] for sp in pos_sp]
    pred_neg = [pred_scores[sp] for sp in neg_sp]

    report = {
        "n_labeled": len(labels), "n_positive": len(pos_sp), "n_negative": len(neg_sp),
        "positive_ranks": {sp: sorted(base_scores.values(), reverse=True).index(base_scores[sp]) + 1
                           for sp in pos_sp},
        "recall_high_tier": round(sum(1 for s in pos if s >= 60) / len(pos), 3) if pos else None,
        "auc_full": round(auc(pos, neg), 3) if pos and neg else None,
        "auc_predictive": round(auc(pred_pos, pred_neg), 3) if pred_pos and pred_neg else None,
        "tier_counts": {
            "High": sum(1 for sp in base_scores if base_scores[sp] >= 60),
            "Medium": sum(1 for sp in base_scores if 30 <= base_scores[sp] < 60),
            "Low": sum(1 for sp in base_scores if base_scores[sp] < 30),
        },
    }
    # 敏感性：各组件权重 ±10%
    sens = {}
    for comp, delta in [("host", 0.1), ("phylogeny", 0.1), ("literature", 0.1), ("human", 0.1)]:
        for direction, sign in [("up", 1), ("down", -1)]:
            altered = []
            for sp in base_scores:
                sr = species_risk[sp]
                extra = 0.0
                if comp == "host":
                    extra = sign * delta * (sr["host_type_score"] + sr["host_breadth_score"])
                elif comp == "phylogeny":
                    extra = sign * delta * sr["phylogeny_score"]
                elif comp == "literature":
                    extra = sign * delta * sr["literature_score"]
                else:
                    extra = sign * delta * sr["human_infection_score"]
                altered.append(min(100.0, max(0.0, sr["species_risk_score"] + extra)))
            sens[f"{comp}_{direction}"] = round(spearman(list(base_scores.values()), altered), 3)
    report["sensitivity_spearman"] = sens
    return report


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="5.3 风险打分：检出可信度 × 物种风险分（高中低三档 + 标签验证）")
    p.add_argument("--config", type=str, default="",
                   help="skills_V2/config.yaml 路径；不给则用 --input-tsv 单表模式")
    p.add_argument("--input-tsv", type=Path, default=None, help="带树列的 vmr_final_classified.tsv")
    p.add_argument("--outdir", type=Path, default=None, help="输出目录；默认输入同目录")
    p.add_argument("--relations", type=Path, default=None, help="extracted_relations.json 路径（可选）")
    p.add_argument("--queries", type=Path, default=None, help="retrieval_queries.tsv 路径（可选）")
    p.add_argument("--labels", type=Path,
                   default=Path(__file__).resolve().parent / "risk_labels.tsv",
                   help="金标准标签 TSV")
    p.add_argument("--lucaprot", type=Path, default=None, help="LucaProt 输出 TSV（序列id + RdRp 概率，可选）")
    p.add_argument("--qc-dir", type=Path, default=None, help="fastp 质控报告目录（默认 results_dir/01_fastp）")
    p.add_argument("--offline", action="store_true", help="跳过 PubTator3 实时检索（传染人项只认 VMR 人类分离株）")
    return p.parse_args()


def load_relations(path: Path | None) -> list[dict] | None:
    if not path or not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def load_query_pmids(path: Path | None) -> dict[str, int]:
    if not path or not path.is_file():
        return {}
    out: dict[str, int] = defaultdict(int)
    with path.open(encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            v = (r.get("virus_name") or "").strip()
            try:
                out[v] += int(r.get("pmid_count") or 0)
            except ValueError:
                pass
    return dict(out)


def merge_relations(paths: list[Path]) -> list[dict]:
    """合并多个样本的 extracted_relations.json，record_id 加样本名前缀避免跨样本串号。"""
    out: list[dict] = []
    for p in paths:
        if not p.is_file():
            continue
        try:
            data = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            continue
        if not isinstance(data, list):
            continue
        sample = p.parent.name
        for rel in data:
            if not isinstance(rel, dict):
                continue
            rel = dict(rel)
            rid = str(rel.get("record_id") or "")
            if rid:
                rel["record_id"] = f"{sample}_{rid}"
            out.append(rel)
    return out


def merge_query_pmids(paths: list[Path]) -> dict[str, int]:
    """合并多个样本的 retrieval_queries.tsv，按 virus_name 累加 pmid_count。"""
    out: dict[str, int] = defaultdict(int)
    for p in paths:
        if not p.is_file():
            continue
        with p.open(encoding="utf-8-sig", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                v = (r.get("virus_name") or "").strip()
                try:
                    out[v] += int(r.get("pmid_count") or 0)
                except ValueError:
                    pass
    return dict(out)


def main() -> int:
    args = parse_args()

    if args.config:
        import yaml
        cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
        results_dir = Path(cfg["paths"]["results_dir"])
        input_tsv = results_dir / "00_results" / "all_vmr_final_classified.tsv"
        outdir = results_dir / "00_results"
        # 默认读 LucaProt 全量输出（l.fasta.csv，含全部 L 片段与 prob 列）
        lucaprot_tsv = args.lucaprot or (results_dir / "00_results" / "l.fasta.csv")
        # 质控目录：paths.fastp_dir 优先，否则 results_dir/01_fastp（与 1-qc.sh 的 FASTP_DIR 规则一致）
        qc_dir = args.qc_dir or Path(cfg["paths"].get("fastp_dir") or (results_dir / "01_fastp"))
        if not input_tsv.is_file():
            print(f"未找到合并总表：{input_tsv}", file=sys.stderr)
            return 2
        ew_dir = results_dir / "08_early_warning"
        relations = merge_relations(sorted(ew_dir.glob("*/extracted_relations.json")))
        query_pmids = merge_query_pmids(sorted(ew_dir.glob("*/retrieval_queries.tsv")))
        targets = [(input_tsv, outdir)]
    else:
        if not args.input_tsv or not args.input_tsv.is_file():
            print("需要 --input-tsv（或 --config 走批处理）", file=sys.stderr)
            return 2
        outdir = args.outdir or args.input_tsv.parent
        lucaprot_tsv = args.lucaprot
        qc_dir = args.qc_dir or (args.input_tsv.parent.parent / "01_fastp")
        relations = load_relations(args.relations)
        query_pmids = load_query_pmids(args.queries)
        targets = [(args.input_tsv, outdir)]

    if not relations:
        print("[警告] 未读到文献证据，文献性状与传染人按无处理", file=sys.stderr)

    labels = load_labels(args.labels)
    lucaprot = load_lucaprot(lucaprot_tsv)
    summaries = []
    for tsv, outdir in targets:
        result = score_table(tsv, outdir, relations, labels, lucaprot, args.offline, query_pmids)
        # 5.3 只负责打分，HTML 报告由 5.4-report.py 读本中间文件生成
        bundle = {
            "rows": result["rows"],
            "rows_full": result.get("rows_full", []),
            "species": result["species"],
            "relations": relations,
            "qc_dir": str(qc_dir) if qc_dir else None,
        }
        bundle_path = outdir / "report_inputs.json"
        bundle_path.write_text(
            json.dumps(bundle, ensure_ascii=False, default=str, indent=2), encoding="utf-8")
        print(json.dumps({"input": str(tsv), "outdir": str(outdir),
                          "row_count": len(result["rows"]),
                          "species_count": len(result["species"]),
                          "bundle": str(bundle_path)}, ensure_ascii=False))
        summaries.append(str(outdir))
    print(json.dumps({"done": len(summaries), "outdirs": summaries}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
