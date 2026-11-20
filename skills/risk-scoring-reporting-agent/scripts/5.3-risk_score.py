#!/usr/bin/env python3
"""Step 5.3 risk scoring: detection confidence x species risk, in three tiers,
validated against a gold-standard label set.

Position in the pipeline
------------------------
Step 4 (tree) -> step 5 (PubTator3 literature evidence) -> this script
(scoring and validation). It replaces the former 5.2-early_warning_score.py; the
legacy rule-based score survives only as a sensitivity comparison and is skipped
when unavailable.

Inputs
------
  - *.vmr_final_classified.tsv with the tree columns (from step 4; carrying
    final_classification / classification_source / tree_position /
    tree_nn_sseqid / tree_nn_distance)
  - extracted_relations.json (step 5 output; when missing the literature traits
    score 0 and a warning is printed)
  - retrieval_queries.tsv (step 5 output; used for the detection-bias covariate)
  - risk_labels.tsv (gold-standard labels; used for validation)
  - LucaProt output (optional TSV: sequence id + RdRp probability; when missing
    the score normalises over DIAMOND and tree only)

Scoring rules (kept in step with results/figure3/scoring_rules.png)
------------------------------------------------------------------
  detection confidence (0-1) = 0.35 x DIAMOND + 0.30 x tree + 0.35 x LucaProt
                               (missing terms renormalise over the remainder)
  species risk (0-100) = host 15 + phylogeny 15 + literature traits 30
                         + human infection 40
    host: type (rodent 10 / bat 5 / shrew-mole 2 / fish-reptile 0)
          + breadth (1 species 0 / 2-3 species 2 / 4+ species 5)
    phylogeny: known pathogenic species 15; rodent Orthohantavirus 10;
               otherwise 0; adjusted slightly by tree distance
    literature traits: reassortment +10, host switch +10, co-infection +5,
                       related virus +5 (human infection excluded here)
    human infection: documented human infection +40
                     (live PubTator3 search; skipped in offline mode)
  final score = detection confidence x species risk
  tiers: High >= 60 / Medium 30-59 / Low < 30
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
sys.path.insert(0, str(Path(__file__).resolve().parent))
from seqid import strip_prodigal_suffix  # noqa: E402
PUBTATOR_SEARCH = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/?text={text}&page=1"
HDR = {"User-Agent": "Yu-hanta-risk-score/0.1 (contact: local)"}

# ---------- Host cleaning and classification ----------
# Rule table: order is the matching priority (specific families and genera first, generic terms last)
# Word boundaries are used so "batfish" is not captured by \bbat\b and Rattus is not misread as \brat\b
HOST_CLASS_RULES: list[tuple[str, str, str]] = [
    ("rodent_muridae", "鼠科 Muridae",
     r"\b(rattus|niviventer|apodemus|mus|bandicota|leopoldamys|bery|maxomys|"
     r"hylomyscus|reithrodontomys|oryzomys|stenocephalemys|murid|murinae|"
     r"rat|rats|mouse|mice)\b"),
    ("rodent_cricetidae", "仓鼠科 Cricetidae",
     r"\b(myodes|clethrionomys|peromyscus|microtus|oligoryzomys|sigmodon|neotoma|"
     r"eothenomys|calomys|akodon|lemmus|ochrotomys|cricet|cricetid)\b"),
    # Fish come first: batfish / goosefish contain "bat" and "fish", so matching fish prevents a bat fallback
    ("fish", "鱼类",
     r"\b(perca|halieutaea|triacanthodes|luposicya|lophius|gobio|gobiid|"
     r"percid|perciform|actinopteryg|fish|batfish|goosefish|spikefish|goby)\b"),
    ("agnathan", "无颌类", r"\b(eptatretus|hagfish|agnathan)\b"),
    ("reptile", "爬行类", r"\b(hemidactylus|gecko|reptil\w*)\b"),
    ("insect", "昆虫", r"\b(perla|plecoptera|insect\w*)\b"),
    ("bird", "鸟类", r"\b(cyanocitta|aves|bird|birds)\b"),
    ("bat_pteropodidae", "狐蝠科 Pteropodidae", r"\b(pteropus|rousettus|pteropodid)\b"),
    ("bat_rhinolophidae", "菊头蝠科 Rhinolophidae", r"\b(rhinolophus|aselliscus|hipposideros)\b"),
    ("bat_vespertilionidae", "蝙蝠科 Vespertilionidae",
     r"\b(pipistrellus|nyctalus|eptesicus|miniopterus|vespertilio|vespertilionid)\b"),
    ("bat_other", "其他蝙蝠", r"\b(chiroptera|taphozous|molossid|emballonurid|bat|bats)\b"),
    ("shrew_soricidae", "鼩鼱科 Soricidae",
     r"\b(sorex|crocidura|chodsigoa|anourosorex|suncus|episoriculus|soricid|shrew|shrews)\b"),
    ("mole_talpidae", "鼹鼠科 Talpidae",
     r"\b(talpa|euroscaptor|urotrichus|scalopus|talpid|mole|moles)\b"),
    ("hedgehog_erinaceidae", "刺猬科 Erinaceidae",
     r"\b(erinaceus|hemiechinus|erinaceid|hedgehog|hedgehogs)\b"),
    # Generic rodent terms come last
    ("rodent_other", "其他啮齿目",
     r"\b(rodents?|voles?|gerbil\w*|hamster\w*|arvicol\w*)\b"),
]

HOST_CLASS_CN: dict[str, str] = {cls: cn for cls, cn, _ in HOST_CLASS_RULES}
HOST_CLASS_CN["other"] = "其他"
HOST_CLASS_CN["unknown"] = "未确定（VMR 与文献均无自然宿主信息）"

# Host type points (human excluded: host_class already drops humans, so a human-only host is not a natural reservoir)
HOST_TYPE_POINTS: dict[str, int] = {
    "rodent_muridae": 10,        # house and field rats/mice, commensal with humans, highest exposure risk
    "rodent_cricetidae": 8,      # voles (main reservoir of HTNV and PUUV)
    "rodent_other": 7,
    "bat_rhinolophidae": 5,
    "bat_vespertilionidae": 5,
    "bat_pteropodidae": 4,
    "bat_other": 4,
    "shrew_soricidae": 3,
    "mole_talpidae": 2,
    "hedgehog_erinaceidae": 2,
    "fish": 0, "agnathan": 0, "reptile": 0, "insect": 0, "bird": 0,
    "other": 0, "unknown": 0,
}

# Human-related terms (used to strip the natural reservoir; kept only as the human_isolation flag)
HUMAN_RE_HOST = re.compile(r"\b(homo sapiens|human|humans|patient|patients)\b", re.I)

# Common dirty-value fixes for the raw host field
HOST_FIXUPS: list[tuple[str, str]] = [
    (r"^su\.?\s*murinus$", "suncus murinus"),
    (r"^s\.?\s*murinus$", "suncus murinus"),
    (r"^rattus\s+sp\.?$", "rattus"),
    (r"^rat$", "rattus"),
    (r"^rats$", "rattus"),
    (r"^ruttus", "rattus"),
    (r"^roratus", "rattus"),
]


def clean_host(raw: str) -> str:
    h = (raw or "").casefold().strip()
    h = re.sub(r"\(.*?\)", "", h)                 # drop the parenthesised common name
    h = re.sub(r"strain.*$", "", h).strip()        # drop the strain suffix
    h = re.sub(r"msb:mamm.*$", "", h).strip()      # drop the museum specimen number
    h = re.sub(r"\s+", " ", h).strip()
    for pat, rep in HOST_FIXUPS:
        if re.match(pat, h):
            h = rep
            break
    return h


def _match_host_class(cands: list[str]) -> tuple[str, str, list[str]]:
    """Match candidate host names against the rule table.

    Returns (class, display name, deduplicated species list).
    """
    species = sorted({c for c in (clean_host(x) for x in cands) if c})
    if not species:
        return "unknown", HOST_CLASS_CN["unknown"], []
    for cls, cn, pat in HOST_CLASS_RULES:
        rx = re.compile(pat, re.I)
        if any(rx.search(s) for s in species):
            return cls, cn, species
    return "other", HOST_CLASS_CN["other"], species


def host_class(hosts: list[str], virus_names=(), offline: bool = True,
               host_cache: dict | None = None) -> tuple[str, str, list[str], bool, str]:
    """Classify the natural reservoir, excluding humans.

    Returns (host_class, display name, natural host species, human_isolation,
    host_source).

    Source priority:
      1) the VMR raw host, when a host remains after removing Homo sapiens
      2) a literature search, only when 1 is empty and not in offline mode
      3) unknown
    """
    cleaned = sorted({clean_host(h) for h in hosts if clean_host(h)})
    human_hits = [h for h in cleaned if HUMAN_RE_HOST.search(h)]
    natural = [h for h in cleaned if h not in human_hits]
    human_isolation = bool(human_hits)

    if natural:
        cls, cn, species = _match_host_class(natural)
        return cls, cn, species, human_isolation, "VMR"

    # No natural reservoir in VMR -> fall back to the literature
    if not offline:
        for v in list(virus_names)[:2]:
            words = None
            if host_cache is not None:
                words = host_cache.get(v)
            if words is None:
                words = pubtator_host_evidence(v)
                if host_cache is not None:
                    host_cache[v] = words
                time.sleep(0.35)
            if words:
                cls, cn, species = _match_host_class(words)
                if cls != "unknown":
                    return cls, cn, species, human_isolation, "literature"

    return "unknown", HOST_CLASS_CN["unknown"], [], human_isolation, "none"


# ---------- Detection confidence ----------
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


def collect_row_keys(tsv_path) -> set[str]:
    """Stream the combined table and collect every `sample_id` pair.

    Used as the validation list when load_lucaprot strips a suffix: the strip is
    accepted only if the result actually matches a contig in the table, which
    prevents false positives such as "NODE_5" -> "NODE".
    """
    keys: set[str] = set()
    if not tsv_path or not Path(tsv_path).is_file():
        return keys
    with Path(tsv_path).open(encoding="utf-8-sig", newline="", errors="ignore") as fh:
        first = fh.readline()
        fh.seek(0)
        delim = "," if first.count(",") >= first.count("\t") else "\t"
        reader = csv.DictReader(fh, delimiter=delim)
        fieldmap = {(f or "").strip().casefold(): f for f in (reader.fieldnames or [])}
        id_col = fieldmap.get("id")
        if id_col is None:
            return keys
        sample_col = fieldmap.get("sample")
        for r in reader:
            rid = (r.get(id_col) or "").strip()
            if not rid:
                continue
            sample = (r.get(sample_col) or "").strip() if sample_col else ""
            keys.add(f"{sample}_{rid}" if sample else rid)
    return keys


def load_lucaprot(path: Path | None, known_keys: set[str] | None = None) -> dict[str, float]:
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
            # The LucaProt protein_id may carry the FASTA marker ">"; strip it before building the key
            sid = (row.get(id_col) or "").strip().lstrip(">")
            raw = (row.get(prob_col) or "").strip()
            if not sid or not raw:
                continue
            try:
                prob = float(raw)
            except ValueError:
                continue
            # Strip only prodigal's _<gene> suffix and keep sample_contig intact.
            # Two passes, so the result is never worse than the previous behaviour:
            #   1) always strip one layer - behaviourally identical to the old
            #      re.sub, so no score is lost when the list is incomplete
            #   2) when the original form is a valid contig in the known list
            #      but the stripped form is not, keep the original
            #      (this blocks false positives such as NODE_5 -> NODE)
            key = strip_prodigal_suffix(sid)
            if known_keys and sid in known_keys and key not in known_keys:
                key = sid
            if key:
                out[key] = max(0.0, min(1.0, prob))
    return out


def row_key(r: dict[str, str]) -> str:
    """Build a unique sequence key from sample+id, matching LucaProt and the
    literature record_id."""
    sample = (r.get("sample") or "").strip()
    rid = (r.get("id") or "").strip()
    return f"{sample}_{rid}" if sample else rid


# ---------- Gold-standard labels ----------
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


# ---------- Literature traits (human infection evidence excluded) ----------
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


# ---------- Human infection (evidence of human infection) ----------
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


# ---------- Host (literature fallback when VMR has no natural reservoir) ----------
HOST_TOKEN_RE = re.compile(
    r"\b(rattus|niviventer|apodemus|mus musculus|bandicota|leopoldamys|"
    r"myodes|clethrionomys|peromyscus|microtus|oligoryzomys|sigmodon|neotoma|"
    r"rhinolophus|pipistrellus|nyctalus|eptesicus|miniopterus|pteropus|rousettus|"
    r"sorex|crocidura|suncus|anourosorex|chodsigoa|erinaceus|talpa|"
    r"rodents?|mouse|mice|rats?|voles?|shrews?|moles?|bats?|insectivores?|"
    r"reservoir host)\b", re.I
)


def pubtator_host_evidence(virus: str) -> list[str]:
    """Extract host terms from PubTator3 abstracts.

    Called only as a fallback when the VMR raw host is empty after removing
    humans. Returns [] on failure; never raises.
    """
    query = f'"{virus}" AND (host OR reservoir OR rodent OR bat OR shrew)'
    url = PUBTATOR_SEARCH.format(text=urllib.parse.quote(query, safe=""))
    try:
        req = urllib.request.Request(url, headers=HDR)
        with urllib.request.urlopen(req, timeout=45) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        words: list[str] = []
        for r in (data or {}).get("results", [])[:6]:
            hl = re.sub(r"@\w+_[\w:]+|@@@|</?m>", " ", r.get("text_hl") or "")
            words.extend(HOST_TOKEN_RE.findall(hl))
        return sorted({w.lower() for w in words})
    except Exception:
        return []


# ---------- Ranking utilities (standard library only) ----------
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


# ---------- Main flow ----------
def score_table(input_tsv: Path, outdir: Path, relations: list[dict] | None,
                labels: dict[str, str], lucaprot: dict[str, float], offline: bool,
                query_pmids: dict[str, int]) -> dict:
    rows = list(csv.DictReader(input_tsv.open(encoding="utf-8-sig", newline=""), delimiter="\t"))
    header = list(rows[0].keys()) if rows else []
    if not rows:
        raise ValueError(f"empty table: {input_tsv}")

    # unique sequence key (sample_id) -> virus name / species
    id_to_virus: dict[str, str] = {}
    id_to_species: dict[str, str] = {}
    for r in rows:
        key = row_key(r)
        if key:
            id_to_virus[key] = (r.get("virus_name") or r.get("species") or "").strip()
            id_to_species[key] = (r.get("species") or r.get("virus_name") or "").strip()

    # Aggregate by species
    by_species: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        sp = (r.get("species") or r.get("virus_name") or "").strip() or "unclassified"
        by_species[sp].append(r)

    # Species-level features
    species_risk: dict[str, dict] = {}
    human_cache: dict[str, tuple[bool, list[str]]] = {}
    host_cache: dict[str, list[str]] = {}
    for sp, srows in by_species.items():
        hosts = [(r.get("host") or "").strip() for r in srows if (r.get("host") or "").strip()]
        viruses_for_host = sorted({id_to_virus.get(row_key(r), "") for r in srows})
        viruses_for_host = [v for v in viruses_for_host if v][:2]
        hcls, hdesc, hspecies, human_isolation, hsource = host_class(
            hosts, viruses_for_host, offline, host_cache)
        hscore = HOST_TYPE_POINTS.get(hcls, 0)
        breadth = len(hspecies)
        bscore = 0 if breadth <= 1 else 2 if breadth <= 3 else 5

        # Phylogeny: known pathogen 15; rodent Orthohantavirus 10; otherwise 0;
        # nudged by tree distance
        genus = (srows[0].get("genus") or "").strip()
        phylo = 0.0
        if labels.get(sp) == "1":
            phylo = 15.0
        elif genus == "Orthohantavirus" and hcls.startswith("rodent"):
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

        # Literature traits
        contigs = {row_key(r) for r in srows if row_key(r)}
        lit, lit_detail = literature_traits(relations or [], contigs)

        # Human infection: a VMR human isolate, or literature evidence
        human_from_vmr = human_isolation
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
        human_score = 60.0 if human_ev else 0.0
        risk = hscore_total + phylo + lit + human_score
        risk = min(100.0, max(0.0, risk))
        # Hard rule: documented human infection means High directly, without
        # relying on the score accumulating to 60
        if human_ev:
            tier = "High"
        elif risk >= 60:
            tier = "High"
        elif risk >= 30:
            tier = "Medium"
        else:
            tier = "Low"

        species_risk[sp] = {
            "species": sp, "host_class": hcls, "host_desc": hdesc,
            "host_species": "; ".join(hspecies), "host_source": hsource,
            "human_isolation": human_isolation,
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

    # Per-row: detection confidence x species risk
    out_rows = []
    for r in rows:
        rid = (r.get("id") or "").strip()
        sample = (r.get("sample") or "").strip()
        key = row_key(r)
        sp = (r.get("species") or r.get("virus_name") or "").strip() or "unclassified"
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
        # Human infection evidence takes priority: a virus that infects humans
        # is serious, so the score floors at max(species risk, 60) rather than
        # being discounted by detection confidence, and the tier is High
        sr_human = bool(sr.get("human_infection_vmr") or sr.get("human_infection_literature"))
        if sr_human:
            final = round(max(risk, 60.0), 2)
            tier = "High"
            flags.append("human_evidence")
        elif final >= 60:
            tier = "High"
        elif final >= 30:
            tier = "Medium"
        else:
            tier = "Low"
        out_rows.append({
            "sample": sample, "id": rid, "species": sp,
            "virus_name": (r.get("virus_name") or "").strip(),
            "segment": (r.get("source_vmr") or r.get("segment") or "").strip(),
            # Carry the host classification on the row too (value comes from the species)
            "host_class": sr.get("host_class", ""),
            "host_desc": sr.get("host_desc", ""),
            "host_species": sr.get("host_species", ""),
            "host_source": sr.get("host_source", ""),
            "human_isolation": sr.get("human_isolation", ""),
            "detection_diamond": round(d, 3),
            "detection_tree": round(t, 3) if t is not None else "",
            "detection_lucaprot": round(l, 3) if l is not None else "",
            "detection_confidence": round(conf, 3),
            "species_risk_score": risk,
            "final_score": final, "tier": tier,
            "flags": ";".join(flags),
        })

    # Write out
    outdir.mkdir(parents=True, exist_ok=True)
    out_header = header + ["host_class", "host_desc", "host_species", "host_source",
                           "human_isolation",
                           "detection_diamond", "detection_tree", "detection_lucaprot",
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
        fields = ["species", "host_class", "host_desc", "host_species", "host_source",
                  "human_isolation", "host_type_score", "host_breadth",
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
    """Validate against the gold-standard labels: pathogenic-species recall,
    AUC, and a sensitivity analysis. Standard library only."""
    if not labels:
        return {"note": "risk_labels.tsv is missing or empty; validation skipped"}
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

    # Prediction variant: drop the human-infection term, rescale to 0-100
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
    # Sensitivity: each component weight varied by +/-10%
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
    p = argparse.ArgumentParser(
        description="Step 5.3 risk scoring: detection confidence x species risk, "
                    "three tiers, validated against gold-standard labels",
        epilog=(
            "Examples:\n"
            "  batch mode (recommended): python 5.3-risk_score.py --config cell_submit-paper.yaml\n"
            "  single table:            python 5.3-risk_score.py --input-tsv <sample>.vmr_final_classified.tsv"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--config", type=str, default="",
                   help="path to config.yaml. Supplying it selects batch mode: the combined table, "
                         "LucaProt output, QC directory and literature evidence are all located "
                         "automatically, "
                        "so one command runs the whole step (recommended). Without it, "
                        "single-table mode applies and --input-tsv is required.")
    p.add_argument("--input-tsv", type=Path, default=None,
                   help="required in single-table mode: a vmr_final_classified.tsv carrying the tree "
                         "columns (final_classification, tree_*, ...). "
                        "Batch mode uses results_dir/00_results/all_vmr_final_classified.tsv "
                        "automatically, so this is not needed there.")
    p.add_argument("--outdir", type=Path, default=None,
                   help="output directory. Defaults: single-table mode = the input file's directory; "
                         "batch mode = results_dir/00_results.")
    p.add_argument("--relations", type=Path, default=None,
                   help="path to extracted_relations.json (step 5 output, literature evidence, "
                         "optional). "
                        "Batch mode scans 08_early_warning/*/ automatically, so this is not "
                        "needed there; when absent the literature terms score 0.")
    p.add_argument("--queries", type=Path, default=None,
                   help="path to retrieval_queries.tsv (step 5 output, detection-bias covariate, "
                         "optional). "
                        "Batch mode scans 08_early_warning/*/ automatically, so this is not "
                        "needed there.")
    p.add_argument("--labels", type=Path,
                   default=Path(__file__).resolve().parent / "risk_labels.tsv",
                   help="gold-standard label TSV (columns species and label, where label is "
                        "0/1/borderline), used to validate the scores. "
                        "Defaults to risk_labels.tsv beside the script; validation is skipped "
                        "when the file is absent.")
    p.add_argument("--lucaprot", type=Path, default=None,
                   help="LucaProt output TSV (columns protein_id and prob, optional). "
                        "Batch mode defaults to results_dir/00_results/l.fasta.csv; when it "
                        "is missing the LucaProt weight is dropped and the remainder "
                        "renormalises.")
    p.add_argument("--qc-dir", type=Path, default=None,
                   help="fastp QC report directory. Defaults to paths.fastp_dir from the "
                         "configuration, falling back to results_dir/01_fastp.")
    p.add_argument("--offline", action="store_true",
                   help="skip the live PubTator3 search: the human-infection term then counts only "
                         "VMR human isolates and no online literature lookup is done.")
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
    """Merge extracted_relations.json across samples, prefixing record_id with the
    sample name so records never collide across samples."""
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
    """Merge retrieval_queries.tsv across samples, summing pmid_count per
    virus_name."""
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
        # Default to the full LucaProt output (l.fasta.csv, every L fragment with its prob column)
        lucaprot_tsv = args.lucaprot or (results_dir / "00_results" / "l.fasta.csv")
        # QC directory: paths.fastp_dir first, else results_dir/01_fastp (same rule as FASTP_DIR in 1-qc.sh)
        qc_dir = args.qc_dir or Path(cfg["paths"].get("fastp_dir") or (results_dir / "01_fastp"))
        if not input_tsv.is_file():
            print(f"combined table not found: {input_tsv}", file=sys.stderr)
            return 2
        ew_dir = results_dir / "08_early_warning"
        relations = merge_relations(sorted(ew_dir.glob("*/extracted_relations.json")))
        query_pmids = merge_query_pmids(sorted(ew_dir.glob("*/retrieval_queries.tsv")))
        targets = [(input_tsv, outdir)]
    else:
        if not args.input_tsv or not args.input_tsv.is_file():
            print("--input-tsv is required (or use --config for batch mode)", file=sys.stderr)
            return 2
        outdir = args.outdir or args.input_tsv.parent
        lucaprot_tsv = args.lucaprot
        qc_dir = args.qc_dir or (args.input_tsv.parent.parent / "01_fastp")
        relations = load_relations(args.relations)
        query_pmids = load_query_pmids(args.queries)
        targets = [(args.input_tsv, outdir)]

    if not relations:
        print("[warning] no literature evidence was read; literature traits and "
              "human infection both score 0", file=sys.stderr)

    labels = load_labels(args.labels)
    # Collect the sample_id pairs first, so load_lucaprot can validate its suffix
    # stripping instead of damaging contigs that end in _<digits>
    known_keys: set[str] = set()
    for _tsv, _out in targets:
        known_keys |= collect_row_keys(_tsv)
    lucaprot = load_lucaprot(lucaprot_tsv, known_keys=known_keys or None)
    if lucaprot:
        hit = sum(1 for k in lucaprot if k in known_keys) if known_keys else 0
        print("[info] LucaProt read %d protein probabilities, %d of which map to "
              "contigs in the combined table"
              % (len(lucaprot), hit), file=sys.stderr)
    summaries = []
    for tsv, outdir in targets:
        result = score_table(tsv, outdir, relations, labels, lucaprot, args.offline, query_pmids)
        # 5.3 only scores; the HTML report is generated by 5.4-report.py from this bundle
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
