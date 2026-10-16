#!/usr/bin/env python3
"""Read two VMR records and build a test literature-evidence bundle."""

from __future__ import annotations

import csv
import json
import re
import time
import yaml
import argparse
from datetime import date
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

INPUT: Path | None = None
OUT: Path | None = None

# ====================== 命令行参数解析 ======================
def parse_args():
    parser = argparse.ArgumentParser(description="PubTator3 文献证据提取脚本")
    parser.add_argument(
        "--config",
        type=str,
        default=str(Path(__file__).resolve().parents[2] / "skills_V2" / "config.yaml"),
        help="配置文件 config.yaml 的路径，默认与原路径一致"
    )
    return parser.parse_args()

# ====================== 加载配置文件 ======================
def load_config(config_path: str) -> dict:
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    return config


def sample_output_dir(input_tsv: Path, input_root: Path, output_root: Path) -> Path:
    relative = input_tsv.relative_to(input_root)
    sample_name = relative.name.removesuffix(".vmr_final_classified.tsv")
    return output_root / relative.parent / sample_name
PUBTATOR_SEARCH = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/?text={text}&page=1"
PUBTATOR_BIOC = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson?pmids={pmids}"
EUTILS_SEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=10&term={term}"
KEYWORDS = re.compile(r"\b(detected|isolated|host|reservoir|spillover|mutation|variant|substitution|nucleocapsid|segment\s+S|reassortment|coinfection)\b", re.I)

# Explicit triggers only. These rules deliberately do not infer causality from co-occurrence.
RELATION_RULES = {
    "virus_detected_in_host": re.compile(
        r"\b(detected in|detection of|identified in|found in|harbored by|isolated from|screened.*\bpositive\b|RNA detected)\b", re.I
    ),
    "virus_associated_with_host": re.compile(
        r"\b(natural host|reservoir host|host species|hosted by|carried by)\b", re.I
    ),
    "virus_has_gene_or_segment": re.compile(
        r"\b(nucleocapsid protein|nucleoprotein|segment\s+[SML]|genomic\s+[SML]\s+and\s+[SML]\s+segments?)\b", re.I
    ),
    "virus_coinfects_with": re.compile(r"\b(co[- ]?infect(?:ions?|ed)?|coinfect(?:ions?|ed)?)\b", re.I),
    "virus_related_to": re.compile(r"\b(closest relatives?|closely related|common evolutionary origin|distantly related)\b", re.I),
    "virus_reassortment": re.compile(r"\b(reassortment|reassorted)\b", re.I),
    "virus_host_switching": re.compile(r"\b(host[- ]switching|cross[- ]species transmission)\b", re.I),
    "variant_reported": re.compile(r"\b(mutation|variant|substitution|amino acid (?:change|exchange))\b", re.I),
}
NEGATION = re.compile(r"\b(no|not|without|lack of|absence of|failed to)\b.{0,55}\b(evidence|detect|identify|find|reassortment|mutation|variant)\b", re.I)
UNCERTAINTY = re.compile(r"\b(may|might|could|possible|possibly|suggest(?:s|ed)?|appears?|putative|potential)\b", re.I)
VARIANT_PATTERN = re.compile(r"\b(?:[A-Z][a-z]{2}|[A-Z])\d{1,5}(?:[A-Z][a-z]{2}|[A-Z]|\*)\b")


def get_json(url: str) -> object:
    request = Request(url, headers={"User-Agent": "Yu-hanta-literature-test/0.1 (contact: local)"})
    with urlopen(request, timeout=45) as response:
        return json.loads(response.read().decode("utf-8"))


def pubtator_pmids(query: str) -> list[str]:
    try:
        data = get_json(PUBTATOR_SEARCH.format(text=quote(query, safe="")))
        return [str(item["pmid"]) for item in data.get("results", []) if item.get("pmid")]
    except Exception:
        return []


def pubmed_pmids(query: str) -> list[str]:
    try:
        data = get_json(EUTILS_SEARCH.format(term=quote(query, safe="")))
        return [str(pmid) for pmid in data["esearchresult"].get("idlist", [])]
    except Exception:
        return []


def extract_passages(pmid: str, virus: str, host: str, gene_hint: str) -> list[dict[str, str]]:
    try:
        docs = get_json(PUBTATOR_BIOC.format(pmids=quote(pmid, safe=",")))
    except Exception as error:
        return [{"pmid": pmid, "section": "", "passage": "", "matched_terms": "", "error": str(error)}]

    if isinstance(docs, dict):
        docs = docs.get("PubTator3", [])
    if not isinstance(docs, list):
        return [{"pmid": pmid, "section": "", "passage": "", "matched_terms": "", "entity_ids": "", "error": "Unexpected BioCJSON root"}]

    results = []
    for doc in docs:
        for passage in doc.get("passages", []):
            text = passage.get("text", "").strip()
            if not text:
                continue
            lowered = text.casefold()
            matched = [term for term in (virus, host, gene_hint) if term and term.casefold() in lowered]
            if matched or KEYWORDS.search(text):
                infons = passage.get("infons", {})
                section = str(infons.get("section") or infons.get("section_type") or infons.get("type") or "")
                entities = []
                for annotation in passage.get("annotations", []):
                    entity_infons = annotation.get("infons", {})
                    entity_id = entity_infons.get("accession") or entity_infons.get("normalized_id") or entity_infons.get("identifier")
                    entity_type = entity_infons.get("type", "")
                    if entity_id:
                        entities.append(f"{annotation.get('text', '')}|{entity_type}|{entity_id}")
                results.append({
                    "pmid": pmid,
                    "section": section,
                    "passage": text,
                    "matched_terms": "; ".join(matched),
                    "entity_ids": "; ".join(entities),
                    "error": "",
                })
    return results


def clean(value: str | None) -> str:
    return (value or "").replace("\r", "").strip()


def entity_matches(passage: dict[str, str], name: str) -> bool:
    return bool(name) and name.casefold() in passage["passage"].casefold()


def entities_from_passage(entity_ids: str) -> list[dict[str, str]]:
    entities = []
    for raw in filter(None, entity_ids.split("; ")):
        parts = raw.split("|", 2)
        if len(parts) == 3:
            entities.append({"text": parts[0], "type": parts[1], "id": parts[2]})
    return entities


GENERIC_VIRUS_NAMES = {
    "hantavirus", "hantaviruses", "orthohantavirus", "orthohantaviruses",
    "coronavirus", "coronaviruses", "bunyavirus", "arenavirus",
}


def is_virus_name(name: str) -> bool:
    """Heuristic: a Species entity is a specific virus if it is '<X> virus' form.

    Requires ' virus' (space + virus) so genus/family names such as
    'Mobatvirus' or 'Hantavirus' (a single concatenated word) are excluded.
    """
    return " virus" in name.casefold()


def related_virus_in_sentence(sentence: str, entities: list[dict[str, str]], exclude: str) -> dict[str, str] | None:
    """Pick another specific virus named in the sentence, excluding the subject virus.

    Uses PubTator3 normalized ids to resolve abbreviations (e.g. 'ACDV' -> the
    full 'Academ virus' entity sharing the same id).
    """
    lowered = sentence.casefold()
    exclude_id = next((entity["id"] for entity in entities if entity["text"].casefold() == exclude.casefold()), None)
    id_to_fullname = {}
    for entity in entities:
        if entity["type"] == "Species" and is_virus_name(entity["text"]) and entity["id"]:
            id_to_fullname.setdefault(entity["id"], entity["text"])

    for entity in entities:
        if entity["type"] != "Species":
            continue
        if entity["text"].casefold() not in lowered:
            continue
        if entity["id"] and entity["id"] == exclude_id:
            continue
        if is_virus_name(entity["text"]):
            return entity
        if entity["id"] in id_to_fullname:
            return {"text": id_to_fullname[entity["id"]], "type": "Species", "id": entity["id"]}
    return None


def trigger_sentence(text: str, start: int, end: int) -> str:
    """Return the sentence containing a matched relation trigger.

    Splits on sentence-final punctuation followed by whitespace and an
    uppercase letter, so abbreviations like 'spp.' or 'e.g.' don't cut
    the sentence short.
    """
    boundaries = [0] + [match.end() for match in re.finditer(r"(?<=[.!?])\s+(?=[A-Z])", text)] + [len(text)]
    left = max((position for position in boundaries if position <= start), default=0)
    right = min((position for position in boundaries if position >= end), default=len(text))
    return text[left:right].strip()


def evidence_level(relation: str, polarity: str, certainty: str, direct: bool) -> str:
    if polarity == "negative":
        return "negative_evidence"
    if not direct:
        return "background_or_related"
    if certainty == "suggested":
        return "hypothesis"
    if relation in {"virus_detected_in_host", "virus_associated_with_host"}:
        return "direct_observational"
    if relation in {"virus_has_gene_or_segment", "virus_coinfects_with"}:
        return "direct_descriptive"
    return "reported_context"


def extract_relations(passage: dict[str, str]) -> list[dict[str, object]]:
    text = passage["passage"]
    virus = passage["virus_name"]
    host = passage["host"]
    gene_hint = passage["gene_hint"]
    has_virus = entity_matches(passage, virus)
    has_host = entity_matches(passage, host)
    has_gene = entity_matches(passage, gene_hint)
    entities = entities_from_passage(passage.get("entity_ids", ""))
    known_virus = next((entity for entity in entities if entity["text"].casefold() == virus.casefold()), None)
    known_host = next((entity for entity in entities if entity["text"].casefold() == host.casefold()), None)
    results = []

    for relation, pattern in RELATION_RULES.items():
        for trigger in pattern.finditer(text):
            evidence_sentence = trigger_sentence(text, trigger.start(), trigger.end())
            sentence_lower = evidence_sentence.casefold()
            virus_in_sentence = virus.casefold() in sentence_lower
            host_in_sentence = host.casefold() in sentence_lower
            gene_in_sentence = gene_hint.casefold() in sentence_lower
            polarity = "negative" if NEGATION.search(evidence_sentence) else "positive"
            certainty = "suggested" if UNCERTAINTY.search(evidence_sentence) else "reported"
            direct = False
            object_name = ""
            object_type = ""
            object_id = ""
            manual_review = False

            if relation in {"virus_detected_in_host", "virus_associated_with_host"}:
                direct = virus_in_sentence and host_in_sentence
                object_name, object_type = host, "Species"
                object_id = known_host["id"] if known_host else ""
                manual_review = not direct
            elif relation == "virus_has_gene_or_segment":
                direct = virus_in_sentence and gene_in_sentence
                object_name, object_type = gene_hint, "GeneOrSegment"
                manual_review = not direct
            elif relation == "virus_coinfects_with":
                direct = virus_in_sentence
                related = related_virus_in_sentence(evidence_sentence, entities, virus)
                if related:
                    object_name, object_type, object_id = related["text"], "Species", related["id"]
                else:
                    object_name, object_type = "co-infecting virus", "Species"
                manual_review = not bool(related)
            elif relation == "virus_related_to":
                direct = virus_in_sentence
                related = related_virus_in_sentence(evidence_sentence, entities, virus)
                if related:
                    object_name, object_type, object_id = related["text"], "Species", related["id"]
                else:
                    object_name, object_type = "related virus", "Species"
                manual_review = not bool(related)
            elif relation == "variant_reported":
                direct = virus_in_sentence
                variants = VARIANT_PATTERN.findall(evidence_sentence)
                object_name, object_type = (variants[0], "Variant") if variants else ("unspecified variant", "Variant")
                manual_review = not bool(variants)
            else:
                direct = virus_in_sentence
                object_name, object_type = relation.removeprefix("virus_"), "Concept"
                manual_review = True

            results.append({
                "record_id": passage["record_id"],
                "pmid": passage["pmid"],
                "section": passage["section"],
                "subject": virus,
                "subject_type": "Species",
                "subject_id": known_virus["id"] if known_virus else "",
                "relation": relation,
                "object": object_name,
                "object_type": object_type,
                "object_id": object_id,
                "polarity": polarity,
                "certainty": certainty,
                "evidence_level": evidence_level(relation, polarity, certainty, direct),
                "direct_target_match": direct,
                "trigger": trigger.group(0),
                "evidence_sentence": evidence_sentence,
                "evidence_text": text,
                "pubtator_entities": entities,
                "manual_review": manual_review,
            })
    return results


def process_sample(input_tsv: Path, output_dir: Path) -> dict[str, object]:
    """Retrieve literature evidence for one classified sample TSV."""
    global INPUT, OUT
    INPUT, OUT = input_tsv, output_dir
    OUT.mkdir(parents=True, exist_ok=True)
    with INPUT.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    queries_out = []
    literature_out = []
    passages_out = []
    relations_out = []
    for row in rows:
        virus = clean(row.get("virus_name"))
        host = clean(row.get("host"))
        description = clean(row.get("description"))
        gene_hint = "nucleocapsid protein" if "nucleocapsid" in description.casefold() else "segment S"
        query_set = [
            ("virus", f'"{virus}"'),
            ("virus_host", f'"{virus}" AND "{host}"'),
            ("virus_host_gene", f'"{virus}" AND "{host}" AND "{gene_hint}"'),
        ]
        all_pmids: list[str] = []
        for query_type, query in query_set:
            source = "PubTator3"
            pmids = pubtator_pmids(query)
            if not pmids:
                source = "PubMed E-utilities fallback"
                pmids = pubmed_pmids(query)
            pmids = list(dict.fromkeys(pmids))
            all_pmids.extend(pmids)
            queries_out.append({
                "record_id": clean(row.get("id")), "virus_name": virus, "host": host,
                "gene_hint": gene_hint, "query_type": query_type, "query": query,
                "source": source, "pmid_count": len(pmids), "pmids": ";".join(pmids),
            })
            time.sleep(0.4)

        for pmid in list(dict.fromkeys(all_pmids))[:10]:
            literature_out.append({
                "record_id": clean(row.get("id")), "virus_name": virus, "host": host,
                "gene_hint": gene_hint, "pmid": pmid,
            })
            for result in extract_passages(pmid, virus, host, gene_hint):
                result.update({
                    "record_id": clean(row.get("id")), "virus_name": virus,
                    "host": host, "gene_hint": gene_hint,
                })
                passages_out.append(result)
                relations_out.extend(extract_relations(result))
            time.sleep(0.4)

    def save(name: str, rows_to_write: list[dict[str, object]]) -> None:
        keys = list(dict.fromkeys(key for row in rows_to_write for key in row))
        with (OUT / name).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows_to_write)

    save("retrieval_queries.tsv", queries_out)
    save("literature_pmids.tsv", literature_out)
    save("evidence_passages.tsv", passages_out)
    save("extracted_relations.tsv", [{**row, "pubtator_entities": json.dumps(row["pubtator_entities"], ensure_ascii=False)} for row in relations_out])
    (OUT / "extracted_relations.json").write_text(
        json.dumps(relations_out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    summary = {
        "retrieval_date": str(date.today()), "input": str(INPUT), "record_count": len(rows),
        "query_count": len(queries_out), "unique_pmids": len({r["pmid"] for r in literature_out}),
        "passage_count": len(passages_out), "relation_count": len(relations_out),
        "manual_review_count": sum(bool(row["manual_review"]) for row in relations_out),
        "output_dir": str(OUT),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    
    # 从配置中提取路径
    project_dir = Path(config["project_dir"])
    diamond_vmr_dir = Path(config["paths"]["diamond_vmr_dir"])
    results_dir = Path(config["paths"]["results_dir"])
    
    input_root = diamond_vmr_dir
    output_root = results_dir / "08_early_warning"

    sample_files = sorted(input_root.rglob("*.vmr_final_classified.tsv"))
    if not sample_files:
        raise FileNotFoundError(f"No *.vmr_final_classified.tsv files under {input_root}")
    summaries = [process_sample(path, sample_output_dir(path, input_root, output_root)) for path in sample_files]
    print(json.dumps({"sample_count": len(summaries), "samples": summaries}, ensure_ascii=False))


if __name__ == "__main__":
    main()
