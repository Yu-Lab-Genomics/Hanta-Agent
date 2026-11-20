# Outputs and contract — Risk Scoring and Reporting Agent

## `results.tsv` — per-contig scores

`<results_dir>/00_results/results.tsv`. One row per contig that reached scoring.

### Identity and annotation

| Column | Meaning |
|---|---|
| `sample` | Sample the contig came from |
| `id` | Contig identifier |
| `species` | Species assigned by the phylogenetic agent |
| `virus_name` | Virus name from the manifest |
| `segment` | Segment, taken from `source_vmr` |
| `host_class` | Broad host class, e.g. `Rodent / Muridae` |
| `host_desc` | Host description in Latin only |
| `host_species` | Host species |
| `host_source` | Where the host annotation came from |
| `human_isolation` | Whether the reference was isolated from a human |

### Detection signals

| Column | Range | Meaning |
|---|---|---|
| `detection_diamond` | 0–1 | E-value and alignment-coverage composite |
| `detection_tree` | 0–1 or empty | Position, support and nearest-neighbour distance composite. Empty when the contig had no tree placement. |
| `detection_lucaprot` | 0–1 or empty | LucaProt RdRp probability. Empty when no prediction was available. |
| `detection_confidence` | 0–1 | Weighted mean of whichever signals were present |

An empty `detection_tree` or `detection_lucaprot` is informative: it means that
signal did not contribute and the confidence is an average over fewer terms. Two
rows with the same `detection_confidence` are not necessarily equally well
supported.

### Scores and tier

| Column | Range | Meaning |
|---|---|---|
| `species_risk_score` | 0–100 | How serious the species is, independent of this detection |
| `final_score` | 0–100 | `detection_confidence × species_risk_score`, except where the human-infection floor applies |
| `tier` | `High` / `Medium` / `Low` | Thresholds at 60 and 30 |
| `flags` | `;`-separated | Currently `human_evidence` when the floor was applied |

## `species_risk.tsv` — per-species attributes

`<results_dir>/00_results/species_risk.tsv`. One row per species.

| Column group | Examples | Meaning |
|---|---|---|
| Identity | `species`, `genus` | Taxonomy |
| Host | `host_class`, `host_desc`, `host_species` | **Latin only.** A mixed Chinese/Latin value must be cleaned. |
| Human | `human_isolation`, `human_infection_vmr`, `human_infection_literature`, `pmid`, `pmid_count` | The evidence that drives the floor |
| Score | species risk value | Feeds every contig of that species |

`host_desc` must not contain CJK characters. If it does, the source manifest was
not cleaned — fix it at the manifest, not here.

## `report_inputs.json` — the bundle

Written by `5.3-risk_score.py`, read by `5.4-report.py`.

```json
{
  "rows": [ … ],           // one entry per scored contig
  "rows_full": [ … ],      // rows before any column subsetting
  "species": { … },        // species_risk content, keyed by species
  "relations": [ … ],      // literature evidence, may be null
  "qc_dir": "…"            // fastp directory, may be null
}
```

| Key | Consumed for |
|---|---|
| `rows` | The main results table |
| `rows_full` | Any figure that needs columns dropped from `rows` |
| `species` | Species-level panels |
| `relations` | Literature-evidence panels |
| `qc_dir` | Sequencing-depth covariates |

`report_inputs.json` is the contract between the two scripts. Do not regenerate
the HTML from `results.tsv` — the bundle carries data the TSV does not.

## `report.html`

A single self-contained file. No external CSS, JS or images, so it can be attached
to a review or a supplementary submission without a working directory.

Requirements the HTML must satisfy:

1. **English only.** Host taxonomy appears as `Muridae`, never as `鼠科 Muridae`.
2. **Neutral wording.** Differences are described as "both", "reference only",
   "this study only". No "missed", "missing", "added", "wrong".
3. **No invented statistics.** Every number comes from the bundle.

## Optional input the report reads

`species_pmids.json`, in the same directory as the bundle. **No script in this
workflow produces it.** It is an external artefact holding per-species PubMed
identifiers. When absent, the report prints a warning and omits the PubMed
section — that is expected, not an error.

## Data flow

```
00_results/all_vmr_final_classified.tsv
00_results/l.fasta.csv                   (optional)
08_early_warning/*/extracted_relations.json   (optional)
                ↓  5.3-risk_score.py
00_results/results.tsv
00_results/species_risk.tsv
00_results/report_inputs.json
                ↓  5.4-report.py
00_results/report.html
```

`results.tsv` and `species_risk.tsv` are the machine-readable deliverables.
`report.html` is the human-readable one. All three come from the same bundle, so
a number quoted in the report can always be traced back to a row in the TSV.
