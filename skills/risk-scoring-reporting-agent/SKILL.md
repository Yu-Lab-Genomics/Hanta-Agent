---
name: risk-scoring-reporting-agent
description: Combine detection confidence with species-level risk into a tiered per-contig score, validate the ranking against a gold-standard label set, and render a self-contained HTML report. Agent 5 of the six-agent Hanta-Agent workflow, comprising two scripts that run in order.
---

# Risk Scoring and Reporting Agent

Convert detections into a ranked, tiered table and a readable report.

This agent has two scripts and they run in order:

1. **`5.3-risk_score.py`** — combine three independent detection signals into a
   confidence value, multiply by the species risk score, assign a tier, and
   validate the ranking with Spearman correlation and AUROC against a label set.
2. **`5.4-report.py`** — render `report_inputs.json` from the previous step into
   a standalone `report.html`.

## When to use

- The phylogenetic agent has produced `00_results/all_vmr_final_classified.tsv`.
- A ranked list of contigs with an explicit high / medium / low tier is needed.
- The scores must be validated against known labels before drawing conclusions.

## When not to use

- The tree step has not run. Tree-derived columns are one of the three signals;
  scoring without them silently degrades confidence to a two-signal average.
- The reference manifest or label set changed but the upstream chain did not.
  Re-run from the sequence-identification agent instead of patching scores here.
- The goal is to filter contigs. **This agent assigns a tier; it does not decide
  what is "real".** See the invariant below.

## Contract

| Item | Value |
|---|---|
| **Input** | `00_results/all_vmr_final_classified.tsv` (batch) or any classified TSV (single-table mode) |
| **Input** | `00_results/l.fasta.csv` — LucaProt probabilities; optional, the weight renormalises when absent |
| **Input** | fastp QC reports, for the sequencing-depth covariates |
| **Input** | `08_early_warning/*/extracted_relations.json` and `retrieval_queries.tsv`; optional literature evidence |
| **Input** | `risk_labels.tsv` beside the script — two columns `species`, `label` (`0` / `1` / `borderline`) |
| **Config** | `--config` selects batch mode; without it, `--input-tsv` selects single-table mode |
| **Output** | `00_results/results.tsv` — per-contig scores and tiers |
| **Output** | `00_results/species_risk.tsv` — per-species risk attributes |
| **Output** | `00_results/report_inputs.json` — the bundle the report script consumes |
| **Output** | `00_results/report.html` |
| **Failure behaviour** | Exit 2 for missing input, missing configuration or an unreadable bundle; a missing optional input is recorded and the remaining weights renormalise |

## Workflow

1. **Pick the mode.** Batch mode (`--config`) auto-locates the combined table,
   LucaProt output, QC directory and literature. Single-table mode
   (`--input-tsv`) is for one sample. Batch mode is the normal path.
2. **Confirm the three signals are present.** `detection_diamond`,
   `detection_tree`, `detection_lucaprot`. A missing signal is not an error; its
   weight is dropped and the rest renormalise. Check the row to see which.
3. **Run scoring.** `python 5.3-risk_score.py --config config.yaml`.
4. **Read the validation block.** Spearman correlation and AUROC against
   `risk_labels.tsv`. A high score that does not separate positives from
   negatives is a broken feature, not a success.
5. **Render.** `python 5.4-report.py --config config.yaml`.

## Invariants

1. **No detection threshold is applied anywhere.** A contig that matched the
   reference database is a related sequence, full stop. Low percent identity
   means the contig covers a fragment, not that the detection is doubtful. Do not
   introduce a `pident` cutoff, and do not label anything "uncertain" on identity
   grounds alone.
2. **`final_score = detection_confidence × species_risk_score`.** The two factors
   are independent: how sure the detection is, times how serious the species is.
   Changing either changes every tier boundary.
3. **Human-infection evidence overrides the tier.** When the species has human
   isolation or human-infection literature, the score is floored at 60 and the
   tier is `High` regardless of detection confidence. This is deliberate: a
   fragment of a known human pathogen still matters.
4. **Tier boundaries are fixed at 60 and 30.** `≥60` High, `≥30` Medium,
   otherwise Low. Do not re-tune them per dataset; the labels were validated
   against these boundaries.
5. **Validation uses the full score for correlation but the risk component alone
   for AUROC.** The AUROC computation drops the human-infection component, because
   keeping it would make the label circular — human infection is part of the
   label. Do not "fix" this.
6. **The report is English-only.** Host taxonomy must appear in Latin
   (`Muridae`), never as a mixed Chinese/Latin string.
7. **Wording is neutral.** Describe differences as "both", "reference only" or
   "this study only". Do not use "missed", "missing", "added" or "wrong".

## Decision points

| Situation | Adjustment | Risk |
|---|---|---|
| LucaProt output absent | Nothing; weights renormalise | Confidence is then a two-signal average, not comparable with runs that had three |
| Literature evidence absent | Nothing; those components score 0 | Genuine risk is understated for species with sparse literature |
| AUROC below 0.6 | Investigate the features, do not move the label set | Adjusting labels to improve AUROC is circular |
| Labels file absent | Validation is skipped | The scores are still produced but unvalidated; say so explicitly |
| Offline environment | Pass `--offline` | Human-infection evidence is then limited to VMR human isolates |

Full parameter reference: `references/parameters.md`.

## Failure handling

| Exit code | Meaning | First action |
|---|---|---|
| 2 | Missing or unreadable input, missing bundle | Read the stderr message; it names the path |
| other | Unexpected error inside scoring | Read the traceback; column-name drift in the upstream table is the usual cause |

Per-case diagnosis: `references/troubleshooting.md`.

## References

- `references/parameters.md` — command line, weights, tier boundaries, configuration.
- `references/outputs-and-contract.md` — every output column and the JSON bundle.
- `references/troubleshooting.md` — failure branches and validation checks.

## Provenance

`scripts/5.3-risk_score.py` and `scripts/5.4-report.py` are copies of
`../../5.3-risk_score.py` and `../../5.4-report.py`. Agent 5 corresponds to step 5
of the pipeline; steps 5.1 and 5.2 belong to the upstream sequence work and are
not part of this agent.

### Bundled helper modules

`scripts/seqid.py` (110 lines) is copied into this directory because
`5.3-risk_score.py` imports it. It exposes FASTA identifier helpers only and
reads no external data. The script resolves its own directory on `sys.path`, so
the module must sit beside it. The same module is bundled with three other
agents; the duplication is intentional so each skill works standalone.

`5.4-report.py` additionally imports `yaml`, which is a third-party package
(PyYAML) and must be installed in the running environment.
