# Parameters — Risk Scoring and Reporting Agent

## 5.3-risk_score.py command line

| Flag | Type | Default | Purpose |
|---|---|---|---|
| `--config` | path | *(empty)* | Selects **batch mode**. Auto-locates the combined table, LucaProt output, QC directory and literature. This is the normal path. |
| `--input-tsv` | path | *(none)* | Selects **single-table mode**. Required when `--config` is absent. |
| `--outdir` | path | see below | Output directory. Single-table mode defaults to the input file's directory; batch mode defaults to `results_dir/00_results`. |
| `--relations` | path | auto-discovered | `extracted_relations.json` literature evidence. Batch mode scans `08_early_warning/*/`. |
| `--queries` | path | auto-discovered | `retrieval_queries.tsv` detection-bias covariates. Batch mode scans `08_early_warning/*/`. |
| `--labels` | path | `<script dir>/risk_labels.tsv` | Gold-standard labels, columns `species` and `label` (`0` / `1` / `borderline`). Absent → validation skipped. |
| `--lucaprot` | path | `results_dir/00_results/l.fasta.csv` | LucaProt probabilities. Absent → that weight renormalises away. |
| `--qc-dir` | path | `paths.fastp_dir`, else `results_dir/01_fastp` | fastp reports for the sequencing-depth covariates. |
| `--offline` | flag | off | Skip live PubTator3 queries. Human-infection evidence is then limited to VMR human isolates. |

Mode selection is exclusive: pass `--config`, or pass `--input-tsv`. Passing
neither exits 2.

## 5.4-report.py command line

| Flag | Type | Default | Purpose |
|---|---|---|---|
| `--config` | path | *(empty)* | Reads `results_dir/00_results/report_inputs.json` and writes `report.html` beside it |
| `--bundle` | path | *(none)* | Explicit path to `report_inputs.json`; the HTML lands in the same directory |
| `--output` | path | `<bundle dir>/report.html` | Override the HTML location |

`--config` takes precedence over `--bundle`. Passing neither exits 2.

## Detection confidence

Each signal produces a value in `[0, 1]`.

### DIAMOND signal

$$s_{\text{diamond}} = 0.6 \cdot s_{\text{evalue}} + 0.4 \cdot c_{\text{cover}}$$

| Term | Definition |
|---|---|
| $s_{\text{evalue}}$ | $\min(1, \max(0, -\log_{10} E / 100))$ |
| $c_{\text{cover}}$ | `length / (abs(send - sstart) + 1)`, capped at 1. Falls back to 0.5 when the coordinates are unreadable. |

### Tree signal

$$s_{\text{tree}} = 0.5 \cdot s_{\text{position}} + 0.3 \cdot s_{\text{support}} + 0.2 \cdot s_{\text{distance}}$$

| Term | Definition |
|---|---|
| $s_{\text{position}}$ | `nested` → 1.0, `sister` → 0.8, anything else → 0.5 |
| $s_{\text{support}}$ | SH-aLRT / 100 when present. Otherwise falls back on `classification_source`: `tree` 0.9, `blast` 0.6, `blast_weak` 0.4, `low_confidence` 0.3, default 0.5. |
| $s_{\text{distance}}$ | $\min(1, 0.05 / \max(d, 0.005))$ where $d$ is the nearest-neighbour distance. 0.5 when absent. |

Returns `None` when neither `tree_position` nor `tree_classification` is present,
which removes the signal and renormalises.

### LucaProt signal

The `prob` column, used directly. Absent → weight dropped.

### Combination

$$c = \frac{\sum_i w_i v_i}{\sum_i w_i}$$

Only signals that produced a value participate, so the denominator adapts.

## Final score and tiers

$$F = c \times R$$

| $F$ | Tier |
|---|---|
| Human-infection evidence present | **High**, and $F$ is floored at $\max(R, 60)$ |
| $\ge 60$ | High |
| $\ge 30$ | Medium |
| otherwise | Low |

$R$ is the species risk score, built from host breadth, human isolation and
human-infection literature.

**The floor is not a bug.** A contig from a species with documented human
infection scores at least 60 even if the detection confidence is low, because a
fragment of a known human pathogen is still a public-health signal.

## Validation

| Statistic | Uses | Why |
|---|---|---|
| Spearman correlation | The **full** `final_score` | The ranking as actually reported |
| AUROC | `species_risk − human_infection` at species level, and `detection_confidence × (species_risk − human_infection)` at sequence level | Removing the human-infection component avoids circularity: human infection is part of the label definition |

Do not change the AUROC feature to the full score. That inflates it to near 1.0
for the wrong reason.

## Configuration keys read

| Key | Used for |
|---|---|
| `paths.results_dir` | Locating every default path |
| `paths.fastp_dir` | QC reports, overriding the `01_fastp` default |
| `threads` | Not used by these two scripts |

There is no `risk_score:` section — the weights and tier boundaries are constants
in the script on purpose, so a scoring run is reproducible from the code alone.
Changing them is a code change and requires re-running validation.
