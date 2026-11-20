# Troubleshooting — Risk Scoring and Reporting Agent

## Failure branches

| Symptom | Cause | Fix |
|---|---|---|
| Exit 2, no input given | Neither `--config` nor `--input-tsv` was passed | Pass one; they are mutually exclusive |
| Exit 2, combined table not found | `00_results/all_vmr_final_classified.tsv` is absent | Run the phylogenetic agent |
| Exit 2, `Requires --bundle (or --config)` | `5.4-report.py` has nothing to read | Run `5.3-risk_score.py` first |
| Exit 2, `Step-5.3 intermediate file not found` | `report_inputs.json` is missing | Same |
| `KeyError` on a column name | The upstream table's header changed | Compare `header = next(csv.reader(open(tsv)))` against the column list in `outputs-and-contract.md` |
| Validation silently skipped | `risk_labels.tsv` is not beside the script | Pass `--labels`, or accept unvalidated scores and say so |

## Diagnostic checks

### Every `detection_tree` is empty

The tree step did not contribute. Causes, in order of likelihood:

1. The tree agent has not run — run it before scoring.
2. The combined table was built without the tree columns.
3. `tree_position` and `tree_classification` are both absent from the header.

An empty tree signal is not an error, but it means the confidence is an average
over two signals instead of three. Results from such a run are not directly
comparable with a run that had all three.

### Every `detection_lucaprot` is empty

`00_results/l.fasta.csv` was not found or not passed. Check:

```bash
ls -la <results_dir>/00_results/l.fasta.csv
head -2 <results_dir>/00_results/l.fasta.csv
```

The header must contain `protein_id` and `prob`. A CSV from a different tool will
be rejected or silently ignored.

### All scores are identical

`species_risk_score` is the same for every contig of a species — that is by
design. If every contig across every species also has the same `final_score`,
then `detection_confidence` is constant, which usually means the three signal
columns are all empty and the confidence fell back to a default.

### Tiers are all High

Either every species carries human-infection evidence, or the floor is being
applied too broadly. Check the `flags` column: rows with `human_evidence` were
floored. If a species that should not be flagged has the flag, the
`human_infection_vmr` or `human_infection_literature` value in
`species_risk.tsv` is wrong.

### AUROC is suspiciously high

Check which feature was used. The correct AUROC feature drops the
human-infection component:

- Species level: `species_risk − human_infection`
- Sequence level: `detection_confidence × (species_risk − human_infection)`

Using the full score makes the AUROC circular and pushes it toward 1.0. That is
not a better model, it is a leak.

## The neutral-wording rule

When describing a difference between two result sets, the permitted phrasing is:

| Permitted | Not permitted |
|---|---|
| both identified | neither missed it |
| reference only / original only | missed / omitted |
| this study only | new / added |
| differs | wrong / erroneous |

This applies to the HTML report, to figure captions and to anything quoted from
them. The reason is factual: a contig absent from one set was not necessarily
missed — the two sets may simply have used different inputs.

## Chinese host names in the HTML

If `Muridae` appears as `鼠科 Muridae`, the source manifest still carries the
bilingual string. The report must show the Latin name only. Fix it upstream in
the manifest and regenerate, rather than patching the HTML — the same string
feeds `species_risk.tsv` and every figure.

To check an existing output:

```bash
grep -c '[一-龥]' <results_dir>/00_results/species_risk.tsv
```

A non-zero count means CJK survived somewhere in the table.

## Safe re-runs

Both scripts overwrite their outputs. Before re-running after a parameter change,
copy the previous `results.tsv` and `species_risk.tsv` aside, because the tier
boundaries and weights are constants in the code — a re-run after an edit is not
comparable with the previous output.
