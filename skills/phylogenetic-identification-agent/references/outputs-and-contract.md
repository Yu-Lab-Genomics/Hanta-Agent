# Outputs and contract — Phylogenetic Identification Agent

## 1. Shared tree (`4-vmr-tree.py`)

`<results_dir>/07_tree/_all/`

| File | Content |
|---|---|
| `tree_input.fa` | Reference protein alignment concatenated with the combined query proteins. This is exactly what MAFFT receives. |
| `tree_trim.fa` | The trimmed alignment produced by trimAl |
| `tree_true_result.treefile` | The Newick tree from IQ-TREE, with SH-aLRT and ultrafast bootstrap values on the nodes |

`tree_true_result.treefile` is the only artefact the classification step reads.
The other two exist so the alignment can be inspected when a topology looks
wrong.

## 2. Per-sample intermediates

`<results_dir>/07_tree/<SAMPLE>/` holds the per-sample query alignment and tree
inputs. These are diagnostics. **Do not build conclusions from a per-sample
tree** — the branch lengths are not comparable with the shared tree and the
classification step does not use them.

## 3. Combined tables (`00_results/`)

| File | Content | Consumed by |
|---|---|---|
| **`all_vmr_final_classified.tsv`** | Union of every sample's classified table, with the tree columns appended | Agent 5 (risk scoring) |
| **`l.fasta`** | Combined L-segment proteins, each header prefixed with its sample name | `4.2-lucaprot.sh` |
| `l.fasta.csv` | LucaProt output, written by the second script | Agent 5, as the `--lucaprot` input |

## 4. Tree columns added to the classified table

The classification step appends columns to the per-contig rows. The scoring agent
reads them by name.

| Column | Meaning |
|---|---|
| `final_classification` | The assigned taxon, or the genus-level fallback when support is insufficient |
| `tree_*` | Placement and support values: the clade the contig joined and its SH-aLRT support |
| `support` / SH-aLRT value | Compared against `phylogeny.sh_alrt_min` to decide species versus genus level |

Exact column names are the ones present in `00_results/all_vmr_final_classified.tsv`;
read the header rather than assuming. The scoring agent already does this.

## 5. Result JSON

`<results_dir>/json/<SAMPLE>.tree.json`:

```json
{
  "status": "success",
  "message": "…",
  "outputs": {
    "tree_input": "…/07_tree/_all/tree_input.fa",
    "tree_trimmed": "…/07_tree/_all/tree_trim.fa",
    "treefile": "…/07_tree/_all/tree_true_result.treefile",
    "classified_tsv": "…/00_results/all_vmr_final_classified.tsv",
    "l_fasta": "…/00_results/l.fasta"
  }
}
```

`<results_dir>/json/lucaprot.json`, written by `4.2-lucaprot.sh`:

```json
{
  "status": "success",
  "input_fasta": "…/00_results/l.fasta",
  "output_csv": "…/00_results/l.fasta.csv",
  "sequence_count": 0,
  "conda_env": "lucaprot",
  "elapsed_seconds": 0
}
```

| Field | Meaning |
|---|---|
| `status` | `success`, `skipped` (input absent) or `failed` |
| `output_csv` | Where the probability table landed. Agent 5 reads this. |
| `sequence_count` | How many proteins were scored. Zero means the FASTA was empty. |

## 6. Reading the LucaProt CSV

`l.fasta.csv` has at least two columns: `protein_id` and `prob`.

| Column | Meaning | Caution |
|---|---|---|
| `protein_id` | The FASTA header, **including the sample prefix** | The prefix is the only link back to a sample; do not strip it when joining |
| `prob` | RdRp probability in `[0, 1]` | This is an independent signal, not a verdict. Agent 5 combines it as one weighted component. |

A `prob` near 0.5 with `threshold: 0.5` is a coin flip. Do not report it as
"not an RdRp".

## Data flow

```
06_diamond_vmr/<S>.vmr_final_classified.tsv   (per sample)
              ↓  4-vmr-tree.py
07_tree/_all/tree_true_result.treefile        (one tree, all samples)
00_results/all_vmr_final_classified.tsv       (with tree columns)
00_results/l.fasta                            (combined L proteins)
              ↓  4.2-lucaprot.sh
00_results/l.fasta.csv                        (RdRp probability)
              ↓
Agent 5 (Risk Scoring and Reporting)
```

If `all_vmr_final_classified.tsv` is missing, the tree step did not finish.
If `l.fasta.csv` is missing but `l.fasta` exists, only the second script failed —
the tree and classification are still valid.
