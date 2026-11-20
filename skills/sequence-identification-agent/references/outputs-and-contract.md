# Outputs and contract — Sequence Identification Agent

## Files produced per sample

All six live in `<results_dir>/06_diamond_vmr/`.

| File | Content | Consumed by |
|---|---|---|
| `<SAMPLE>.vmr_merged.tsv` | One row per contig across all three hit tables, with `source_vmr` recording which database won | Diagnostic |
| `<SAMPLE>.vmr_ids.txt` | The matched contig identifiers, one per line | Diagnostic |
| `<SAMPLE>.vmr_hits.fa` | The contig sequences for those identifiers | Diagnostic |
| `<SAMPLE>.vmr_clustered.fa` | CD-HIT-EST representatives | Diagnostic |
| `<SAMPLE>.vmr_clustered.fa.clstr` | Cluster membership, with alignment coverage because `-p 1` | Diagnostic |
| `<SAMPLE>.vmr_final.tsv` | Representative contigs with their hit columns | Diagnostic |
| **`<SAMPLE>.vmr_final_classified.tsv`** | **Representative contigs joined to the reference manifest** | **Agent 4 (tree) and Agent 5 (scoring)** |

The JSON in `outputs` names every one of these; downstream code must read
`outputs.classified_tsv` rather than assuming the filename.

## Result JSON

Written to `<results_dir>/json/<SAMPLE>.vmr_cluster.json`.

```json
{
  "status": "success",
  "message": "…",
  "outputs": {
    "merged_tsv": "…/<SAMPLE>.vmr_merged.tsv",
    "ids_txt": "…/<SAMPLE>.vmr_ids.txt",
    "hits_fasta": "…/<SAMPLE>.vmr_hits.fa",
    "clustered_fasta": "…/<SAMPLE>.vmr_clustered.fa",
    "final_tsv": "…/<SAMPLE>.vmr_final.tsv",
    "classified_tsv": "…/<SAMPLE>.vmr_final_classified.tsv"
  }
}
```

| Field | Meaning |
|---|---|
| `status` | `success` when the agent finished, including the "no hits" case. `failed` on a handled error. |
| `outputs.*` | **Absolute paths.** Never reconstruct a filename from the sample name; read these keys. |

## The classified table

`vmr_final_classified.tsv` is the contract with the rest of the workflow. One
row per CD-HIT-EST representative.

| Column group | Examples | Origin |
|---|---|---|
| Contig identity | contig id, length | MEGAHIT |
| Hit metrics | `pident`, `length`, `evalue`, `bitscore` | DIAMOND |
| Database provenance | `vmr_s`, `vmr_m`, `vmr_l`, `source_vmr` | This agent |
| Reference identity | `acc_base`, accession, definition | Manifest join |
| Reference taxonomy | species, genus, host columns | Manifest join |

### Reading `vmr_s` / `vmr_m` / `vmr_l` correctly

These three columns record **which reference databases the contig hit**, not
which segment it ends up assigned to. A contig can hit all three simultaneously —
in the production dataset 119 contigs do exactly that. The assigned segment is
`source_vmr` and nothing else.

Do not compute segment proportions from the three hit columns; that inflates the
counts. Use `source_vmr`.

## Upstream and downstream

```
Agent 2 → 06_diamond_vmr/<S>.vmr_{l,m,s}.tsv
              ↓
        Agent 3 (this one)
              ↓
        06_diamond_vmr/<S>.vmr_final_classified.tsv
              ↓
        Agent 4 → 00_results/all_vmr_final_classified.tsv
```

The tree agent concatenates the per-sample classified tables. If a sample is
missing from that union, check whether this agent wrote a JSON for it.
