# Troubleshooting — Phylogenetic Identification Agent

## The single most common failure

```
Alphabet 'e' is unknown
```

MAFFT is refusing the input because it was given nucleotides where it expected
amino acids. The `E` it does not recognise is glutamic acid.

**Cause.** `DB/tree/vmr_L_aln.fasta` contains nucleotide sequences. The tree
script concatenates that file directly with the query proteins, so a nucleotide
reference produces a mixed input that MAFFT cannot parse at all.

**Fix.** Rebuild the reference alignment with the database-update agent. That
agent's `build_tree_ref.py` translates the cleaned nucleotide library, takes the
longest protein per reference, and aligns the proteins. Verify before rebuilding:

```bash
grep -v '^>' DB/tree/vmr_L_aln.fasta | head -c 200
```

If the characters are dominated by `A C G T N`, the file is nucleotide and is
wrong. Amino-acid alignments contain `E`, `L`, `K`, `-`, and usually `X`.

## Failure branches

| Symptom | Cause | Fix |
|---|---|---|
| Exit 2, `no classified TSV found` | `06_diamond_vmr/*.vmr_final_classified.tsv` is absent | Run the sequence-identification agent |
| Exit 2, configuration error | `-c` was omitted; this script has no auto-discovery | Always pass `-c` |
| Exit 2, reference file missing | `DB/tree/vmr_L_aln.fasta` or `DB/diamond/vmr_l.dmnd` absent | Run the database-update agent |
| IQ-TREE runs for hours | `-m MFP` testing many models on a large alignment | Reduce `iqtree3.bootstrap` / `alrt`, keep the model search |
| Tree is a star phylogeny | Alignment too divergent, or the trim removed all informative sites | Inspect `tree_trim.fa`; consider `trimal.method: gappyout` |
| Support values all below the threshold | Genuinely ambiguous placement | Do not lower `sh_alrt_min` to force a species call |
| Every contig collapses into one clade | The reference set contains near-identical representatives, or duplicate sequences crept into `vmr_L_aln.fasta` | Deduplicate the reference set and rebuild |

## LucaProt-specific failures

### LucaProt cannot be located

```
lucaprot.code_dir is empty in the configuration and no LucaProt checkout was found automatically.
Clone it first:
    git clone https://github.com/alibaba/LucaProt.git <repo>/scripts/LucaProt
Looked for predict_many_samples.py under:
    <repo>/scripts/LucaProt/src
    <repo>/LucaProt/src
    <script dir>/../scripts/LucaProt/src
```

The three paths are printed so you can see exactly what was tried. Fixes, in
order of preference:

1. **Clone it where the script expects** — this keeps `config.yaml` untouched:
   ```bash
   cd <repo>
   git clone https://github.com/alibaba/LucaProt.git scripts/LucaProt
   ```
2. **Point at an existing checkout** — when LucaProt is already somewhere else:
   ```yaml
   lucaprot:
     code_dir: /absolute/path/to/LucaProt/src
   ```

A candidate only counts as found when `predict_many_samples.py` exists inside it.
A directory that exists but holds no such file is treated as not found, so a
half-finished clone reports the same message rather than failing later.

### Other LucaProt failures

| Symptom | Cause | Fix |
|---|---|---|
| Exit 2, `l.fasta` not found | The tree step has not run | Run `4-vmr-tree.py` first |
| Exit 2, output exists | `l.fasta.csv` is already present | Back it up, then pass `--overwrite` |
| `conda run -n lucaprot` fails | The environment does not exist or lacks the model dependencies | Check `conda env list`; the value must be the `lucaprot` env, not the pipeline env |
| Out of disk space during prediction | `emb/` is a cache that grows with every protein analysed, and reaches tens of gigabytes on large runs | Delete `scripts/LucaProt/emb` and rerun; it is recomputed |
| Prediction is very slow, GPU idle | `gpu_id` points at a device that is not visible, so the model silently falls back to CPU | Set `gpu_id: 0` on a single-GPU machine, `-1` to select CPU deliberately |
| CUDA error, no device | `gpu_id` points at a device that is not visible | Set `gpu_id: 0` |
| All probabilities exactly 0.5 | Input proteins are truncated to nothing, or the checkpoint path is wrong | Check `truncation_seq_length` against the protein lengths in `l.fasta` |

## Verifying the tree step without running it

`--dry-run` composes and logs the commands but does not execute them. Useful to
confirm the reference paths before committing to a multi-hour IQ-TREE run:

```bash
python 4-vmr-tree.py -c config.yaml --dry-run
```

Then read the log. If the paths in the composed MAFFT command point at a
nucleotide reference, stop and fix the database instead of running the tree.

## Before trusting a classification

1. Check `tree_true_result.treefile` in a tree viewer. The reference species
   should be where you expect them; a sample clade nested inside the wrong genus
   is a data problem, not a threshold problem.
2. Check the support distribution. If most nodes sit just under
   `sh_alrt_min`, the alignment — not the threshold — is the thing to fix.
3. Check that `l.fasta` headers still carry the sample prefix. Losing it breaks
   the join back to samples and every downstream per-sample statistic.
