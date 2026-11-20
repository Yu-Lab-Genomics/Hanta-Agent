---
name: phylogenetic-identification-agent
description: Build the shared RdRp reference tree from all samples, place each contig on it, assign species-level classification using SH-aLRT support, and run LucaProt to obtain an independent deep-learning RdRp probability for the combined L segment proteins. Agent 4 of the six-agent Hanta-Agent workflow, comprising two scripts that must run in order.
---

# Phylogenetic Identification Agent

Decide which virus a contig belongs to, using a tree rather than a BLAST hit.

This agent has two scripts and they must run in the order below:

1. **`4-vmr-tree.py`** — align all representative contigs against the reference
   alignment, build one shared tree across every sample, cut it into
   species-level groups at the bootstrap support threshold, and write the
   combined L-segment protein FASTA.
2. **`4.2-lucaprot.sh`** — run LucaProt on that FASTA to get an independent
   RdRp probability per protein. It consumes `00_results/l.fasta`, so it cannot
   run first.

## When to use

- The sequence-identification agent has produced `06_diamond_vmr/<S>.vmr_final_classified.tsv`.
- Contigs need species assignment from phylogenetic placement rather than from
  the best hit.
- The tree must be rebuilt after adding samples or changing the reference set.

## When not to use

- No clustered contigs exist yet — run the sequence-identification agent first.
- The reference alignment or DIAMOND index is missing. This agent does not build
  them; use the database-update agent.
- Only a LucaProt probability is wanted. Then run `4.2-lucaprot.sh` alone, but
  `00_results/l.fasta` must already exist.

## Contract

| Item | Value |
|---|---|
| **Input** | `06_diamond_vmr/<SAMPLE>.vmr_final_classified.tsv` (all samples, scanned) |
| **Input** | `DB/tree/vmr_L_aln.fasta` — **protein** alignment of one representative per reference species |
| **Input** | `DB/diamond/vmr_l.dmnd` — DIAMOND index for the L segment |
| **Input** | `scripts/LucaProt/src` — external LucaProt checkout, auto-discovered; optional |
| **Config** | `-c <path>` is **required**; there is no auto-discovery |
| **Output** | `07_tree/_all/{tree_input.fa, tree_trim.fa, tree_true_result.treefile}` — the shared tree |
| **Output** | `07_tree/<SAMPLE>/` — per-sample intermediates |
| **Output** | `00_results/all_vmr_final_classified.tsv` — union of all samples with tree columns |
| **Output** | `00_results/l.fasta` — combined L-segment proteins, sample-prefixed |
| **Output** | `json/<SAMPLE>.tree.json`; `json/lucaprot.json` from the second script |
| **Downstream interface** | `00_results/all_vmr_final_classified.tsv` for scoring; `00_results/l.fasta.csv` for the LucaProt probability |
| **Failure behaviour** | Exit 2 for missing configuration, missing classified tables or missing reference files; MAFFT / trimAl / IQ-TREE failures propagate |

## Workflow

1. **Verify the reference alignment is protein, not nucleotide.** `vmr_L_aln.fasta`
   is concatenated directly with the query proteins before MAFFT runs. A
   nucleotide file there makes MAFFT fail with `Alphabet 'e' is unknown` — the
   nucleotide alphabet has no amino acid `E`.
2. **Preview.** `--dry-run` composes the commands and writes them to the logs
   without running the tools.
3. **Build the tree.** Run `4-vmr-tree.py -c <config>`. One shared tree is built
   across all samples by design: per-sample trees are not comparable, and the
   classification step needs a single topology.
4. **Read the support threshold.** A clade is accepted as a species group when
   its SH-aLRT support reaches `phylogeny.sh_alrt_min` (default 80). Contigs that
   fall below it stay at genus level.
5. **Run LucaProt.** `4.2-lucaprot.sh -c <config>` reads `00_results/l.fasta`,
   calls the `lucaprot` conda environment, and writes `00_results/l.fasta.csv`.
   It refuses to overwrite an existing CSV unless `--overwrite` is given.
   LucaProt itself is an external project and is **not** bundled here. Clone it
   into the repository first and the step configures itself:
   ```bash
   git clone https://github.com/alibaba/LucaProt.git scripts/LucaProt
   ```
   With `lucaprot.code_dir` left empty, the script searches
   `<repo>/scripts/LucaProt/src`, `<repo>/LucaProt/src` and
   `<script dir>/../scripts/LucaProt/src`, in that order, and only accepts a
   candidate that actually contains `predict_many_samples.py`. It exits 2 and
   lists every path it tried when none matches.

## External dependency: LucaProt

LucaProt is about 870 MB to clone (1048 files, including 95 MB of model weights)
and creates an `emb/` cache beside it that grows with every protein analysed.
Both are covered by `.gitignore`, so neither can accidentally enter this
repository.

| Situation | What to do |
|---|---|
| Standard setup | Clone into `<repo>/scripts/LucaProt`, leave `code_dir` empty |
| LucaProt lives elsewhere | Set `lucaprot.code_dir` to its `src` directory |
| No GPU available | Set `lucaprot.gpu_id: -1`; it runs on CPU, slowly |
| Running out of disk | Delete `scripts/LucaProt/emb`; it is a cache and is recomputed |
| Skipping LucaProt entirely | Do not run this script. The scoring agent drops that weight and renormalises the remaining two |

## Invariants

1. **One tree for all samples.** `07_tree/_all/` is the only tree that matters.
   Do not build per-sample trees and merge them; the branch lengths would not be
   comparable.
2. **The tree is built on proteins.** Both the reference alignment and the query
   sequences must be amino acids. This is the most common failure mode.
3. **Classification depends on SH-aLRT, not on bootstrap alone.** The threshold is
   `phylogeny.sh_alrt_min`; changing it changes how many contigs reach species
   level, which changes every downstream risk score.
4. **`l.fasta` carries a sample prefix in each header.** That prefix is how
   LucaProt results are mapped back to samples. Do not strip it.
5. **`4.2-lucaprot.sh` is not read-only by default.** Without `--overwrite` it
   exits 2 when the CSV exists; with `--overwrite` it replaces it. Back up first.

## Decision points

| Situation | Adjustment | Risk |
|---|---|---|
| Too few contigs reach species level | Lower `phylogeny.sh_alrt_min` below 80 | Accepts weakly supported clades; unsupported species calls enter the risk score |
| Species groups look over-split | Raise `sh_alrt_min`, or re-check the reference alignment for mislabelled representatives | Fewer, broader groups |
| Tree building exceeds the time budget | Lower `iqtree3.bootstrap` / `alrt`, never `-m MFP` | Reduced support resolution |
| MAFFT is slow on a large set | Keep `mafft.preset: auto`; do not force `--genafpair` | `auto` already selects per input size |
| LucaProt needed without a tree | Run `4.2-lucaprot.sh` alone against an existing `l.fasta` | It will still fail if `l.fasta` is absent |

Full parameter reference: `references/parameters.md`.

## Failure handling

| Exit code | Meaning | First action |
|---|---|---|
| 2 | Missing config, missing classified tables, or missing reference files | Read the stderr message; it names the exact path |
| other | MAFFT, trimAl or IQ-TREE failed | Read `log/` for the tool output |
| `Alphabet 'e' is unknown` | The reference alignment is nucleotide | Rebuild `DB/tree/vmr_L_aln.fasta` with the database-update agent |
| `4.2-lucaprot.sh` exit 2 | `l.fasta` absent, LucaProt not found, or CSV exists without `--overwrite` | Run the tree step first; clone LucaProt; or pass `--overwrite` after backing up |

Per-case diagnosis: `references/troubleshooting.md`.

## References

- `references/parameters.md` — configuration keys and every tool argument.
- `references/outputs-and-contract.md` — output trees, the classification columns, and the JSON schemas.
- `references/troubleshooting.md` — failure branches and what to check.

## Provenance

`scripts/4-vmr-tree.py` is a copy of `../../4-vmr-tree.py`, unmodified.
`scripts/4.2-lucaprot.sh` is a copy of `../../4.2-lucaprot.sh`, rewritten on
2026-11-20: shebang and error handling added, all three hard-coded absolute paths
removed, the GPU id made configurable, and JSON output introduced. Also on
2026-11-20, `lucaprot.code_dir` and `lucaprot.python` became optional: an empty
`code_dir` triggers auto-discovery of `<repo>/scripts/LucaProt/src`, so a fresh
clone needs no configuration.

### Bundled helper modules

| Module | Lines | Purpose |
|---|---|---|
| `scripts/seqid.py` | 110 | FASTA identifier handling: `read_fasta_ids()`, `strip_prodigal_suffix()`, `prodigal_orig_id()` |
| `scripts/species_norm.py` | 84 | `normalize_species()` — canonical species-name form |

Both are pure functions with no external data files, and both are copied into
this directory because the scripts put their own directory on `sys.path`. The
same modules are bundled with `database-update-agent`,
`risk-scoring-reporting-agent` and `sequence-identification-agent`. The
duplication is intentional: each skill must work when copied out on its own.
