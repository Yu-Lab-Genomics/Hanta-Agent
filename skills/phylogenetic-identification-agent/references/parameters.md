# Parameters — Phylogenetic Identification Agent

Two scripts, two configuration surfaces.

- `4-vmr-tree.py` reads `mafft:`, `trimal:`, `iqtree3:` and `phylogeny:`.
- `4.2-lucaprot.sh` reads `lucaprot:`.

## 4-vmr-tree.py command line

| Flag | Type | Default | Purpose |
|---|---|---|---|
| `-s` / `--sample` | string | all samples | Restrict to one sample. The shared tree in `07_tree/_all/` is still built from whatever samples are selected. |
| `-c` / `--config` | path | **required** | Configuration file. There is no auto-discovery; the script exits 2 without it. |
| `--manifest` | path | `DB/manifest.csv` | Reference manifest used to label tree leaves. |
| `--dry-run` | flag | off | Compose the commands and log them. Does not run MAFFT, trimAl or IQ-TREE. |

## Alignment configuration

| Key | Default | Maps to | Effect |
|---|---|---|---|
| `mafft.preset` | `auto` | `mafft --auto` | MAFFT picks the algorithm from the input size. Leave as is; forcing a preset on an unknown input size either wastes hours or degrades the alignment. |
| `phylogeny.mafft_preset` | `auto` | legacy duplicate | Superseded by `mafft.preset`; leave empty. |
| `trimal.method` | `automated1` | `trimal -automated1` | Heuristic that chooses gap and similarity thresholds automatically. The usual choice for phylogenetics. |
| `trimal.out_format` | `fasta` | `-fasta` | Output format. Must be a format IQ-TREE can read. |

## Tree configuration

| Key | Default | Maps to | Effect |
|---|---|---|---|
| `iqtree3.model` | `MFP` | `-m MFP` | ModelFinder Plus: tests substitution models and picks the best. **Do not replace with a fixed model** — a wrong model silently distorts branch lengths. |
| `iqtree3.bootstrap` | `1000` | `-B 1000` | Ultrafast bootstrap replicates. |
| `iqtree3.alrt` | `1000` | `-alrt 1000` | SH-aLRT replicates. These are the values the classification threshold reads. |
| `phylogeny.iqtree_bootstrap` | `1000` | legacy duplicate | Superseded by `iqtree3.bootstrap`. |
| `phylogeny.iqtree_alrt` | `1000` | legacy duplicate | Superseded by `iqtree3.alrt`. |

## Classification configuration

| Key | Default | Meaning |
|---|---|---|
| `phylogeny.sh_alrt_min` | `80` | Minimum SH-aLRT support for a clade to be accepted as a species group. **Changing this changes the number of species-level calls and therefore every risk score.** |

## Software paths

| Key | Default | Note |
|---|---|---|
| `software.mafft` | `mafft` | Bare name on `PATH`, or an absolute path |
| `software.trimal` | `trimal` | Same |
| `software.iqtree3` | `iqtree3` | Note the `3` — IQ-TREE 3, not 2 |
| `software.seqkit` | `seqkit` | Used for sequence extraction |
| `software.diamond` | `diamond` | Used for the reference check |
| `threads` | `64` | Passed to IQ-TREE and MAFFT |

## 4.2-lucaprot.sh command line

| Flag | Type | Default | Purpose |
|---|---|---|---|
| `-c` / `--config` | path | auto-discovered | Configuration file. Discovery walks up four directory levels from the script and accepts `config.yaml`, `config-cell-paper.yaml` or `cell_submit-paper.yaml`. |
| `--dry-run` | flag | off | Resolve the command and write the log, then stop. Nothing is deleted and no prediction runs. |
| `--overwrite` | flag | off | Replace an existing `l.fasta.csv`. Without it the script exits 2 rather than overwriting. |
| `-h` / `--help` | flag | — | Usage message |

The input FASTA, output directory and embedding directory are **not** command
line flags. They all derive from the configuration, so the script never depends
on the current working directory.

## lucaprot configuration

### How code_dir is resolved

`4.2-lucaprot.sh` needs the directory holding `predict_many_samples.py`. It
resolves it in this order:

1. `lucaprot.code_dir`, when set — used verbatim.
2. **Automatic discovery**, when `code_dir` is empty. The script tries, relative
to its own location:

   | Candidate | Layout it covers |
   |---|---|
   | `<repo>/scripts/LucaProt/src` | The documented layout: clone LucaProt into `<repo>/scripts/LucaProt` |
   | `<repo>/LucaProt/src` | LucaProt cloned at the repository root instead |
   | `<script dir>/../scripts/LucaProt/src` | Flat layout |

3. Nothing found — exit 2, printing every path it tried.

A candidate only counts as found when `predict_many_samples.py` exists inside it,
so an empty or half-cloned directory is not silently accepted.

**The recommended setup is the zero-configuration one:** clone LucaProt into
`<repo>/scripts/LucaProt`, leave `lucaprot.code_dir` empty, and the step works on
any machine without editing the configuration.

```bash
git clone https://github.com/alibaba/LucaProt.git <repo>/scripts/LucaProt
```

### Keys

| Key | Default | Meaning |
|---|---|---|
| `lucaprot.code_dir` | *(empty)* | Directory containing `predict_many_samples.py`. **Empty means auto-discover** — see above. Set it only when LucaProt lives outside the repository |
| `lucaprot.conda_env` | `lucaprot` | Conda environment holding the LucaProt dependencies. Invoked as `conda run -n <env>` |
| `lucaprot.python` | *(empty)* | Command run **inside** that environment. Empty falls back to `python3`. The full invocation is `conda run -n <conda_env> <python> predict_many_samples.py`, so an absolute interpreter path is only needed when the environment name is not enough |
| `lucaprot.gpu_id` | `2` | GPU index. Was hard-coded to `2` before the rewrite; set it to `0` when only one device is visible, `-1` for CPU |
| `lucaprot.threshold` | `0.5` | RdRp probability threshold for the binary call |
| `lucaprot.truncation_seq_length` | `4096` | Maximum protein length fed to the model; longer proteins are truncated |
| `lucaprot.dataset_name` | `rdrp_40_extend` | Checkpoint dataset name |
| `lucaprot.dataset_type` | `protein` | — |
| `lucaprot.task_type` | `binary_class` | — |
| `lucaprot.model_type` | `sefn` | Model architecture |
| `lucaprot.time_str` | `20230201140320` | Checkpoint timestamp |
| `lucaprot.step` | `100000` | Checkpoint global step |
| `lucaprot.print_per_number` | `10` | Progress printing interval |
| `lucaprot.emb_dir` | *(empty)* | Empty means `<code_dir>/../emb`, that is `scripts/LucaProt/emb`. This is a **cache**: it grows with every new protein and may be deleted to reclaim space; the next run recomputes it |
| `lucaprot.extra_args` | *(empty)* | Appended verbatim to the LucaProt call |

## Input and output paths derived by the script

| Path | Derivation |
|---|---|
| Input FASTA | `<paths.results_dir>/00_results/l.fasta` |
| Output CSV | `<paths.results_dir>/00_results/l.fasta.csv` |
| Result JSON | `<paths.results_dir>/json/lucaprot.json` |
| Log | `<paths.results_dir>/log/lucaprot.log` |

## Ordering rule

`4.2-lucaprot.sh` must run **after** `4-vmr-tree.py`, because `l.fasta` is an
output of the tree step. Running them the other way round exits 2 with a missing
input message. This ordering is a property of the data flow, not a convention.
