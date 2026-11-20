# Configuration reference — Sequence Handle Agent

All keys come from the shared `../config.yaml`. This agent reads only the
`project_dir`, `conda_env`, `threads`, `paths`, `diamond` and `software`
sections; every other section belongs to a different agent.

---

## 1. How the configuration is located

The script resolves the default configuration by walking up from its own
directory and returning the first of these that exists:

```
<dir>/config.yaml
<dir>/config-cell-paper.yaml
<dir>/cell_submit-paper.yaml
```

walked over `<script dir>`, its parent, grandparent and great-grandparent.

Resulting behaviour:

| Where the script lives | Resolved configuration |
|---|---|
| `skills_new/2-species_class.py` | `skills_new/config-cell-paper.yaml` |
| `final/sequence-handle-agent/scripts/2-species_class.py` | `final/config.yaml` |

An explicit `-c /abs/path/config.yaml` always takes precedence. **Prefer it in
production**: when a directory holds both `config-cell-paper.yaml` and
`cell_submit-paper.yaml`, the search order decides which one wins, and that is
easy to get wrong silently.

The previous version of this script hard-coded
`script_dir.parent.parent / "skills_V2" / "config.yaml"`, which pointed at a
stale configuration and broke as soon as the script was moved.

---

## 2. Global keys

| Key | Shipped value | Meaning |
|---|---|---|
| `project_dir` | `/workspace/liangpeiyu/project/hantaKG_Agent` | Required. The database directory is resolved as `<project_dir>/DB/diamond` |
| `conda_env` | `lpy-hanta-v2` | Environment used for every `conda run -n` invocation |
| `threads` | `64` | Passed to DIAMOND as `--threads` — the full budget, not divided |

There is no per-sample parallelism here: samples are processed one after another
and each one gets all `threads` cores. This differs from the quality-control
agent, which splits the budget across concurrent samples.

---

## 3. `paths` — where things are read and written

| Key | Resolution | Used as |
|---|---|---|
| `results_dir` | `paths.results_dir`, else `<project_dir>/results` | Root of `06_diamond_vmr/` |
| `json_dir` | `paths.json_dir`, else `<results_dir>/json` | Scanned for upstream records; also receives the result JSON |
| `log_dir` | `paths.log_dir`, else `<results_dir>/log` | Per-sample log directory |

None of these are created by hand: the agent creates `06_diamond_vmr/` and the
per-sample log directory as needed.

---

## 4. `diamond` — alignment parameters

| Key | Shipped value | Passed as | Effect |
|---|---|---|---|
| `evalue` | `0.00001` | `--evalue` | Sequence-level significance threshold. Raising it admits weak hits that later inflate detections |
| `max_target_seqs` | `1` | `--max-target-seqs` | Best hit per database per query. This is a **per-database** limit, so one contig can still appear in both `vmr_l` and `vmr_m` |

The remaining DIAMOND arguments are fixed in the script:

```
--threads <threads>
-d DB/diamond/vmr_{l,m,s}.dmnd
-q <outputs.contigs_fasta from the upstream JSON>
-o 06_diamond_vmr/<SAMPLE>.vmr_{l,m,s}.tsv
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
--header
```

`blastx` is used, so the query is translated in all six frames. The query must
therefore be **nucleotide** contigs, not proteins.

---

## 5. `software` — executable overrides

The shipped configuration carries a `software` block used by later agents
(`seqkit`, `diamond`, `mafft`, `trimal`, `iqtree3`). **This agent ignores it** —
DIAMOND is invoked by name inside the conda environment. If DIAMOND ever needs an
absolute path, the change has to be made in the script, not in the YAML.

---

## 6. Upstream discovery rule

The agent does not take a sample list. It walks `json_dir/*.json`, parses each
file, and for every top-level key whose value contains a non-empty
`outputs.contigs_fasta`, registers that sample.

Consequences worth knowing:

- A record written by the quality-control agent qualifies; a record written by
  this agent does not, because its `outputs` block holds `vmr_l`/`vmr_m`/`vmr_s`
  instead of `contigs_fasta`. There is no risk of the agent consuming its own
  output.
- A record left behind by a `--dry-run` **will** be picked up. If a dry-run JSON
  claims `status: success` for a sample whose contigs were never assembled, the
  agent will register the sample and then fail on the missing file.
- An unparsable JSON raises `Failed to read JSON` and stops the whole run, rather than
  silently skipping the sample.

Passing `-s <SAMPLE>` restricts the run to that key, but discovery still happens
from the JSON directory — the sample must already have a qualifying record.

---

## 7. Command-line flags

| Flag | Meaning |
|---|---|
| `-s`, `--sample` | Process only this sample key |
| `-c`, `--config` | Configuration path; overrides the upward search |
| `--dry-run` | Write the command into the log and skip the DIAMOND call; the JSON is still written as `success` |
| `--overwrite` | Permit replacing existing hit tables |

`--overwrite` is not a background option: without it, any pre-existing
`06_diamond_vmr/<SAMPLE>.vmr_*.tsv` makes that sample fail with exit code 2 and
the message `output already exists, pass --overwrite to replace it`.

---

## 8. Sections this agent ignores

`input`, `cleanup`, `fastp`, `prinseq`, `bowtie2`, `megahit`, `prodigal`,
`mafft`, `trimal`, `iqtree3`, `phylogeny`, `vmr_merge_cluster`, `tree`, `report`
and any other section belong to other agents. Changing them here has no effect.
