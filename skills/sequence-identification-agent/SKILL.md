---
name: sequence-identification-agent
description: Merge the three per-database DIAMOND hit tables of every sample into one record per contig, extract the matched contig sequences, dereplicate them with CD-HIT-EST, and join the reference manifest so that host and taxonomy columns survive into the final classification table. Agent 3 of the six-agent Hanta-Agent workflow.
---

# Sequence Identification Agent

Turn raw DIAMOND hits into a non-redundant, manifest-annotated contig set.

The previous agent writes one hit table per reference database (`vmr_l`, `vmr_m`,
`vmr_s`). A single contig can hit more than one database, so the three tables
overlap. This agent resolves that overlap, drops near-identical contigs, and
attaches the reference metadata needed by every later step.

## When to use

- The sequence-handle agent has produced `06_diamond_vmr/<SAMPLE>.vmr_{l,m,s}.tsv`.
- Contigs must be dereplicated before tree building.
- The reference manifest changed and the annotated tables must be regenerated.

## When not to use

- The DIAMOND hit tables do not exist yet — run the sequence-handle agent first.
- There are no hits for a sample. The agent still records a result but produces
  no clustered FASTA; do not interpret that as a failure.
- The task is tree building, risk scoring or reporting.

## Contract

| Item | Value |
|---|---|
| **Input** | `06_diamond_vmr/<SAMPLE>.vmr_{l,m,s}.tsv` — discovered by scanning the directory |
| **Input** | `DB/manifest.csv` — reference manifest, second column `acc_base` |
| **Input** | `<results_dir>/04_megahit/<SAMPLE>/final.contigs.fa` — the sequences the hits refer to |
| **Config** | `-c <path>`; the default is resolved by walking up from the script and lands on `skills/config.yaml` |
| **Output** | Six files per sample in `06_diamond_vmr/`, plus `json/<SAMPLE>.vmr_cluster.json` and per-tool logs |
| **Downstream interface** | `outputs.classified_tsv` — the manifest-annotated table the tree step consumes |
| **Failure behaviour** | Exit code 2 for missing configuration, missing hit tables or missing contigs; tool failures propagate |

## Workflow

1. **Confirm the upstream state.** Each sample needs all three hit tables and a
   contig FASTA. A missing table is reported as `missing VMR tables` and stops that
   sample.
2. **Confirm the manifest.** `DB/manifest.csv` must have `acc_base` as its second
   column; the reader raises if it does not. Species and host columns are taken
   from the remaining columns.
3. **Preview.** `--dry-run` writes the commands into the logs and skips the tools.
   It does not delete anything.
4. **Run.** One sample with `-s`, or every sample found by omitting it.
5. **Read the JSON.** `status: success` with the message `no VMR hit identifiers` means
   the sample has no hits at all — legitimate, but downstream will find nothing
   to score.

## Invariants

1. **Every contig appears once in the output**, under the database given by
   `source_vmr`. The merged table is the single source of truth for segment
   assignment; do not re-derive the segment from the per-database columns.
2. **Dereplication is by nucleotide identity, not by name.** CD-HIT-EST at
   `-c 0.99 -aS 0.8` keeps one representative per cluster; the representative
   keeps the cluster's span of hits.
3. **The manifest join must not drop rows silently.** A contig whose accession is
   absent from the manifest keeps the classification but loses the extra columns;
   check the `manifest matched x/y` message in the JSON.
4. **`--dry-run` writes a `status: success` JSON.** Delete it before the real run.
5. **Downstream reads `outputs.classified_tsv`, not the individual intermediate
   files.** Treat the other five as diagnostic.

## Decision points

| Situation | Adjustment | Risk |
|---|---|---|
| Near-identical contigs survive as separate records | Tighten `vmr_merge_cluster.cd_hit_est_c` above 0.99 | Fewer records, but genuinely distinct variants can be merged |
| A legitimate variant is being collapsed | Raise `cd_hit_est_aS` coverage requirement, or lower `cd_hit_est_c` | More records, more near-duplicates in the tree |
| Many contigs are missing from the manifest join | Refresh `DB/manifest.csv` from the current VMR release | A stale manifest silently drops host information |
| Clustering is too slow | Lower `cd_hit_est_n` from 10 (it must match `-c`) | Word size and identity threshold must stay compatible |

Full parameter reference: `references/parameters.md`.

## Failure handling

| Exit code | Meaning | First action |
|---|---|---|
| 2 | Missing config, missing hit tables, missing contigs, or manifest format error | Read `message` in `json/<SAMPLE>.vmr_cluster.json` |
| other | `seqkit` or `cd-hit-est` failed | Read `log/<SAMPLE>/vmr_seqkit.log` or `vmr_cdhit.log` |

Per-case diagnosis: `references/troubleshooting.md`.

## References

- `references/parameters.md` — configuration keys and every tool argument.
- `references/outputs-and-contract.md` — the six output files and the JSON schema.
- `references/troubleshooting.md` — failure branches and what to check.

## Provenance

`scripts/3-vmr-merge-cluster.py` is a copy of `../../3-vmr-merge-cluster.py`.
Both copies share the 2026-11-20 fix: the default configuration is resolved by
walking up from the script instead of being hard-coded to
`skills_V2/config.yaml`.

### Bundled helper modules

`scripts/species_norm.py` (84 lines) is copied into this directory because
`3-vmr-merge-cluster.py` imports it. It exposes a single function,
`normalize_species()`, and reads no external data. The script resolves its own
directory on `sys.path`, so the module must sit beside it — do not move it to a
shared location without also changing that import path.

The same module is also bundled with `phylogenetic-identification-agent`, which
needs it as well. The duplication is intentional: each skill must work when
copied out on its own.
