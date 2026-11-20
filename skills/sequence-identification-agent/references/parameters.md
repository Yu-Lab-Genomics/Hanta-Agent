# Parameters — Sequence Identification Agent

All parameters live in the `vmr_merge_cluster:` section of `config.yaml`.
Every key is optional; a missing key falls back to the default shown below, and
an empty value means "do not pass the flag at all".

## Command line

| Flag | Type | Default | Purpose |
|---|---|---|---|
| `-s` / `--sample` | string | all samples found | Process one sample. Omitting it scans `06_diamond_vmr/` and processes every sample group. |
| `-c` / `--config` | path | auto-discovered | Configuration file. Discovery walks up four directory levels from the script and accepts `config.yaml`, `config-cell-paper.yaml` or `cell_submit-paper.yaml`. |
| `--manifest` | path | `DB/manifest.csv` | Reference manifest. Override only when testing against a trimmed manifest. |
| `--dry-run` | flag | off | Compose the commands and write them to the logs. Does not execute `seqkit` or `cd-hit-est`, and does not delete anything. |

## Configuration keys

| Key | Default | Maps to | Effect |
|---|---|---|---|
| `seqkit_grep_extra_args` | *(empty)* | appended to `seqkit grep` | Escape hatch; leave empty. |
| `cd_hit_est_c` | `0.99` | `-c` | Nucleotide identity threshold. **The single most important value.** Lower it and genuinely distinct variants merge. |
| `cd_hit_est_aS` | `0.8` | `-aS` | Fraction of the shorter sequence that must align. Raise it to reject partial overlaps. |
| `cd_hit_est_n` | `10` | `-n` | Word size. CD-HIT requires word size 5 for `-c` 0.7–1.0, 4 for 0.6–0.7, 3 for 0.5–0.6, 2 for 0.4–0.5. **Must be changed together with `-c`.** |
| `cd_hit_est_g` | `1` | `-g` | `1` = accurate mode, slowest, best cluster quality. `0` = fast mode. |
| `cd_hit_est_G` | `0` | `-G` | `0` = do not use global sequence identity. |
| `cd_hit_est_d` | `0` | `-d` | Length of the description written into the `.clstr` file. `0` keeps the full definition line. |
| `cd_hit_est_p` | `1` | `-p` | `1` = print alignment coverage in the `.clstr` file. |
| `cd_hit_est_M` | `0` | `-M` | Memory limit in MB. `0` = unlimited. |
| `cd_hit_est_T` | *(empty)* | `-T` | Threads. Empty = CD-HIT default (1). Set this on multi-core machines. |
| `cd_hit_est_extra_args` | *(empty)* | appended verbatim | Escape hatch; leave empty. |

## Paths used by this agent

These come from the `paths:` section and from the results layout, not from
`vmr_merge_cluster:`.

| Path | Source | Used for |
|---|---|---|
| `<results_dir>/06_diamond_vmr/` | `paths.results_dir` | Input hit tables and all six outputs |
| `<results_dir>/04_megahit/<SAMPLE>/final.contigs.fa` | results layout | Sequence source for `seqkit grep` |
| `<results_dir>/json/` | results layout | `<SAMPLE>.vmr_cluster.json` |
| `<results_dir>/log/<SAMPLE>/` | results layout | `vmr_seqkit.log`, `vmr_cdhit.log` |

## Tuning guidance

**Identity threshold.** 0.99 is deliberately strict: the three databases contain
overlapping reference genomes, so the same contig often hits several of them.
Dropping below 0.95 starts merging distinct virus variants, which changes the
tree topology and therefore every downstream risk score. Change `-c` only
together with `-n`, and re-run the whole downstream chain afterwards.

**Threads.** This step is CPU-bound on `cd-hit-est`. Leaving `cd_hit_est_T`
empty serialises it. On a shared cluster, set it to your allocation rather than
to the number of physical cores.

**Accuracy mode.** `-g 1` is the default for a reason: `-g 0` clusters by the
first sequence above threshold, so the representative depends on input order and
the result is not reproducible. Do not change it for published analyses.
