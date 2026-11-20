# Configuration reference — Quality Control Agent

All keys live in `../config.yaml` (a copy of `skills_new/config-cell-paper.yaml`).
The file is shared by all six agents of the Hanta-Agent workflow; this agent reads
only the sections listed here and ignores the rest.

The script is silent: it resolves every value at start-up, exports it as a shell
variable, and records the exact command in `log/<SAMPLE>/<step>.log`.

---

## 1. Global keys

| Key | Default in shipped config | Meaning | Change impact |
|---|---|---|---|
| `project_dir` | `/workspace/liangpeiyu/project/hantaKG_Agent` | Project root; the rRNA index is resolved relative to it | Wrong value makes the bowtie2 index unfindable |
| `conda_env` | `lpy-hanta-v2` | Conda environment used for every tool call (`conda run -n`) | Tools silently missing if wrong |
| `threads` | `64` | Total core budget | Samples run in parallel at `threads / 16`, each using 16 threads; the remainder stays idle rather than being oversubscribed |

`threads` is never passed directly to a tool. Each tool receives 16 threads and
the batch scheduler decides how many samples run concurrently:

```
PER_PROC_THREADS = 16
PARALLEL         = threads // 16   (minimum 1)
```

With `threads: 64` four libraries run at once. With `threads: 192` twelve run at
once. If `threads` is below 16, one sample runs and every tool receives the whole
budget.

---

## 2. `paths` — input and output roots

| Key | Shipped value | Meaning |
|---|---|---|
| `data_dir` | `.../data_cell_2/sra` | Directory scanned for paired-end reads |
| `results_dir` | `.../results_cell_paper` | Root of every output directory |

Optional overrides, resolved as `paths.<key>` and otherwise falling back to
`<results_dir>/<default>`:

| Key | Default |
|---|---|
| `fastp_dir` | `01_fastp` |
| `prinseq_dir` | `02_prinseq` |
| `clean_dir` | `03_clean_non_rRNA` |
| `megahit_dir` | `04_megahit` |
| `prodigal_dir` | `05_prodigal` |
| `log_dir` | `log/<SAMPLE>` |
| `json_dir` | `json` |

Leave these unset unless a project needs a different directory layout. Downstream
agents locate their inputs through the JSON, not through these names.

---

## 3. `input` — read discovery

| Key | Meaning |
|---|---|
| `r1_suffixes` | Comma-separated suffixes that identify the forward read |
| `r2_suffixes` | Comma-separated suffixes that identify the reverse read |

Shipped values cover the common SRA, ENA and in-house naming conventions:

```
_r1.fq.gz, _R1.fq.gz, .r1.fq.gz, .R1.fq.gz, _1.fq.gz, _1.fastq.gz, _R1.fastq.gz, _r1.fastq.gz
_r2.fq.gz, _R2.fq.gz, .r2.fq.gz, .R2.fq.gz, _2.fq.gz, _2.fastq.gz, _R2.fastq.gz, _r2.fastq.gz
```

The sample name is the filename with the matched suffix removed. When a library
is given explicitly with `-s`, the suffixes are still used to build the file
names, so a mismatch produces "input reads not found".

---

## 4. `cleanup` — intermediate removal

| Key | Shipped value | Meaning |
|---|---|---|
| `intermediates` | `true` | When `true`, the fastp, PRINSEQ++ and non-rRNA FASTQ files are deleted immediately after MEGAHIT succeeds |

**This is irreversible.** The fastp HTML and JSON reports are always kept. The
JSON `outputs` block is blanked for deleted files, so downstream code never sees
a dangling path.

Set it to `false` while developing or whenever the intermediate reads may be
needed for a second assembly attempt.

---

## 5. `fastp` — adapter and quality trimming

| Key | Shipped value | Meaning |
|---|---|---|
| `disable_adapter_trimming` | `false` | Turn the automatic adapter trimming off |
| `detect_adapter_for_pe` | `true` | Detect the adapter sequence for paired-end data |
| `cut_front` | `true` | Sliding-window trimming from the 5' end |
| `cut_tail` | `true` | Sliding-window trimming from the 3' end |
| `cut_right` | `true` | Sliding window moving right; everything after a low-quality window is dropped |
| `cut_front_window_size` | `5` | Window size in bp |
| `cut_tail_window_size` | `5` | Window size in bp |
| `cut_right_window_size` | `5` | Window size in bp |
| `cut_front_mean_quality` | `20` | Mean Phred threshold for the front window |
| `cut_tail_mean_quality` | `20` | Mean Phred threshold for the tail window |
| `cut_right_mean_quality` | `20` | Mean Phred threshold for the right-sliding window |
| `qualified_quality_phred` | `20` | A base counts as qualified at or above this Phred score |
| `unqualified_percent_limit` | `20` | Drop the read when more than this percentage of bases are unqualified |
| `length_required` | `50` | Drop reads shorter than this after trimming |
| `n_base_limit` | `0` | Drop reads containing more than this many `N` bases (`0` = none tolerated) |
| `extra_args` | empty | Appended verbatim to the fastp command line |

`fastp` writes `<SAMPLE>_fastp.html` and `<SAMPLE>_fastp.json` next to the
trimmed reads. These are the only per-library quality reports the workflow keeps
and they are what Figure 4b links to.

---

## 6. `prinseq` — low-complexity filtering

| Key | Shipped value | Meaning |
|---|---|---|
| `lc_entropy` | `0.5` | Discard sequences whose entropy is below this value |
| `lc_dust` | `0.5` | Discard sequences whose DUST score is above this value |
| `extra_args` | empty | Appended verbatim to the PRINSEQ++ command line |

PRINSEQ++ writes uncompressed FASTQ, which is the largest intermediate in the
pipeline. Output is split into `_good_out` (paired), `_single_out` (one mate
survived) and `_bad_out`; the agent resolves and uses the paired output.

---

## 7. `bowtie2` — rRNA depletion

| Key | Shipped value | Meaning |
|---|---|---|
| `index_prefix` | `rrna_index` | Index basename inside `<project_dir>/DB/sortmerna` |
| `preset` | `very-sensitive` | One of `very-fast`, `fast`, `sensitive`, `very-sensitive` |
| `ref_dir` | `<project_dir>/DB/sortmerna` | Directory holding the index |
| `extra_args` | empty | Appended verbatim to the bowtie2 command line |

The command uses `--un-conc-gz`, so **unmapped** reads — the ones that are not
rRNA — are written to disk; aligned reads go to `/dev/null` via `--no-unal` and
`-S /dev/null`. If the index has the wrong prefix, bowtie2 fails immediately with
a missing-index error rather than producing empty output.

---

## 8. `megahit` — assembly

| Key | Shipped value | Meaning |
|---|---|---|
| `min_contig_len` | `500` | Only contigs at least this long are written to `final.contigs.fa` |
| `overwrite` | `true` | Remove an existing output directory before assembling |
| `extra_args` | empty | Appended verbatim to the megahit command line |

With `overwrite: false`, a pre-existing `04_megahit/<SAMPLE>` aborts the run with
exit code 2 instead of silently reusing stale contigs.

Raising `min_contig_len` reduces spurious downstream hits but discards partial
genomes; lowering it increases recall at the cost of more fragmented ORFs. The
value of 500 matches the assembly settings reported in the manuscript.

---

## 9. `prodigal` — ORF prediction

| Key | Shipped value | Meaning |
|---|---|---|
| `mode` | `meta` | `meta` for metagenomes, `single` for an isolated genome |
| `min_protein_aa` | `200` | Proteins shorter than this are removed into a separate file |
| `extra_args` | empty | Appended verbatim to the prodigal command line |

Prodigal runs with `-q`, so no banner reaches the log. Three outputs are written
(`.faa`, `.ffn`, `.gff`), then a length filter produces
`<SAMPLE>_final.contigs_proteins.min200aa.faa`. Only this filtered file is handed
to the classification agent.

---

## 10. Sections this agent ignores

`software`, `diamond`, `mafft`, `trimal`, `iqtree3`, `phylogeny`,
`vmr_merge_cluster`, `tree`, `report` and any other section in the shared
configuration belong to downstream agents. Editing them has no effect here.
