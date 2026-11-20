# Troubleshooting — Quality Control Agent

The agent prints nothing to stdout. Every diagnosis starts from
`json/<SAMPLE>.qc.json`: read `status`, then `failed_step`, then open
`log/<SAMPLE>/<failed_step>.log`. The log begins with the exact command that was
executed, which is normally enough to identify the problem.

Exit codes: `1` = a tool stage failed, `2` = invocation or configuration problem,
`127` = an executable is missing.

---

## 1. Preflight failures (`failed_step: preflight`, exit 127)

| Message | Cause | Fix |
|---|---|---|
| `required tool not found: <tool>` | A tool is not on `PATH` | Run everything through `conda run -n <conda_env>`; confirm `conda_env` in `config.yaml` |
| `conda command not found` | `conda` itself is not on `PATH` | Source the conda profile, or point `PATH` at the installation |
| `conda environment <env> is not usable` | The named environment does not exist | `conda env list`; correct `conda_env` |
| `python3 is required to parse config.yaml` | No `python3` on `PATH` | The YAML reader is Python-based; set `PYTHON=/abs/path/python3` |

These checks run before any data is touched, so no partial output is left behind.

---

## 2. Configuration failures (exit 2)

| Message | Cause | Fix |
|---|---|---|
| `config file not found: <path>` | The default search failed and no `-c` was given | The agent expects `config.yaml` two directories above `scripts/`. Confirm `final/config.yaml` exists, or pass `-c` explicitly |
| `failed to parse config.yaml` | Malformed YAML, or a value that cannot be read | Validate the YAML; check for tabs, unbalanced quotes, or a missing `project_dir` |
| `unknown option: <flag>` | Unsupported argument | Only `-s`, `-c`, `--dry-run` and `-h` are accepted |
| `unexpected argument: <arg>` | A second positional argument | Exactly one sample name may be positional |
| `r1_suffixes and r2_suffixes must have the same number of entries` | The two comma-separated lists differ in length | Keep the lists parallel, one R2 suffix per R1 suffix |

---

## 3. Input discovery failures (exit 2)

| Message | Cause | Fix |
|---|---|---|
| `no paired-end FASTQ files found for sample <SAMPLE> in <data_dir>` | The sample name does not match the file on disk, or the suffix list misses the naming convention | Check the actual filename and either rename or extend `input.r1_suffixes` / `input.r2_suffixes` |
| `no sample could be auto-discovered in <data_dir> and -s was not given` | `data_dir` is empty, wrong, or the suffixes match nothing | Verify `paths.data_dir`; list a few filenames and compare against the suffix list |
| `cannot parse paired-end output files: <prefix>` | PRINSEQ++ produced no `_good_out` pair, typically because every read was rejected | Inspect `02_prinseq/`; the `_bad_out` and `_single_out` files show what survived |

The most common cause of a "sample not found" error is a case mismatch: `_R1` vs
`_r1`, or `_1.fastq.gz` vs `_1.fq.gz`. Both are covered by the shipped suffix
list, but a non-standard naming scheme is not.

---

## 4. Stage failures

### 4.1 fastp

Symptom: exit 1, `failed_step: fastp`, log ends with a fastp error. Typical
causes are a truncated `.gz` file or a read pair where R1 and R2 have different
read counts.

Check in the log:

```bash
grep -E "ERROR|WARNING|reads" <results_dir>/log/<SAMPLE>/fastp.log
```

If `after_filtering.total_reads` in `<SAMPLE>_fastp.json` is zero, the quality
thresholds are too strict for this library — raise `length_required` down or
relax `qualified_quality_phred` only after confirming the run is genuinely low
quality.

### 4.2 PRINSEQ++

Symptom: exit 1, or `cannot parse paired-end output files: <prefix>`. Entropy and DUST thresholds of `0.5`
are permissive; a library that loses everything here is usually mis-specified
(single-end data supplied as paired) rather than genuinely low-complexity.

### 4.3 bowtie2

| Message | Cause | Fix |
|---|---|---|
| `bowtie2 produced no non-rRNA read output` | `--un-conc-gz` wrote nothing, usually a missing or mismatched index | Confirm `<project_dir>/DB/sortmerna/<index_prefix>.*.bt2` exists; check `index_prefix` |
| `no usable reads remain after rRNA removal` | Every read mapped to rRNA, or the input files are empty | Check that the upstream PRINSEQ++ output is non-empty; a genuinely rRNA-dominated library must be reported as failed rather than pushed downstream |

Both are treated as fatal on purpose: an empty depleted library cannot produce a
meaningful assembly.

### 4.4 MEGAHIT

| Message | Cause | Fix |
|---|---|---|
| `MEGAHIT output already exists; set megahit.overwrite: true or delete <dir> manually` | `overwrite` is `false` and the directory exists | Either set `overwrite: true` or remove the directory deliberately |
| `refusing to delete an unexpected MEGAHIT directory: <path>` | Safety guard triggered because the computed path does not match `<megahit_dir>/<SAMPLE>` | Never rewrite the path variables. This guard exists so a stray configuration cannot trigger a recursive delete |
| `MEGAHIT produced no final.contigs.fa` | Assembly finished but no contig passed `--min-contig-len` | Lower `megahit.min_contig_len`, or accept that the library has no assembable content |

`megahit failed` with exit 1 usually means memory exhaustion. MEGAHIT scales
with the number of concurrent samples; reduce `threads` so that fewer samples run
in parallel.

### 4.5 Prodigal

| Message | Cause | Fix |
|---|---|---|
| `Prodigal produced no protein sequences` | The contig FASTA is empty or unreadable | Verify `04_megahit/<SAMPLE>/final.contigs.fa` is non-empty |
| `prodigal failed` | `-p meta` vs `-p single` mismatch, or a corrupted contig file | For isolate data set `prodigal.mode: single` |

A successful Prodigal run can still yield an empty `min200aa.faa` when every
predicted protein is shorter than the threshold. That is not an error — the
filtered file simply has no records — but the sample should then be reported as
producing no usable protein for downstream classification.

---

## 5. Silent problems that do not raise an error

| Symptom | Explanation | Action |
|---|---|---|
| `outputs.fastp_r1` is missing from the JSON although the stage succeeded | `cleanup.intermediates: true` deleted it after MEGAHIT | Expected. The fastp HTML/JSON reports remain |
| **A whole `04_megahit/<SAMPLE>` directory vanished after a `--dry-run`** | Fixed on 2026-11-20: the overwrite guard used to run before the dry-run check, so `rm -rf` executed even in a dry run. Scripts older than that patch still carry the bug | Restore from a backup, or re-run the sample for real. Never dry-run in a production `results_dir` |
| A `--dry-run` left a `status: success` JSON behind | The dry run reports success for skipped stages — this is still true after the patch | Delete `json/<SAMPLE>.qc.json` before the real run, or use a throwaway `results_dir` |
| A `--dry-run` overwrote the sample's logs | The dry run still writes `log/<SAMPLE>/*.log`, replacing earlier runs | Use a throwaway `results_dir` when the log content matters |
| Sample appears successful but downstream finds no hits | Assembly produced only very short contigs, or the library genuinely contains no target | Inspect `final.contigs.fa` length distribution; compare against `min_contig_len` |
| Run takes far longer than expected | Parallelism is `threads / 16`; a low `threads` value serialises the batch | Raising `threads` also raises memory pressure on MEGAHIT |
| Disk fills during the run | PRINSEQ++ writes uncompressed FASTQ, the largest intermediate | Set `cleanup.intermediates: true`, or clean `02_prinseq/` between batches |

---

## 6. Recovery patterns

**Retry a single sample.** Correct the cause, then re-invoke with `-s <SAMPLE>`.
With `megahit.overwrite: true` the previous partial assembly is removed first.
Stages before MEGAHIT re-run from scratch, because each stage overwrites its own
output.

**Retry only the assembly.** Not supported inside this agent. Re-running the
whole sample is the supported path; the cost is dominated by PRINSEQ++ and
MEGAHIT, not by fastp.

**Investigate a batch-wide failure.** If every sample fails at the same stage,
the cause is environmental (missing environment, missing index, wrong
`project_dir`) rather than per-library. Confirm with one sample before touching
the batch.

**Never** delete `json/<SAMPLE>.qc.json` to "reset" a sample. Re-running
overwrites it; deleting it only removes the only record of what happened.
