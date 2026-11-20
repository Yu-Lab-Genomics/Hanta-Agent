# Troubleshooting — Sequence Handle Agent

Read `json/<SAMPLE>.species_class.json` first: `status`, then `message`, then the
affected `logs.<db>` file. Every log begins with the full DIAMOND command.

Exit codes: `2` for anything the script detects and handles itself, the DIAMOND
return code for anything it does not.

---

## 1. Configuration failures (exit 2)

| Message | Cause | Fix |
|---|---|---|
| `config file not found: <path>` | Neither `-c` nor the upward search found a configuration | Pass `-c /abs/path/config.yaml`. The search looks for `config.yaml`, then `config-cell-paper.yaml`, then `cell_submit-paper.yaml`, walking up four levels |
| `project_dir is missing from config.yaml` | The key is absent or empty | `project_dir` is mandatory; every database path is derived from it |
| `JSON directory not found: <path>` | `paths.json_dir` (or `<results_dir>/json`) does not exist | Run the quality-control agent first; the directory is created by that agent |
| `DIAMOND DB directory not found: <path>` | `<project_dir>/DB/diamond` does not exist | Check `project_dir`; this directory must hold `vmr_l.dmnd`, `vmr_m.dmnd`, `vmr_s.dmnd` |
| `no processable contigs_fasta found` | No JSON in the directory exposes `outputs.contigs_fasta` | The quality-control JSONs are missing, or all of them lack that key. Confirm the upstream agent ran into the same `results_dir` |

A wrong `results_dir` is the most common root cause of the last two: the agent
reads JSONs from one project and databases from another, and reports an empty
sample set.

---

## 2. Per-sample failures (exit 2)

| Message | Cause | Fix |
|---|---|---|
| `contigs_fasta does not exist: <path>` | The upstream JSON names a contig file that was deleted or never written | Check the upstream record's `status` and `outputs.contigs_fasta`; re-run the quality-control agent for that sample |
| `database does not exist: <path>` | One of the three `.dmnd` files is missing | Rebuild the affected database, or restore the file. The sample stops at the first missing database, so check all three |
| `output already exists, pass --overwrite to replace it: <path>` | A hit table from an earlier run is still on disk | Add `--overwrite`, or move the old tables aside |
| `<db> alignment failed, log: <log>` | DIAMOND returned non-zero for that database | Read the log. A malformed query, a corrupted `.dmnd`, or a disk-full condition are the usual causes |

`--overwrite` is checked per database, in the order `vmr_l`, `vmr_m`, `vmr_s`.
The first existing table aborts the sample, so a partially completed previous run
must be cleaned up entirely before a retry.

---

## 3. Silent problems that do not raise an error

| Symptom | Explanation | Action |
|---|---|---|
| `status: success` but a table contains only the three header lines | The contigs had no hit in that database | Legitimate. If all three are empty, the sample has no detection |
| `status: success` and the JSON was written by a `--dry-run` | The dry run skips DIAMOND but still records success | Delete the JSON and re-run for real |
| A sample that was never assembled gets picked up | Discovery is by presence of `outputs.contigs_fasta`, not by upstream `status` | It fails on the missing file; either ignore it or remove the stale upstream record |
| Contig appears in both `vmr_l` and `vmr_m` | `max_target_seqs` is per database | Expected. The dereplication agent resolves the priority; do not filter here |
| Hit rate drops after a database update | A `.dmnd` was rebuilt from a different reference set | Compare `DB/diamond/vmr_*.fasta` against the manifest before rerunning a batch |
| The run stops with no per-sample JSON | An unparsable JSON in the directory raised `Failed to read JSON` | Find and repair the malformed file; the agent aborts the whole run rather than skip it |

---

## 4. Recovery patterns

**Retry one sample.** Fix the cause and re-run with `-s <SAMPLE> --overwrite`.
Because the three databases are processed in sequence, a failure part-way through
leaves the earlier tables on disk; `--overwrite` clears them.

**Retry the whole batch.** Use `--overwrite`, otherwise every sample that already
has tables fails immediately and the run reports nothing useful.

**A sample keeps failing on the same database.** Check whether that `.dmnd` is
the only one built from a different reference release; a partially updated
database set produces inconsistent hit counts across the three segments.

**Never delete the species_class JSON to "reset" a sample.** Re-running
overwrites it, and deleting it removes the record of which query was aligned.

---

## 5. Interaction with the quality-control agent

This agent has no way to verify that an assembly is complete. If the upstream
JSON says `success` but the contig FASTA is truncated or empty, DIAMOND returns
zero with an empty output and this agent records `success` with empty tables. A
batch-level sanity check on total hit counts catches that faster than inspecting
individual samples.

The reverse is not a risk: this agent's JSON never contains
`outputs.contigs_fasta`, so it cannot be mistaken for an upstream record by a
later run.
