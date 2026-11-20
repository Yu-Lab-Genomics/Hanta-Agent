---
name: sequence-handle-agent
description: Align assembled contigs against the three hantavirus reference databases (vmr_l, vmr_m, vmr_s) with DIAMOND blastx and record which segment database each contig hits. Consumes the per-sample JSON written by the quality-control agent and writes one hit table per database plus a structured JSON result. Agent 2 of the six-agent Hanta-Agent workflow.
---

# Sequence Handle Agent

Take the contigs produced by the quality-control agent and align them, at the
protein level, against the three hantavirus reference databases. The output is
the raw hit evidence that dereplication and phylogeny build on.

The agent discovers its own input. It scans the JSON directory for any record
that exposes `outputs.contigs_fasta` and aligns each one against every reference
database. **There is no need to pass a contig path** — doing so is not supported.

## When to use

- The quality-control agent has finished for one or more samples and their JSON
  records show `status: success`.
- Contigs must be screened against the hantavirus L, M and S reference databases
  before dereplication and tree building.
- A reference database has been rebuilt and the alignment must be refreshed.

## When not to use

- Contigs do not exist yet — run the quality-control agent first.
- The task is dereplication or clustering of the hits — that is the next agent.
- The task is tree building or species assignment — later agents.
- The input is protein FASTA rather than assembled nucleotide contigs. DIAMOND is
  invoked with `blastx`, which translates the query itself.

## Contract

| Item | Value |
|---|---|
| **Input** | Every `json/*.json` containing `outputs.contigs_fasta`; the contig FASTA itself is read from that path |
| **Upstream requirement** | The producing record must have `status: success`, otherwise the contig file may be missing or empty |
| **Config** | `-c <path>`. The default is resolved by walking up from the script and lands on `final/config.yaml`; an explicit `-c` is still preferred for production runs |
| **Prerequisite data** | `DB/diamond/vmr_l.dmnd`, `vmr_m.dmnd`, `vmr_s.dmnd` |
| **Prerequisite environment** | conda env named by `conda_env`, containing DIAMOND |
| **Output** | `06_diamond_vmr/<SAMPLE>.vmr_{l,m,s}.tsv`, `json/<SAMPLE>.species_class.json`, `log/<SAMPLE>/species_class_vmr_{l,m,s}.log` |
| **Downstream interface** | `outputs.vmr_l`, `outputs.vmr_m`, `outputs.vmr_s` inside the JSON |
| **Failure behaviour** | Exit code 2 for configuration, database or overwrite problems; the DIAMOND exit code propagates otherwise. A failed sample still gets a JSON with `status: failed` |

## Workflow

1. **Confirm the upstream state.** Pick the samples whose `<SAMPLE>.qc.json` has
   `status: success` and a non-empty `outputs.contigs_fasta`. Samples still in
   progress will be discovered but will fail on the missing file.
2. **Confirm the databases.** All three `.dmnd` files must exist in
   `DB/diamond/`. One missing database aborts that sample immediately with
   "database does not exist".
3. **Preview.** Run with `--dry-run` and read
   `log/<SAMPLE>/species_class_vmr_*.log` to check the command, the E-value
   threshold and the target-count limit. Unlike the quality-control agent, this
   dry run deletes nothing — but it does still write a `status: success` JSON.
4. **Run.** One sample with `-s`, or every discovered sample by omitting it.
   **Add `--overwrite` whenever the hit tables may already exist**, otherwise
   each such sample stops with "output already exists".
5. **Read the JSON.** `status: success` and three populated `outputs` keys mean
   the sample is ready for dereplication. Anything else must be diagnosed from
   `message`.

## Invariants

1. **Always pass `-c` explicitly in production.** The fallback search works, but
   it may pick `config-cell-paper.yaml` or `cell_submit-paper.yaml` in a
   directory that holds more than one configuration. Being explicit removes the
   ambiguity.
2. **Re-running over existing output requires `--overwrite`.** Without it the
   sample fails with exit code 2 before DIAMOND starts.
3. **Never hand a contig path to the agent.** Discovery goes through the JSON so
   that a sample which failed upstream cannot be silently aligned.
4. **`--dry-run` also writes `status: success`.** Delete that JSON before the
   real run, or the orchestrator will believe the sample is finished.
5. **Downstream relies on `outputs.vmr_l`, `outputs.vmr_m` and `outputs.vmr_s`.**
   The file names encode the sample and the database; do not rename them.
6. **Do not raise `diamond.evalue` to make marginal hits appear.** The threshold
   of 1e-5 is shared with every downstream claim about detection.

## Decision points

| Situation | Adjustment | Risk |
|---|---|---|
| Fewer hits than expected | Check that all three `.dmnd` files were rebuilt after the last database update | Comparing against a stale database silently loses species |
| Too many weak hits | Lower `diamond.evalue` below 1e-5 | Fewer distant-but-real hits; document any change |
| A contig hits two databases | Expected — `max_target_seqs` is per database, so a contig can appear in both `vmr_l` and `vmr_m` | The next agent resolves the priority; do not filter here |
| Long run time on a large batch | Raise `threads` | DIAMOND scales well but competes with other agents for the same cores |

Full parameter reference: `references/parameters.md`.

## Failure handling

| Exit code | Meaning | First action |
|---|---|---|
| 2 | Missing config, missing project_dir, missing JSON directory, missing DB directory, no samples found, missing contig FASTA, missing `.dmnd`, or pre-existing output without `--overwrite` | Read `message` in the sample JSON |
| other non-zero | DIAMOND itself failed | Read `log/<SAMPLE>/species_class_<db>.log` |

Per-case diagnosis: `references/troubleshooting.md`.

## References

- `references/parameters.md` — configuration keys, the upstream JSON discovery
  rule, and every DIAMOND argument.
- `references/outputs-and-contract.md` — hit-table layout and the complete JSON
  schema.
- `references/troubleshooting.md` — all nine failure branches.

## Provenance

The executable in `scripts/2-species_class.py` is a copy of
`../../2-species_class.py`. Both copies share the same fix: the default
configuration is now located by walking up from the script instead of being
hard-coded to `skills_V2/config.yaml`. The authoritative pipeline source is
`skills_new/2-species_class.py`; if either copy changes, both must be updated
together.
