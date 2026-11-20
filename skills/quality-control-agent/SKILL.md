---
name: quality-control-agent
description: Raw-read quality control and assembly for metavirome libraries (fastp adapter/quality trimming, PRINSEQ++ low-complexity filtering, bowtie2 rRNA depletion, MEGAHIT assembly, Prodigal ORF prediction with a >=200 aa filter). Produces the contigs and protein sequences that every downstream step consumes, plus a structured per-sample JSON result for the orchestrator. Agent 1 of the six-agent Hanta-Agent workflow.
---

# Quality Control Agent

Turn paired-end raw reads into assembled contigs and predicted proteins. This is
the first agent of the Hanta-Agent workflow: nothing downstream can start until
it has produced `final.contigs.fa` and a protein FASTA.

The agent is a silent bash pipeline. It writes **no progress output to stdout**;
success and failure are reported exclusively through `json/<SAMPLE>.qc.json`.
Reading that file is the only reliable way to know what happened.

## When to use

- Raw paired-end FASTQ for one library, a subset of libraries, or a whole batch
  must be converted into contigs and proteins.
- A library failed downstream and needs to be re-processed from reads.
- A new batch of libraries has been deposited and the classification agent has
  nothing to consume yet.

## When not to use

- The inputs are already contigs or proteins — start at the classification agent.
- Only one step of this agent needs re-running (for example re-assembly only).
  Do not re-run the whole agent to fix one stage; adjust the configuration or use
  the parent workflow's per-step entry point instead.
- The task is alignment, phylogeny, risk scoring, or reporting — those are the
  downstream agents.

## Contract

| Item | Value |
|---|---|
| **Input** | `paths.data_dir/<SAMPLE>_1.fastq.gz` and `<SAMPLE>_2.fastq.gz` (suffixes configurable under `input:`) |
| **Config** | `../config.yaml`, two directories above `scripts/` — the agent locates it automatically |
| **Prerequisite environment** | conda env named by `conda_env` (default `lpy-hanta-v2`) containing fastp, prinseq++, bowtie2, megahit, prodigal |
| **Prerequisite data** | bowtie2 rRNA index at `<project_dir>/DB/sortmerna/<index_prefix>` (default `rrna_index`) |
| **Output** | Five stage directories under `paths.results_dir` plus `json/<SAMPLE>.qc.json` and `log/<SAMPLE>/` — see `references/outputs-and-contract.md` |
| **Downstream interface** | `outputs.contigs_fasta` and `outputs.prodigal_min_aa_faa` inside the JSON |
| **Failure behaviour** | Non-zero exit code, and a JSON with `status: failed` plus `failed_step` and `message` |

## Workflow

1. **Confirm the configuration.** Read `../config.yaml` and check that
   `paths.data_dir` holds the intended reads and `paths.results_dir` is the
   intended output root. Both are absolute paths shared by all six agents, so a
   mistake here sends results to the wrong place.
2. **Dry run first — in a scratch directory.** Copy `config.yaml`, change
   `paths.results_dir` to a throwaway path, and only then invoke with
   `--dry-run`. The dry run prints the resolved commands, but it also deletes an
   existing MEGAHIT directory and overwrites the sample's JSON and logs (see
   Invariant 2). Inspect the fastp thresholds, the MEGAHIT `--min-contig-len`,
   and the Prodigal `-p` mode, then discard the scratch directory.
3. **Run.** Either one sample (`-s <SAMPLE>`) or the whole batch (omit `-s`; the
   agent discovers every R1 file in `data_dir` and runs samples in parallel).
   Parallelism is derived automatically: `threads / 16` samples at a time, each
   using 16 threads, so the configured core budget is never oversubscribed.
4. **Read the JSON, not the screen.** Open
   `<results_dir>/json/<SAMPLE>.qc.json`. `status: success` means every stage
   finished; anything else must be diagnosed from `failed_step` and then from
   `log/<SAMPLE>/<failed_step>.log`.
5. **Hand off.** Pass `outputs.contigs_fasta` and `outputs.prodigal_min_aa_faa`
   to the classification agent. Do not hard-code these paths downstream.

## Invariants

1. **Never judge success from stdout.** The pipeline prints nothing on purpose.
   Only `json/<SAMPLE>.qc.json` carries the verdict.
2. **`--dry-run` is not read-only.** Since the 2026-11-20 patch it no longer
   deletes an existing MEGAHIT directory, but it still writes a
   `status: success` record into `json/<SAMPLE>.qc.json` and overwrites every
   file in `log/<SAMPLE>/`. An earlier version executed `rm -rf` on
   `04_megahit/<SAMPLE>` even in a dry run, which permanently destroyed an
   assembly. Prefer a scratch `paths.results_dir` for any dry run.
3. **`cleanup.intermediates` is destructive.** When `true`, the fastp, PRINSEQ++
   and non-rRNA FASTQ files are deleted as soon as MEGAHIT finishes, and they
   cannot be recovered. Confirm the current value before every run; the shipped
   configuration has it set to `true`.
4. **`megahit.overwrite: true` removes `results_dir/04_megahit/<SAMPLE>`.** The
   script refuses to delete any path that does not match that exact pattern, so
   never rewrite the variable to point elsewhere.
5. **Do not edit the JSON by hand.** It is the sole machine-readable interface
   between this agent and the orchestrator.
6. **Downstream depends on two output keys only** — `outputs.contigs_fasta` and
   `outputs.prodigal_min_aa_faa`. Treat every other path as internal.

## Decision points

| Situation | Adjustment | Risk |
|---|---|---|
| Reads carry obvious adapter contamination | Set `fastp.disable_adapter_trimming: false` and `detect_adapter_for_pe: true` | None; this is the intended setting |
| Short-insert libraries lose too many reads | Lower `fastp.length_required` below 50 | Shorter reads can create chimeric assemblies |
| Assemblies are fragmented | Lower `megahit.min_contig_len` below 500 | More short contigs, more false downstream hits |
| RdRp domains are being missed | Lower `prodigal.min_protein_aa` below 200 | Partial ORFs enter the tree with weak support |
| Disk pressure on intermediate FASTQ | Set `cleanup.intermediates: true` | Irreversible loss of intermediate reads |
| rRNA depletion is too weak | Raise `bowtie2.preset` to `very-sensitive` (already the default) | Slower run |

Full parameter reference: `references/parameters.md`.

## Failure handling

| Exit code | Meaning | First action |
|---|---|---|
| 1 | A pipeline stage failed | Read `log/<SAMPLE>/<failed_step>.log` |
| 2 | Bad arguments, missing config, or an existing MEGAHIT directory with overwrite off | Fix the invocation or the directory, then retry |
| 127 | A required executable is missing | Activate or repair the conda environment |

Per-stage diagnosis: `references/troubleshooting.md`.

## References

- `references/parameters.md` — every configuration key, its default, and the
  effect of changing it.
- `references/outputs-and-contract.md` — output file layout and the complete
  JSON schema.
- `references/troubleshooting.md` — failure branches and how to resolve them.

## Provenance

The executable in `scripts/1-qc.sh` is a copy of `../../1-qc.sh` at the time
this skill was built. The authoritative pipeline source is
`skills_new/1-qc.sh`; if either copy changes, both must be updated together.
Configuration is a copy of `skills_new/config-cell-paper.yaml`.
