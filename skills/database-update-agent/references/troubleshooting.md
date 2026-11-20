# Troubleshooting — Database Update Agent

## Symptom: a rebuild succeeded but nothing changed downstream

**The most likely cause.** The files were written somewhere the consumers do not
read.

This is not hypothetical. Before 2026-11-20, `build_clean_db.py` wrote its
DIAMOND indexes to `DB/clean/vmr_*.dmnd` while agents 2 and 4 read
`DB/diamond/vmr_*.dmnd`. The script reported success, the index existed, and the
pipeline kept using the old database. Nothing in the output said anything was
wrong.

Check the actual paths:

```bash
ls -la <project_dir>/DB/diamond/vmr_*.dmnd
ls -la <project_dir>/DB/tree/vmr_L_aln.fasta
```

If either is missing or older than the rebuild, the rebuild went elsewhere. Run
the script with `--dry-run` and read the `DIAMOND output` / `output:` lines — they
name the directory that will be used.

Then confirm `project_dir` in `config.yaml` is the directory that contains `DB/`.
That single key decides where all six agents look.

## `Alphabet 'e' is unknown`

**Where it appears.** `4-vmr-tree.py`, when MAFFT parses `tree_input.fa`.

**What it means.** MAFFT was handed nucleotides and found the letter `E`, which is
glutamic acid. The nucleotide alphabet has no `E`, so parsing stops. The rejected
file is the concatenation of the backbone `DB/tree/vmr_L_aln.fasta` and the query
proteins. The query side is protein by construction, so the backbone is the
nucleotide one.

**Why it happens.** An earlier version of the backbone builder filtered nucleotide
records and wrote them out unaligned under an `*_aln.fasta` name. The name
suggested an alignment; the contents were raw nucleotides.

**Fix.**

1. Confirm the diagnosis:
   ```bash
   grep -v '^>' <project_dir>/DB/tree/vmr_L_aln.fasta | head -c 200
   ```
   If the output is essentially `ACGTN`, the file is nucleotide.
2. Rebuild it:
   ```bash
   python build_clean_db.py
   python build_tree_ref.py
   ```
3. Re-run the same `grep` and confirm `E`, `L`, `K` and `-` appear.
4. Re-run the tree agent.

## Database directory cannot be resolved

| Message | Cause | Fix |
|---|---|---|
| `Cannot determine the database directory` | No `--db-dir`, and no configuration found upward from the script | Pass `--db-dir <path>`, or `--config <path>` |
| `project_dir is missing from <config>` | The configuration exists but has no `project_dir` key | Add `project_dir: /path/to/project` at the top level |
| `no configuration file found` (bash script) | Same as above, from `fetch_rrna_index.sh` | Pass `--db-dir` or `-c` |

## Input problems

| Message | Cause | Fix |
|---|---|---|
| `cleaned library directory not found: <db>/clean` | The curated nucleotide library was never installed, or `project_dir` is wrong | Install `DB/clean/` with the three `*_clean.fasta` files |
| `input missing: <db>/clean/vmr_L_clean.fasta` | The directory exists but a segment file does not | Supply all three segments, or accept that the library is incomplete |
| `manifest.csv not found` | `install_manifest.py` has not run | Run it, or pass `--source` |
| `protein library vmr_l.fa not found` | `build_clean_db.py` has not run, or wrote elsewhere | Run it, then re-check `DB/diamond/` |
| `source manifest not found` | `clean/manifest_clean.csv` is absent | Supply it with the cleaned library, or pass `--source` |
| `source manifest is missing required column(s)` | The manifest is malformed | Columns `acc_base`, `species` and `segment` are mandatory |

## Build failures

| Symptom | Cause | Fix |
|---|---|---|
| Prodigal produces no proteins | Input is amino acids, or is empty | Check with `grep -v '^>' <seg>_clean.fasta \| head -c 100` |
| Translation coverage far below 1.0 | The library is mostly fragments | Verify this is intended; short records legitimately yield no protein |
| DIAMOND `makedb` fails | `vmr_<seg>.fa` is empty | Inspect the Prodigal output and the coverage line |
| MAFFT takes hours | `--auto` picked a slow strategy for a large alignment | Raise `--threads`; do not force a preset |
| `fewer than 2 usable references, cannot align` | Almost no representative has a protein in `vmr_l.fa` | Check `protein library hits` and whether Prodigal ran on the right input |
| Tree has one tip per species but a species is missing | No representative passed `is_complete()`, or the species name is inconsistent in the manifest | Check the manifest for name variants |

## rRNA index problems

| Symptom | Cause | Fix |
|---|---|---|
| `--check` reports files MISSING | Never downloaded, or downloaded elsewhere | Rerun without `--check`, with `--url` |
| `checksum verification failed` | A truncated or corrupted download | Delete the named file and rerun; the downloader retries three times before giving up |
| Download fails part way | Network interruption on a 1.2 GB transfer | Rerun; already-complete files are skipped |
| `index already installed and complete` when you wanted to re-fetch | The guard is doing its job | Add `--force` |
| bowtie2 reports it cannot find the index | The prefix is right but the files are not | Confirm `ls <project_dir>/DB/sortmerna/rrna_index*.bt2` and that `project_dir` is correct |

## Verifying a database without running the pipeline

All four scripts have a safe mode. Use them before committing to a long build:

```bash
python build_clean_db.py --dry-run      # composed prodigal and diamond commands
python build_tree_ref.py --dry-run      # composed MAFFT command
python install_manifest.py --check      # column validation, writes nothing
bash fetch_rrna_index.sh --check        # index completeness, writes nothing
```

Then read the printed paths. If they do not sit under the `DB/` that
`config.yaml` points at, stop and fix `project_dir` instead of running the build.

## Safety rules

- **Never rebuild against a live results directory while an analysis is running.**
  The DIAMOND indexes are read by agent 2; replacing them mid-run mixes two
  database versions in one result set.
- **Back up `DB/clean/` before rebuilding.** The Prodigal step overwrites
  `*_proteins.faa` and `*_genes.ffn` in place.
- **Treat a successful exit as the start of verification, not the end.** Both
  Python scripts report success once the tools return 0; whether the outputs are
  correct requires the checks above.
- **`fetch_rrna_index.sh` never deletes.** A failed download leaves a `.partial`
  file, so an interrupted transfer is never mistaken for a complete index.
