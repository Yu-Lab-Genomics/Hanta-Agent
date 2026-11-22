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
| `manifest.csv not found` | The raw metadata table is absent | Restore it from the repository. It is a required input, not something this agent produces |
| `nucleotide input missing for segment(s): L` | One of `diamond/vmr_{L,M,S}_99.fasta` is absent | Restore it from the repository. All three are needed before curation can run |
| `cleaned library directory not found: <db>/clean` | `clean_db.py` has not run, or `project_dir` is wrong | Run `clean_db.py`. It needs `DB/manifest.csv` and `DB/diamond/vmr_{L,M,S}_99.fasta` in place |
| `input missing: <db>/clean/vmr_L_clean.fasta` | `clean_db.py` ran but this segment produced nothing | Check the `vmr_L_clean.fasta: N sequences` line it printed; zero means no record of that segment survived curation |
| `protein library vmr_l.fa not found` | `build_clean_db.py` has not run, or wrote elsewhere | Run it, then re-check `DB/diamond/` |
| `source manifest not found` | `clean/manifest_clean.csv` is absent | Run `clean_db.py` first, or pass `--source` |
| `source manifest is missing required column(s)` | The manifest is malformed | Columns `acc_base`, `species` and `segment` are mandatory |

### clean_db.py keeps almost nothing

If `kept records` is far below `input records`, read the `dropped ...` lines
printed between them. Each names one rule and counts what it removed:

| Line | What triggered it |
|---|---|
| `dropped placeholder species` | The species name ends in `sp.` or is a bare genus name |
| `dropped vague name` | The name carries no place name and cannot be tied to a location |
| `dropped neither complete nor long enough` | NCBI does not call the record complete, and it is under the segment floor (L 5,000 / M 3,000 / S 1,400 nt) |

If the counts look wrong, the input manifest is probably not the one you expect.
Check its row count against the 2,504 of the released database first.

### clean_db.py wrote fewer sequences than manifest rows

This is expected, not a bug. Curation writes a sequence only when its
`acc_base` appears in the curated manifest, and `diamond/vmr_{L,M,S}_99.fasta`
holds far more records than the manifest does:

| Segment | In `vmr_<seg>_99.fasta` | In `manifest.csv` | Written to `clean/` |
|---|---|---|---|
| L | 1,895 | 570 | 570 |
| M | 1,682 | 694 | 694 |
| S | 2,438 | 1,240 | 1,240 |

The nucleotide file is the CD-HIT pool; the manifest decides what the database
contains. If a sequence you need is missing from `clean/`, add its record to the
manifest rather than to the FASTA.

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
| `--check` reports files MISSING | Never downloaded, or downloaded elsewhere | Rerun without `--check`; `--url` defaults to the release archive |
| `download failed: ... (partial data left at ...)` | Network interruption on a 1.2 GB transfer | Rerun. An archive already on disk is reused unless `--force` is given |
| `archive checksum mismatch` | The archive is truncated, or the expected hash is stale | Delete `<db>/sortmerna/rrna_index.tar.gz` and rerun. If it keeps failing after a re-download, re-upload the asset and refresh `ARCHIVE_SHA256` |
| `archive does not contain: ...` | The tarball is missing one or more `.bt2` files | Repack the archive. The extracted tree is kept under `.extract-*` so you can see what was actually inside |
| `extraction failed` | Corrupt gzip stream, or a full disk | Check `df -h <db>/sortmerna`; the archive is left in place so it can be retested |
| `no SHA256 configured for the archive` | `ARCHIVE_SHA256` ships empty | Not an error. Fill it in after uploading, or pass `--sha256` per run |
| `checksum verification failed` (files mode) | A truncated or corrupted `.bt2` | Delete the named file and rerun; the downloader retries three times before giving up |
| `index already installed and complete` when you wanted to re-fetch | The guard is doing its job | Add `--force` |
| bowtie2 reports it cannot find the index | The prefix is right but the files are not | Confirm `ls <project_dir>/DB/sortmerna/rrna_index*.bt2` and that `project_dir` is correct |

## Verifying a database without running the pipeline

All five scripts have a safe mode. Use them before committing to a long build:

```bash
python clean_db.py --dry-run           # curation counts, writes nothing
python build_clean_db.py --dry-run     # composed prodigal and diamond commands
python build_tree_ref.py --dry-run     # composed MAFFT command
python install_manifest.py --check     # column validation, writes nothing
bash fetch_rrna_index.sh --check       # index completeness, writes nothing
```

Then read the printed paths. If they do not sit under the `DB/` that
`config.yaml` points at, stop and fix `project_dir` instead of running the build.

## Safety rules

- **Never rebuild against a live results directory while an analysis is running.**
  The DIAMOND indexes are read by agent 2; replacing them mid-run mixes two
  database versions in one result set.
- **`DB/clean/` is safe to delete; the two raw inputs are not.** `clean_db.py`
  regenerates `clean/` from `DB/manifest.csv` and `diamond/vmr_{L,M,S}_99.fasta`
  in seconds. Those two inputs come from NCBI and nothing here can recreate them.
- **Back up `DB/clean/` before rerunning `build_clean_db.py`** if you care about
  the previous Prodigal output: the step overwrites `*_proteins.faa` and
  `*_genes.ffn` in place.
- **Treat a successful exit as the start of verification, not the end.** Every
  Python script reports success once the tools return 0; whether the outputs are
  correct requires the checks above.
- **`fetch_rrna_index.sh` never deletes on failure.** A failed download leaves a
  `.partial` file and a failed extraction leaves `.extract-*`, so an interrupted
  run is never mistaken for a complete index. A *successful* archive install does
  delete the downloaded `.tar.gz` to reclaim about 1 GB; pass `--keep-archive` to
  keep it.
