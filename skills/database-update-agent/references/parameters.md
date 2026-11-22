# Parameters — Database Update Agent

Five scripts. None of them reads a `config.yaml` section of its own; they only
use `project_dir` to derive the database directory.

## Shared path resolution

`clean_db.py`, `build_clean_db.py`, `build_tree_ref.py`, `install_manifest.py` and
`fetch_rrna_index.sh` all resolve the database directory the same three ways, in
this priority order:

1. `--db-dir <path>` — used verbatim.
2. `--config <path>` / `-c <path>` — read `project_dir`, use `<project_dir>/DB`.
3. Automatic discovery — walk up from the script, at most four levels, looking
   for `config.yaml`, `config-cell-paper.yaml` or `cell_submit-paper.yaml`, then
   apply rule 2.

If none of these yields a directory, the script exits 2 and says so.

| Setting | Where | Effect |
|---|---|---|
| `project_dir` | top level of `config.yaml` | The one key this agent reads. `<project_dir>/DB` becomes the database root |
| `threads` | top level | **Not read here.** `build_tree_ref.py` uses `--threads` instead |

## clean_db.py

Applies the curation rules and writes the library the next two steps consume.

| Flag | Type | Default | Purpose |
|---|---|---|---|
| `--db-dir` | path | *(resolved)* | Database root |
| `--config` | path | *(discovered)* | Used only to derive `<project_dir>/DB` |
| `--dry-run` | flag | off | Report the curation outcome and stop. Writes nothing |

### Reads

| Path | Content |
|---|---|
| `<db>/manifest.csv` | Full metadata table, 2,504 rows |
| `<db>/diamond/vmr_L_99.fasta` | Dereplicated L-segment nucleotides |
| `<db>/diamond/vmr_M_99.fasta` | Dereplicated M-segment nucleotides |
| `<db>/diamond/vmr_S_99.fasta` | Dereplicated S-segment nucleotides |

All four are required. A missing manifest exits 2 before any work starts. A
missing nucleotide file is reported *after* the manifest has been read, so the
`kept records` count still prints first — useful, because it tells you how much
would survive regardless of which segment file is absent.

### Writes

| Path | Content |
|---|---|
| `<db>/clean/manifest_clean.csv` | Curated metadata, same columns as the input |
| `<db>/clean/vmr_{L,M,S}_clean.fasta` | Only the records that survived curation |

### Curation rules

A record is dropped when any of these holds:

| Rule | Reported as |
|---|---|
| Species name contains `sp.`, or is a bare genus name | `placeholder species` |
| Species name is in the hand-maintained vague list | `vague name` |
| Neither marked complete by NCBI nor above the segment length floor | `neither complete nor long enough` |

The bare genus names are `Orthohantavirus`, `Mobatvirus`, `Thottimvirus`,
`Loanvirus`, `Actinovirus`, `Percilovirus`, `Agnathovirus` and `Reptillovirus`.
The vague list currently holds one entry, `Skink hantavirus`, which names no
place and therefore cannot be tied to a location.

"Marked complete" means the `description` or `title` column contains
`complete genome`, `complete cds` or `complete sequence`, case-insensitively.

| Segment | Length floor |
|---|---|
| L | 5,000 nt |
| M | 3,000 nt |
| S | 1,400 nt |

A record that NCBI calls complete is kept regardless of length, so there is no
absolute lower bound beyond the floor.

### Expected output

```
input records : 2504
kept records  : 2504
species       : 113 (ICTV 57 + traditional 56)
  vmr_L_clean.fasta: 570 sequences
  vmr_M_clean.fasta: 694 sequences
  vmr_S_clean.fasta: 1240 sequences
```

Two numbers deserve an explanation, because both look wrong at first glance.

**`input records` equals `kept records`.** The shipped `manifest.csv` was already
curated before release, so no rule fires and this step changes nothing. That is
the expected result. The step still matters for two reasons: it regenerates
`clean/` from scratch if those files are lost, and it is the place where a newly
downloaded manifest would actually be filtered.

**The sequence counts are far below the input file counts.** The nucleotide
inputs hold more records than the manifest does:

| Segment | Records in `vmr_<seg>_99.fasta` | Of those, present in `manifest.csv` | Written to `clean/` |
|---|---|---|---|
| L | 1,895 | 570 | 570 |
| M | 1,682 | 694 | 694 |
| S | 2,438 | 1,240 | 1,240 |

The filter is by `acc_base`: a sequence is written only if its accession appears
in the curated manifest. CD-HIT at 99% keeps every near-identical variant, while
the manifest lists one entry per selected record, so most CD-HIT output never
survives. This is by design — the manifest decides what the database contains,
and `vmr_{L,M,S}_99.fasta` is just the pool it draws from.

## build_clean_db.py

| Flag | Type | Default | Purpose |
|---|---|---|---|
| `--db-dir` | path | *(resolved)* | Database root. Inputs come from `<db>/clean/`, outputs go to `<db>/diamond/` |
| `--config` | path | *(discovered)* | Used only to derive `<project_dir>/DB` |
| `--prodigal` | string | `prodigal` | Prodigal executable |
| `--diamond` | string | `diamond` | DIAMOND executable |
| `--dry-run` | flag | off | Print every composed command and stop. Nothing is created or written |

### Reads

| Path | Content |
|---|---|
| `<db>/clean/vmr_L_clean.fasta` | Cleaned L-segment nucleotides |
| `<db>/clean/vmr_M_clean.fasta` | Cleaned M-segment nucleotides |
| `<db>/clean/vmr_S_clean.fasta` | Cleaned S-segment nucleotides |

All three are required. A missing one exits 2 with its path.

### Writes

| Path | Content |
|---|---|
| `<db>/diamond/vmr_{l,m,s}.fa` | Longest protein per nucleotide record |
| `<db>/diamond/vmr_{l,m,s}.dmnd` | DIAMOND index |
| `<db>/clean/vmr_{L,M,S}_clean_proteins.faa` | All Prodigal proteins, before the filter |
| `<db>/clean/vmr_{L,M,S}_clean_genes.ffn` | Predicted CDS nucleotides |

The `.faa` file is deliberately kept: comparing it against `vmr_<seg>.fa` is how
you confirm the longest-protein filter did what you expect.

### Prodigal invocation

```bash
prodigal -i <seg>_clean.fasta -a <seg>_clean_proteins.faa \
         -d <seg>_clean_genes.ffn -o /dev/null -p meta -q
```

| Flag | Meaning |
|---|---|
| `-a` | Protein FASTA output |
| `-d` | Nucleotide CDS output |
| `-o /dev/null` | Discard the GenBank-style summary; it is large and unused |
| `-p meta` | Metagenomic mode — no assumption that the sequence starts at a gene boundary |
| `-q` | Quiet |

`-p meta` matters: reference sequences may be partial, and single mode would
mispredict the start.

### DIAMOND invocation

```bash
diamond makedb --in <db>/diamond/vmr_<seg>.fa -d <db>/diamond/vmr_<seg>
```

No parameters are exposed. Default sensitivity and block size are adequate for a
few thousand reference proteins.

## build_tree_ref.py

| Flag | Type | Default | Purpose |
|---|---|---|---|
| `--db-dir` | path | *(resolved)* | Database root |
| `--config` | path | *(discovered)* | Used only to derive `<project_dir>/DB` |
| `--out` | path | `<db>/tree/vmr_L_aln.fasta` | Alignment output path |
| `--mafft` | string | `mafft` | MAFFT executable |
| `--threads` | int | `8` | Passed to MAFFT |
| `--mafft-extra` | string | `--auto` | MAFFT arguments. Advanced — leave as is |
| `--dry-run` | flag | off | Print the composed command without running MAFFT |

### Reads

| Path | Content |
|---|---|
| `<db>/manifest.csv` | Taxonomy and length metadata, used to pick one representative per species |
| `<db>/diamond/vmr_l.fa` | L-segment proteins from `build_clean_db.py` |

Both required. A missing one exits 2 with its path and names the script to run
first.

### Writes

| Path | Content |
|---|---|
| `<db>/tree/vmr_L_aln.fasta` | The protein alignment the tree agent uses as backbone |
| `<db>/tree/vmr_L_ref_prot.fa` | The unaligned proteins, kept as a record and used as MAFFT input |

### Representative selection

Per species, in priority order:

1. Complete genome (`is_complete()`), then
2. Longest sequence (`length_nt`)

One tip per species. More than one makes near-identical sequences crowd the tree
and pulls support values down across the whole topology.

`is_complete()` checks the `description` column, plus `title` when present.
`manifest.csv` has no `title` column, and the lookup already tolerates its
absence, so nothing breaks.

### MAFFT invocation

```bash
mafft --thread <threads> --auto <db>/tree/vmr_L_ref_prot.fa > <out>
```

`--auto` selects the algorithm from the input size. Do not replace it with a
fixed strategy unless the alignment size is stable and known.

## install_manifest.py

| Flag | Type | Default | Purpose |
|---|---|---|---|
| `--db-dir` | path | *(resolved)* | Database root |
| `--config` | path | *(discovered)* | Used only to derive `<project_dir>/DB` |
| `--source` | path | `<db>/clean/manifest_clean.csv` | Explicit source file |
| `--force` | flag | off | Replace an existing `manifest.csv` |
| `--check` | flag | off | Validate only; write nothing |

| Behaviour | Result |
|---|---|
| Source present, target absent | Copies and reports rows and columns |
| Source present, target present, no `--force` | Reports `already installed` and exits 0 |
| Source present, target present, `--force` | Replaces |
| Source absent, target present | Reports `source not present` and exits 0 |
| Source absent, target absent | Exits 2 |
| Source or target missing `acc_base`, `species` or `segment` | Exits 2 and names the column |

Idempotent by design, so a setup script can call it unconditionally.

## fetch_rrna_index.sh

Fetches the bowtie2 index. Two modes, chosen by the shape of the URL:

| Mode | Trigger | Behaviour |
|---|---|---|
| archive | URL ends in `.tar.gz` or `.tgz` | One download, then the six `.bt2` files are extracted into `<db>/sortmerna/` |
| files | any other URL | Six separate downloads against a base URL |

Archive mode is the default and is what a fresh clone uses. Files mode is kept
because it lets you replace one corrupt `.bt2` without pulling the whole 1.2 GB
again.

| Flag | Type | Default | Purpose |
|---|---|---|---|
| `--url` | string | the bundled release archive | Source to fetch. `.tar.gz` or `.tgz` selects archive mode, anything else files mode |
| `--db-dir` | path | *(resolved)* | Database root; files land in `<db>/sortmerna/` |
| `-c` / `--config` | path | *(discovered)* | Used to derive `<project_dir>/DB` |
| `--sha256` | hex | *(empty)* | Expected SHA256 of the archive. Archive mode only; overrides the built-in value |
| `--keep-archive` | flag | off | Keep the `.tar.gz` after extraction. Archive mode only |
| `--check` | flag | off | Verify what is installed and exit. Writes nothing |
| `--make-checksums` | flag | off | Write `SHA256SUMS` for the installed files and print the SHA256 of any archive present |
| `--with-fasta` | flag | off | Files mode only: also fetch `rrna_merged.fasta` |
| `--force` | flag | off | Re-download even when the index looks intact |
| `--dry-run` | flag | off | Print the plan and exit |

### Archive mode, step by step

1. If all six `.bt2` files are present and `--force` is absent, exit 0 without
   touching the network.
2. Download to `<db>/sortmerna/<name>.partial`, then rename. An archive already
   on disk is reused unless `--force` is given.
3. Verify the archive against `--sha256`, or against the value built into the
   script. With neither set, the check is skipped and a warning is printed.
4. Warn if free space under `<db>/sortmerna/` is below twice the archive size.
5. Extract into `<db>/sortmerna/.extract-<random>/`.
6. Locate each of the six files by name anywhere in the extracted tree, so an
   archive that wraps everything in a directory works unchanged.
7. Move the six files into `<db>/sortmerna/` and delete the scratch directory.
8. Delete the archive unless `--keep-archive` was given.
9. Run the same completeness check that `--check` performs.

A failure at step 2 leaves `<name>.partial` behind; a failure at step 5 or 6
leaves the `.extract-*` directory. Neither is cleaned up, so the cause can be
inspected before retrying.

### Where the URL and checksum come from

Both values are resolved with the same precedence:

```
command line (--url / --sha256)  >  config.yaml  >  the constants in the script
```

The constants at the top of the script are:

| Constant | Ships as |
|---|---|
| `DEFAULT_URL` | `.../releases/download/sortmerna-v1.0.0/rrna_index.tar.gz` |
| `BUILTIN_ARCHIVE_SHA256` | empty, which means the download is not verified |

`config.yaml` supplies them through `bowtie2.fetch_url` and
`bowtie2.fetch_sha256`. `1-qc.sh` reads the same two keys, so one edit covers
both the automatic preflight download and this script.

After uploading a new archive, obtain the hash with `--make-checksums` and paste
it into `config.yaml`, into `BUILTIN_ARCHIVE_SHA256`, or pass it per run with
`--sha256`.

### Files mode

| File | Approximate size |
|---|---|
| `rrna_index.1.bt2` | 353 MB |
| `rrna_index.2.bt2` | 181 MB |
| `rrna_index.3.bt2` | 7.6 MB |
| `rrna_index.4.bt2` | 181 MB |
| `rrna_index.rev.1.bt2` | 353 MB |
| `rrna_index.rev.2.bt2` | 181 MB |
| **total** | **≈ 1.26 GB** |
| `rrna_merged.fasta` | 1028 MB, only with `--with-fasta` |

`rrna_merged.fasta` is **not** needed at run time. `1-qc.sh` invokes bowtie2 with
`-x <index-prefix>`, which reads the `.bt2` files and never opens the FASTA. Fetch
it only if you plan to rebuild the index later.

### Integrity

| Behaviour | Result |
|---|---|
| `--sha256` or `ARCHIVE_SHA256` set, archive mode | The archive is hashed and compared; a mismatch exits 1 |
| Neither set, archive mode | Installation proceeds and prints that integrity was not verified |
| `SHA256SUMS` present, files mode | Every downloaded file is checked with `sha256sum -c --ignore-missing` |
| `SHA256SUMS` absent, files mode | Installation proceeds and prints that integrity was not verified |
| A checksum fails | Exit 1 with the failing file or archive named |
| A download fails, files mode | That file is counted, the rest continue, exit 1 at the end |

### Uploading the assets

For the archive layout a single file is enough:

```bash
bash fetch_rrna_index.sh --db-dir <path> --make-checksums   # prints the archive SHA256
```

Upload `rrna_index.tar.gz`, then paste the printed hash into `ARCHIVE_SHA256`.
For the per-file layout, `--make-checksums` writes `<db>/sortmerna/SHA256SUMS`
and prints the file list to upload alongside it.

## GitHub release size limits

The two limits are often confused:

| Limit | Value | Applies to |
|---|---|---|
| Repository file | 100 MB | Files added with `git add` |
| Release asset | 2048 MB | Files attached to a release |

The largest file here is 1028 MB, comfortably inside the release limit. None of
them may enter the repository itself, so `DB/sortmerna/` must be gitignored.
