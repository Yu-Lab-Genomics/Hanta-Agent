# Parameters — Database Update Agent

Four scripts. None of them reads a `config.yaml` section of its own; they only
use `project_dir` to derive the database directory.

## Shared path resolution

`build_clean_db.py`, `build_tree_ref.py`, `install_manifest.py` and
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

| Flag | Type | Default | Purpose |
|---|---|---|---|
| `--url` | string | *(none)* | Base URL of the release assets. Required unless `--check` or `--make-checksums` |
| `--db-dir` | path | *(resolved)* | Database root; files land in `<db>/sortmerna/` |
| `-c` / `--config` | path | *(discovered)* | Used to derive `<project_dir>/DB` |
| `--check` | flag | off | Verify what is installed and exit. Writes nothing |
| `--make-checksums` | flag | off | Write `SHA256SUMS` for the files present, then list what to upload |
| `--with-fasta` | flag | off | Also fetch `rrna_merged.fasta` |
| `--force` | flag | off | Re-download even when the index looks intact |
| `--dry-run` | flag | off | Print the URLs that would be fetched |

### Files fetched

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
| `SHA256SUMS` present | Every downloaded file is verified with `sha256sum -c --ignore-missing` |
| `SHA256SUMS` absent | Installation proceeds and prints that integrity was not verified |
| A checksum fails | Exit 1 with the failing file named |
| A download fails | That file is counted, the rest continue, exit 1 at the end |

Partially downloaded files are kept as `<name>.partial` rather than deleted, so a
failure can be inspected or resumed.

### Uploading the assets

```bash
bash fetch_rrna_index.sh --db-dir <path> --make-checksums
```

This writes `<db>/sortmerna/SHA256SUMS` and prints the exact file list to upload.
Attach those seven files to the release, then pass the release's download base URL
to `--url`.

## GitHub release size limits

The two limits are often confused:

| Limit | Value | Applies to |
|---|---|---|
| Repository file | 100 MB | Files added with `git add` |
| Release asset | 2048 MB | Files attached to a release |

The largest file here is 1028 MB, comfortably inside the release limit. None of
them may enter the repository itself, so `DB/sortmerna/` must be gitignored.
