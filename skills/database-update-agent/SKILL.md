---
name: database-update-agent
description: Rebuild and install the reference database that the other five agents read — fetch the rRNA index, promote the manifest, translate the cleaned nucleotide library with Prodigal, build the DIAMOND indexes, and rebuild the protein alignment used as tree backbone. Agent 6 of the six-agent Hanta-Agent workflow. Run only when setting up a new machine or when the reference library changes.
---

# Database Update Agent

Put a working database under `<project_dir>/DB/` and make sure the other agents
find it there.

This agent has five scripts. They are independent, so run only the ones you need:

| Script | Purpose | Run when |
|---|---|---|
| `clean_db.py` | Curate `DB/manifest.csv` + `diamond/vmr_{L,M,S}_99.fasta` into `DB/clean/` | First setup, or the raw library changed |
| `install_manifest.py` | Promote `clean/manifest_clean.csv` to `DB/manifest.csv` | After `clean_db.py` |
| `build_clean_db.py` | Prodigal + longest protein + `diamond makedb` | After `install_manifest.py` |
| `build_tree_ref.py` | One representative per species, MAFFT protein alignment | After `build_clean_db.py` |
| `fetch_rrna_index.sh` | Download the rRNA index: one release archive, or six files | First setup, or the index is missing |

## The one thing to get right

Every other agent computes its database path as `<project_dir>/DB/...`, where
`project_dir` comes from `config.yaml`. This agent writes into that same
`<project_dir>/DB/`. If `project_dir` points somewhere else than the directory
holding the database, the other agents will not see anything this agent built.

The scripts therefore take `--db-dir` (explicit) or `--config` (derive
`<project_dir>/DB`). When neither is given, they walk up from their own location
looking for `config.yaml`, `config-cell-paper.yaml` or `cell_submit-paper.yaml`.

## Layout this agent produces

```
<project_dir>/DB/
├── manifest.csv                 install_manifest.py
├── clean/                       clean_db.py
│   ├── vmr_{L,M,S}_clean.fasta
│   ├── manifest_clean.csv
│   └── *_clean_{proteins.faa,genes.ffn}   build_clean_db.py, inspection only
├── diamond/                     build_clean_db.py
│   ├── vmr_{L,M,S}_99.fasta     raw input, CD-HIT output
│   ├── vmr_{l,m,s}.fa
│   └── vmr_{l,m,s}.dmnd
├── sortmerna/                   fetch_rrna_index.sh
│   ├── rrna_index.*.bt2
│   └── SHA256SUMS
└── tree/                        build_tree_ref.py
    ├── vmr_L_ref_prot.fa
    └── vmr_L_aln.fasta
```

Three tiers are worth keeping apart. `DB/manifest.csv` and
`diamond/vmr_{L,M,S}_99.fasta` are **raw material**; `clean/` is an
**intermediate**; `diamond/*.fa`, `diamond/*.dmnd` and `tree/vmr_L_aln.fasta`
are the **finished products** the other agents open. The released repository
ships all three tiers, so a fresh clone can skip straight to running samples.

## When to use

- Setting up the pipeline on a new machine.
- The reference library gained or lost sequences.
- The taxonomic manifest was updated and the tree backbone must follow.
- `4-vmr-tree.py` fails with `Alphabet 'e' is unknown`, which means the backbone
  alignment is nucleotide and must be rebuilt.

## When not to use

- Nothing in the reference library changed and the database is already in place.
  `build_clean_db.py` and `build_tree_ref.py` are slow, and rebuilding
  invalidates every downstream comparison.
- Only sample data changed. Re-run the sample pipeline instead.
- `DB/manifest.csv` or `DB/diamond/vmr_{L,M,S}_99.fasta` is absent. Those two
  are the only things this agent cannot create: they come from NCBI. Everything
  downstream of them, including `DB/clean/`, is produced here.

## Contract

| Item | Value |
|---|---|
| **Input** | `DB/manifest.csv` — taxonomy and accession map, the single required metadata file |
| **Input** | `DB/diamond/vmr_{L,M,S}_99.fasta` — raw nucleotide library, one file per segment |
| **Input** | `DB/diamond/vmr_l.fa` — required by `build_tree_ref.py`, written by `build_clean_db.py` |
| **Config** | `--config <path>` or `--db-dir <path>`; auto-discovered upward from the script when omitted |
| **Output** | `DB/clean/vmr_{L,M,S}_clean.fasta` and `DB/clean/manifest_clean.csv` |
| **Output** | `DB/manifest.csv` |
| **Output** | `DB/diamond/vmr_{l,m,s}.fa` and `.dmnd` |
| **Output** | `DB/sortmerna/rrna_index.*.bt2` |
| **Output** | `DB/tree/vmr_L_aln.fasta` |
| **Downstream interface** | `DB/diamond/*.dmnd` (agents 2, 4), `DB/manifest.csv` (agents 3, 4), `DB/tree/vmr_L_aln.fasta` (agent 4), `DB/sortmerna/rrna_index` (agent 1) |
| **Failure behaviour** | Exit 2 for a missing or unresolvable database directory, a missing input, or a manifest lacking required columns; tool failures propagate; `fetch_rrna_index.sh` exits 1 on a download or checksum failure |
| **Idempotency** | Every script may be re-run. `clean_db.py` overwrites `DB/clean/`; `install_manifest.py` needs `--force` to replace an existing manifest; `fetch_rrna_index.sh` re-downloads only with `--force` |

## Workflow

1. **Resolve the database directory.** Confirm `project_dir` in `config.yaml`
   points at the directory whose `DB/` subdirectory holds (or will hold) the
   database. All five scripts honour this.

2. **Fetch the rRNA index.** `bash fetch_rrna_index.sh --db-dir <DB>`. About
   1.2 GB across six `.bt2` files. `--url` defaults to the published archive, so
   nothing else is needed; the six files arrive in one download and are
extracted in place. Verify afterwards with `--check`. To prepare the assets
yourself, run `--make-checksums` first and upload what it lists.

3. **Curate the library.** `python clean_db.py`. It reads `DB/manifest.csv` and
   `DB/diamond/vmr_{L,M,S}_99.fasta`, drops records whose species name is too
   vague to classify and records shorter than the per-segment length floor, and
   writes `DB/clean/`. It prints `input records` and `kept records`; if those
   two numbers are equal the manifest was already curated and the step is a
   no-op, which is the normal case for the released database.

4. **Install the manifest.** `python install_manifest.py --force`. Idempotent:
   without `--force` it reports an existing `DB/manifest.csv` and exits 0. Use
   `--check` to validate columns without writing.

5. **Build the protein databases.** `python build_clean_db.py`. It writes
   `.fa` and `.dmnd` straight into `DB/diamond/`, which is where agents 2 and 4
   read them. Use `--dry-run` first to see the composed commands.

6. **Check the translation coverage.** The script prints
   `vmr_l.fa: N proteins (prodigal covered n/m nucleotide records)`. Records
   without a predicted protein are short fragments and are expected to be
   skipped; the printed hint says so. Coverage far below expectation means the
   input is not what you think it is.

7. **Build the tree backbone.** `python build_tree_ref.py`. It reads
   `DB/diamond/vmr_l.fa` and `DB/manifest.csv`, and writes
   `DB/tree/vmr_L_aln.fasta`.

8. **Verify the backbone is protein** before running any tree:
   `grep -v '^>' DB/tree/vmr_L_aln.fasta | head -c 200` must show amino-acid
   characters, not just `ACGTN`.

9. **Re-run the sample pipeline** from the tree agent onward. Every score
   downstream of the tree is now stale.

## Invariants

1. **Outputs land where the consumers read them.** The DIAMOND indexes go into
   `DB/diamond/` and the backbone into `DB/tree/`, not into `DB/clean/`. Writing
   them anywhere else means someone has to copy files by hand, and a forgotten
   copy is indistinguishable from a failed build.
2. **The tree backbone must be protein.** `4-vmr-tree.py` concatenates the
   backbone directly with query proteins and hands the mixture to MAFFT. A
   nucleotide backbone makes MAFFT reject the whole file with
   `Alphabet 'e' is unknown`. This is the single reason `build_tree_ref.py`
   exists in its current form: an earlier version produced an unaligned
   nucleotide file and every tree run failed.
3. **One representative per species.** `build_tree_ref.py` picks, per species,
   complete genomes first and longest sequence second. More than one
   representative per species makes near-identical tips crowd the tree and
   collapses support values.
4. **`DB/clean/` is reproducible from the two raw inputs.** `clean_db.py`
   rebuilds it from `DB/manifest.csv` plus `diamond/vmr_{L,M,S}_99.fasta` in a
   few seconds. It is therefore safe to delete when space is tight, but do not
   delete the two inputs: nothing in this agent can regenerate those.
5. **Longest protein per record, not all proteins.** Prodigal emits several
   proteins per nucleotide record; the script keeps the longest. Keeping all of
   them would make one record dominate the DIAMOND index.
6. **`install_manifest.py` is idempotent; `fetch_rrna_index.sh` checks integrity.**
   Both are safe to call from a setup script that may run twice. Re-downloading
   needs `--force`; replacing the manifest needs `--force`.
7. **Rebuilding invalidates prior results.** After a database rebuild, every
   classification, score and tier from earlier runs is no longer comparable.

## Decision points

| Situation | Adjustment | Risk |
|---|---|---|
| Database lives outside the project directory | Pass `--db-dir` to each script | The other agents still compute `<project_dir>/DB`, so they will not see it — you would also need a symlink |
| A species is over-represented in the backbone | Tighten representative selection; the script already prefers complete genomes | More representatives crowd the tree and depress support |
| Translation coverage is low | Check that the input is nucleotide and complete, not a fragment set | Low coverage is often legitimate for short contigs |
| MAFFT is slow | Raise `--threads`, or pass `--mafft-extra` with a faster preset | `--auto` already adapts; forcing a preset can degrade the alignment |
| Release assets are unavailable | Build the index from SILVA with `bowtie2-build` into `DB/sortmerna/rrna_index` | Takes 20-40 minutes and needs the SILVA download |
| Different segment needed | `build_clean_db.py` already does all three; `build_tree_ref.py` is L-only by design | The tree is built on L (the RdRp), the standard hantavirus marker |

Full parameter reference: `references/parameters.md`.

## Failure handling

| Exit code | Meaning | First action |
|---|---|---|
| 2 | Database directory unresolved, input missing, or manifest columns absent | Read stderr; it names the path |
| 1 | (`fetch_rrna_index.sh` only) download, checksum or extraction failure | Read stderr: it names the `.partial` archive or the `.extract-*` directory left behind |
| other | Prodigal, DIAMOND or MAFFT failed | Read the tool output; the scripts print the composed command |

Per-case diagnosis: `references/troubleshooting.md`.

## References

- `references/parameters.md` — every flag and every path each script touches.
- `references/outputs-and-contract.md` — output files, formats and consumers.
- `references/troubleshooting.md` — failure branches, including the alphabet error.

## Provenance

`scripts/` holds copies of `../../build_clean_db.py`, `../../build_tree_ref.py`,
`../../install_manifest.py`, `../../clean_db.py` and `../../fetch_rrna_index.sh`.

On 2026-11-20 the first two were changed so that they no longer hard-code the
project root and no longer write into `DB/clean/`. Before that change,
`build_clean_db.py` wrote its DIAMOND indexes to `DB/clean/vmr_*.dmnd` while
agents 2 and 4 read `DB/diamond/vmr_*.dmnd` — a silent mismatch that made a
successful rebuild have no effect. `install_manifest.py` and
`fetch_rrna_index.sh` are new in the same change.

On 2026-11-22 `clean_db.py` was added and adapted to `--db-dir`. Until then the
curation step ran outside the agent, which meant the documented workflow started
from a directory the reader had no instructions to produce. The whole chain
`manifest.csv` → `clean/` → `diamond/` → `tree/` was then run end to end against
the released inputs; the resulting `vmr_{l,m,s}.fa` and `.dmnd` files matched the
published database byte for byte, and the tree backbone matched in sequence
count (84) with a 7-column alignment-width difference attributable to the MAFFT
version.

### Bundled helper modules

`scripts/seqid.py` (110 lines) is copied into this directory because
`build_clean_db.py` imports `prodigal_orig_id()`, `read_fasta_ids()` and
`strip_prodigal_suffix()` from it. It reads no external data. The scripts resolve
their own directory on `sys.path`, so the module must sit beside them. The same
module is bundled with three other agents; the duplication is intentional so each
skill works standalone.
