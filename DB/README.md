# Reference Database

Everything the pipeline reads at run time lives here.

Two parts:

| Part | Size | How you get it |
|---|---|---|
| Metadata, protein libraries, DIAMOND indexes, tree backbone | 25.7 MB | ships with this repository |
| rRNA index (6 bowtie2 files) | 1.2 GB | downloaded from the release page |

The reference set itself is **2,504 nucleotide records from 113 hantavirus species**,
built from the whole family *Hantaviridae* (`txid1980413`) rather than only the ICTV
type sequences, so that recently described viruses are represented too.

---

## 1. Layout

```
DB/
├── README.md                      this file
├── manifest.csv                   2,504 records x 20 columns, 113 species
├── diamond/
│   ├── vmr_l.fa                   570 proteins  (L segment, RdRp)
│   ├── vmr_l.dmnd                 DIAMOND index built from the file above
│   ├── vmr_m.fa                   694 proteins  (M segment, glycoprotein)
│   ├── vmr_m.dmnd
│   ├── vmr_s.fa                   1,237 proteins (S segment, nucleoprotein)
│   ├── vmr_s.dmnd
│   ├── vmr_L_99.fasta             1,895 raw nucleotides, CD-HIT 99% (L)
│   ├── vmr_M_99.fasta             1,682 raw nucleotides (M)
│   └── vmr_S_99.fasta             2,438 raw nucleotides (S)
├── clean/
│   ├── manifest_clean.csv         curated manifest, 2,504 rows
│   ├── vmr_L_clean.fasta          570 curated nucleotides
│   ├── vmr_M_clean.fasta          694
│   └── vmr_S_clean.fasta          1,240
├── sortmerna/                     DOWNLOADED, not in git
│   ├── rrna_index.1.bt2
│   ├── rrna_index.2.bt2
│   ├── rrna_index.3.bt2
│   ├── rrna_index.4.bt2
│   ├── rrna_index.rev.1.bt2
│   └── rrna_index.rev.2.bt2
└── tree/
    └── vmr_L_aln.fasta            84 sequences, one representative per species
```

Three tiers, and it helps to know which is which:

- **Raw material** — `diamond/vmr_{L,M,S}_99.fasta`. CD-HIT output straight from the
  download. Nothing reads these at run time; they exist so `clean/` can be rebuilt.
- **Intermediates** — `clean/`. Produced by `clean_db.py`, consumed by
  `build_clean_db.py`. Keeping them in the repository lets a fresh clone skip the
  curation step.
- **Finished products** — `diamond/vmr_{l,m,s}.fa`, `diamond/*.dmnd`,
  `tree/vmr_L_aln.fasta`. These are what the pipeline actually opens.

The L-segment count differs on purpose between tiers: 1,895 raw nucleotides become 570
after curation, because most records are fragments of species already represented.

---

## 2. Who reads what

| File | Read by | Purpose |
|---|---|---|
| `diamond/vmr_l.dmnd` | `2-species_class.py`, `4-vmr-tree.py` | L-segment protein search, genus and species assignment |
| `diamond/vmr_m.dmnd` | `2-species_class.py` | M-segment protein search |
| `diamond/vmr_s.dmnd` | `2-species_class.py` | S-segment protein search |
| `manifest.csv` | `2-species_class.py`, `3-vmr-merge-cluster.py`, `4-vmr-tree.py` | accession to species map, host, country, collection date |
| `sortmerna/rrna_index.*.bt2` | `1-qc.sh` | bowtie2 rRNA read removal |
| `tree/vmr_L_aln.fasta` | `4-vmr-tree.py` | reference backbone for the phylogenetic tree |
| `diamond/vmr_{l,m,s}.fa` | nobody at run time | source for rebuilding `.dmnd` |
| `clean/*` | `build_clean_db.py` | input for the protein extraction step |

`1-qc.sh` uses **bowtie2**, not SortMeRNA, even though the directory is called
`sortmerna/`. The name is historical. Only the six `.bt2` files are opened; the
960 MB `rrna_merged.fasta` the index was built from is not needed and is not shipped.

---

## 3. Fetch the rRNA index

**Normally you do not need to do this.** `1-qc.sh` checks for the index during
preflight and downloads it before processing any sample, so a fresh clone runs
as-is. Fetch it ahead of time when compute nodes have no internet access, or when
you want to confirm the download on its own:

```bash
bash skills/database-update-agent/scripts/fetch_rrna_index.sh --db-dir DB
```

That fetches a single archive, `rrna_index.tar.gz`, extracts the six `.bt2` files
into `DB/sortmerna/`, and deletes the archive to reclaim the space.

Useful flags:

```bash
--check           # report what is present and correct, download nothing
--dry-run         # show the plan
--force           # re-download even if the files already look right
--keep-archive    # keep rrna_index.tar.gz after extraction
--sha256 <hex>    # verify the archive against a known hash
```

To fetch the six files separately instead, pass a base URL that does not end in
`.tar.gz`; the script switches to per-file mode and honours an optional
`SHA256SUMS` from the same location.

To turn the automatic download off, set `bowtie2.auto_fetch: false` in
`skills/config.yaml`. Agent 1 then stops with the manual command printed instead.

---

## 4. Rebuild from raw material

The four steps below regenerate every file in `clean/`, `diamond/` and `tree/` from
`manifest.csv` plus `diamond/vmr_{L,M,S}_99.fasta`. They are also what
`database-update-agent` runs.

```bash
DBDIR=/path/to/Hanta-Agent/DB
S=skills/database-update-agent/scripts

# 1. Curate: drop vague names and fragments below the length floor
python "$S/clean_db.py" --db-dir "$DBDIR"

# 2. Promote the curated manifest to the one the pipeline reads
python "$S/install_manifest.py" --db-dir "$DBDIR" --force

# 3. Predict proteins with Prodigal, then index them with DIAMOND
python "$S/build_clean_db.py" --db-dir "$DBDIR"

# 4. Pick one representative per species and align them with MAFFT
python "$S/build_tree_ref.py" --db-dir "$DBDIR"
```

Each script accepts `--dry-run`. Add `--prodigal`, `--diamond`, `--mafft` if those
binaries are not on `PATH`. To refresh the source data instead of rebuilding from it,
see `skills/database-update-agent/SKILL.md`.

---

## 5. Verification counts

If a rebuild produces different numbers, something upstream changed. These are the
values for the released database.

| Item | Value |
|---|---|
| `manifest.csv` data rows | 2,504 |
| `manifest.csv` columns | 20 |
| Species | 113 (57 ICTV + 56 additional) |
| L-segment proteins | 570 / 1,197,088 aa |
| M-segment proteins | 694 / 777,743 aa |
| S-segment proteins | 1,237 / 579,584 aa |
| Tree backbone sequences | 84 |

Per-segment totals agree with `diamond dbinfo`. Note that the S segment has 1,240
curated nucleotides but only 1,237 proteins: three short records produce no open
reading frame, which is expected.

---

## 6. Not tracked in this repository

Excluded by `.gitignore` for size or because they are reproducible.

| Path | Size | Why |
|---|---|---|
| `sortmerna/*.bt2` | 1.2 GB | release asset, fetch with the script above |
| `sortmerna/*.tar.gz` | ~1.2 GB | the release archive itself, deleted after a successful extraction |
| `sortmerna/rrna_merged.fasta` | 960 MB | only needed to rebuild the bowtie2 index, never read at run time |
| `clean/*_clean_proteins.faa` | 2.8 MB | Prodigal output, regenerated by step 3 |
| `clean/*_clean_genes.ffn` | 7.8 MB | same |

---

## 7. Provenance

- **Sequences**: NCBI Nucleotide, `txid1980413[Organism:exp]`, which recurses into
  unclassified nodes so that unassigned hantavirus records are included.
- **Species assignment**: ICTV VMR plus manual curation for taxa not yet in VMR.
- **Redundancy removal**: CD-HIT-EST at 99% identity per segment.
- **Curation floor**: L 5,000 nt, M 3,000 nt, S 1,400 nt; records described only as
  "hantavirus" with no species attribution are dropped.
- **Date**: curated 2026-10-22, released database assembled 2026-10-27.

Filtering the reference set for similarity is deliberately *not* part of the pipeline.
This is a hantavirus-specific database, so any DIAMOND hit counts as a related
sequence; a low identity score means the query is a fragment, not that the hit is
spurious.
