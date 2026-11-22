# Outputs and contract — Database Update Agent

## The layout the other five agents expect

Every consumer computes its path as `<project_dir>/DB/...`, where `project_dir`
comes from `config.yaml`. This agent writes into that same `<project_dir>/DB/`.

```
<project_dir>/DB/
├── manifest.csv                        <- shipped, reinstalled by install_manifest.py
├── clean/                              <- clean_db.py
│   ├── vmr_L_clean.fasta
│   ├── vmr_M_clean.fasta
│   ├── vmr_S_clean.fasta
│   ├── manifest_clean.csv
│   ├── vmr_{L,M,S}_clean_proteins.faa  <- prodigal intermediate
│   └── vmr_{L,M,S}_clean_genes.ffn     <- prodigal intermediate
├── diamond/
│   ├── vmr_{L,M,S}_99.fasta            <- raw material, shipped with the repository
│   ├── vmr_l.fa                        <- build_clean_db.py
│   ├── vmr_l.dmnd
│   ├── vmr_m.fa
│   ├── vmr_m.dmnd
│   ├── vmr_s.fa
│   └── vmr_s.dmnd
├── sortmerna/                          <- fetch_rrna_index.sh
│   ├── rrna_index.1.bt2
│   ├── rrna_index.2.bt2
│   ├── rrna_index.3.bt2
│   ├── rrna_index.4.bt2
│   ├── rrna_index.rev.1.bt2
│   ├── rrna_index.rev.2.bt2
│   └── SHA256SUMS                      <- files mode only, optional
└── tree/                               <- build_tree_ref.py
    ├── vmr_L_ref_prot.fa
    └── vmr_L_aln.fasta
```

Two of these are input to the whole chain and cannot be generated here:
`manifest.csv` and `diamond/vmr_{L,M,S}_99.fasta`. Everything else is produced by
this agent.

## Who reads what

| File | Read by | Purpose |
|---|---|---|
| `diamond/vmr_l.dmnd` | `2-species_class.py`, `4-vmr-tree.py` | L-segment (RdRp) protein alignment, genus and species assignment |
| `diamond/vmr_m.dmnd` | `2-species_class.py` | M-segment (glycoprotein) alignment |
| `diamond/vmr_s.dmnd` | `2-species_class.py` | S-segment (nucleoprotein) alignment |
| `manifest.csv` | `3-vmr-merge-cluster.py`, `4-vmr-tree.py` | acc_base to species mapping, host and location metadata |
| `sortmerna/rrna_index.*.bt2` | `1-qc.sh` | bowtie2 rRNA depletion |
| `tree/vmr_L_aln.fasta` | `4-vmr-tree.py` | Reference backbone for the MAFFT tree |
| `diamond/vmr_{L,M,S}_99.fasta` | `clean_db.py` | Raw nucleotide pool that curation selects from |
| `diamond/vmr_{l,m,s}.fa` | Nobody at run time | Source for rebuilding `.dmnd` |
| `clean/vmr_{L,M,S}_clean.fasta` | `build_clean_db.py` | Prodigal input |
| `clean/manifest_clean.csv` | `install_manifest.py` | Source for the installed manifest |

The consumer side is verified: a `config.yaml` whose `project_dir` points at a
directory with this `DB/` layout makes `1-qc.sh` invoke
`bowtie2 -x <project_dir>/DB/sortmerna/rrna_index` and `2-species_class.py`
invoke `diamond blastx -d <project_dir>/DB/diamond/vmr_l.dmnd`. No path is
hard-coded on their side.

## The format requirement that breaks trees

`tree/vmr_L_aln.fasta` must be amino acids. The tree agent concatenates it with
the query proteins:

```python
tree_input = ref.read_text() + combined_query.read_text()
mafft --auto tree_input
```

Mixing alphabets makes MAFFT abort with `Alphabet 'e' is unknown`. Verify:

```bash
grep -v '^>' <db>/tree/vmr_L_aln.fasta | head -c 200
```

| Observed characters | Verdict |
|---|---|
| `ACGTN` almost exclusively | Nucleotide — wrong, rebuild |
| `E L K D R - X` present | Protein — correct |

## Console output worth reading

`clean_db.py`:

```
input records : 2504
kept records  : 2504
species       : 113 (ICTV 57 + traditional 56)
  vmr_L_clean.fasta: 570 sequences
  vmr_M_clean.fasta: 694 sequences
  vmr_S_clean.fasta: 1240 sequences
```

| Field | Meaning |
|---|---|
| `input records` | Rows in `manifest.csv` |
| `kept records` | Rows that passed curation. Equal to the input for the shipped database, which was curated before release |
| Per-segment `sequences` | Records written to `clean/`, selected by `acc_base` |

When `kept records` is much below `input records`, look at the `dropped ...`
lines printed between the two: each names a rule and counts the records it
removed.

`build_clean_db.py`, one line per segment:

```
  vmr_l.fa: 570 proteins (prodigal covered 570/570 nucleotide records)
  vmr_l.dmnd built
  vmr_m.fa: 694 proteins (prodigal covered 694/694 nucleotide records)
  vmr_m.dmnd built
  vmr_s.fa: 1237 proteins (prodigal covered 1237/1240 nucleotide records)
  note: 3 nucleotide records yielded no protein (short fragments, expected)
  vmr_s.dmnd built
```

| Field | Meaning |
|---|---|
| Protein count | Records written to `vmr_<seg>.fa` |
| `prodigal covered n/m` | Records with at least one predicted protein, out of total input records |
| The note line | Uncovered records are short fragments and are normally expected |

A coverage ratio near 1.0 means full-length references. A ratio near 0.5 means
the library contains many fragments — check that this matches the curation intent
before proceeding. The shipped database sits at 570/570, 694/694 and 1237/1240.

`build_tree_ref.py`:

```
84 species -> 84 representatives
protein library hits: 84
protein length: min 178 / median 2155 / max 2737 aa
running: mafft --thread 8 --auto <db>/tree/vmr_L_ref_prot.fa
output: <db>/tree/vmr_L_aln.fasta
  84 records / alignment length [3461]
  checks passed: protein alphabet, equal-length alignment
```

| Field | Meaning |
|---|---|
| `species -> representatives` | Species count in the manifest versus usable L-segment representatives |
| `protein library hits` | Representatives with a protein in `vmr_l.fa`. Fewer than the previous number means Prodigal predicted no ORF for some |
| Alignment length | Larger than the longest protein above, because gaps inflate the column count |
| `checks passed` | The two self-checks both succeeded |

The alignment length is not stable across MAFFT versions: the same 84 sequences
produced 3,454 columns under one release and 3,461 under another. This does not
matter downstream, because `4-vmr-tree.py` concatenates the backbone with the
query proteins as plain text and re-aligns the whole mixture, so the existing gap
positions are recomputed anyway.

## Manifest, before and after

| Source | Installed | Difference |
|---|---|---|
| `clean/manifest_clean.csv` | `manifest.csv` | None. The two files are byte-identical apart from the name |

`install_manifest.py` validates that `acc_base`, `species` and `segment` are
present before copying, so a malformed manifest fails here rather than during
tree building.

## Downstream impact

```
DB/manifest.csv                           (raw material, shipped with the repo)
DB/diamond/vmr_{L,M,S}_99.fasta           (raw material, shipped with the repo)
                ↓  clean_db.py
DB/clean/vmr_{L,M,S}_clean.fasta
DB/clean/manifest_clean.csv
                ↓  install_manifest.py
DB/manifest.csv                           (installed copy)
                ↓  build_clean_db.py
DB/diamond/vmr_{l,m,s}.fa  +  .dmnd
                ↓  build_tree_ref.py
DB/tree/vmr_L_aln.fasta
                ↓
Agents 1 -> 2 -> 3 -> 4 -> 5   (all downstream results are now stale)
```

After any database rebuild, the following are invalid until the pipeline is
re-run end to end:

- Every per-sample hit table (`06_diamond_vmr/`)
- The shared tree and all classifications (`07_tree/`, `00_results/`)
- Every score, tier and report (`results.tsv`, `species_risk.tsv`, `report.html`)

Do not mix results computed against different database versions in one figure.
The identity values would not be comparable, and the difference would be
attributable to the database rather than to the samples.
