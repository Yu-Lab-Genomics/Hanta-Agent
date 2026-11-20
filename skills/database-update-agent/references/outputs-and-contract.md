# Outputs and contract — Database Update Agent

## The layout the other five agents expect

Every consumer computes its path as `<project_dir>/DB/...`, where `project_dir`
comes from `config.yaml`. This agent writes into that same `<project_dir>/DB/`.

```
<project_dir>/DB/
├── manifest.csv                        <- install_manifest.py
├── clean/                              <- INPUT, supplied externally
│   ├── vmr_L_clean.fasta
│   ├── vmr_M_clean.fasta
│   ├── vmr_S_clean.fasta
│   ├── manifest_clean.csv
│   ├── vmr_{L,M,S}_clean_proteins.faa  <- prodigal intermediate
│   └── vmr_{L,M,S}_clean_genes.ffn     <- prodigal intermediate
├── diamond/                            <- build_clean_db.py
│   ├── vmr_l.fa
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
│   └── SHA256SUMS
└── tree/                               <- build_tree_ref.py
    ├── vmr_L_ref_prot.fa
    └── vmr_L_aln.fasta
```

## Who reads what

| File | Read by | Purpose |
|---|---|---|
| `diamond/vmr_l.dmnd` | `2-species_class.py`, `4-vmr-tree.py` | L-segment (RdRp) protein alignment, genus and species assignment |
| `diamond/vmr_m.dmnd` | `2-species_class.py` | M-segment (glycoprotein) alignment |
| `diamond/vmr_s.dmnd` | `2-species_class.py` | S-segment (nucleoprotein) alignment |
| `manifest.csv` | `3-vmr-merge-cluster.py`, `4-vmr-tree.py` | acc_base to species mapping, host and location metadata |
| `sortmerna/rrna_index.*.bt2` | `1-qc.sh` | bowtie2 rRNA depletion |
| `tree/vmr_L_aln.fasta` | `4-vmr-tree.py` | Reference backbone for the MAFFT tree |
| `diamond/vmr_{l,m,s}.fa` | Nobody at run time | Source for rebuilding `.dmnd` |
| `clean/*` | `build_clean_db.py` | Input only |

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

`build_clean_db.py`, one line per segment:

```
  vmr_l.fa: 570 proteins (prodigal covered 560/570 nucleotide records)
  note: 10 nucleotide records yielded no protein (short fragments, expected)
  vmr_l.dmnd built
```

| Field | Meaning |
|---|---|
| Protein count | Records written to `vmr_<seg>.fa` |
| `prodigal covered n/m` | Records with at least one predicted protein, out of total input records |
| The note line | Uncovered records are short fragments and are normally expected |

A coverage ratio near 1.0 means full-length references. A ratio near 0.5 means
the library contains many fragments — check that this matches the curation intent
before proceeding.

`build_tree_ref.py`:

```
95 species -> 84 representatives
protein library hits: 84
protein length: min 2102 / median 6431 / max 6563 aa
running: mafft --thread 8 --auto <db>/tree/vmr_L_ref_prot.fa
output: <db>/tree/vmr_L_aln.fasta
  84 records / alignment length 6588
  checks passed: protein alphabet, equal-length alignment
```

| Field | Meaning |
|---|---|
| `species -> representatives` | Species count in the manifest versus usable L-segment representatives |
| `protein library hits` | Representatives with a protein in `vmr_l.fa`. Fewer than the previous number means Prodigal predicted no ORF for some |
| Alignment length | Should be close to the maximum protein length above |
| `checks passed` | The two self-checks in step 5 both succeeded |

## Manifest, before and after

| Source | Installed | Difference |
|---|---|---|
| `clean/manifest_clean.csv` | `manifest.csv` | None. The two files are byte-identical apart from the name |

`install_manifest.py` validates that `acc_base`, `species` and `segment` are
present before copying, so a malformed manifest fails here rather than during
tree building.

## Downstream impact

```
DB/clean/vmr_{L,M,S}_clean.fasta          (input, supplied externally)
DB/clean/manifest_clean.csv
                ↓  install_manifest.py
DB/manifest.csv
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
