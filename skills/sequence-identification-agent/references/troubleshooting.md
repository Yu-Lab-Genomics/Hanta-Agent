# Troubleshooting — Sequence Identification Agent

## Read the JSON first

`<results_dir>/json/<SAMPLE>.vmr_cluster.json` carries a `status` and a `message`
for every handled failure. Start there, then follow the branch below.

## Failure branches

| Symptom in `message` | Cause | Fix |
|---|---|---|
| `missing VMR tables` | The sample has none of `vmr_{l,m,s}.tsv` in `06_diamond_vmr/` | Run the sequence-handle agent for this sample |
| `no VMR hit identifiers` (status `success`) | Every hit table is empty | Legitimate: the sample has no hantavirus hits. Downstream will score nothing. |
| `the second column of manifest.csv must be acc_base` | `DB/manifest.csv` was replaced or reordered | Restore the column order; the manifest reader depends on positional column 2 |
| Contig FASTA not found | `04_megahit/<SAMPLE>/final.contigs.fa` is missing | Re-run the QC agent — the assembly was never produced or was removed |
| `seqkit` exits non-zero | The ID list contains a header that is not in the contig FASTA | Compare `vmr_ids.txt` against the `> ` headers in `final.contigs.fa`; a DIAMOND database built from a different reference set is the usual cause |
| `cd-hit-est` exits non-zero | `cd_hit_est_n` is incompatible with `cd_hit_est_c` | See the word-size table in `parameters.md` |

## Symptom-driven checks

### Fewer representatives than expected

1. Check the `manifest matched x/y` ratio in the JSON. A low ratio means the
   manifest is stale relative to the DIAMOND database, not that clustering
   over-merged.
2. Inspect the `.clstr` file. A cluster with many members and 99 % identity is
   expected; a cluster spanning two different segments means `cd_hit_est_c` is
   too low.
3. Confirm `-g` is still `1`. With `-g 0` the representative depends on input
   order, so results are not reproducible between runs.

### All contigs collapse into one cluster

Almost always `cd_hit_est_c` was lowered, or `-aS` was set to `0`. Restore
`0.99` / `0.8` and re-run. If the collapse persists, the contig FASTA probably
contains a single very long concatenated sequence.

### Segment proportions do not match the publication

You are reading `vmr_s` / `vmr_m` / `vmr_l` as if they were mutually exclusive.
They are not: they record database hits. Use `source_vmr`.

### A dry run left a success JSON

`--dry-run` writes a `status: success` JSON without producing any real output.
Delete `<results_dir>/json/<SAMPLE>.vmr_cluster.json` before the real run,
otherwise the next reader will believe the sample is complete.

## Safe-testing rule

Never point a test run at a production `results_dir`. Copy one sample's inputs
into a throwaway directory under `/tmp` and run there. `--dry-run` in this agent
is genuinely non-destructive, but it still writes JSON and logs into the results
directory you point it at.
