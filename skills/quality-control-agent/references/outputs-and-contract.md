# Outputs and JSON contract — Quality Control Agent

The agent writes five stage directories, one JSON result per sample, and one log
directory per sample. Everything is rooted at `paths.results_dir`.

```
<results_dir>/
├── 01_fastp/            <SAMPLE>_fastp_R1.fq.gz  <SAMPLE>_fastp_R2.fq.gz
│                        <SAMPLE>_fastp.html      <SAMPLE>_fastp.json
├── 02_prinseq/          <SAMPLE>_prinseq_good_out_fwd.fq.gz
│                        <SAMPLE>_prinseq_good_out_rev.fq.gz
│                        <SAMPLE>_prinseq_single_out_fwd/rev.fq.gz
│                        <SAMPLE>_prinseq_bad_out_fwd/rev.fq.gz
├── 03_clean_non_rRNA/   <SAMPLE>_clean_non_rRNA_R1.fq.gz  ..._R2.fq.gz
├── 04_megahit/<SAMPLE>/ final.contigs.fa  + MEGAHIT intermediates
├── 05_prodigal/         <SAMPLE>_final.contigs_proteins.faa
│                        <SAMPLE>_final.contigs_proteins.min200aa.faa
│                        <SAMPLE>_final.contigs_genes.ffn
│                        <SAMPLE>_final.contigs.gff
├── json/
│   └── <SAMPLE>.qc.json          <-- the only interface the orchestrator reads
└── log/
    └── <SAMPLE>/
        ├── pipeline.log          start/stop record for the whole sample
        ├── fastp.log             one file per stage; the command is echoed
        ├── prinseq.log
        ├── bowtie2.log
        ├── megahit.log
        └── prodigal.log
```

---

## 1. Stage-by-stage outputs

### 1.1 fastp

| File | Content |
|---|---|
| `<SAMPLE>_fastp_R1.fq.gz`, `..._R2.fq.gz` | Trimmed paired reads |
| `<SAMPLE>_fastp.html` | Human-readable quality report |
| `<SAMPLE>_fastp.json` | Machine-readable report; `summary.before_filtering.total_reads` and `summary.after_filtering.total_reads` give the per-library read counts used in the manuscript |

The HTML report is what the result report links to for every positive library.
**Never delete this pair of report files**, even when intermediates are cleaned.

### 1.2 PRINSEQ++

PRINSEQ++ names its outputs with the `-out_name` prefix plus a suffix, and the
agent resolves whichever pair actually exists, in this order:

```
<SAMPLE>_prinseq_good_out_fwd.fq.gz  + _rev.fq.gz     <- the usual case
<SAMPLE>_prinseq_good_out_1.fq.gz    + _2.fq.gz
<SAMPLE>_prinseq_good_out.1.fq.gz    + .2.fq.gz
<SAMPLE>_prinseq_good_out_R1.fq.gz   + _R2.fq.gz
```

`_good_out` holds pairs where both mates survived; `_single_out` holds runs where
only one mate survived; `_bad_out` holds rejected reads. Only `_good_out` is used
downstream. These files are uncompressed after decompression and are the largest
intermediate in the pipeline.

> During `--dry-run` the resolver cannot stat the files and simply reports the
> first candidate shape, so the paths printed in a dry run are indicative, not
> guaranteed.

### 1.3 bowtie2

`--un-conc-gz` writes unmapped reads as `<prefix>.1` and `<prefix>.2`; the agent
renames them to `<SAMPLE>_clean_non_rRNA_R1.fq.gz` and `_R2.fq.gz`. If either
renamed file is missing or empty the stage fails with exit code 1 — an empty
rRNA-depleted library means the whole sample is unusable.

### 1.4 MEGAHIT

`04_megahit/<SAMPLE>/final.contigs.fa` is the assembly. It is checked for
non-zero size immediately; an empty assembly aborts the run.

### 1.5 Prodigal

| File | Content |
|---|---|
| `..._proteins.faa` | All predicted proteins |
| `..._proteins.min200aa.faa` | Proteins of at least `min_protein_aa` (default 200) |
| `..._genes.ffn` | Predicted coding sequences (nucleotide) |
| `..._final.contigs.gff` | Coordinates in GFF format |

Only the length-filtered `.faa` is handed downstream.

---

## 2. The JSON result

Path: `<results_dir>/json/<SAMPLE>.qc.json`. One top-level key per sample.

### 2.1 Shape

```json
{
  "<SAMPLE>": {
    "sample":        "SRR24112974",
    "status":        "success",
    "failed_step":   "",
    "message":       "pipeline finished",
    "exit_code":     0,
    "threads":       64,
    "log_dir":       "/abs/path/results_cell_paper/log/SRR24112974",
    "result_json":   "/abs/path/results_cell_paper/json/SRR24112974.qc.json",
    "steps": {
      "fastp":    {"status": "ok", "log": "/abs/.../log/SRR24112974/fastp.log"},
      "prinseq":  {"status": "ok", "log": "/abs/.../log/SRR24112974/prinseq.log"},
      "bowtie2":  {"status": "ok", "log": "/abs/.../log/SRR24112974/bowtie2.log"},
      "megahit":  {"status": "ok", "log": "/abs/.../log/SRR24112974/megahit.log"},
      "prodigal": {"status": "ok", "log": "/abs/.../log/SRR24112974/prodigal.log"}
    },
    "outputs": {
      "fastp_r1":            "/abs/.../01_fastp/SRR24112974_fastp_R1.fq.gz",
      "fastp_r2":            "/abs/.../01_fastp/SRR24112974_fastp_R2.fq.gz",
      "prinseq_r1":          "/abs/.../02_prinseq/SRR24112974_prinseq_good_out_fwd.fq.gz",
      "prinseq_r2":          "/abs/.../02_prinseq/SRR24112974_prinseq_good_out_rev.fq.gz",
      "clean_non_rrna_r1":   "/abs/.../03_clean_non_rRNA/SRR24112974_clean_non_rRNA_R1.fq.gz",
      "clean_non_rrna_r2":   "/abs/.../03_clean_non_rRNA/SRR24112974_clean_non_rRNA_R2.fq.gz",
      "contigs_fasta":       "/abs/.../04_megahit/SRR24112974/final.contigs.fa",
      "prodigal_faa":        "/abs/.../05_prodigal/SRR24112974_final.contigs_proteins.faa",
      "prodigal_min_aa_faa": "/abs/.../05_prodigal/SRR24112974_final.contigs_proteins.min200aa.faa"
    }
  }
}
```

### 2.2 Field meanings

| Field | Type | Meaning |
|---|---|---|
| `status` | `success` \| `failed` | **The only reliable verdict.** Ignore stdout |
| `failed_step` | string | `fastp`, `prinseq`, `bowtie2`, `megahit`, `prodigal`, `preflight`, or empty |
| `message` | string | Human-readable reason, in the language of the script's author |
| `exit_code` | int | 0 on success, otherwise the code passed to `exit` |
| `threads` | int | Effective core budget used for the run |
| `steps.<stage>.status` | `ok` \| `failed` | Per-stage outcome; only stages that actually ran appear |
| `outputs.<key>` | path | Absolute path; **omitted entirely** when the file was cleaned up or never produced |

### 2.3 The two keys downstream may rely on

```
outputs.contigs_fasta
outputs.prodigal_min_aa_faa
```

Everything else is internal. A downstream agent that finds either key missing
must treat the sample as not ready rather than guess a path.

### 2.4 Failure JSON

On failure the same file is written with `status: failed` and the offending stage:

```json
{
  "SRR24112974": {
    "sample": "SRR24112974",
    "status": "failed",
    "failed_step": "megahit",
    "message": "MEGAHIT did not produce final.contigs.fa",
    "exit_code": 1,
    "steps": {"megahit": {"status": "failed", "log": "/abs/.../megahit.log"}}
  }
}
```

Argument and preflight errors are written by a lighter path that emits only
`sample`, `status`, `message` and `exit_code`.

---

## 3. Verifying a batch

For a completed batch, every sample should satisfy:

1. `json/<SAMPLE>.qc.json` exists and parses.
2. `status == "success"`.
3. `outputs.contigs_fasta` exists and is non-empty.
4. `outputs.prodigal_min_aa_faa` exists and is non-empty.

A minimal check, run from the results directory:

```bash
python3 - <<'PY'
import json, glob, os
for f in sorted(glob.glob("json/*.qc.json")):
    d = json.load(open(f, encoding="utf-8"))
    for sample, rec in d.items():
        ctg = rec.get("outputs", {}).get("contigs_fasta", "")
        faa = rec.get("outputs", {}).get("prodigal_min_aa_faa", "")
        ok = (rec.get("status") == "success"
              and ctg and os.path.getsize(ctg) > 0
              and faa and os.path.getsize(faa) > 0)
        print(f"{sample}\t{'OK' if ok else 'CHECK'}\t{rec.get('failed_step','')}")
PY
```

---

## 4. What this agent does *not* produce

- No taxon assignment — that is the classification agent.
- No coverage, identity or clustering statistics beyond the raw tool logs.
- No aggregated batch summary — the orchestrator reads the per-sample JSON files
  and builds its own view.
