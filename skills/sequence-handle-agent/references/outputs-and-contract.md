# Outputs and JSON contract — Sequence Handle Agent

Everything is rooted at `paths.results_dir`.

```
<results_dir>/
├── 06_diamond_vmr/
│   ├── <SAMPLE>.vmr_l.tsv        DIAMOND hits against the L-segment database
│   ├── <SAMPLE>.vmr_m.tsv        DIAMOND hits against the M-segment database
│   ├── <SAMPLE>.vmr_s.tsv        DIAMOND hits against the S-segment database
│   └── ... one triple per processed sample
├── json/
│   └── <SAMPLE>.species_class.json    <-- the interface the orchestrator reads
└── log/
    └── <SAMPLE>/
        ├── species_class_vmr_l.log
        ├── species_class_vmr_m.log
        └── species_class_vmr_s.log
```

---

## 1. Hit tables

One table per sample per database. The file has a three-line DIAMOND header
followed by tab-separated records:

```
# DIAMOND v2.2.5. http://github.com/bbuchfink/diamond
# Invocation: diamond blastx --threads ... -d .../vmr_l.dmnd -q .../final.contigs.fa ...
# Fields: Query Seq - id, Subject Seq - id, Percentage of identical matches, ...
```

Every non-`#` line has these twelve columns:

| # | Column | Meaning |
|---|---|---|
| 1 | qseqid | Contig identifier, as in `final.contigs.fa` |
| 2 | sseqid | Reference accession |
| 3 | pident | Amino-acid identity, % |
| 4 | length | Alignment length, amino acids |
| 5 | mismatch | Mismatches |
| 6 | gapopen | Gap openings |
| 7 | qstart | Alignment start on the query |
| 8 | qend | Alignment end on the query |
| 9 | sstart | Alignment start on the subject |
| 10 | send | Alignment end on the subject |
| 11 | evalue | Expect value |
| 12 | bitscore | Bit score |

A table with no records beyond the header means the sample has no hit in that
database. **That is not an error** — the run still reports `success`. Absence of
hits in all three databases is what the next agent interprets as "no detection".

When counting hits, skip the three `#` lines. A quick check:

```bash
grep -vc '^#' <results_dir>/06_diamond_vmr/<SAMPLE>.vmr_l.tsv
```

---

## 2. The JSON result

Path: `<results_dir>/json/<SAMPLE>.species_class.json`, one top-level key per
sample. This file has a **different shape** from the quality-control agent's
JSON: there is no `steps` block, `query` is recorded explicitly, and the log
paths live under `logs`.

### 2.1 Shape

```json
{
  "<SAMPLE>": {
    "sample":    "SRR24112974",
    "status":    "success",
    "message":   "success",
    "exit_code": 0,
    "query":     "/abs/.../04_megahit/SRR24112974/final.contigs.fa",
    "outputs": {
      "vmr_l": "/abs/.../06_diamond_vmr/SRR24112974.vmr_l.tsv",
      "vmr_m": "/abs/.../06_diamond_vmr/SRR24112974.vmr_m.tsv",
      "vmr_s": "/abs/.../06_diamond_vmr/SRR24112974.vmr_s.tsv"
    },
    "logs": {
      "vmr_l": "/abs/.../log/SRR24112974/species_class_vmr_l.log",
      "vmr_m": "/abs/.../log/SRR24112974/species_class_vmr_m.log",
      "vmr_s": "/abs/.../log/SRR24112974/species_class_vmr_s.log"
    }
  }
}
```

### 2.2 Field meanings

| Field | Type | Meaning |
|---|---|---|
| `status` | `success` \| `failed` | The verdict for this sample |
| `message` | string | `success`, or the reason for failure |
| `exit_code` | int | 0 on success, 2 for handled failures, the DIAMOND code otherwise |
| `query` | path | The contig FASTA that was aligned — useful for auditing which assembly was used |
| `outputs.<db>` | path | Written even when the table is empty, as long as DIAMOND itself returned 0 |
| `logs.<db>` | path | The log holds the full command on its first line |

### 2.3 The three keys downstream relies on

```
outputs.vmr_l
outputs.vmr_m
outputs.vmr_s
```

The dereplication agent reads these three paths and merges them into one record
per contig. Everything else in this JSON is metadata.

### 2.4 Failure JSON

```json
{
  "SRR24112974": {
    "sample": "SRR24112974",
    "status": "failed",
    "message": "output already exists, pass --overwrite to replace it: /abs/.../SRR24112974.vmr_l.tsv",
    "exit_code": 2,
    "query": "/abs/.../final.contigs.fa",
    "outputs": {"vmr_l": "..."},
    "logs":    {"vmr_l": "..."}
  }
}
```

`outputs` and `logs` are filled incrementally, so on failure they contain only
the databases that were attempted. Downstream must check `status` before reading
them.

---

## 3. Verifying a batch

For each sample:

1. `json/<SAMPLE>.species_class.json` exists and parses.
2. `status == "success"`.
3. All three `outputs` keys resolve to existing files.
4. Each table has at least the three header lines.

A minimal check, run from the results directory:

```bash
python3 - <<'PY'
import json, glob, os
for f in sorted(glob.glob("json/*.species_class.json")):
    d = json.load(open(f, encoding="utf-8"))
    for sample, rec in d.items():
        outs = rec.get("outputs", {})
        paths = [outs.get(k, "") for k in ("vmr_l", "vmr_m", "vmr_s")]
        ok = (rec.get("status") == "success"
              and all(p and os.path.isfile(p) for p in paths))
        empty = [k for k in ("vmr_l", "vmr_m", "vmr_s")
                 if outs.get(k) and os.path.isfile(outs[k])
                 and sum(1 for ln in open(outs[k]) if not ln.startswith("#")) == 0]
        print(f"{sample}\t{'OK' if ok else 'CHECK'}\tno-hits-in:{','.join(empty) or '-'}")
PY
```

`no-hits-in: vmr_l,vmr_m,vmr_s` for a sample is a legitimate outcome, not a
failure.

---

## 4. What this agent does *not* produce

- No contig sequence extraction — that is the dereplication agent.
- No taxon assignment — hits must first be merged and clustered.
- No aggregate counts across samples — the orchestrator sums the per-sample JSON
  files.
- No coverage or identity statistics beyond the raw DIAMOND columns.
