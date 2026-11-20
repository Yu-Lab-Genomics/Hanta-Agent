#!/usr/bin/env python3
import argparse
import csv
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from species_norm import normalize_species
"""
Merge the three per-database DIAMOND hit tables into one record per contig.

Steps
-----
1. Scan <sample>.vmr_s.tsv / vmr_m.tsv / vmr_l.tsv and group them by sample id.
2. Skip DIAMOND comment lines starting with '#'.
3. Take the identifier from the first column.
4. Merge the hit status of the three tables.
5. Deduplicate by identifier.
6. Write:
   - <sample>.vmr_merged.tsv
   - <sample>.vmr_ids.txt
   - the dry-run status JSON
7. Outside dry-run mode, continue with:
   - seqkit grep -f id.txt final.contigs.fa > hits.fa
   - cd-hit-est ...
   - parse the .clstr file
   - write <sample>.vmr_final.tsv
"""
DBS = ["s", "m", "l"]
PICK_ORDER = ["l", "m", "s"]
FIELDS = [
    "id", "vmr_s", "vmr_m", "vmr_l",
    "source_vmr", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]


def read_yaml(path):
    cfg = {}
    section = None
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        indent = len(raw) - len(raw.lstrip(" "))
        text = raw.strip()
        if ":" not in text:
            continue
        key, val = text.split(":", 1)
        key = key.strip()
        val = val.split(" #", 1)[0].strip().strip("'\"")
        if indent == 0:
            if val == "":
                cfg[key] = {}
                section = key
            else:
                cfg[key] = val
                section = None
        elif section:
            cfg.setdefault(section, {})[key] = val
    return cfg


def cfg_path(cfg, project_dir, key, default):
    paths = cfg.get("paths", {})
    if isinstance(paths, dict) and paths.get(key):
        return Path(paths[key])
    return Path(default)


def cfg_section(cfg, name):
    section = cfg.get(name, {})
    return section if isinstance(section, dict) else {}


def cfg_value(section, key, default):
    val = section.get(key, "")
    return str(val if val not in ("", None) else default)


def cfg_extra_args(section, key):
    val = section.get(key, "")
    if not val:
        return []
    return shlex.split(str(val))


def build_seqkit_cmd(section, ids, contigs):
    return ["seqkit", "grep"] + cfg_extra_args(section, "seqkit_grep_extra_args") + ["-f", ids, contigs]


def build_cdhit_cmd(section, hit_fa, clustered, threads):
    return [
        "cd-hit-est",
        "-i", hit_fa,
        "-o", clustered,
        "-c", cfg_value(section, "cd_hit_est_c", "0.99"),
        "-aS", cfg_value(section, "cd_hit_est_aS", "0.8"),
        "-n", cfg_value(section, "cd_hit_est_n", "10"),
        "-g", cfg_value(section, "cd_hit_est_g", "1"),
        "-G", cfg_value(section, "cd_hit_est_G", "0"),
        "-d", cfg_value(section, "cd_hit_est_d", "0"),
        "-p", cfg_value(section, "cd_hit_est_p", "1"),
        "-M", cfg_value(section, "cd_hit_est_M", "0"),
        "-T", cfg_value(section, "cd_hit_est_T", threads),
    ] + cfg_extra_args(section, "cd_hit_est_extra_args")


def run(cmd, log_file, dry_run):
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with log_file.open("w", encoding="utf-8") as log:
        log.write("CMD: " + " ".join(map(str, cmd)) + "\n")
        if dry_run:
            log.write("DRY_RUN: skip\n")
            return 0
        p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
        return p.returncode


def run_to_file(cmd, out_file, log_file, dry_run):
    out_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with log_file.open("w", encoding="utf-8") as log:
        log.write("CMD: " + " ".join(map(str, cmd)) + " > " + str(out_file) + "\n")
        if dry_run:
            log.write("DRY_RUN: skip\n")
            return 0
        with out_file.open("w") as out:
            p = subprocess.run(cmd, stdout=out, stderr=log)
        return p.returncode


def samples_from_dir(vmr_dir):
    samples = {}
    for p in sorted(vmr_dir.glob("*.vmr_*.tsv")):
        name = p.name
        if ".vmr_" not in name:
            continue
        sample, rest = name.split(".vmr_", 1)
        db = rest.split(".", 1)[0]
        if db in DBS:
            samples.setdefault(sample, {})[db] = p
    return samples


def read_hits(path):
    rows = []
    seen = set()
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            sid = parts[0]
            if sid and sid not in seen:
                seen.add(sid)
                rows.append((sid, parts[1:]))
    return rows


def merge_tables(files, merged_tsv, ids_txt):
    table = {}
    order = []
    for db in DBS:
        for sid, rest in read_hits(files[db]):
            if sid not in table:
                table[sid] = {"s": "", "m": "", "l": "", "rows": {}}
                order.append(sid)
            table[sid][db] = db
            if db not in table[sid]["rows"]:
                table[sid]["rows"][db] = rest
    merged_tsv.parent.mkdir(parents=True, exist_ok=True)
    with merged_tsv.open("w", encoding="utf-8") as out:
        out.write("\t".join(FIELDS) + "\n")
        for sid in order:
            row = table[sid]
            source = ""
            detail = []
            for db in PICK_ORDER:
                if db in row["rows"]:
                    source = db
                    detail = row["rows"][db]
                    break
            detail = detail + [""] * (len(FIELDS) - 5 - len(detail))
            values = [sid, row["s"], row["m"], row["l"], source] + detail[:len(FIELDS) - 5]
            out.write("\t".join(values) + "\n")
    ids_txt.write_text("\n".join(order) + ("\n" if order else ""), encoding="utf-8")
    return order


def find_contigs(results_dir, sample):
    p = results_dir / "04_megahit" / sample / "final.contigs.fa"
    if p.is_file():
        return p
    hits = sorted((results_dir / "04_megahit").glob(sample + "*/final.contigs.fa"))
    if hits:
        return hits[0]
    return p


def parse_cluster_reps(clstr):
    reps = []
    if not clstr.is_file():
        return reps
    for line in clstr.read_text(encoding="utf-8", errors="ignore").splitlines():
        if "*" not in line or ">" not in line:
            continue
        sid = line.split(">", 1)[1].split("...", 1)[0]
        if sid:
            reps.append(sid)
    return reps


def filter_final(merged_tsv, reps, final_tsv):
    keep = set(reps)
    kept = 0
    with merged_tsv.open("r", encoding="utf-8") as src, final_tsv.open("w", encoding="utf-8") as out:
        header = src.readline()
        out.write(header)
        for line in src:
            sid = line.split("\t", 1)[0]
            if sid in keep:
                out.write(line)
                kept += 1
    return kept


def append_manifest(final_tsv, classified_tsv, manifest_csv):
    with manifest_csv.open("r", encoding="utf-8-sig", newline="") as src:
        manifest_reader = csv.reader(src)
        manifest_header = next(manifest_reader)
        if len(manifest_header) < 2 or manifest_header[1] != "acc_base":
            raise ValueError("the second column of manifest.csv must be acc_base")
        extra_header = manifest_header[2:]
        species_extra_i = manifest_header.index("species") - 2 if "species" in manifest_header else None
        manifest = {}
        for row in manifest_reader:
            if len(row) < 2:
                continue
            manifest[row[1]] = row[2:]

    classified_tsv.parent.mkdir(parents=True, exist_ok=True)
    matched = 0
    total = 0
    with final_tsv.open("r", encoding="utf-8", newline="") as src, \
            classified_tsv.open("w", encoding="utf-8", newline="") as out:
        reader = csv.reader(src, delimiter="\t")
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        header = next(reader)
        try:
            sseqid_index = header.index("sseqid")
        except ValueError as exc:
            raise ValueError("vmr_final.tsv has no sseqid column") from exc
        writer.writerow(header + extra_header)
        for row in reader:
            total += 1
            key = row[sseqid_index].split("_", 1)[0] if len(row) > sseqid_index else ""
            extra = manifest.get(key)
            if extra is not None:
                matched += 1
                extra = list(extra)
                if species_extra_i is not None and species_extra_i < len(extra):
                    extra[species_extra_i] = normalize_species(extra[species_extra_i])
                writer.writerow(row + extra)
            else:
                writer.writerow(row)
    return total, matched


def write_json(path, sample, status, message, outputs):
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {sample: {"sample": sample, "status": status, "message": message, "outputs": outputs}}
    path.write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")


def process_sample(sample, files, cfg, args):
    project_dir = Path(cfg["project_dir"])
    conda_env = cfg.get("conda_env", "") or "lpy-hanta-v2"
    results_dir = cfg_path(cfg, project_dir, "results_dir", project_dir / "results")
    vmr_dir = results_dir / "06_diamond_vmr"
    log_dir = results_dir / "log" / sample
    threads = str(cfg.get("threads", "32") or "32")
    prefix = vmr_dir / sample
    merged = Path(str(prefix) + ".vmr_merged.tsv")
    ids = Path(str(prefix) + ".vmr_ids.txt")
    hit_fa = Path(str(prefix) + ".vmr_hits.fa")
    clustered = Path(str(prefix) + ".vmr_clustered.fa")
    final = Path(str(prefix) + ".vmr_final.tsv")
    classified = Path(str(prefix) + ".vmr_final_classified.tsv")
    manifest_csv = Path(args.manifest)
    result_json = results_dir / "json" / (sample + ".vmr_cluster.json")
    outputs = {"merged_tsv": str(merged), "ids_txt": str(ids), "hits_fasta": str(hit_fa), "clustered_fasta": str(clustered), "final_tsv": str(final), "classified_tsv": str(classified)}

    missing = [d for d in DBS if d not in files]
    if missing:
        write_json(result_json, sample, "failed", "missing VMR tables: " + ",".join(missing), outputs)
        return 2
    contigs = find_contigs(results_dir, sample)
    if not contigs.is_file():
        write_json(result_json, sample, "failed", "contigs not found: %s" % contigs, outputs)
        return 2

    ids_list = merge_tables(files, merged, ids)
    if not ids_list:
        write_json(result_json, sample, "success", "no VMR hit identifiers", outputs)
        return 0

    vmr_cfg = cfg_section(cfg, "vmr_merge_cluster")
    code = run_to_file(["conda", "run", "-n", conda_env] + build_seqkit_cmd(vmr_cfg, ids, contigs), hit_fa, log_dir / "vmr_seqkit.log", args.dry_run)
    if code != 0:
        write_json(result_json, sample, "failed", "seqkit failed", outputs)
        return code

    code = run(["conda", "run", "-n", conda_env] + build_cdhit_cmd(vmr_cfg, hit_fa, clustered, threads), log_dir / "vmr_cdhit.log", args.dry_run)
    if code != 0:
        write_json(result_json, sample, "failed", "cd-hit-est failed", outputs)
        return code
    if args.dry_run:
        write_json(result_json, sample, "success", "dry-run complete", outputs)
        return 0

    reps = parse_cluster_reps(Path(str(clustered) + ".clstr"))
    kept = filter_final(merged, reps, final)
    try:
        total, matched = append_manifest(final, classified, manifest_csv)
    except (OSError, ValueError) as exc:
        write_json(result_json, sample, "failed", "manifest append failed: %s" % exc, outputs)
        return 2
    write_json(result_json, sample, "success",
               "done, %s representatives kept, manifest matched %s/%s"
               % (kept, matched, total), outputs)
    return 0


def find_default_config(script_dir):
    """Locate the shared configuration by walking up from this script.

    Handles both supported layouts:
      * skill layout  <root>/<skill>/scripts/<script>.py  -> <root>/config.yaml
      * flat layout   <root>/<script>.py                  -> <root>/config.yaml
        (or a named project config in the same directory)

    An explicit ``-c`` always wins.
    """
    for d in (script_dir, script_dir.parent,
              script_dir.parent.parent, script_dir.parent.parent.parent):
        for name in ("config.yaml", "config-cell-paper.yaml",
                     "cell_submit-paper.yaml"):
            cand = d / name
            if cand.is_file():
                return cand
    return script_dir.parent.parent / "skills_V2" / "config.yaml"


def main():
    script_dir = Path(__file__).resolve().parent
    default_config = find_default_config(script_dir)
    ap = argparse.ArgumentParser()
    ap.add_argument("-s", "--sample", default="")
    ap.add_argument("-c", "--config", default=str(default_config))
    ap.add_argument("--manifest", default=None,
                    help="path to manifest.csv; defaults to project_dir/DB/manifest.csv "
                         "from the configuration")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    cfg = read_yaml(args.config)
    if not cfg.get("project_dir"):
        print("project_dir is missing from config.yaml", file=sys.stderr)
        return 2
    project_dir = Path(cfg["project_dir"])
    if not args.manifest:
        args.manifest = str(project_dir / "DB" / "manifest.csv")
    results_dir = cfg_path(cfg, project_dir, "results_dir", project_dir / "results")
    vmr_dir = results_dir / "06_diamond_vmr"
    groups = samples_from_dir(vmr_dir)
    if args.sample:
        groups = {args.sample: groups.get(args.sample, {})}
    if not groups:
        print("no VMR TSV found", file=sys.stderr)
        return 2
    final_code = 0
    for sample, files in groups.items():
        code = process_sample(sample, files, cfg, args)
        if code != 0:
            final_code = code
    return final_code


if __name__ == "__main__":
    sys.exit(main())