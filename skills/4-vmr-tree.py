#!/usr/bin/env python3
import argparse
import csv
import json
import re
import shlex
import subprocess
import sys
from pathlib import Path

"""

基于 06_diamond_vmr 目录下的 *.vmr_final_classified.tsv 文件，使用 IQ-TREE 构建系统发育树，并对每个 query 序列进行分类。
输出： *.vmr_final_classified.tsv 
提取所有比对的l蛋白，到00_results目录下，用于lucaprot的分析。

排错机制，在json/*.vmr_cluster.json
输出	位置
总表（新）	<results_dir>/00_results/all_vmr_final_classified.tsv
l.fasta	<results_dir>/00_results/l.fasta
每样本 TSV（原地更新，加 6 列）	<diamond_vmr_dir>/<样本名>.vmr_final_classified.tsv
共享树	<results_dir>/07_tree/_all/
每样本中间文件	<results_dir>/07_tree/<样本名>/

[重新运行，从2代码开始运行]
"""

def read_yaml(path):
    cfg, section = {}, None
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        indent = len(raw) - len(raw.lstrip(" "))
        text = raw.strip()
        if ":" not in text:
            continue
        key, val = text.split(":", 1)
        val = val.split(" #", 1)[0].strip().strip("'\"")
        if indent == 0:
            section = key if not val else None
            cfg[key] = {} if not val else val
        elif section:
            cfg.setdefault(section, {})[key.strip()] = val
    return cfg


def value(cfg, section, key, default):
    obj = cfg.get(section, {})
    val = obj.get(key, "") if isinstance(obj, dict) else ""
    return str(val) if val not in ("", None) else str(default)


def extra(cfg, section, key):
    val = value(cfg, section, key, "")
    return shlex.split(val) if val else []


def software(cfg, name, default):
    return value(cfg, "software", name, default)


def run(cmd, log, dry=False, stdout_path=None):
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as fh:
        fh.write("CMD: " + shlex.join(map(str, cmd)) + "\n")
        if dry:
            fh.write("DRY_RUN: skip\n")
            return 0
        if stdout_path:
            with stdout_path.open("w", encoding="utf-8") as out:
                p = subprocess.run(cmd, stdout=out, stderr=fh)
        else:
            p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
        return p.returncode


def parse_table(path):
    with path.open(encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        raise ValueError("empty TSV: %s" % path)
    return rows[0], rows[1:]


def write_table(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def fasta_ids(path):
    return [line[1:].strip().split()[0] for line in path.open(encoding="utf-8", errors="ignore") if line.startswith(">")]


def parse_diamond(path):
    hits = {}
    if not path.is_file():
        return hits
    with path.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            try:
                score = float(p[11])
                pident = float(p[2])
            except ValueError:
                continue
            if p[0] not in hits or score > hits[p[0]]["bitscore"]:
                hits[p[0]] = {"sseqid": p[1], "bitscore": score, "pident": pident}
    return hits


class Node:
    def __init__(self):
        self.name = ""
        self.support = None
        self.length = 0.0
        self.children = []
        self.parent = None
        self.depth_nodes = 0
        self.depth_len = 0.0


def parse_newick(text):
    """解析 Newick（保留枝长与自展支持度）。"""
    s = text.strip()
    if s.endswith(";"):
        s = s[:-1]
    s = s.strip()
    n = len(s)
    i = 0

    def _node():
        nonlocal i
        node = Node()
        if s[i] == "(":
            i += 1
            while True:
                child = _node()
                child.parent = node
                node.children.append(child)
                if s[i] == ",":
                    i += 1
                    continue
                if s[i] == ")":
                    i += 1
                    break
                raise ValueError("bad newick at %d" % i)
            j = i
            while j < n and s[j] not in ":,()":
                j += 1
            label = s[i:j]
            i = j
            if label.strip():
                try:
                    node.support = float(label)
                except ValueError:
                    node.support = None
        else:
            j = i
            while j < n and s[j] not in ":,()":
                j += 1
            node.name = s[i:j]
            i = j
        if i < n and s[i] == ":":
            i += 1
            j = i
            while j < n and s[j] not in ",()":
                j += 1
            bl = s[i:j]
            node.length = float(bl) if bl.strip() else 0.0
            i = j
        return node

    root = _node()
    _annotate_tree(root)
    return root


def _annotate_tree(node, d_nodes=0, d_len=0.0):
    node.depth_nodes = d_nodes
    node.depth_len = d_len
    for c in node.children:
        _annotate_tree(c, d_nodes + 1, d_len + c.length)


def _tree_lca(a, b):
    while a.depth_nodes > b.depth_nodes:
        a = a.parent
    while b.depth_nodes > a.depth_nodes:
        b = b.parent
    while a is not b:
        a = a.parent
        b = b.parent
    return a


def _patristic(a, b):
    l = _tree_lca(a, b)
    return a.depth_len + b.depth_len - 2.0 * l.depth_len


def _subtree_leaves(node):
    out = []

    def walk(x):
        if not x.children:
            out.append(x)
        else:
            for c in x.children:
                walk(c)

    walk(node)
    return out


def _is_ancestor(anc, node):
    x = node
    while x is not None:
        if x is anc:
            return True
        x = x.parent
    return False


def classify_tree_nearest(tree_file, query_ids, species_map):
    """按“最近姐妹支”分类：取 patristic 距离最近的参考叶子的物种。

    返回 dict: qid(叶子名) -> dict(species, position, nn, dist, support)
      position: nested（落在该种单系支内）/ sister（贴支外）
    """
    root = parse_newick(tree_file.read_text(encoding="utf-8", errors="ignore"))
    all_leaves = _subtree_leaves(root)
    qpat = re.compile(r"_\d+$")
    query_leaves = [l for l in all_leaves if qpat.search(l.name)]
    ref_leaves = [l for l in all_leaves if not qpat.search(l.name)]

    ref_species = {}
    for l in ref_leaves:
        sp = species_map.get(l.name, species_map.get(l.name.split("_", 1)[0], ""))
        if sp:
            ref_species[l.name] = sp

    by_species = {}
    for l in ref_leaves:
        sp = ref_species.get(l.name)
        if sp:
            by_species.setdefault(sp, []).append(l)

    def _mrca(leaves):
        m = leaves[0]
        for x in leaves[1:]:
            m = _tree_lca(m, x)
        return m

    species_mrca = {sp: _mrca(ls) for sp, ls in by_species.items()}

    result = {}
    for q in query_leaves:
        support = None
        p = q.parent
        while p is not None and support is None:
            support = p.support
            p = p.parent
        cand_node = q.parent
        cands = []
        while cand_node is not None:
            leaves = _subtree_leaves(cand_node)
            cands = [l for l in leaves if l in ref_leaves]
            if cands:
                break
            cand_node = cand_node.parent
        if not cands:
            result[q.name] = {"species": "", "position": "", "nn": "", "dist": None, "support": support}
            continue
        nearest = min(cands, key=lambda c: _patristic(q, c))
        dist = _patristic(q, nearest)
        sp = ref_species[nearest.name]
        pos = "nested" if (sp in species_mrca and _is_ancestor(species_mrca[sp], q)) else "sister"
        result[q.name] = {"species": sp, "position": pos, "nn": nearest.name, "dist": dist, "support": support}
    return result


def build_species_map(manifest_csv, classified_tsv=None):
    """从 manifest.csv 构建 acc_base→species 映射；manifest 缺失时回退到 classified TSV。"""
    species_map = {}
    if manifest_csv is not None and Path(manifest_csv).is_file():
        with Path(manifest_csv).open("r", encoding="utf-8-sig", newline="") as fh:
            mreader = csv.reader(fh)
            mheader = next(mreader)
            macc_i = mheader.index("acc_base")
            msp_i = mheader.index("species")
            for mrow in mreader:
                if len(mrow) > max(macc_i, msp_i):
                    acc = mrow[macc_i].strip()
                    if acc:
                        species_map[acc] = mrow[msp_i].strip()
        return species_map
    if classified_tsv is None or not Path(classified_tsv).is_file():
        return species_map
    header, rows = parse_table(classified_tsv)
    idx = {name: i for i, name in enumerate(header)}
    species_i = idx.get("species")
    sseq_i = idx.get("sseqid")
    if species_i is None or sseq_i is None:
        return species_map
    for r in rows:
        if len(r) > max(species_i, sseq_i):
            ref = r[sseq_i].split("_", 1)[0]
            if r[species_i].strip():
                species_map[ref] = r[species_i].strip()
    return species_map


def discover(vmr_dir, sample=""):
    files = {}
    for p in sorted(vmr_dir.glob("*.vmr_final_classified.tsv")):
        s = p.name.split(".vmr_final_classified.tsv", 1)[0]
        if not sample or sample == s:
            files[s] = p
    return files


def extract_sample(sample, classified, cfg, args):
    """提取单个样本的 L 蛋白（每 contig 取最佳 RdRp），写带样本名前缀的 query fasta。"""
    project = Path(cfg.get("project_dir", "."))
    paths = cfg.get("paths", {})
    results = Path(paths.get("results_dir", project / "results"))
    prodigal = results / "05_prodigal"
    diamond_db = project / "DB" / "diamond" / "vmr_l.dmnd"
    out_dir = results / "07_tree" / sample
    logs = out_dir / "logs"
    out_dir.mkdir(parents=True, exist_ok=True)

    header, rows = parse_table(classified)
    if len(header) < 4:
        raise ValueError("TSV 至少需要 4 列")
    id_i, l_i = 0, 3
    targets = [r for r in rows if len(r) > l_i and r[l_i].strip()]

    if not targets:
        return {"sample": sample, "classified": classified, "header": header,
                "rows": rows, "selected": {}, "prefixed_fa": None}

    protein_faa = prodigal / f"{sample}_final.contigs_proteins.faa"
    if not protein_faa.is_file():
        raise FileNotFoundError(str(protein_faa))

    ids = [r[id_i].strip() for r in targets if r[id_i].strip()]
    ids_file = out_dir / f"{sample}.l_ids.txt"
    ids_file.write_text("\n".join(ids) + "\n", encoding="utf-8")
    query_faa = out_dir / f"{sample}.l_proteins.faa"
    pattern = "(?:" + "|".join(re.escape(seq_id) for seq_id in ids) + ")_[0-9]+"
    seqkit = [software(cfg, "seqkit", "seqkit"), "grep", "-n", "-r", "-p", pattern, str(protein_faa)] + extra(cfg, "vmr_merge_cluster", "seqkit_grep_extra_args")
    code = run(seqkit, logs / "seqkit.log", args.dry_run, query_faa)
    if code:
        raise RuntimeError("seqkit failed")
    diamond_tsv = out_dir / f"{sample}.vmr_l.diamond.tsv"
    dcmd = [software(cfg, "diamond", "diamond"), "blastp", "--threads", str(cfg.get("threads", 32)), "-d", str(diamond_db), "-q", str(query_faa), "-o", str(diamond_tsv), "--evalue", value(cfg, "diamond", "evalue", "0.00001"), "--max-target-seqs", value(cfg, "diamond", "max_target_seqs", "1"), "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "--header"]
    code = run(dcmd, logs / "diamond.log", args.dry_run)
    if code:
        raise RuntimeError("diamond failed")
    if args.dry_run:
        return {"sample": sample, "classified": classified, "header": header,
                "rows": rows, "selected": {}, "prefixed_fa": None}

    hits = parse_diamond(diamond_tsv)
    selected = {}
    for qid, hit in hits.items():
        selected[qid] = hit
    # 按 contig 分组，取每个 contig 最高分的一条 RdRp
    contig_groups = {}
    for prot_id, hit in selected.items():
        mat = re.match(r"(.*)_\d+$", prot_id)
        ctid = mat.group(1) if mat else prot_id
        contig_groups.setdefault(ctid, []).append((prot_id, hit))
    final_selected = {}
    for ctid, plist in contig_groups.items():
        plist.sort(key=lambda x: x[1]["bitscore"], reverse=True)
        best_pid, best_hit = plist[0]
        final_selected[best_pid] = best_hit
    selected = final_selected

    seq_ids = set(selected)
    selected_faa = out_dir / f"{sample}.rdrp.faa"          # 原 id（每样本一份）
    prefixed_fa = out_dir / f"{sample}.rdrp.prefixed.fa"   # 带样本名前缀（用于合并建树）
    with query_faa.open(encoding="utf-8", errors="ignore") as src, \
         selected_faa.open("w", encoding="utf-8") as dst, \
         prefixed_fa.open("w", encoding="utf-8") as pdst:
        keep = False
        for line in src:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                keep = name in seq_ids
                if keep:
                    dst.write(line)
                    pdst.write(">" + sample + "_" + name + "\n")
            elif keep:
                dst.write(line)
                pdst.write(line)

    return {"sample": sample, "classified": classified, "header": header,
            "rows": rows, "selected": selected, "prefixed_fa": prefixed_fa}


def build_combined_tree(cfg, args, combined_query, all_dir, ref):
    """用 73 参考 + 合并 query 建一棵树，返回 treefile 路径。"""
    logs = all_dir / "logs"
    combined = all_dir / "tree_input.fa"
    combined.write_text(ref.read_text(encoding="utf-8") + combined_query.read_text(encoding="utf-8"), encoding="utf-8")
    aln = all_dir / "tree_aln.fa"
    code = run([software(cfg, "mafft", "mafft"), "--thread", str(cfg.get("threads", 32)), "--" + value(cfg, "mafft", "preset", "auto"), str(combined)], logs / "mafft.log", False, aln)
    if code: raise RuntimeError("mafft failed")
    trimmed = all_dir / "tree_trim.fa"
    code = run([software(cfg, "trimal", "trimal"), "-in", str(aln), "-out", str(trimmed), "-" + value(cfg, "trimal", "method", "automated1"), "-" + value(cfg, "trimal", "out_format", "fasta")], logs / "trimal.log")
    if code: raise RuntimeError("trimal failed")
    prefix = all_dir / "tree_true_result"
    iq = [software(cfg, "iqtree3", "iqtree3"), "-s", str(trimmed), "-m", value(cfg, "iqtree3", "model", "MFP"), "-B", value(cfg, "iqtree3", "bootstrap", "1000"), "-T", str(cfg.get("threads", 32)), "-pre", str(prefix)]
    code = run(iq, logs / "iqtree.log")
    if code: raise RuntimeError("iqtree3 failed")
    return Path(str(prefix) + ".treefile")


def classify_sample(data, classifications):
    """用共享树对单个样本分类，写回 TSV（去掉样本名前缀映射回原 id）。"""
    sample = data["sample"]
    classified = data["classified"]
    header = data["header"]
    rows = data["rows"]
    selected = data["selected"]

    out_header = header + ["final_classification", "classification_source",
                            "tree_position", "tree_bootstrap", "tree_nn_sseqid", "tree_nn_distance"]
    if not selected:
        for r in rows:
            r.extend(["", "", "", "", "", ""])
        write_table(classified, out_header, rows)
        return

    for r in rows:
        rid = r[0]
        qid = next((k for k in selected if re.sub(r"_\d+$", "", k) == rid), None)
        prefixed = f"{sample}_{qid}" if qid else None
        tr = classifications.get(prefixed) if prefixed else None
        if tr and tr.get("species"):
            support = tr.get("support")
            val = [tr["species"], "tree", tr.get("position", ""),
                   ("" if support is None else "%.1f" % support),
                   tr.get("nn", ""),
                   ("" if tr.get("dist") is None else "%.6f" % tr["dist"])]
        else:
            val = ["", "", "", "", "", ""]
        r.extend(val)
    write_table(classified, out_header, rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-s", "--sample", default="")
    ap.add_argument("-c", "--config", required=True)
    ap.add_argument("--manifest", default=None,
                    help="manifest.csv 路径；默认取 config 里 project_dir/DB/manifest.csv")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    cfg = read_yaml(args.config)
    paths = cfg.get("paths", {})
    project = Path(cfg.get("project_dir", "."))
    if not args.manifest:
        args.manifest = str(project / "DB" / "manifest.csv")
    results = Path(paths.get("results_dir", project / "results"))
    vmr = Path(paths.get("diamond_vmr_dir", results / "06_diamond_vmr"))
    files = discover(vmr, args.sample)
    if not files:
        print("no classified TSV found", file=sys.stderr); return 2

    species_map = build_species_map(Path(args.manifest), None)

    # Phase 1：逐样本提取 L 蛋白（写带样本名前缀的 query）
    samples_data = []
    for sample, path in files.items():
        try:
            samples_data.append(extract_sample(sample, path, cfg, args))
        except Exception as exc:
            print(f"{sample}: {exc}", file=sys.stderr); return 2

    if args.dry_run:
        print(json.dumps({"status": "dry-run", "samples": [d["sample"] for d in samples_data]}, ensure_ascii=False))
        return 0

    # Phase 2：合并所有样本 query，建一棵树
    all_dir = results / "07_tree" / "_all"
    all_dir.mkdir(parents=True, exist_ok=True)
    combined_query = all_dir / "combined_query.fa"
    with combined_query.open("w", encoding="utf-8") as out:
        for d in samples_data:
            if d["prefixed_fa"] is not None and d["prefixed_fa"].is_file():
                out.write(d["prefixed_fa"].read_text(encoding="utf-8"))
    # 合并一份 l.fasta（带样本名前缀），放 results/00_results/l.fasta
    l_fasta = results / "00_results" / "l.fasta"
    l_fasta.parent.mkdir(parents=True, exist_ok=True)
    l_fasta.write_text(combined_query.read_text(encoding="utf-8"), encoding="utf-8")

    ref = project / "DB" / "tree" / "vmr_L_aln.fasta"
    if combined_query.stat().st_size == 0:
        classifications = {}
    else:
        treefile = build_combined_tree(cfg, args, combined_query, all_dir, ref)
        classifications = classify_tree_nearest(treefile, set(), species_map)

    # Phase 3：逐样本用同一棵树分类写回
    for d in samples_data:
        try:
            classify_sample(d, classifications)
            print(json.dumps({"sample": d["sample"], "status": "success",
                              "classified": len(d["selected"]), "output": str(d["classified"])}, ensure_ascii=False))
        except Exception as exc:
            print(f"{d['sample']}: {exc}", file=sys.stderr); return 2

    # 汇总所有样本成一个总表，输出到 00_results
    tree_cols = ["final_classification", "classification_source",
                 "tree_position", "tree_bootstrap", "tree_nn_sseqid", "tree_nn_distance"]
    combined_header = None
    combined_rows = []
    for d in samples_data:
        if combined_header is None:
            combined_header = ["sample"] + d["header"] + tree_cols
        for r in d["rows"]:
            combined_rows.append([d["sample"]] + r)
    if combined_header is not None:
        combined_path = results / "00_results" / "all_vmr_final_classified.tsv"
        combined_path.parent.mkdir(parents=True, exist_ok=True)
        write_table(combined_path, combined_header, combined_rows)
        print(json.dumps({"combined_table": str(combined_path), "rows": len(combined_rows)}, ensure_ascii=False))

    return 0

if __name__ == "__main__":
    sys.exit(main())