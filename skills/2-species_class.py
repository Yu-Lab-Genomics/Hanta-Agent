#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
"""
对于megahit组装后的contigs.fasta文件，使用DIAMOND进行物种分类，输出到06_diamond_vmr目录下

text                                                                                                                                                                                          
读取 skills_V2/config.yaml                                                                                                                                                                    
读取 <results_dir>/json/.json                                                                                                                                                         
提取 outputs.contigs_fasta                                                                                                                                                           
对 vmr_l / vmr_m / vmr_s 三个 DIAMOND 数据库生成并执行 blastx                                                                                                                                 
输出到 <results_dir>/06_diamond_vmr/                                                                                                                                                                   
日志写到 <results_dir>/log/<sample>/                                                                                                                                                                
结构化 JSON 写到 <results_dir>/json/<sample>.species_class.json 
"""

DB_NAMES = ["vmr_l", "vmr_m", "vmr_s"]
OUTFMT = [
    "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]


def read_simple_yaml(path):
    cfg = {}
    section = None
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            indent = len(line) - len(line.lstrip(" "))
            text = line.strip()
            if ":" not in text:
                continue
            key, value = text.split(":", 1)
            key = key.strip()
            value = value.split(" #", 1)[0].strip()
            value = value.strip('"\'')
            if indent == 0:
                if value == "":
                    cfg[key] = {}
                    section = key
                else:
                    cfg[key] = value
                    section = None
            elif section:
                if not isinstance(cfg.get(section), dict):
                    cfg[section] = {}
                cfg[section][key] = value
    return cfg


def get_config_value(cfg, section, key, default):
    data = cfg.get(section, {})
    if isinstance(data, dict):
        return data.get(key, default) or default
    return default


def get_path_value(cfg, project_dir, key, default_name):
    paths = cfg.get("paths", {})
    if isinstance(paths, dict):
        value = paths.get(key, "")
        if value:
            return Path(value)
    return project_dir / default_name


def find_samples(json_dir, wanted_sample):
    samples = {}
    for json_file in sorted(json_dir.glob("*.json")):
        try:
            data = json.loads(json_file.read_text(encoding="utf-8"))
        except Exception as e:
            raise RuntimeError("读取 JSON 失败: %s: %s" % (json_file, e))
        if not isinstance(data, dict):
            continue
        for sample, item in data.items():
            if wanted_sample and sample != wanted_sample:
                continue
            if not isinstance(item, dict):
                continue
            outputs = item.get("outputs", {})
            if not isinstance(outputs, dict):
                continue
            fasta = outputs.get("contigs_fasta")
            if fasta:
                samples[sample] = {
                    "sample": sample,
                    "source_json": str(json_file),
                    "contigs_fasta": fasta,
                }
    return samples


def write_result(result_json, sample, status, message, exit_code, query, outputs, logs):
    payload = {
        sample: {
            "sample": sample,
            "status": status,
            "message": message,
            "exit_code": exit_code,
            "query": query,
            "outputs": outputs,
            "logs": logs,
        }
    }
    result_json.parent.mkdir(parents=True, exist_ok=True)
    result_json.write_text(json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")


def run_cmd(cmd, log_file, dry_run):
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with open(log_file, "w", encoding="utf-8") as log:
        log.write("CMD: " + " ".join(cmd) + "\n")
        if dry_run:
            log.write("DRY_RUN: skip execution\n")
            return 0
        p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
        return p.returncode


def main():
    script_dir = Path(__file__).resolve().parent
    default_config = script_dir.parent.parent / "skills_V2" / "config.yaml"

    parser = argparse.ArgumentParser(description="DIAMOND species classification for hantavirus contigs")
    parser.add_argument("-s", "--sample", default="", help="只处理指定样本")
    parser.add_argument("-c", "--config", default=str(default_config), help="配置文件路径")
    parser.add_argument("--dry-run", action="store_true", help="只生成命令和结果 JSON，不执行 DIAMOND")
    parser.add_argument("--overwrite", action="store_true", help="允许覆盖已有 tab 输出")
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.is_file():
        print("找不到配置文件: %s" % config_path, file=sys.stderr)
        return 2

    cfg = read_simple_yaml(config_path)
    project_dir = cfg.get("project_dir")
    if not project_dir:
        print("config.yaml 中缺少 project_dir", file=sys.stderr)
        return 2

    project_dir = Path(project_dir)
    conda_env = cfg.get("conda_env", "") or "lpy-hanta-v2"
    threads = str(cfg.get("threads", "32") or "32")
    evalue = str(get_config_value(cfg, "diamond", "evalue", "0.00001"))
    max_target = str(get_config_value(cfg, "diamond", "max_target_seqs", "1"))

    results_dir = get_path_value(cfg, project_dir, "results_dir", "results")
    json_dir = get_path_value(cfg, results_dir, "json_dir", "json")
    db_dir = project_dir / "DB" / "diamond"
    out_dir = get_path_value(cfg, results_dir, "diamond_vmr_dir", "06_diamond_vmr")

    if not json_dir.is_dir():
        print("找不到 JSON 目录: %s" % json_dir, file=sys.stderr)
        return 2
    if not db_dir.is_dir():
        print("找不到 DIAMOND DB 目录: %s" % db_dir, file=sys.stderr)
        return 2

    samples = find_samples(json_dir, args.sample)
    if not samples:
        print("没有找到可处理的 contigs_fasta", file=sys.stderr)
        return 2

    out_dir.mkdir(parents=True, exist_ok=True)
    final_code = 0

    for sample, info in samples.items():
        query = Path(info["contigs_fasta"])
        log_root = get_path_value(cfg, results_dir, "log_dir", "log")
        log_dir = log_root / sample
        result_json = json_dir / (sample + ".species_class.json")
        outputs = {}
        logs = {}

        if not query.is_file():
            write_result(result_json, sample, "failed", "contigs_fasta 不存在: %s" % query, 2, str(query), outputs, logs)
            final_code = 2
            continue

        sample_failed = False
        message = "success"

        for db_name in DB_NAMES:
            db_file = db_dir / (db_name + ".dmnd")
            out_file = out_dir / (sample + "." + db_name + ".tsv")
            log_file = log_dir / ("species_class_" + db_name + ".log")
            outputs[db_name] = str(out_file)
            logs[db_name] = str(log_file)

            if not db_file.is_file():
                sample_failed = True
                message = "数据库不存在: %s" % db_file
                final_code = 2
                break

            if out_file.exists() and not args.overwrite:
                sample_failed = True
                message = "输出已存在，如需覆盖请加 --overwrite: %s" % out_file
                final_code = 2
                break

            cmd = [
                "conda", "run", "-n", conda_env, "diamond", "blastx",
                "--threads", threads,
                "-d", str(db_file),
                "-q", str(query),
                "-o", str(out_file),
                "--evalue", evalue,
                "--max-target-seqs", max_target,
                "--outfmt",
            ] + OUTFMT + ["--header"]

            code = run_cmd(cmd, log_file, args.dry_run)
            if code != 0:
                sample_failed = True
                message = "%s 比对失败，日志: %s" % (db_name, log_file)
                final_code = code
                break

        if sample_failed:
            write_result(result_json, sample, "failed", message, final_code, str(query), outputs, logs)
        else:
            write_result(result_json, sample, "success", message, 0, str(query), outputs, logs)

    return final_code


if __name__ == "__main__":
    sys.exit(main())
