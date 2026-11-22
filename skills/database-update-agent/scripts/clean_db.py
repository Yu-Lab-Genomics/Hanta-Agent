#!/usr/bin/env python3
"""Build the cleaned reference library from the dereplicated download.

This is the first step of the database build. It applies the curation rules to
the full NCBI download and writes the library that the later steps consume.

Curation rules
--------------
1. Drop genus-level placeholders ("Orthohantavirus sp.", "Hantaviridae sp.", ...).
2. Keep a record when NCBI marks it complete, or when it reaches the segment
   length floor (L >= 5000, M >= 3000, S >= 1400 nt).
3. Drop vague names that carry no place name ("Skink hantavirus").
4. No absolute lower length bound beyond rule 2.

Where things live (all under <db_dir>, which defaults to <project_dir>/DB):

    <db>/manifest.csv                     input  - full metadata table
    <db>/diamond/vmr_{L,M,S}_99.fasta     input  - dereplicated nucleotides
    <db>/clean/manifest_clean.csv         OUTPUT - curated metadata
    <db>/clean/vmr_{L,M,S}_clean.fasta    OUTPUT - curated nucleotides

Usage
-----
    python clean_db.py                       # config.yaml -> project_dir/DB
    python clean_db.py --db-dir /data/DB     # explicit database directory
    python clean_db.py --dry-run             # report what would be kept
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

GENERA = ["Orthohantavirus", "Mobatvirus", "Thottimvirus", "Loanvirus",
          "Actinovirus", "Percilovirus", "Agnathovirus", "Reptillovirus"]
THR = {"L": 5000, "M": 3000, "S": 1400}
VAGUE = {"Skink hantavirus"}  # no place name, deliberately excluded
SEGMENTS = [("L", "L"), ("M", "M"), ("S", "S")]
CONFIG_NAMES = ("config.yaml", "config-cell-paper.yaml", "cell_submit-paper.yaml")


def find_default_config(script_dir):
    """Walk up from the script looking for a configuration file."""
    for d in (script_dir, script_dir.parent,
              script_dir.parent.parent, script_dir.parent.parent.parent):
        for name in CONFIG_NAMES:
            cand = d / name
            if cand.is_file():
                return cand
    return None


def read_project_dir(config_path):
    """Read project_dir from a simple two-level YAML file."""
    for raw in config_path.read_text(encoding="utf-8").splitlines():
        line = raw.rstrip()
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if line[:1].strip() and ":" in line:
            key, value = line.split(":", 1)
            if key.strip() == "project_dir":
                return value.split(" #", 1)[0].strip().strip("'\"")
    return ""


def resolve_db_dir(args):
    """Resolve the database directory from --db-dir, --config, or discovery."""
    if args.db_dir:
        return Path(args.db_dir).expanduser().resolve()
    cfg = Path(args.config).expanduser() if args.config else find_default_config(
        Path(__file__).resolve().parent)
    if cfg and cfg.is_file():
        project_dir = read_project_dir(cfg)
        if project_dir:
            return (Path(project_dir).expanduser() / "DB").resolve()
        print("project_dir is missing from %s" % cfg, file=sys.stderr)
        return None
    return None


def kind(species):
    """Classify a species name as a placeholder, an ICTV binomial or a traditional name."""
    sp = (species or "").strip()
    if re.search(r"\bsp\.", sp) or sp in GENERA:
        return "placeholder"
    for g in GENERA:
        if sp.startswith(g + " ") and " sp" not in sp:
            return "ictv"
    return "traditional"


def is_complete(row):
    t = ((row.get("description") or "") + " " + (row.get("title") or "")).lower()
    return ("complete genome" in t) or ("complete cds" in t) or ("complete sequence" in t)


def length_nt(row):
    try:
        return int(float(row.get("length_nt") or 0))
    except ValueError:
        return 0


def segment_of(row):
    s = (row.get("segment") or "").upper()
    for x in ("L", "M", "S"):
        if s == x or s.startswith(x):
            return x
    return "?"


def main():
    ap = argparse.ArgumentParser(
        description="Build the cleaned reference library from the dereplicated "
                    "download (curation rules applied to manifest + sequence)")
    ap.add_argument("--db-dir", type=Path, default=None,
                    help="database directory; defaults to <project_dir>/DB from the "
                         "configuration")
    ap.add_argument("--config", type=str, default="",
                    help="path to config.yaml; used to derive <project_dir>/DB")
    ap.add_argument("--dry-run", action="store_true",
                    help="report the curation outcome without writing anything")
    args = ap.parse_args()

    db_dir = resolve_db_dir(args)
    if db_dir is None:
        print("Cannot determine the database directory. Pass --db-dir, or make "
              "sure project_dir is set in config.yaml.", file=sys.stderr)
        return 2

    src_manifest = db_dir / "manifest.csv"
    diamond = db_dir / "diamond"
    out_dir = db_dir / "clean"

    if not src_manifest.is_file():
        print("manifest.csv not found: %s" % src_manifest, file=sys.stderr)
        print("This is the full metadata table that ships with the repository.",
              file=sys.stderr)
        return 2

    print("database directory : %s" % db_dir)
    print("metadata           : %s" % src_manifest)
    print("nucleotides        : %s" % diamond)
    print("cleaned output     : %s" % out_dir)
    print("")

    with src_manifest.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        fields = reader.fieldnames
        rows = [dict(r) for r in reader]

    keep = []
    drop = Counter()
    for r in rows:
        sp = (r.get("species") or "").strip()
        if kind(sp) == "placeholder":
            drop["placeholder species"] += 1
            continue
        if sp in VAGUE:
            drop["vague name"] += 1
            continue
        seg = segment_of(r)
        if not (is_complete(r) or (seg in THR and length_nt(r) >= THR[seg])):
            drop["neither complete nor long enough"] += 1
            continue
        keep.append(r)

    print("input records : %d" % len(rows))
    for k, v in drop.most_common():
        print("  dropped %-34s %d" % (k, v))
    print("kept records  : %d" % len(keep))
    species = {r.get("species", "").strip() for r in keep}
    print("species       : %d (ICTV %d + traditional %d)" % (
        len(species),
        len({s for s in species if kind(s) == "ictv"}),
        len({s for s in species if kind(s) == "traditional"})))

    missing = [seg for seg, _ in SEGMENTS
               if not (diamond / ("vmr_%s_99.fasta" % seg)).is_file()]
    if missing:
        print("nucleotide input missing for segment(s): %s"
              % ", ".join(missing), file=sys.stderr)
        print("Expected <db>/diamond/vmr_{L,M,S}_99.fasta", file=sys.stderr)
        return 2

    if args.dry_run:
        print("")
        print("dry run: nothing was written")
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)

    # Metadata
    with (out_dir / "manifest_clean.csv").open("w", encoding="utf-8",
                                               newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter=",",
                                lineterminator="\n")
        writer.writeheader()
        for r in keep:
            writer.writerow({k: r.get(k, "") for k in fields})

    # Sequences, filtered to the records that survived curation
    keep_acc = {(r.get("acc_base") or "").strip() for r in keep}
    for seg, _ in SEGMENTS:
        src = diamond / ("vmr_%s_99.fasta" % seg)
        dst = out_dir / ("vmr_%s_clean.fasta" % seg)
        n = 0
        writing = False
        with src.open(encoding="utf-8", errors="ignore") as fh, \
                dst.open("w", encoding="utf-8") as out:
            for line in fh:
                if line.startswith(">"):
                    acc = line[1:].split()[0]
                    writing = acc in keep_acc
                    if writing:
                        out.write(line)
                        n += 1
                elif writing:
                    out.write(line)
        print("  vmr_%s_clean.fasta: %d sequences" % (seg, n))

    print("")
    print("next: python install_manifest.py --db-dir %s" % db_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
