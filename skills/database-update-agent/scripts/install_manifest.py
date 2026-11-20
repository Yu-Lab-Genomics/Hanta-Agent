#!/usr/bin/env python3
"""Install the reference manifest as <db>/manifest.csv.

The pipeline scripts (3-vmr-merge-cluster.py, 4-vmr-tree.py) read
<db>/manifest.csv. The curated library ships the same table under
<db>/clean/manifest_clean.csv. This script performs that one-step promotion and
refuses to run if the required columns are absent, so a malformed manifest is
caught here rather than halfway through tree building.

Where things live (all under <db_dir>, which defaults to <project_dir>/DB):

    <db>/clean/manifest_clean.csv   input, supplied with the cleaned library
    <db>/manifest.csv               OUTPUT, read by agents 3 and 4

Usage
-----
    python install_manifest.py                      # config.yaml -> project_dir/DB
    python install_manifest.py --db-dir /data/DB
    python install_manifest.py --force              # replace an existing file
    python install_manifest.py --check              # validate only, write nothing

If <db>/manifest.csv already exists and --force is not given, the script reports
that it is already installed and exits 0. That makes it safe to call in a setup
script that may run more than once.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

CONFIG_NAMES = ("config.yaml", "config-cell-paper.yaml", "cell_submit-paper.yaml")
REQUIRED_COLUMNS = ("acc_base", "species", "segment")


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


def inspect(path):
    """Return (row count, header list) for a CSV, or raise with a clear message."""
    with path.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("the file is empty")
        rows = sum(1 for _ in reader)
    return rows, header


def main():
    ap = argparse.ArgumentParser(
        description="Install <db>/clean/manifest_clean.csv as <db>/manifest.csv")
    ap.add_argument("--db-dir", type=Path, default=None,
                    help="database directory; defaults to <project_dir>/DB from the "
                         "configuration")
    ap.add_argument("--config", type=str, default="",
                    help="path to config.yaml; used to derive <project_dir>/DB")
    ap.add_argument("--source", type=Path, default=None,
                    help="explicit source file; defaults to "
                         "<db>/clean/manifest_clean.csv")
    ap.add_argument("--force", action="store_true",
                    help="replace an existing manifest.csv")
    ap.add_argument("--check", action="store_true",
                    help="validate the source and the installed file, write nothing")
    args = ap.parse_args()

    db_dir = resolve_db_dir(args)
    if db_dir is None:
        print("Cannot determine the database directory. Pass --db-dir, or make "
              "sure project_dir is set in config.yaml.", file=sys.stderr)
        return 2

    source = args.source or (db_dir / "clean" / "manifest_clean.csv")
    target = db_dir / "manifest.csv"

    print("database directory : %s" % db_dir)
    print("source             : %s" % source)
    print("target             : %s" % target)
    print("")

    # Validate whatever is available. --check accepts the target on its own.
    checked = False
    if source.is_file():
        try:
            rows, header = inspect(source)
        except (OSError, ValueError) as exc:
            print("cannot read the source manifest: %s" % exc, file=sys.stderr)
            return 2
        absent = [c for c in REQUIRED_COLUMNS if c not in header]
        if absent:
            print("source manifest is missing required column(s): %s"
                  % ", ".join(absent), file=sys.stderr)
            return 2
        print("source: %d rows, %d columns" % (rows, len(header)))
        checked = True
    elif target.is_file():
        print("source not present; an installed manifest already exists")
    else:
        print("source manifest not found: %s" % source, file=sys.stderr)
        print("This file ships with the cleaned library. Supply it, or pass "
              "--source.", file=sys.stderr)
        return 2

    if args.check:
        if target.is_file():
            try:
                rows, header = inspect(target)
            except (OSError, ValueError) as exc:
                print("cannot read %s: %s" % (target, exc), file=sys.stderr)
                return 2
            absent = [c for c in REQUIRED_COLUMNS if c not in header]
            if absent:
                print("installed manifest is missing required column(s): %s"
                      % ", ".join(absent), file=sys.stderr)
                return 2
            print("installed: %d rows, %d columns" % (rows, len(header)))
        print("")
        print("check only: nothing was written")
        return 0

    if target.is_file() and not args.force:
        rows, header = inspect(target)
        print("already installed: %d rows, %d columns" % (rows, len(header)))
        print("pass --force to replace it")
        return 0

    if not checked:
        print("cannot install: the source manifest is not available", file=sys.stderr)
        return 2

    try:
        shutil.copy2(source, target)
    except OSError as exc:
        print("copy failed: %s" % exc, file=sys.stderr)
        return 2

    rows, header = inspect(target)
    print("installed: %d rows, %d columns -> %s" % (rows, len(header), target))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
