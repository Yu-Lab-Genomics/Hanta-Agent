#!/usr/bin/env python3
"""Build the tree reference set from the cleaned library (one per species),
output as a **protein** alignment.

Why this matters
----------------
4-vmr-tree.py concatenates the reference and the query proteins verbatim and
hands the mixture to MAFFT::

    tree_input = ref.read_text() + combined_query.read_text()
    mafft --auto tree_input

So <db>/tree/vmr_L_aln.fasta **must be protein**. An earlier version of this
script only filtered nucleotide records; the resulting *_aln.fasta was in fact
unaligned nucleotide and made MAFFT abort with `Alphabet 'e' is unknown` - the
nucleotide alphabet has no amino acid E.

Where things live (all under <db_dir>, which defaults to <project_dir>/DB):

    <db>/manifest.csv          taxonomy, one row per reference record
    <db>/diamond/vmr_l.fa      L-segment proteins, written by build_clean_db.py
    <db>/tree/vmr_L_ref_prot.fa   unaligned proteins, kept on disk
    <db>/tree/vmr_L_aln.fasta     OUTPUT, the protein backbone 4-vmr-tree.py reads

The alignment goes straight into tree/ because that is where 4-vmr-tree.py looks
for it. It is read from diamond/ rather than from clean/ because build_clean_db.py
now writes the protein libraries there.

Usage
-----
    python build_tree_ref.py                       # config.yaml -> project_dir/DB
    python build_tree_ref.py --db-dir /data/DB     # explicit database directory
    python build_tree_ref.py --out FILE --threads 16
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

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


def is_complete(r):
    t = ((r.get("description") or "") + " " + (r.get("title") or "")).lower()
    return ("complete genome" in t) or ("complete cds" in t) or ("complete sequence" in t)


def nt(r):
    try:
        return int(float(r.get("length_nt") or 0))
    except ValueError:
        return 0


def read_fasta_dict(path):
    """Read a FASTA into {id: seq}; keep the longest when an id repeats."""
    out = {}
    cur, buf = None, []

    def flush():
        if cur is None:
            return
        s = "".join(buf)
        if cur not in out or len(s) > len(out[cur]):
            out[cur] = s

    with path.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                head = line[1:].strip()
                cur = head.split(None, 1)[0] if head else None
                buf = []
            elif cur is not None:
                buf.append(line.strip())
    flush()
    return out


def main():
    ap = argparse.ArgumentParser(
        description="Build the tree reference set (protein alignment)")
    ap.add_argument("--db-dir", type=Path, default=None,
                    help="database directory; defaults to <project_dir>/DB from the "
                         "configuration")
    ap.add_argument("--config", type=str, default="",
                    help="path to config.yaml; used to derive <project_dir>/DB")
    ap.add_argument("--out", type=Path, default=None,
                    help="output file; defaults to <db>/tree/vmr_L_aln.fasta")
    ap.add_argument("--mafft", default="mafft",
                    help="MAFFT executable (default: found on PATH)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--mafft-extra", default="--auto",
                    help="MAFFT arguments (default --auto)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the composed command without running anything")
    args = ap.parse_args()

    db_dir = resolve_db_dir(args)
    if db_dir is None:
        print("Cannot determine the database directory. Pass --db-dir, or make "
              "sure project_dir is set in config.yaml.", file=sys.stderr)
        return 2

    manifest = db_dir / "manifest.csv"
    prot_fa = db_dir / "diamond" / "vmr_l.fa"
    tree_dir = db_dir / "tree"

    if not manifest.is_file():
        print("manifest.csv not found: %s" % manifest, file=sys.stderr)
        print("Run install_manifest.py, or copy it from the cleaned library.",
              file=sys.stderr)
        return 2
    if not prot_fa.is_file():
        print("protein library vmr_l.fa not found: %s" % prot_fa, file=sys.stderr)
        print("Run build_clean_db.py first.", file=sys.stderr)
        return 2

    out_path = args.out or (tree_dir / "vmr_L_aln.fasta")

    # 1) One representative per species, L segment only
    with manifest.open(encoding="utf-8-sig", newline="") as fh:
        rows = [r for r in csv.DictReader(fh)
                if (r.get("segment") or "").upper().startswith("L")]
    best = {}
    for r in rows:
        sp = (r.get("species") or "").strip()
        if not sp:
            continue
        cur = best.get(sp)
        if cur is None or (is_complete(r), nt(r)) > (is_complete(cur), nt(cur)):
            best[sp] = r

    ref_ids = [r["acc_base"].strip() for r in best.values()
               if (r.get("acc_base") or "").strip()]
    print("%d species -> %d representatives" % (len(best), len(ref_ids)))

    # 2) Pull the matching proteins. The reference must be protein.
    prots = read_fasta_dict(prot_fa)
    picked = []
    missing = []
    for acc in ref_ids:
        seq = prots.get(acc)
        if seq:
            picked.append((acc, seq))
        else:
            missing.append(acc)
    print("protein library hits: %d" % len(picked))
    if missing:
        print("[warning] %d representatives have no protein in the library "
              "(prodigal predicted no ORF), skipped: %s%s"
              % (len(missing), ", ".join(missing[:8]),
                 " ..." if len(missing) > 8 else ""), file=sys.stderr)
    if len(picked) < 2:
        print("fewer than 2 usable references, cannot align", file=sys.stderr)
        return 2

    # 3) Write the unaligned proteins (kept as a record and used as MAFFT input)
    raw_path = tree_dir / "vmr_L_ref_prot.fa"
    if not args.dry_run:
        tree_dir.mkdir(parents=True, exist_ok=True)
        with raw_path.open("w", encoding="utf-8") as fh:
            for acc, seq in picked:
                fh.write(">%s\n%s\n" % (acc, seq))
    lens = sorted(len(s) for _, s in picked)
    print("protein length: min %d / median %d / max %d aa"
          % (lens[0], lens[len(lens) // 2], lens[-1]))

    # 4) Align the proteins with MAFFT
    cmd = [args.mafft, "--thread", str(args.threads)] + args.mafft_extra.split() + [str(raw_path)]
    print("running: " + " ".join(cmd))
    if args.dry_run:
        print("")
        print("dry run: MAFFT was not executed and %s was not written" % out_path)
        return 0

    with out_path.open("w", encoding="utf-8") as out:
        code = subprocess.run(cmd, stdout=out).returncode
    if code:
        print("mafft failed (rc=%d)" % code, file=sys.stderr)
        return 1

    # 5) Self-checks: record count, equal length (a real alignment), protein alphabet
    got = read_fasta_dict(out_path)
    lengths = {len(s) for s in got.values()}
    print("output: %s" % out_path)
    print("  %d records / alignment length %s" % (len(got), sorted(lengths)))
    if len(got) != len(picked):
        print("[warning] output has %d records but %d were written"
              % (len(got), len(picked)), file=sys.stderr)
    if len(lengths) != 1:
        print("[warning] sequences differ in length, this is not a valid alignment",
              file=sys.stderr)
        return 1
    bad = set()
    for seq in got.values():
        bad |= {c for c in set(seq.upper()) if c not in "ACDEFGHIKLMNPQRSTVWY-XBZUO*"}
    if bad:
        print("[warning] non-protein characters found: %s" % "".join(sorted(bad)),
              file=sys.stderr)
        return 1
    print("  checks passed: protein alphabet, equal-length alignment")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
