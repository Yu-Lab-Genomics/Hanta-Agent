#!/usr/bin/env python3
"""Rebuild the DIAMOND protein databases from the cleaned nucleotide library.

Pipeline: prodigal translation -> keep the longest protein per record -> makedb.

Where things live (all under <db_dir>, which defaults to <project_dir>/DB):

    <db>/clean/vmr_{L,M,S}_clean.fasta          input, supplied externally
    <db>/clean/vmr_{L,M,S}_clean_proteins.faa   prodigal output, kept on disk
    <db>/clean/vmr_{L,M,S}_clean_genes.ffn
    <db>/diamond/vmr_{l,m,s}.fa                 OUTPUT, consumed by the pipeline
    <db>/diamond/vmr_{l,m,s}.dmnd               OUTPUT, consumed by the pipeline

The DIAMOND indexes are written straight into diamond/ because that is where
2-species_class.py and 4-vmr-tree.py look for them. Writing them anywhere else
would mean copying files by hand, which is exactly the mismatch this script now
avoids.

Usage
-----
    python build_clean_db.py                       # config.yaml -> project_dir/DB
    python build_clean_db.py --db-dir /data/DB     # explicit database directory
    python build_clean_db.py --config my.yaml      # explicit configuration
    python build_clean_db.py --dry-run             # print what would run
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from seqid import prodigal_orig_id, read_fasta_ids, strip_prodigal_suffix  # noqa: E402

SEGMENTS = [("L", "l"), ("M", "m"), ("S", "s")]
CONFIG_NAMES = ("config.yaml", "config-cell-paper.yaml", "cell_submit-paper.yaml")


def find_default_config(script_dir: Path):
    """Walk up from the script looking for a configuration file."""
    for d in (script_dir, script_dir.parent,
              script_dir.parent.parent, script_dir.parent.parent.parent):
        for name in CONFIG_NAMES:
            cand = d / name
            if cand.is_file():
                return cand
    return None


def read_project_dir(config_path: Path) -> str:
    """Read project_dir from a simple two-level YAML file.

    Only the subset this pipeline uses is supported, so PyYAML is not required.
    """
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


def run(cmd, dry_run):
    print("  " + " ".join(str(c) for c in cmd[:6]) + " ...")
    if dry_run:
        return
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def longest_protein(nuc_fa, proteins, out_fa):
    """Prodigal emits several proteins per nucleotide record; keep the longest.

    The original identifier is recovered in this priority order:
      1) `ID=<index>_<gene>` in the prodigal header - the index refers to the
         line order of the input FASTA, so this is guess-free and the most
         reliable route;
      2) fall back to strip_prodigal_suffix, validated against the known
         nucleotide identifiers. Only one suffix layer is stripped, so
         NC_055171 is not truncated to NC.

    Returns (proteins written, nucleotide records seen by prodigal, total input
    nucleotide records).
    """
    orig_ids = read_fasta_ids(nuc_fa)
    known = set(orig_ids)
    best = {}
    seen = set()
    cur_key, cur_seq = None, []

    def flush():
        if cur_key is None:
            return
        seq = "".join(cur_seq)
        seen.add(cur_key)
        if cur_key not in best or len(seq) > best[cur_key][0]:
            best[cur_key] = (len(seq), seq)

    with proteins.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                head = line[1:].strip()
                key = prodigal_orig_id(head, orig_ids)
                if not key:
                    token = head.split(None, 1)[0] if head else ""
                    key = strip_prodigal_suffix(token, known)
                cur_key = key or None
                cur_seq = []
            elif cur_key is not None:
                cur_seq.append(line.strip())
    flush()

    with out_fa.open("w", encoding="utf-8") as fh:
        for key, (_, seq) in best.items():
            fh.write(">" + key + "\n" + seq + "\n")
    return len(best), len(seen), len(orig_ids)


def main():
    ap = argparse.ArgumentParser(
        description="Rebuild the DIAMOND protein databases (prodigal + longest "
                    "protein + makedb) and write them to <db>/diamond/")
    ap.add_argument("--db-dir", type=Path, default=None,
                    help="database directory; defaults to <project_dir>/DB from the "
                         "configuration")
    ap.add_argument("--config", type=str, default="",
                    help="path to config.yaml; used to derive <project_dir>/DB")
    ap.add_argument("--prodigal", default="prodigal",
                    help="prodigal executable (default: found on PATH)")
    ap.add_argument("--diamond", default="diamond",
                    help="diamond executable (default: found on PATH)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the composed commands without running anything")
    args = ap.parse_args()

    db_dir = resolve_db_dir(args)
    if db_dir is None:
        print("Cannot determine the database directory. Pass --db-dir, or make "
              "sure project_dir is set in config.yaml.", file=sys.stderr)
        return 2

    clean = db_dir / "clean"
    diamond = db_dir / "diamond"
    if not clean.is_dir():
        print("cleaned library directory not found: %s" % clean, file=sys.stderr)
        print("This directory holds vmr_{L,M,S}_clean.fasta and must be supplied "
              "before the databases can be rebuilt.", file=sys.stderr)
        return 2

    if not args.dry_run:
        diamond.mkdir(parents=True, exist_ok=True)

    print("database directory : %s" % db_dir)
    print("nucleotide library : %s" % clean)
    print("DIAMOND output     : %s" % diamond)
    print("")

    for seg, low in SEGMENTS:
        nuc = clean / ("vmr_%s_clean.fasta" % seg)
        if not nuc.is_file():
            print("input missing: %s" % nuc, file=sys.stderr)
            return 2
        prot = clean / ("vmr_%s_clean_proteins.faa" % seg)
        run([args.prodigal, "-i", str(nuc), "-a", str(prot),
             "-d", str(clean / ("vmr_%s_clean_genes.ffn" % seg)),
             "-o", "/dev/null", "-p", "meta", "-q"], args.dry_run)

        out_fa = diamond / ("vmr_%s.fa" % low)
        if args.dry_run:
            print("  vmr_%s.fa would be written to %s" % (low, out_fa))
        else:
            n_best, n_seen, n_total = longest_protein(nuc, prot, out_fa)
            print("  vmr_%s.fa: %d proteins (prodigal covered %d/%d nucleotide records)"
                  % (low, n_best, n_seen, n_total))
            if n_seen < n_total:
                print("  note: %d nucleotide records yielded no protein "
                      "(short fragments, expected)" % (n_total - n_seen))

        run([args.diamond, "makedb", "--in", str(out_fa),
             "-d", str(diamond / ("vmr_%s" % low))], args.dry_run)
        print("  vmr_%s.dmnd built" % low)

    if args.dry_run:
        print("")
        print("dry run: nothing was executed and no file was written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
