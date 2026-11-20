#!/usr/bin/env python3
"""Shared helpers for sequence identifier handling.

Background
----------
Prodigal writes protein headers in a fixed shape::

    ><original id>_<gene> # <start> # <end> # <strand> # ID=<index>_<gene>;partial=...
    >PX111623.2_1 # 1 # 6453 # 1 # ID=1_1;partial=11;start_type=Edge;...

The **first number after `ID=` is the position of that sequence in the input
FASTA, counting from 1**; the second is the gene number within that sequence.
That index recovers the original identifier with no guessing at all.

What the old code got wrong::

    re.sub(r"_\\d+$", "", "NC_055171_1")  -> "NC_055171"   # correct
    re.sub(r"_\\d+$", "", "NC_055171")    -> "NC"          # wrong: it ate the accession digits

So "does this carry a prodigal suffix" cannot be decided by string inspection.
It must be decided from a known identifier list, or from the index prodigal
writes into the header. This module provides three functions:
"""
from __future__ import annotations

import re
from pathlib import Path

# Strip one layer only: the greedy (.*) anchors the match on the rightmost _<digits>
_ONE_SUFFIX = re.compile(r"^(.*)_(\d+)$")
# prodigal's ID=<index>_<gene>
_PRODIGAL_ID = re.compile(r"(?:^|[;\s])ID=(\d+)_(\d+)(?:[;\s]|$)")


def read_fasta_ids(path) -> list[str]:
    """Read FASTA identifiers in file order (the first token after `>`)."""
    ids: list[str] = []
    with Path(path).open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                parts = line[1:].strip().split(None, 1)
                if parts:
                    ids.append(parts[0])
    return ids


def strip_prodigal_suffix(sid: str, known_ids=None) -> str:
    """Strip prodigal's `_<gene>` suffix, at most one layer.

    Parameters
    ----------
    sid : the identifier to process.
    known_ids : the known set of original identifiers (a set, the keys of a
        dict, or a list all work). When given, the decision order is:

        1. `sid` is itself in the list   -> return it unchanged (blocks NC_055171 -> NC)
        2. stripping one layer hits      -> return the stripped form (normal prodigal suffix)
        3. neither is in the list        -> return it unchanged (conservative, never mangles)

        When None, one layer is stripped unconditionally: the same result as the
        old code, but only one layer.

    Examples
    --------
    >>> strip_prodigal_suffix("NC_055171", {"NC_055171"})
    'NC_055171'
    >>> strip_prodigal_suffix("NC_055171_1", {"NC_055171"})
    'NC_055171'
    >>> strip_prodigal_suffix("k141_1118_1", {"k141_1118"})
    'k141_1118'
    """
    s = (sid or "").strip()
    if not s:
        return ""

    if known_ids is not None:
        try:
            hit_self = s in known_ids
        except TypeError:                      # an unhashable container was passed (a list, say)
            known = set(known_ids)
            hit_self = s in known
        else:
            known = known_ids
        if hit_self:
            return s
        m = _ONE_SUFFIX.match(s)
        if m and m.group(1) in known:
            return m.group(1)
        return s

    m = _ONE_SUFFIX.match(s)
    return m.group(1) if m else s


def prodigal_orig_id(header: str, orig_ids: list[str]) -> str:
    """Recover the original identifier from prodigal's `ID=<index>_<gene>`.

    Parameters
    ----------
    header : a prodigal header line, with or without the leading `>`.
    orig_ids : the identifiers of the FASTA given to prodigal, **in input order**.

    Returns
    -------
    The original identifier, or an empty string when the header cannot be parsed
    or the index is out of range. The caller decides the fallback.
    """
    m = _PRODIGAL_ID.search(header or "")
    if not m:
        return ""
    idx = int(m.group(1))
    if 1 <= idx <= len(orig_ids):
        return orig_ids[idx - 1]
    return ""
