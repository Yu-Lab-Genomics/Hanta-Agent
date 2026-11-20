#!/usr/bin/env python3
"""Normalise hantavirus species names to one canonical spelling.

Rules (option A: traditional names are kept as separate species):
  1) Quoted candidate name -> traditional name
     (e.g. Orthohantavirus sp. 'seewisense' -> Seewis virus)
  2) Legacy genus corrected: Hantavirus sp. / Hantavirus sp. XXX -> Orthohantavirus sp.
  3) Genus without the sp. suffix: Orthohantavirus -> Orthohantavirus sp.
  4) Hantavirus + strain name (an unclassified strain) -> Orthohantavirus sp.
  5) Everything else (binomials, traditional names, family-level sp.) is kept as is

Usage:
  from species_norm import normalize_species
  norm = normalize_species("Orthohantavirus")  # -> "Orthohantavirus sp."
"""
from __future__ import annotations

import re

# The eight official ICTV genera of Hantaviridae
GENERA = [
    "Orthohantavirus", "Mobatvirus", "Thottimvirus", "Loanvirus",
    "Actinovirus", "Percilovirus", "Agnathovirus", "Reptillovirus",
]

# Quoted candidate name -> traditional name (NCBI taxonomy candidate spelling
# mapped to the traditional name)
CANDIDATE_ALIAS = {
    "Orthohantavirus sp. 'seewisense'": "Seewis virus",
    "Orthohantavirus sp. 'yakeshiense'": "Yakeshi virus",
    "Orthohantavirus sp. 'fusongense'": "Fusong virus",
    "Orthohantavirus sp. 'necocliense'": "Necocli virus",
    "Orthohantavirus sp. 'oxbowense'": "Oxbow virus",
    # chicoanense has no traditional name; its candidate spelling is the
    # canonical one, so it is deliberately absent from this table
}

# Legacy genus (called Hantavirus before 2019, Orthohantavirus now)
_HANTAVIRUS_STRAIN_RE = re.compile(r"^Hantavirus(?:\s+sp\.?\s*.*)?$")


def normalize_species(species: str) -> str:
    """Normalise a species name to its single canonical spelling. Empty input
    is returned unchanged."""
    sp = (species or "").strip()
    if not sp:
        return sp

    # 1) quoted candidate name -> traditional name
    if sp in CANDIDATE_ALIAS:
        return CANDIDATE_ALIAS[sp]

    # 2) legacy genus: Hantavirus sp. / Hantavirus sp. XXX -> Orthohantavirus sp.
    if sp == "Hantavirus sp." or sp.startswith("Hantavirus sp. "):
        return "Orthohantavirus sp."

    # 3) genus given without the sp. suffix
    if sp in GENERA:
        return sp + " sp."

    # 4) Hantavirus + strain name (an unclassified strain) -> Orthohantavirus sp.
    if sp.startswith("Hantavirus ") or sp.startswith("Hantavirus\t"):
        return "Orthohantavirus sp."

    # 5) everything else is kept (binomials, traditional names, family-level sp.)
    return sp


if __name__ == "__main__":
    tests = [
        "Orthohantavirus",                     # -> Orthohantavirus sp.
        "Hantavirus sp.",                      # -> Orthohantavirus sp.
        "Hantavirus sp. Paranoa",              # -> Orthohantavirus sp.
        "Hantavirus A9",                       # -> Orthohantavirus sp.
        "Thottimvirus",                        # -> Thottimvirus sp.
        "Loanvirus",                           # -> Loanvirus sp.
        "Mobatvirus",                          # -> Mobatvirus sp.
        "Orthohantavirus sp.",                 # kept as is
        "Hantaviridae sp.",                    # kept as is
        "Orthohantavirus sp. 'seewisense'",    # -> Seewis virus
        "Orthohantavirus sp. 'yakeshiense'",   # -> Yakeshi virus
        "Orthohantavirus sp. 'chicoanense'",   # kept as is
        "Orthohantavirus seoulense",           # kept as is
        "Huangpi virus",                       # kept as is
        "Orthohantavirus maporalense",         # kept as is
    ]
    for t in tests:
        print("%-40s -> %s" % (t, normalize_species(t)))
