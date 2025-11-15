"""Shim module for backward-compatible `import helix_api`.

This lives under `src/` so it is importable when tests add `src` to `PYTHONPATH`.
It simply re-exports helpers from `helix.api`.
"""
from __future__ import annotations

from helix.api import (
    PROTEIN_AVAILABLE,
    dna_summary,
    fold_rna,
    protein_summary,
    spectrum_leaderboard,
    triage_report,
)

__all__ = [
    "dna_summary",
    "triage_report",
    "fold_rna",
    "spectrum_leaderboard",
    "protein_summary",
    "PROTEIN_AVAILABLE",
]

