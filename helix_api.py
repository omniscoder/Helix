"""Backwards-compatible shim for the legacy `helix_api` helpers.

Tests and examples historically imported a top-level `helix_api` module.
The modern implementations live in `helix.api`, so we re-export them here.
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

