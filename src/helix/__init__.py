"""Helix core package."""

from importlib import metadata

from . import bioinformatics, codon, cyclospectrum, triage, string, seed, motif, crispr, prime
from .api import dna_summary, triage_report, fold_rna, spectrum_leaderboard, protein_summary

try:  # pragma: no cover - metadata only at runtime
    __version__ = metadata.version("veri-helix")
except metadata.PackageNotFoundError:  # pragma: no cover - source tree / editable installs
    __version__ = "0.0.0"

__all__ = [
    "bioinformatics",
    "codon",
    "cyclospectrum",
    "triage",
    "string",
    "seed",
    "motif",
    "crispr",
    "prime",
    "dna_summary",
    "triage_report",
    "fold_rna",
    "spectrum_leaderboard",
    "protein_summary",
    "__version__",
]
