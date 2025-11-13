"""Helix core package."""

from importlib import metadata

from . import bioinformatics, codon, cyclospectrum, triage, string, seed, motif, crispr, prime, pcr
from .api import dna_summary, triage_report, fold_rna, spectrum_leaderboard, protein_summary
from .crispr.model import (
    CasSystem,
    CasSystemType,
    DigitalGenome,
    GuideRNA,
    PAMRule,
    TargetSite,
)
from .crispr.simulator import CutEvent, find_candidate_sites, rank_off_targets, simulate_cuts
from .prime.model import PegRNA, PrimeEditOutcome, PrimeEditor
from .prime.simulator import locate_prime_target_site, simulate_prime_edit

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
    "pcr",
    "CasSystemType",
    "PAMRule",
    "CasSystem",
    "GuideRNA",
    "TargetSite",
    "DigitalGenome",
    "CutEvent",
    "find_candidate_sites",
    "simulate_cuts",
    "rank_off_targets",
    "prime",
    "PegRNA",
    "PrimeEditor",
    "PrimeEditOutcome",
    "locate_prime_target_site",
    "simulate_prime_edit",
    "dna_summary",
    "triage_report",
    "fold_rna",
    "spectrum_leaderboard",
    "protein_summary",
    "__version__",
]
