from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cyclospectrum import (
    cyclic_spectrum,
    leaderboard_cyclopeptide_sequencing,
    linear_spectrum,
    score_peptide,
)


def is_rotation(candidate: str, reference: str) -> bool:
    if len(candidate) != len(reference):
        return False
    doubled = reference * 2
    return candidate in doubled


def test_linear_and_cyclic_spectra_match_rosalind_example():
    peptide = "NQEL"
    assert linear_spectrum(peptide) == [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]
    assert cyclic_spectrum(peptide) == [
        0,
        113,
        114,
        128,
        129,
        227,
        242,
        242,
        257,
        355,
        356,
        370,
        371,
        484,
    ]


def test_score_peptide_counts_shared_masses():
    peptide = "NQEL"
    spectrum = cyclic_spectrum(peptide)
    expected = len(spectrum)
    assert score_peptide(peptide, spectrum) == expected
    truncated = spectrum[:-2]
    assert score_peptide(peptide, truncated) == expected - 2


def test_leaderboard_recovers_peptide_from_perfect_spectrum():
    peptide = "NQEL"
    spectrum = cyclic_spectrum(peptide)
    hits = leaderboard_cyclopeptide_sequencing(spectrum, leaderboard_size=5)
    candidates = [entry[0] for entry in hits]
    assert any(is_rotation(candidate, peptide) for candidate in candidates)
