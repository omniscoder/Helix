"""Peptide spectrum utilities for Helix."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

from amino_acids import mass_index_table


def peptide_mass(peptide: str) -> int:
    """Return the integer mass of a peptide string."""
    try:
        return sum(mass_index_table[aa] for aa in peptide.upper())
    except KeyError as exc:
        raise ValueError(f"Unknown residue '{exc.args[0]}' in peptide '{peptide}'.") from exc


def _prefix_masses(peptide: str) -> List[int]:
    masses = [0]
    running = 0
    for aa in peptide.upper():
        try:
            running += mass_index_table[aa]
        except KeyError as exc:
            raise ValueError(f"Unknown residue '{exc.args[0]}' in peptide '{peptide}'.") from exc
        masses.append(running)
    return masses


def linear_spectrum(peptide: str) -> List[int]:
    """Return the theoretical linear spectrum for a peptide."""
    prefix = _prefix_masses(peptide)
    spectrum = [0]
    for i in range(len(prefix)):
        for j in range(i + 1, len(prefix)):
            spectrum.append(prefix[j] - prefix[i])
    return sorted(spectrum)


def cyclic_spectrum(peptide: str) -> List[int]:
    """Return the theoretical cyclic spectrum for a peptide."""
    if not peptide:
        return [0]
    prefix = _prefix_masses(peptide)
    peptide_total = prefix[-1]
    spectrum = [0, peptide_total]
    for i in range(1, len(prefix) - 1):
        for j in range(i + 1, len(prefix)):
            mass = prefix[j] - prefix[i]
            spectrum.append(mass)
            spectrum.append(peptide_total - mass)
    return sorted(spectrum)


def theoretical_spectrum(peptide: str, *, cyclic: bool = True) -> List[int]:
    """Return either the cyclic or linear theoretical spectrum."""
    return cyclic_spectrum(peptide) if cyclic else linear_spectrum(peptide)


def score_peptide(
    peptide: str,
    experimental_spectrum: Sequence[int],
    *,
    cyclic: bool = True,
) -> int:
    """Score a peptide by the shared masses with an experimental spectrum."""
    theo = Counter(theoretical_spectrum(peptide, cyclic=cyclic))
    exp = Counter(experimental_spectrum)
    score = 0
    for mass, count in theo.items():
        score += min(count, exp.get(mass, 0))
    return score


def _is_consistent(peptide: str, experimental_counter: Counter[int]) -> bool:
    """Check whether the peptide's linear spectrum fits inside the experimental one."""
    linear = Counter(linear_spectrum(peptide))
    for mass, count in linear.items():
        if count > experimental_counter.get(mass, 0):
            return False
    return True


def _trim_leaderboard(entries: List[Tuple[str, int]], limit: int) -> List[Tuple[str, int]]:
    if len(entries) <= limit:
        return entries
    entries.sort(key=lambda item: item[1], reverse=True)
    cutoff_score = entries[limit - 1][1]
    return [entry for entry in entries if entry[1] >= cutoff_score]


def leaderboard_cyclopeptide_sequencing(
    experimental_spectrum: Sequence[int],
    *,
    leaderboard_size: int = 10,
    alphabet: Iterable[str] | None = None,
) -> List[Tuple[str, int]]:
    """Return high-scoring peptides for an experimental spectrum."""
    if leaderboard_size <= 0:
        raise ValueError("leaderboard_size must be > 0")

    letters = list(alphabet) if alphabet else sorted(mass_index_table.keys())
    parent_mass = max(experimental_spectrum)
    experimental_counter = Counter(experimental_spectrum)

    leaderboard: List[Tuple[str, int]] = [("", 0)]
    best: List[Tuple[str, int]] = []
    best_score = 0

    while leaderboard:
        expanded: List[str] = []
        for peptide, _ in leaderboard:
            for letter in letters:
                expanded.append(peptide + letter)

        leaderboard = []
        for peptide in expanded:
            mass = peptide_mass(peptide)
            if mass == parent_mass:
                score = score_peptide(peptide, experimental_spectrum, cyclic=True)
                leaderboard.append((peptide, score))
                if score > best_score:
                    best_score = score
                    best = [(peptide, score)]
                elif score == best_score:
                    best.append((peptide, score))
            elif mass < parent_mass and _is_consistent(peptide, experimental_counter):
                score = score_peptide(peptide, experimental_spectrum, cyclic=False)
                leaderboard.append((peptide, score))

        leaderboard = _trim_leaderboard(leaderboard, leaderboard_size)

    if not best:
        best = leaderboard
    # Deduplicate peptides while preserving score ordering.
    dedup: Dict[str, int] = {}
    for peptide, score in sorted(best, key=lambda item: item[1], reverse=True):
        dedup.setdefault(peptide, score)
    return [(pep, dedup[pep]) for pep in dedup]
