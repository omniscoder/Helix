"""Annotated Nussinov RNA folding implementation for Helix."""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

import numpy as np

import bioinformatics

CANONICAL_PAIRS = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
WOBBLE_PAIRS = {("G", "U"), ("U", "G")}


@dataclass(frozen=True)
class NussinovResult:
    sequence: str
    dp_table: np.ndarray
    pairs: List[Tuple[int, int]]
    structure: str

    def score(self) -> int:
        return int(self.dp_table[0, len(self.sequence) - 1]) if self.sequence else 0


def _prepare_sequence(raw: str) -> str:
    dna = bioinformatics.normalize_sequence(raw)
    return dna.replace("T", "U")


def _can_pair(base_a: str, base_b: str, allow_wobble: bool) -> bool:
    pair = (base_a, base_b)
    if pair in CANONICAL_PAIRS:
        return True
    return allow_wobble and pair in WOBBLE_PAIRS


def _traceback(
    dp: np.ndarray,
    sequence: str,
    i: int,
    j: int,
    min_loop_length: int,
    allow_wobble: bool,
    pairs: List[Tuple[int, int]],
) -> None:
    if i >= j:
        return
    if dp[i, j] == dp[i + 1, j]:
        _traceback(dp, sequence, i + 1, j, min_loop_length, allow_wobble, pairs)
        return
    if dp[i, j] == dp[i, j - 1]:
        _traceback(dp, sequence, i, j - 1, min_loop_length, allow_wobble, pairs)
        return
    if (
        j - i - 1 >= min_loop_length
        and _can_pair(sequence[i], sequence[j], allow_wobble)
        and dp[i, j] == dp[i + 1, j - 1] + 1
    ):
        pairs.append((i, j))
        _traceback(dp, sequence, i + 1, j - 1, min_loop_length, allow_wobble, pairs)
        return

    for k in range(i + 1, j):
        if dp[i, j] == dp[i, k] + dp[k + 1, j]:
            _traceback(dp, sequence, i, k, min_loop_length, allow_wobble, pairs)
            _traceback(dp, sequence, k + 1, j, min_loop_length, allow_wobble, pairs)
            return


def _pairs_to_dot_bracket(length: int, pairs: Sequence[Tuple[int, int]]) -> str:
    structure = ["." for _ in range(length)]
    for i, j in pairs:
        structure[i] = "("
        structure[j] = ")"
    return "".join(structure)


def nussinov(
    sequence: str,
    *,
    min_loop_length: int = 3,
    allow_wobble_pairs: bool = True,
) -> NussinovResult:
    """Compute an RNA secondary structure using the Nussinov recurrence."""
    rna = _prepare_sequence(sequence)
    n = len(rna)
    dp = np.zeros((n, n), dtype=int)
    if n == 0:
        return NussinovResult(rna, dp, [], "")

    for span in range(1, n):
        for i in range(n - span):
            j = i + span
            best = max(dp[i + 1, j], dp[i, j - 1])
            if j - i - 1 >= min_loop_length and _can_pair(rna[i], rna[j], allow_wobble_pairs):
                best = max(best, dp[i + 1, j - 1] + 1)

            split_scores = [dp[i, k] + dp[k + 1, j] for k in range(i + 1, j)]
            if split_scores:
                best = max(best, max(split_scores))
            dp[i, j] = best

    pairs: List[Tuple[int, int]] = []
    _traceback(dp, rna, 0, n - 1, min_loop_length, allow_wobble_pairs, pairs)
    pairs.sort()
    structure = _pairs_to_dot_bracket(n, pairs)
    return NussinovResult(rna, dp, pairs, structure)


def fold_to_dot_bracket(sequence: str, **kwargs) -> str:
    """Convenience helper that returns only the dot-bracket string."""
    return nussinov(sequence, **kwargs).structure
