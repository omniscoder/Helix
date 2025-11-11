"""Simple seed-and-extend alignment with band + x-drop."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict


@dataclass
class SeedMatch:
    ref_pos: int
    read_pos: int
    length: int


def _score(a: str, b: str) -> int:
    return 2 if a == b else -2


def extend_alignment(
    seed: SeedMatch,
    ref: str,
    read: str,
    band: int = 64,
    xdrop: int = 10,
) -> Dict[str, object]:
    """Perform a banded Smith-Waterman extension around the seed."""
    ref = ref.upper()
    read = read.upper()

    window = max(band, seed.length + band)
    ref_start = max(0, seed.ref_pos - window)
    ref_end = min(len(ref), seed.ref_pos + seed.length + window)
    read_start = max(0, seed.read_pos - window)
    read_end = min(len(read), seed.read_pos + seed.length + window)

    ref_slice = ref[ref_start:ref_end]
    read_slice = read[read_start:read_end]
    n, m = len(ref_slice), len(read_slice)

    gap_penalty = -3
    best_score = 0
    best_pos = (0, 0)

    dp = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[0] * (m + 1) for _ in range(n + 1)]  # 1 diag, 2 up, 3 left

    for i in range(1, n + 1):
        j_min = max(1, i - band)
        j_max = min(m, i + band)
        for j in range(j_min, j_max + 1):
            match = dp[i - 1][j - 1] + _score(ref_slice[i - 1], read_slice[j - 1])
            delete = dp[i - 1][j] + gap_penalty
            insert = dp[i][j - 1] + gap_penalty
            best = max(0, match, delete, insert)
            dp[i][j] = best
            if best == 0:
                trace[i][j] = 0
            elif best == match:
                trace[i][j] = 1
            elif best == delete:
                trace[i][j] = 2
            else:
                trace[i][j] = 3

            if best > best_score:
                best_score = best
                best_pos = (i, j)
        row_best = max(dp[i])
        if best_score - row_best > xdrop and i > band:
            break

    i, j = best_pos
    cigar = []
    matches = 0
    while i > 0 and j > 0 and trace[i][j] != 0:
        direction = trace[i][j]
        if direction == 1:
            op = "M"
            i -= 1
            j -= 1
            matches += 1
        elif direction == 2:
            op = "D"
            i -= 1
        else:
            op = "I"
            j -= 1
        if cigar and cigar[-1][0] == op:
            cigar[-1] = (op, cigar[-1][1] + 1)
        else:
            cigar.append((op, 1))
    cigar_str = "".join(f"{length}{op}" for op, length in reversed(cigar))

    return {
        "ref_start": ref_start + i,
        "ref_end": ref_start + best_pos[0],
        "read_start": read_start + j,
        "read_end": read_start + best_pos[1],
        "score": best_score,
        "cigar": cigar_str,
        "matches": matches,
    }
