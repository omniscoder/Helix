"""Bit-parallel <=k mismatch finder."""
from __future__ import annotations

def k_mismatch_positions(text: str, pattern: str, k: int) -> list[int]:
    """Return start indices where pattern occurs in text with <= k mismatches."""
    text = text.upper()
    pattern = pattern.upper()
    n, m = len(text), len(pattern)
    if m == 0 or n < m:
        return []
    if m > 64:
        hits = []
        for i in range(0, n - m + 1):
            window = text[i : i + m]
            mismatches = sum(ch1 != ch2 for ch1, ch2 in zip(window, pattern))
            if mismatches <= k:
                hits.append(i)
        return hits
    masks = {base: 0 for base in "ACGTN"}
    for idx, base in enumerate(pattern):
        for key in masks:
            if base == key or base == "N" or key == "N":
                masks[key] |= 1 << idx
    hits = []
    for i in range(0, n - m + 1):
        window = text[i : i + m]
        mismatches = 0
        for j, ch in enumerate(window):
            mask = masks.get(ch, 0)
            if not (mask & (1 << j)):
                mismatches += 1
                if mismatches > k:
                    break
        if mismatches <= k:
            hits.append(i)
    return hits
