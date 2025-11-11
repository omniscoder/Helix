"""FM-index utilities for exact substring search."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List


def _ensure_sentinel(text: str) -> str:
    text = text.upper()
    if "$" in text:
        raise ValueError("Input text must not contain '$'.")
    return text + "$"


def build_bwt(text: str) -> tuple[str, List[int], str]:
    """Return (BWT string, suffix array, augmented text) for `text`."""
    augmented = _ensure_sentinel(text)
    sa = sorted(range(len(augmented)), key=lambda i: augmented[i:])
    bwt_chars = []
    for idx in sa:
        prev = augmented[idx - 1] if idx > 0 else augmented[-1]
        bwt_chars.append(prev)
    return "".join(bwt_chars), sa, augmented


@dataclass(frozen=True)
class FMIndex:
    text: str
    bwt: str
    suffix_array: List[int]
    c_table: Dict[str, int]
    occ_table: Dict[str, List[int]]


def build_fm(text: str) -> FMIndex:
    """Build an FM-index for `text`."""
    bwt, sa, augmented = build_bwt(text)
    alphabet = sorted(set(bwt))
    counts: Dict[str, int] = {ch: 0 for ch in alphabet}
    c_table: Dict[str, int] = {}
    total = 0
    for ch in alphabet:
        c_table[ch] = total
        total += bwt.count(ch)

    occ_table: Dict[str, List[int]] = {ch: [0] * (len(bwt) + 1) for ch in alphabet}
    for i, ch in enumerate(bwt, start=1):
        for sym in alphabet:
            occ_table[sym][i] = occ_table[sym][i - 1]
        occ_table[ch][i] += 1

    return FMIndex(text=augmented, bwt=bwt, suffix_array=sa, c_table=c_table, occ_table=occ_table)


def rank(index: FMIndex, char: str, pos: int) -> int:
    """Return the number of occurrences of `char` in BWT[:pos]."""
    if char not in index.occ_table:
        return 0
    pos = max(0, min(pos, len(index.bwt)))
    return index.occ_table[char][pos]


def select(index: FMIndex, char: str, nth: int) -> int | None:
    """Return the position of the nth (0-indexed) occurrence of `char` in the BWT."""
    if char not in index.occ_table or nth < 0:
        return None
    target = nth + 1
    counts = index.occ_table[char]
    for i in range(1, len(counts)):
        if counts[i] == target:
            return i - 1
    return None


def search(index: FMIndex, pattern: str) -> Dict[str, object]:
    """Return JSON-friendly search results for `pattern`."""
    pattern = pattern.upper()
    if not pattern:
        return {"pattern": pattern, "count": len(index.text) - 1, "positions": list(range(len(index.text) - 1))}

    l, r = 0, len(index.bwt)
    for ch in reversed(pattern):
        if ch not in index.c_table:
            return {"pattern": pattern, "count": 0, "positions": []}
        l = index.c_table[ch] + rank(index, ch, l)
        r = index.c_table[ch] + rank(index, ch, r)
        if l >= r:
            return {"pattern": pattern, "count": 0, "positions": []}

    positions = [index.suffix_array[i] for i in range(l, r)]
    positions = [pos for pos in positions if pos < len(index.text) - 1]
    positions.sort()
    return {"pattern": pattern, "count": len(positions), "positions": positions}
