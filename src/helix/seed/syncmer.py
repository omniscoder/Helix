"""Syncmer sampling (open syncmers with canonical hashing)."""
from __future__ import annotations

from typing import List, Tuple

from .minimizer import DNA2, _canon_hash


def syncmers(seq: str, k: int, s: int, position: int | None = None) -> List[Tuple[int, str, int]]:
    """Return syncmers as [(pos, kmer, hash)].

    Uses an open-syncmer scheme: choose the k-mers whose minimum s-mer occurs at `position`.
    If `position` is None, use the center ( (k-s)//2 ).
    """
    seq = seq.upper()
    n = len(seq)
    if k <= 0 or s <= 0 or s > k or n < k:
        return []
    if position is None:
        position = (k - s) // 2
    if not 0 <= position <= k - s:
        raise ValueError("position must be between 0 and k-s")

    out: List[Tuple[int, str, int]] = []
    for i in range(n - k + 1):
        kmer = seq[i : i + k]
        if any(c not in DNA2 for c in kmer):
            continue
        s_hashes = []
        for offset in range(k - s + 1):
            s_mer = kmer[offset : offset + s]
            if any(c not in DNA2 for c in s_mer):
                s_hashes.append((offset, None))
            else:
                val, _ = _canon_hash(s_mer)
                s_hashes.append((offset, val))
        valid = [item for item in s_hashes if item[1] is not None]
        if not valid:
            continue
        min_val = min(val for _, val in valid if val is not None)
        if any(offset == position and val == min_val for offset, val in valid):
            k_hash, _ = _canon_hash(kmer)
            out.append((i, kmer, k_hash))
    return out
