"""Deterministic minimizer extraction with canonical hashing."""
from __future__ import annotations

from collections import deque
from typing import List, Tuple

DNA2 = {"A": 0, "C": 1, "G": 2, "T": 3}
RC = {"A": "T", "C": "G", "G": "C", "T": "A"}


def _encode_kmer(kmer: str) -> int:
    v = 0
    for c in kmer:
        v = (v << 2) | DNA2[c]
    return v


def _revcomp(kmer: str) -> str:
    return "".join(RC[b] for b in reversed(kmer))


def _canon_hash(kmer: str) -> Tuple[int, bool]:
    forward = _encode_kmer(kmer)
    reverse = _encode_kmer(_revcomp(kmer))
    if forward <= reverse:
        return forward, False
    return reverse, True


def minimizers(seq: str, k: int, w: int) -> List[Tuple[int, str, int]]:
    """Return [(pos, kmer, hash)] minimizers using canonical k-mers over A/C/G/T."""
    seq = seq.upper()
    n = len(seq)
    if k <= 0 or w <= 0 or n < k or w < 1:
        return []

    km_hash: List[Tuple[int, str, int | None]] = []
    for i in range(n - k + 1):
        kmer = seq[i : i + k]
        if any(c not in DNA2 for c in kmer):
            km_hash.append((i, kmer, None))
        else:
            canon, _ = _canon_hash(kmer)
            km_hash.append((i, kmer, canon))

    out: List[Tuple[int, str, int]] = []
    dq: deque[int] = deque()
    for i, (pos, kmer, h) in enumerate(km_hash):
        if h is not None:
            while dq and km_hash[dq[-1]][2] is not None and km_hash[dq[-1]][2] >= h:
                dq.pop()
            dq.append(i)
        left = i - (w - 1)
        while dq and dq[0] < left:
            dq.popleft()
        if left >= 0 and dq:
            j = dq[0]
            pj, kj, hj = km_hash[j]
            if hj is not None:
                out.append((pj, kj, hj))
    return out
