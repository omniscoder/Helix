"""MinHash sketches for sequences."""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Iterable, List, Sequence, Tuple

DNA = "ACGT"


def _normalize(seq: str) -> str:
    return "".join(base for base in seq.upper() if base in DNA)


def _hash64(kmer: str, seed: int) -> int:
    data = f"{seed}|{kmer}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(data).digest()[:8], "little")


@dataclass(frozen=True)
class MinHashSketch:
    k: int
    size: int
    hashes: Tuple[int, ...]

    def to_dict(self) -> dict:
        return {"k": self.k, "size": self.size, "hashes": list(self.hashes)}


def compute_sketch(seq: str, k: int = 21, sketch_size: int = 1000, seeds: Sequence[int] | None = None) -> MinHashSketch:
    sequence = _normalize(seq)
    if len(sequence) < k:
        raise ValueError("Sequence shorter than k")
    seeds = seeds or list(range(sketch_size))
    hashes = []
    for seed in seeds:
        min_hash = min(
            _hash64(sequence[i : i + k], seed)
            for i in range(len(sequence) - k + 1)
        )
        hashes.append(min_hash)
    return MinHashSketch(k=k, size=len(seeds), hashes=tuple(hashes))
