"""Distances over sketches."""
from __future__ import annotations

import math

from .minhash import MinHashSketch


def mash_distance(sketch_a: MinHashSketch, sketch_b: MinHashSketch) -> float:
    if sketch_a.size != sketch_b.size or sketch_a.k != sketch_b.k:
        raise ValueError("Sketches must have the same size and k")
    matches = sum(1 for a, b in zip(sketch_a.hashes, sketch_b.hashes) if a == b)
    if matches == 0:
        return 1.0
    jaccard = matches / sketch_a.size
    return min(1.0, max(0.0, -1.0 / sketch_a.k * math.log(2 * jaccard / (1 + jaccard))))
