"""Microbenchmark for CRISPR physics kernels."""

from __future__ import annotations

import argparse
import time
from typing import Sequence

import numpy as np

from helix.crispr.model import CasSystem, CasSystemType, GuideRNA, PAMRule
from helix.crispr.physics import score_pairs_encoded_multi, create_crispr_physics, _encode_sequence_to_uint8


_BASES = np.array(list("ACGT"))


def _random_sequences(count: int, length: int, rng: np.random.Generator) -> list[str]:
    encoded = rng.integers(0, 4, size=(count, length))
    return ["".join(_BASES[row]) for row in encoded]


def _bench(func, *args, **kwargs) -> float:
    start = time.perf_counter()
    func(*args, **kwargs)
    return time.perf_counter() - start


def run_shapes(shapes: Sequence[tuple[int, int, int]], backend: str, seed: int) -> None:
    rng = np.random.default_rng(seed)
    cas = CasSystem(
        name="demo",
        system_type=CasSystemType.CAS9,
        pam_rules=[PAMRule(pattern="NGG")],
        cut_offset=3,
    )
    for g, n, l in shapes:
        guide_strs = _random_sequences(g, l, rng)
        guides = np.stack([_encode_sequence_to_uint8(seq) for seq in guide_strs])
        try:
            physics_list = [create_crispr_physics(cas, GuideRNA(sequence=seq), backend=backend) for seq in guide_strs]
        except RuntimeError as exc:
            print(f"[warn] backend '{backend}' unavailable: {exc}")
            return
        window_strs = _random_sequences(n, l, rng)
        windows = np.stack([_encode_sequence_to_uint8(seq) for seq in window_strs])
        pam = np.zeros((g, n), dtype=np.float32)
        elapsed = _bench(
            score_pairs_encoded_multi,
            physics_list,
            guides,
            windows,
            pam_penalties=pam,
        )
        total_pairs = g * n
        print(f"backend={backend} shape=({g},{n},{l}) -> {total_pairs/elapsed/1e6:.2f} MPairs/s")


def main() -> None:
    parser = argparse.ArgumentParser(description="CRISPR engine microbenchmark")
    parser.add_argument("--backend", default="cpu-reference", help="Physics backend to exercise")
    parser.add_argument(
        "--shapes",
        nargs="*",
        default=["1x512x20", "1x4096x20", "1x16384x20"],
        help="List of GxNxL triples",
    )
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    shapes = []
    for token in args.shapes:
        g, n, l = token.split("x")
        shapes.append((int(g), int(n), int(l)))
    run_shapes(shapes, args.backend, args.seed)


if __name__ == "__main__":
    main()
