"""HyperLogLog sketches (Dashing-style)."""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from typing import Sequence

DNA = "ACGT"


def _normalize(seq: str) -> str:
    return "".join(base for base in seq.upper() if base in DNA)


def _hash64(kmer: str) -> int:
    return int.from_bytes(hashlib.sha256(kmer.encode("utf-8")).digest()[:8], "little")


@dataclass(frozen=True)
class HLLSketch:
    p: int
    registers: tuple[int, ...]

    def estimate(self) -> float:
        m = 1 << self.p
        alpha = 0.7213 / (1 + 1.079 / m)
        harmonic = sum(2.0 ** (-register) for register in self.registers)
        raw = alpha * m * m / harmonic
        if raw <= 2.5 * m:
            zeros = self.registers.count(0)
            if zeros:
                return m * math.log(m / zeros)
        elif raw > (1 / 30) * (1 << 32):
            raw = - (1 << 32) * math.log(1 - raw / (1 << 32))
        return raw

    def to_dict(self) -> dict:
        return {"p": self.p, "registers": list(self.registers), "estimate": self.estimate()}


def compute_hll(seq: str, k: int = 21, p: int = 10) -> HLLSketch:
    sequence = _normalize(seq)
    if len(sequence) < k:
        raise ValueError("Sequence shorter than k")
    m = 1 << p
    registers = [0] * m
    bits = 64
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i : i + k]
        h = _hash64(kmer)
        idx = h & (m - 1)
        w = h >> p
        zero_count = bits - p - w.bit_length() + 1 if w else bits - p + 1
        registers[idx] = max(registers[idx], zero_count)

    return HLLSketch(p=p, registers=tuple(registers))


def union_hll(a: HLLSketch, b: HLLSketch) -> HLLSketch:
    if a.p != b.p:
        raise ValueError("HLL sketches must share the same p parameter")
    registers = tuple(max(x, y) for x, y in zip(a.registers, b.registers))
    return HLLSketch(p=a.p, registers=registers)
