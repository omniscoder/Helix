"""Shared utilities for Helix motif discovery."""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Sequence

BASES = ["A", "C", "G", "T"]
EPSILON = 1e-12


def normalize_sequence(seq: str) -> str:
    return "".join(base for base in seq.upper() if base in BASES)


def initial_pwm(width: int) -> List[Dict[str, float]]:
    return [{base: 1.0 / len(BASES) for base in BASES} for _ in range(width)]


def compute_log_prob(kmer: str, pwm: Sequence[Dict[str, float]]) -> float:
    return sum(math.log(max(pwm[i].get(base, EPSILON), EPSILON)) for i, base in enumerate(kmer))


@dataclass
class MotifResult:
    pwm: List[Dict[str, float]]
    log_likelihood: float
    solver: str

    def consensus(self) -> str:
        return "".join(max(column.items(), key=lambda item: item[1])[0] for column in self.pwm)

    def as_json(self) -> Dict[str, object]:
        return {
            "pwm": self.pwm,
            "log_likelihood": self.log_likelihood,
            "consensus": self.consensus(),
            "solver": self.solver,
        }
