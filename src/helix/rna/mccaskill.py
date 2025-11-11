"""Simplified McCaskill partition function + MEA structures."""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List


BASE_PAIRS = {
    ("A", "U"): -1.0,
    ("U", "A"): -1.0,
    ("G", "C"): -2.0,
    ("C", "G"): -2.0,
    ("G", "U"): -0.5,
    ("U", "G"): -0.5,
}


@dataclass
class ThermoParams:

    temperature: float = 310.15  # Kelvin (37C)
    gas_constant: float = 0.001987  # kcal/(mol*K)

    @property
    def beta(self) -> float:
        return 1.0 / (self.gas_constant * self.temperature)


def _can_pair(a: str, b: str) -> bool:
    return (a, b) in BASE_PAIRS


def _pair_energy(a: str, b: str) -> float:
    return BASE_PAIRS.get((a, b), 0.0)


def _normalize_sequence(seq: str) -> str:
    return seq.upper().replace("T", "U")


def partition(seq: str, params: Dict[str, float] | None = None) -> Dict[str, object]:
    """Compute a simplified partition function + base-pair probabilities."""
    sequence = _normalize_sequence(seq)
    n = len(sequence)
    thermo = ThermoParams(**(params or {}))
    beta = thermo.beta

    Q = [[0.0 for _ in range(n)] for _ in range(n)]
    pair_weights = [[0.0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        Q[i][i] = 1.0

    def get_Q(i: int, j: int) -> float:
        if i > j:
            return 1.0
        return Q[i][j]

    for length in range(2, n + 1):
        for i in range(0, n - length + 1):
            j = i + length - 1
            total = get_Q(i, j - 1)
            for k in range(i, j):
                if _can_pair(sequence[k], sequence[j]):
                    left = get_Q(i, k - 1)
                    inside = get_Q(k + 1, j - 1)
                    energy = _pair_energy(sequence[k], sequence[j])
                    boltz = math.exp(-energy * beta)
                    contrib = left * inside * boltz
                    total += contrib
                    pair_weights[k][j] += contrib
            Q[i][j] = total

    Z = get_Q(0, n - 1)
    posterior = [[0.0 for _ in range(n)] for _ in range(n)]
    if Z > 0:
        for i in range(n):
            for j in range(i + 1, n):
                if pair_weights[i][j] > 0:
                    prob = pair_weights[i][j] / Z
                    posterior[i][j] = posterior[j][i] = prob

    free_energy = - (1 / beta) * math.log(Z) if Z > 0 else math.inf
    meta = {
        "sequence": sequence,
        "length": n,
        "temperature": thermo.temperature,
        "gas_constant": thermo.gas_constant,
    }
    return {
        "partition_function": Z,
        "free_energy": free_energy,
        "posterior": posterior,
        "meta": meta,
    }


def mea(seq: str, posterior: List[List[float]], gamma: float = 1.0) -> Dict[str, object]:
    """Return MEA structure (dot-bracket) using posterior probabilities."""
    sequence = _normalize_sequence(seq)
    n = len(sequence)
    if len(posterior) != n:
        raise ValueError("Posterior matrix size must match sequence length.")

    unpaired_weights = [max(0.0, 1.0 - sum(row)) for row in posterior]
    dp = [[0.0 for _ in range(n)] for _ in range(n)]
    trace = [[None for _ in range(n)] for _ in range(n)]

    for length in range(1, n + 1):
        for i in range(0, n - length + 1):
            j = i + length - 1
            best = dp[i][j - 1] + unpaired_weights[j] if j > 0 else unpaired_weights[j]
            trace_choice = ("unpair", j)
            for k in range(i, j):
                candidate = dp[i][k] + dp[k + 1][j]
                if candidate > best:
                    best = candidate
                    trace_choice = ("split", k)
            if _can_pair(sequence[i], sequence[j]):
                pair_score = (dp[i + 1][j - 1] + gamma * posterior[i][j]) if i + 1 <= j - 1 else gamma * posterior[i][j]
                if pair_score > best:
                    best = pair_score
                    trace_choice = ("pair", None)
            dp[i][j] = best
            trace[i][j] = trace_choice

    structure = ["." for _ in range(n)]

    def traceback(i: int, j: int) -> None:
        if i >= j or trace[i][j] is None:
            return
        action, value = trace[i][j]
        if action == "unpair":
            traceback(i, j - 1)
        elif action == "split":
            k = value
            traceback(i, k)
            traceback(k + 1, j)
        elif action == "pair":
            structure[i] = "("
            structure[j] = ")"
            traceback(i + 1, j - 1)

    traceback(0, n - 1)
    dot_bracket = "".join(structure)
    return {
        "sequence": sequence,
        "structure": dot_bracket,
        "score": dp[0][n - 1],
        "gamma": gamma,
    }
