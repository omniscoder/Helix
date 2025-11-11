"""Expectation-maximization motif discovery."""
from __future__ import annotations

import math
from typing import List, Sequence

from .core import BASES, EPSILON, MotifResult, compute_log_prob, initial_pwm, normalize_sequence


def discover_em(
    sequences: Sequence[str],
    width: int,
    iterations: int = 50,
    pseudocount: float = 0.1,
) -> MotifResult:
    clean_sequences = [normalize_sequence(seq) for seq in sequences if len(normalize_sequence(seq)) >= width]
    if not clean_sequences:
        raise ValueError("No sequences long enough for the specified width.")

    pwm = initial_pwm(width)
    log_likelihood = 0.0

    for _ in range(iterations):
        responsibilities: List[List[float]] = []
        log_likelihood = 0.0
        for seq in clean_sequences:
            positions = len(seq) - width + 1
            scores = [compute_log_prob(seq[pos : pos + width], pwm) for pos in range(positions)]
            max_log = max(scores)
            exp_scores = [math.exp(score - max_log) for score in scores]
            total = sum(exp_scores) + EPSILON
            responsibilities.append([score / total for score in exp_scores])
            log_likelihood += math.log(total) + max_log

        new_pwm = [{base: pseudocount for base in BASES} for _ in range(width)]
        total_resp = [pseudocount * len(BASES) for _ in range(width)]
        for seq, resp in zip(clean_sequences, responsibilities):
            for pos, weight in enumerate(resp):
                kmer = seq[pos : pos + width]
                for offset, base in enumerate(kmer):
                    new_pwm[offset][base] += weight
                    total_resp[offset] += weight

        for idx in range(width):
            for base in BASES:
                new_pwm[idx][base] /= total_resp[idx]
        pwm = new_pwm

    return MotifResult(pwm=pwm, log_likelihood=log_likelihood, solver="em")
