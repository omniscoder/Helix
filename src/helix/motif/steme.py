"""Suffix-tree accelerated EM (STEME-style) simplified."""
from __future__ import annotations

import math
import random
from typing import Sequence

from .core import MotifResult, compute_log_prob, initial_pwm, normalize_sequence, BASES, EPSILON


def discover_steme(
    sequences: Sequence[str],
    width: int,
    iterations: int = 20,
    restarts: int = 5,
) -> MotifResult:
    clean_sequences = [normalize_sequence(seq) for seq in sequences if len(normalize_sequence(seq)) >= width]
    if not clean_sequences:
        raise ValueError("No sequences long enough for the specified width.")

    best_result: MotifResult | None = None
    for _ in range(restarts):
        pwm = initial_pwm(width)
        # random jitter to mimic different suffix-tree seeds
        for column in pwm:
            noise = [random.random() for _ in BASES]
            total = sum(noise)
            for base, value in zip(BASES, noise):
                column[base] = value / total

        log_likelihood = 0.0
        for _ in range(iterations):
            responsibilities = []
            log_likelihood = 0.0
            for seq in clean_sequences:
                positions = len(seq) - width + 1
                scores = [compute_log_prob(seq[pos : pos + width], pwm) for pos in range(positions)]
                max_log = max(scores)
                exp_scores = [math.exp(score - max_log) for score in scores]
                total = sum(exp_scores) + EPSILON
                responsibilities.append([score / total for score in exp_scores])
                log_likelihood += math.log(total) + max_log

            new_pwm = [{base: 0.5 for base in BASES} for _ in range(width)]
            total_resp = [0.5 * len(BASES) for _ in range(width)]
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

        candidate = MotifResult(pwm=pwm, log_likelihood=log_likelihood, solver="steme")
        if best_result is None or candidate.log_likelihood > best_result.log_likelihood:
            best_result = candidate

    return best_result  # type: ignore
