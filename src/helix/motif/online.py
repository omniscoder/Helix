"""Online EM (EXTREME-style) motif updates."""
from __future__ import annotations

from typing import Sequence

from .core import MotifResult, initial_pwm, normalize_sequence, BASES


def discover_online(
    sequences: Sequence[str],
    width: int,
    learning_rate: float = 0.3,
    passes: int = 3,
) -> MotifResult:
    clean_sequences = [normalize_sequence(seq) for seq in sequences if len(normalize_sequence(seq)) >= width]
    if not clean_sequences:
        raise ValueError("No sequences long enough for the specified width.")

    pwm = initial_pwm(width)
    for _ in range(passes):
        for seq in clean_sequences:
            positions = len(seq) - width + 1
            if positions <= 0:
                continue
            for pos in range(positions):
                kmer = seq[pos : pos + width]
                for idx, base in enumerate(kmer):
                    for b in BASES:
                        target = 1.0 if b == base else 0.0
                        pwm[idx][b] = (1 - learning_rate) * pwm[idx][b] + learning_rate * target
                for idx in range(width):
                    total = sum(pwm[idx][b] for b in BASES) or 1.0
                    for b in BASES:
                        pwm[idx][b] /= total

    return MotifResult(pwm=pwm, log_likelihood=0.0, solver="online")
