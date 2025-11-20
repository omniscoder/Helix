"""Shared sequence encoding helpers for Helix engines."""

from __future__ import annotations

import numpy as np

ASCII_TO_BASE = np.full(256, 4, dtype=np.uint8)
ASCII_TO_BASE[ord("A")] = 0
ASCII_TO_BASE[ord("C")] = 1
ASCII_TO_BASE[ord("G")] = 2
ASCII_TO_BASE[ord("T")] = 3


def encode_sequence_to_uint8(sequence: str) -> np.ndarray:
    """Return a uint8 array with the Helix DNA encoding (A=0, C=1, G=2, T=3, other=4)."""

    if not sequence:
        return np.zeros((0,), dtype=np.uint8)
    upper = sequence.upper()
    raw = np.frombuffer(upper.encode("ascii"), dtype=np.uint8)
    return ASCII_TO_BASE[raw]


__all__ = ["encode_sequence_to_uint8", "ASCII_TO_BASE"]
