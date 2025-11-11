"""Toy thermodynamic parameters for RNA folding."""
from __future__ import annotations

from typing import Callable, Dict, Tuple

PAIRS: set[Tuple[str, str]] = {
    ("A", "U"),
    ("U", "A"),
    ("G", "C"),
    ("C", "G"),
    ("G", "U"),
    ("U", "G"),
}

STACK: Dict[Tuple[str, str, str, str], float] = {
    ("A", "U", "A", "U"): -1.0,
    ("U", "A", "U", "A"): -1.0,
    ("G", "C", "G", "C"): -2.0,
    ("C", "G", "C", "G"): -2.0,
    ("G", "C", "A", "U"): -1.5,
    ("U", "A", "C", "G"): -1.5,
    ("G", "U", "A", "U"): -0.5,
    ("U", "A", "U", "G"): -0.5,
    ("G", "C", "U", "G"): -1.2,
    ("G", "U", "C", "G"): -0.8,
    ("C", "G", "A", "U"): -1.3,
    ("U", "A", "G", "C"): -1.3,
}


def _hairpin_penalty(loop_len: int) -> float:
    return 1.0 + 0.2 * max(0, loop_len - 3)


def _bulge_penalty(size: int) -> float:
    if size <= 0:
        return 0.0
    return 2.0 + 0.4 * (size - 1)


def _internal_penalty(left: int, right: int) -> float:
    return 2.0 + 0.2 * (left + right)


def _stack_energy(a: str, b: str, c: str, d: str) -> float:
    return STACK.get((a, b, c, d), -0.6)


DEFAULTS = {
    "hairpin_min": 3,
    "hairpin_penalty": _hairpin_penalty,
    "bulge_penalty": _bulge_penalty,
    "internal_penalty": _internal_penalty,
    "stack_energy": _stack_energy,
    "beta": 1.0,
}
