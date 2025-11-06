"""Quick GC skew visualization using Helix's bioinformatics helpers.

Run with:
    python examples/gc_skew_demo.py
"""
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np

import bioinformatics


def load_sequence(source: Path | None = None) -> str:
    """Return an uppercase DNA string with whitespace removed."""
    if source is None:
        raw = bioinformatics.seq
    else:
        raw = source.read_text()
    return "".join(raw.upper().split())


def compute_skew_profile(genome: str) -> Iterable[int]:
    """Delegate to `bioinformatics.skew` and drop the numpy array interface."""
    profile: np.ndarray = bioinformatics.skew(genome)
    return profile.tolist()


def main() -> None:
    genome = load_sequence()
    print(f"Loaded sequence length: {len(genome)} bases")

    skew_profile = compute_skew_profile(genome)
    x = range(len(skew_profile))

    plt.plot(x, skew_profile)
    plt.title("GC Skew (Helix sample fragment)")
    plt.xlabel("Nucleotide position")
    plt.ylabel("Cumulative skew")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
