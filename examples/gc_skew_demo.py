"""Quick GC skew visualization using Helix's bioinformatics helpers."""
from __future__ import annotations

import argparse
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot GC skew for a DNA fragment.")
    parser.add_argument(
        "--input",
        type=Path,
        help="Optional path to a text/FASTA file; defaults to the built-in sample.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    genome = load_sequence(args.input)
    print(f"Loaded sequence length: {len(genome)} bases")

    skew_profile = compute_skew_profile(genome)
    x = range(len(skew_profile))

    plt.plot(x, skew_profile)
    plt.title("GC Skew")
    plt.xlabel("Nucleotide position")
    plt.ylabel("Cumulative skew")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
