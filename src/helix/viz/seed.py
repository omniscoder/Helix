"""Seed density visualization helpers."""
from __future__ import annotations

from pathlib import Path
from typing import List, Tuple


def plot_density(seeds: List[Tuple[int, str, int]], seq_length: int, output: Path, title: str = "Seed density") -> None:
    import matplotlib.pyplot as plt  # type: ignore

    if not seeds:
        fig, ax = plt.subplots(figsize=(8, 2))
        ax.text(0.5, 0.5, "No seeds", ha="center", va="center")
    else:
        positions = [pos for pos, _, _ in seeds]
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.hist(positions, bins=50, range=(0, seq_length))
        ax.set_xlabel("Position")
        ax.set_ylabel("Count")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)
