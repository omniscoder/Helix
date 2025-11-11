"""Motif visualization helpers."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List

BASES = ["A", "C", "G", "T"]


def plot_pwm(pwm: List[Dict[str, float]], output: Path, title: str = "Motif PWM") -> None:
    import matplotlib.pyplot as plt  # type: ignore
    import numpy as np

    matrix = np.array([[position.get(base, 0.0) for base in BASES] for position in pwm]).T
    fig, ax = plt.subplots(figsize=(len(pwm) * 0.5 + 1, 3))
    im = ax.imshow(matrix, cmap="viridis", aspect="auto", vmin=0, vmax=1)
    ax.set_yticks(range(len(BASES)))
    ax.set_yticklabels(BASES)
    ax.set_xticks(range(len(pwm)))
    ax.set_xlabel("Position")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="Probability")
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)
