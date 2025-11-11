"""RNA visualization helpers (dot-plots, entropy)."""
from __future__ import annotations

from pathlib import Path
from typing import List


def plot_dotplot(posterior: List[List[float]], output: Path, title: str = "RNA dot-plot") -> None:
    import matplotlib.pyplot as plt  # type: ignore
    import numpy as np

    matrix = np.array(posterior)
    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(matrix, cmap="viridis", origin="lower", vmin=0, vmax=1)
    fig.colorbar(im, ax=ax, label="P(i,j)")
    ax.set_xlabel("j")
    ax.set_ylabel("i")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


def plot_arc(dotbracket: str, output: Path, title: str = "RNA arc diagram") -> None:
    import matplotlib.pyplot as plt  # type: ignore
    import numpy as np

    stack = []
    pairs = []
    for idx, char in enumerate(dotbracket):
        if char == "(":
            stack.append(idx)
        elif char == ")" and stack:
            start = stack.pop()
            pairs.append((start, idx))

    n = len(dotbracket)
    fig, ax = plt.subplots(figsize=(6, 2))
    ax.axhline(0, color="black", linewidth=0.5)
    for i, j in pairs:
        xs = np.linspace(i, j, 50)
        ys = 0.1 * np.sin(np.linspace(0, np.pi, 50))
        ax.plot(xs, ys, color="tab:blue")
    ax.set_xlim(0, max(n - 1, 1))
    ax.set_ylim(-0.02, 0.35)
    ax.set_title(title)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


def plot_entropy(entropy: List[float], output: Path, title: str = "RNA entropy") -> None:
    import matplotlib.pyplot as plt  # type: ignore

    fig, ax = plt.subplots(figsize=(6, 2))
    ax.plot(range(len(entropy)), entropy, color="tab:orange")
    ax.set_xlabel("Position")
    ax.set_ylabel("Entropy (nats)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)
