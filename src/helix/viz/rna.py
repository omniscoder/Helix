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
