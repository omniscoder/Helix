"""RNA visualization helpers (dot-plots, arcs, entropy)."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np

from ._utils import VizSpec, apply_rc, finalize


def plot_rna_dotplot(
    *,
    posterior: List[List[float]],
    vmin: float = 0.0,
    vmax: float = 1.0,
    save: Optional[str] = None,
    save_viz_spec: Optional[str] = None,
    extra_meta: Optional[Dict[str, Any]] = None,
):
    apply_rc()
    size = len(posterior)
    if any(len(row) != size for row in posterior):
        raise AssertionError("posterior must be square")
    matrix = np.array(posterior, dtype=float)
    assert matrix.ndim == 2 and matrix.shape[0] == matrix.shape[1], "posterior must be square"
    n = matrix.shape[0]
    mask = np.tri(n, n, k=0, dtype=bool)
    masked = np.ma.array(matrix, mask=mask)
    fig, ax = plt.subplots(figsize=(4.5, 4.0))
    im = ax.imshow(masked, origin="lower", interpolation="nearest", vmin=vmin, vmax=vmax)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(-0.5, n - 0.5)
    ax.set_xlabel("j (base index)")
    ax.set_ylabel("i (base index)")
    ax.set_title("RNA pairing posterior (upper triangle)")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("P(i,j)")

    upper = matrix[np.triu_indices(n, k=1)]
    q = np.quantile(upper, [0.0, 0.25, 0.5, 0.75, 0.95, 1.0]) if upper.size else np.zeros(6)

    meta = {"n": int(n), "vmin": float(vmin), "vmax": float(vmax)}
    if extra_meta:
        meta.update(extra_meta)
    spec = VizSpec(
        kind="rna_dotplot",
        meta=meta,
        primitives={
            "nonzero_pairs": int((upper > 0).sum()),
            "quantiles": [float(x) for x in q.tolist()],
        },
    )
    return finalize(fig, spec, save=save, save_viz_spec=save_viz_spec)


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
