"""Distance heatmap visualization for sketch outputs."""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np

from ._utils import VizSpec, apply_rc, finalize


def plot_distance_heatmap(
    *,
    matrix: Sequence[Sequence[float]],
    labels: Sequence[str],
    method: str = "minhash",
    save: Optional[str] = None,
    save_viz_spec: Optional[str] = None,
    extra_meta: Optional[Dict[str, Any]] = None,
):
    """Render a distance matrix heatmap with a colorbar."""
    apply_rc()
    mat = np.array(matrix, dtype=float)
    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        raise AssertionError("matrix must be square")
    if len(labels) != mat.shape[0]:
        raise AssertionError("labels length must match matrix size")
    n = mat.shape[0]

    fig_size = max(4.0, 0.7 * n)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    im = ax.imshow(mat, origin="lower", cmap="magma")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels)
    ax.set_title(f"{method.upper()} distance heatmap")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("distance")

    meta = {"n": n, "method": method}
    if extra_meta:
        meta.update(extra_meta)
    spec = VizSpec(
        kind="distance_heatmap",
        meta=meta,
        primitives={
            "min": float(np.min(mat)),
            "max": float(np.max(mat)),
            "mean": float(np.mean(mat)),
        },
    )
    return finalize(fig, spec, save=save, save_viz_spec=save_viz_spec)
