"""Minimizer density visualization."""
from __future__ import annotations

from typing import Any, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from ._utils import VizSpec, apply_rc, finalize


def _normalize_minimizer_positions(minimizers: Sequence[Any]) -> List[int]:
    """Return a permissive list of start positions extracted from mixed payloads."""
    positions: List[int] = []
    for entry in minimizers:
        candidate: Any | None = None
        if isinstance(entry, int):
            candidate = entry
        elif isinstance(entry, (tuple, list)) and entry:
            candidate = entry[0]
        elif isinstance(entry, dict):
            for key in ("pos", "position", "ref_start"):
                if key in entry:
                    candidate = entry[key]
                    break
        else:
            candidate = entry
        if candidate is None:
            continue
        try:
            positions.append(int(candidate))
        except (TypeError, ValueError):
            continue
    return positions


def _clamp_bin_count(bin_count: int, sequence_length: int) -> int:
    length = max(1, int(sequence_length))
    bins = max(1, int(bin_count))
    return min(bins, length)


def _compute_density(positions: np.ndarray, length: int, bin_count: int) -> Tuple[np.ndarray, np.ndarray]:
    if length <= 0 or bin_count <= 0:
        return np.array([]), np.array([])
    bins = np.linspace(0, length, bin_count + 1, dtype=float)
    hist, edges = np.histogram(positions, bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, hist.astype(int)


def plot_minimizer_density(
    *,
    sequence_length: int,
    minimizers: Sequence[Any],
    bin_count: int = 200,
    save: Optional[str] = None,
    save_viz_spec: Optional[str] = None,
):
    apply_rc()
    normalized = np.array(_normalize_minimizer_positions(minimizers), dtype=int)
    valid_positions = normalized[(normalized >= 0) & (normalized < max(1, sequence_length))]
    bins = _clamp_bin_count(bin_count, sequence_length)
    centers, density = _compute_density(valid_positions, sequence_length, bins)

    fig, ax = plt.subplots(figsize=(8, 2.8))
    ax.plot(centers, density, linewidth=1.5)
    ax.set_xlim(0, sequence_length)
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("Minimizers per bin")
    ax.set_title("Minimizer density")

    spec = VizSpec(
        kind="minimizer_density",
        meta={"sequence_length": int(sequence_length), "bin_count": int(bins)},
        primitives={
            "points": int(centers.size),
            "total_minimizers": int(valid_positions.size),
            "density_sum": int(density.sum()),
            "density_max": int(density.max() if density.size else 0),
            "x_min": 0,
            "x_max": int(sequence_length),
        },
    )
    return finalize(fig, spec, save=save, save_viz_spec=save_viz_spec)
