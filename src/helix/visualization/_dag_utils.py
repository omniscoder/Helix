"""Internal helpers shared across edit DAG visualizations."""
from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import networkx as nx


def format_stage_label(stage: Optional[str]) -> str:
    """Return a human-friendly stage/phase label."""
    if not stage:
        return "Unknown"
    return stage.replace("_", " ").title()


def stage_palette(stages: Sequence[str], cmap_name: str) -> Dict[str, Tuple[float, float, float, float]]:
    """Assign deterministic colors to stage names."""
    unique = []
    for stage in stages:
        if stage not in unique:
            unique.append(stage)
    if not unique:
        return {}
    cmap = plt.get_cmap(cmap_name)
    if len(unique) == 1:
        return {unique[0]: cmap(0.1)}
    palette = {}
    for idx, stage in enumerate(unique):
        palette[stage] = cmap(idx / max(1, len(unique) - 1))
    return palette


def compute_layout(graph: nx.DiGraph, layout: str, seed: int = 0) -> Dict[str, Tuple[float, float]]:
    """Return node positions for the requested layout."""
    layout = layout.lower()
    if layout == "kamada-kawai":
        return nx.kamada_kawai_layout(graph)
    if layout == "spectral":
        return nx.spectral_layout(graph)
    if layout == "circular":
        return nx.circular_layout(graph)
    if layout == "shell":
        return nx.shell_layout(graph)
    if layout == "planar" and nx.check_planarity(graph)[0]:
        return nx.planar_layout(graph)
    if layout == "random":
        return nx.random_layout(graph, seed=seed)
    return nx.spring_layout(graph, seed=seed)
