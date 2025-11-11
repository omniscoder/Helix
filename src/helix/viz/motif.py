"""Motif visualization helpers."""
from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import matplotlib.pyplot as plt

from ._utils import VizSpec, apply_rc, finalize

BASES = ["A", "C", "G", "T"]
COLORS = {
    "A": "#1b9e77",
    "C": "#d95f02",
    "G": "#7570b3",
    "T": "#e7298a",
}


def _normalize_column(column: Dict[str, float]) -> Dict[str, float]:
    values = [max(0.0, float(column.get(base, 0.0))) for base in BASES]
    total = sum(values)
    if total <= 0:
        return {base: 1.0 / len(BASES) for base in BASES}
    return {base: values[idx] / total for idx, base in enumerate(BASES)}


def _column_info(column: Dict[str, float]) -> tuple[float, Dict[str, float]]:
    entropy = 0.0
    for prob in column.values():
        if prob > 0:
            entropy -= prob * math.log2(prob)
    info = max(0.0, 2.0 - entropy)
    heights = {base: column[base] * info for base in BASES}
    return info, heights


def plot_motif_logo(
    *,
    pwm: Sequence[Dict[str, float]],
    title: str = "Motif logo",
    save: Optional[str] = None,
    save_viz_spec: Optional[str] = None,
):
    """Render a stacked information logo for PWM columns."""
    apply_rc()
    width = len(pwm)
    fig_width = max(4.0, 0.5 * width + 1.0)
    fig, ax = plt.subplots(figsize=(fig_width, 3.5))

    infos: List[float] = []
    for idx, raw_column in enumerate(pwm):
        column = _normalize_column(raw_column)
        info, heights = _column_info(column)
        infos.append(info)
        y_pos = 0.0
        for base, height in sorted(heights.items(), key=lambda item: item[1]):
            if height <= 0:
                continue
            ax.bar(
                idx + 0.5,
                height,
                bottom=y_pos,
                width=0.8,
                color=COLORS.get(base, "#666666"),
                edgecolor="white",
                linewidth=0.4,
            )
            ax.text(
                idx + 0.5,
                y_pos + height / 2,
                base,
                ha="center",
                va="center",
                color="white",
                fontsize=max(8, 10 + height * 2),
                fontweight="bold",
            )
            y_pos += height

    ax.set_xlim(0, width)
    ax.set_ylim(0, 2.0)
    ax.set_xticks([i + 0.5 for i in range(width)])
    ax.set_xticklabels(range(1, width + 1))
    ax.set_ylabel("Information (bits)")
    ax.set_xlabel("Position")
    ax.set_title(title)

    spec = VizSpec(
        kind="motif_logo",
        meta={"columns": width},
        primitives={
            "mean_information": float(sum(infos) / width) if width else 0.0,
            "max_information": float(max(infos) if infos else 0.0),
        },
    )
    return finalize(fig, spec, save=save, save_viz_spec=save_viz_spec)


def plot_pwm(pwm: List[Dict[str, float]], output: Path, title: str = "Motif logo") -> None:
    """Backward-compatible wrapper used by CLI --plot flag."""
    plot_motif_logo(pwm=pwm, title=title, save=str(output))
