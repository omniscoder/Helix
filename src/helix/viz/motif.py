"""Motif visualization helpers."""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

import matplotlib.pyplot as plt

from ._utils import VizSpec, apply_rc, finalize

DEFAULT_ALPHABET = ["A", "C", "G", "T"]
COLORS = {
    "A": "#1b9e77",
    "C": "#d95f02",
    "G": "#7570b3",
    "T": "#e7298a",
}

PWMColumn = Union[Dict[str, float], Sequence[float]]


def _normalize_pwm(pwm: Sequence[PWMColumn], alphabet: Sequence[str]) -> List[Dict[str, float]]:
    normalized: List[Dict[str, float]] = []
    for column in pwm:
        if isinstance(column, dict):
            values = {symbol: float(column.get(symbol, 0.0)) for symbol in alphabet}
        else:
            values_list = list(column)
            if len(values_list) != len(alphabet):
                raise ValueError("PWM column length must match alphabet length.")
            values = {alphabet[idx]: float(val) for idx, val in enumerate(values_list)}
        total = sum(max(0.0, v) for v in values.values())
        if total <= 0:
            normalized.append({symbol: 1.0 / len(alphabet) for symbol in alphabet})
        else:
            normalized.append({symbol: max(0.0, values[symbol]) / total for symbol in alphabet})
    return normalized


def _background_map(background: Optional[Sequence[float]], alphabet: Sequence[str]) -> Dict[str, float]:
    if background is None:
        prob = 1.0 / len(alphabet)
        return {symbol: prob for symbol in alphabet}
    values = list(background)
    if len(values) != len(alphabet):
        raise ValueError("Background length must match alphabet.")
    total = sum(max(0.0, v) for v in values)
    if total <= 0:
        prob = 1.0 / len(alphabet)
        return {symbol: prob for symbol in alphabet}
    return {alphabet[idx]: max(0.0, val) / total for idx, val in enumerate(values)}


def _column_info(column: Dict[str, float], background: Dict[str, float], alphabet: Sequence[str]) -> tuple[float, Dict[str, float]]:
    info = 0.0
    heights: Dict[str, float] = {}
    for symbol in alphabet:
        p = column.get(symbol, 0.0)
        q = background.get(symbol, 1.0 / len(alphabet))
        if p > 0 and q > 0:
            info += p * math.log2(p / q)
        heights[symbol] = p  # heights scaled later by info
    info = max(0.0, info)
    heights = {symbol: heights[symbol] * info for symbol in alphabet}
    return info, heights


def plot_motif_logo(
    *,
    pwm: Sequence[PWMColumn],
    title: str = "Motif logo",
    alphabet: Optional[Sequence[str]] = None,
    background: Optional[Sequence[float]] = None,
    save: Optional[str] = None,
    save_viz_spec: Optional[str] = None,
    extra_meta: Optional[Dict[str, Any]] = None,
):
    """Render a stacked information logo for PWM columns."""
    apply_rc()
    symbols = list(alphabet or DEFAULT_ALPHABET)
    normalized_pwm = _normalize_pwm(pwm, symbols)
    bg_map = _background_map(background, symbols)

    width = len(normalized_pwm)
    fig_width = max(4.0, 0.5 * width + 1.0)
    fig, ax = plt.subplots(figsize=(fig_width, 3.5))

    infos: List[float] = []
    for idx, column in enumerate(normalized_pwm):
        info, heights = _column_info(column, bg_map, symbols)
        infos.append(info)
        y_pos = 0.0
        for symbol, height in sorted(heights.items(), key=lambda item: item[1]):
            if height <= 0:
                continue
            ax.bar(
                idx + 0.5,
                height,
                bottom=y_pos,
                width=0.8,
                color=COLORS.get(symbol, "#666666"),
                edgecolor="white",
                linewidth=0.4,
            )
            ax.text(
                idx + 0.5,
                y_pos + height / 2,
                symbol,
                ha="center",
                va="center",
                color="white",
                fontsize=max(8, 10 + height * 2),
                fontweight="bold",
            )
            y_pos += height

    ax.set_xlim(0, width)
    ax.set_ylim(0, 2.5)
    ax.set_xticks([i + 0.5 for i in range(width)])
    ax.set_xticklabels(range(1, width + 1))
    ax.set_ylabel("Information (bits)")
    ax.set_xlabel("Position")
    ax.set_title(title)

    meta = {"columns": width, "alphabet": "".join(symbols)}
    if extra_meta:
        meta.update(extra_meta)
    spec = VizSpec(
        kind="motif_logo",
        meta=meta,
        primitives={
            "mean_information": float(sum(infos) / width) if width else 0.0,
            "max_information": float(max(infos) if infos else 0.0),
        },
    )
    return finalize(fig, spec, save=save, save_viz_spec=save_viz_spec)


def plot_pwm(pwm: List[Dict[str, float]], output: Path, title: str = "Motif logo") -> None:
    """Backward-compatible wrapper used by CLI --plot flag."""
    plot_motif_logo(pwm=pwm, title=title, save=str(output))
