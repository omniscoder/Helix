"""Shared helpers for Helix visualizations."""
from __future__ import annotations

import io
import json
from dataclasses import dataclass, asdict
from importlib import metadata
from typing import Any, Dict, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RC = {
    "figure.dpi": 120,
    "savefig.dpi": 120,
    "font.size": 10,
    "axes.grid": True,
    "axes.facecolor": "white",
}


def apply_rc() -> None:
    for key, value in RC.items():
        plt.rcParams[key] = value


try:  # pragma: no cover - importlib metadata path only runs once
    HELIX_VERSION = metadata.version("helix")
except metadata.PackageNotFoundError:  # pragma: no cover - editable installs
    HELIX_VERSION = "dev"


@dataclass(frozen=True)
class VizSpec:
    kind: str
    meta: Dict[str, Any]
    primitives: Dict[str, Any]
    spec_version: str = "1.0"

    def to_json(self) -> str:
        return json.dumps(asdict(self), sort_keys=True, separators=(",", ":"))


def _footer_text(spec: VizSpec) -> str:
    parts = [f"Helix {HELIX_VERSION}", spec.kind]
    if spec.meta:
        meta_bits = []
        for key in sorted(spec.meta):
            value = spec.meta[key]
            if isinstance(value, (int, float, str)):
                meta_bits.append(f"{key}={value}")
        if meta_bits:
            parts.append(", ".join(meta_bits))
    return " • ".join(parts)


def finalize(
    fig: plt.Figure,
    spec: VizSpec,
    save: Optional[str] = None,
    save_viz_spec: Optional[str] = None,
) -> Tuple[plt.Figure, Dict[str, Any]]:
    footer = _footer_text(spec)
    fig.text(
        0.01,
        0.01,
        footer,
        fontsize=8,
        color="#555555",
        ha="left",
        va="bottom",
        alpha=0.85,
    )
    if save:
        fig.savefig(save, bbox_inches="tight", facecolor="white")
    if save_viz_spec:
        with open(save_viz_spec, "w", encoding="utf-8") as handle:
            handle.write(spec.to_json())
    return fig, asdict(spec)


def fig_bytes_sha256(fig: plt.Figure) -> str:
    import hashlib

    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", facecolor="white")
    return hashlib.sha256(buf.getvalue()).hexdigest()
