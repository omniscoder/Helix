"""CRISPR track visualization (guide + outcomes)."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt


def _guide_interval(guide: Dict[str, object]) -> Tuple[int, int]:
    start = int(guide.get("start", 0))
    end = int(guide.get("end", start))
    if start > end:
        start, end = end, start
    return start, end


def render_crispr_track(
    payload: Dict[str, object],
    *,
    save: str | None = None,
    show: bool = False,
    extra_meta: Dict[str, object] | None = None,
    save_viz_spec: str | None = None,
) -> Tuple[plt.Figure, Dict[str, object]]:
    site = payload.get("site", {})
    guide = payload.get("guide", {})
    outcomes: List[Dict[str, object]] = payload.get("outcomes", [])
    site_length = int(site.get("length", 0)) or len(site.get("sequence", ""))
    guide_start, guide_end = _guide_interval(guide)

    fig, axes = plt.subplots(2, 1, figsize=(8, 5), gridspec_kw={"height_ratios": [1, 1]})

    # Sequence track
    ax = axes[0]
    ax.set_title("Guide + PAM context")
    ax.set_xlim(0, max(site_length, guide_end + 1))
    ax.set_ylim(0, 1)
    ax.axhline(0.5, color="#9ca3af", linewidth=1)
    ax.add_patch(
        plt.Rectangle((guide_start, 0.25), guide_end - guide_start, 0.5, color="#60a5fa", alpha=0.7)
    )
    pam = guide.get("pam_site")
    if isinstance(pam, dict):
        pam_start = int(pam.get("start", guide_end))
        pam_end = int(pam.get("end", pam_start))
        ax.add_patch(
            plt.Rectangle((pam_start, 0.2), max(1, pam_end - pam_start), 0.6, color="#f97316", alpha=0.5)
        )
    ax.text(
        guide_start,
        0.85,
        f"Guide {guide.get('id', '')} ({guide.get('strand', '+')})",
        fontsize=10,
        color="#1f2937",
    )
    ax.set_yticks([])
    ax.set_xlabel("Position (nt)")

    # Outcome distribution
    ax2 = axes[1]
    ax2.set_title("Cut/repair outcomes")
    labels = [entry["label"] for entry in outcomes]
    probs = [entry.get("probability", 0) for entry in outcomes]
    ax2.bar(labels, probs, color="#34d399")
    ax2.set_ylabel("Probability")
    ax2.set_ylim(0, max(probs + [0.1]) * 1.2)
    for label, prob in zip(labels, probs):
        ax2.text(label, prob + 0.01, f"{prob:.2f}", ha="center", va="bottom")

    fig.tight_layout()

    spec = {
        "kind": "viz_crispr_track",
        "site_length": site_length,
        "guide": {
            "id": guide.get("id"),
            "start": guide_start,
            "end": guide_end,
            "strand": guide.get("strand"),
        },
        "outcomes": [{"label": entry["label"], "probability": entry.get("probability") } for entry in outcomes],
    }
    if extra_meta:
        spec["meta"] = extra_meta

    if save:
        fig.savefig(save, dpi=200)
    if save_viz_spec:
        Path(save_viz_spec).write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8")
    if show:  # pragma: no cover - interactive path
        plt.show()
    plt.close(fig)
    return fig, spec
