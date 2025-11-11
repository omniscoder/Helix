"""Seed chaining visualization."""
from __future__ import annotations

from typing import Dict, List, Optional

import matplotlib.pyplot as plt

from ._utils import VizSpec, apply_rc, finalize


def plot_seed_chain(
    *,
    ref_length: int,
    qry_length: int,
    chains: List[List[Dict[str, int]]],
    save: Optional[str] = None,
    save_viz_spec: Optional[str] = None,
):
    apply_rc()
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    line_count = 0
    palette = plt.rcParams["axes.prop_cycle"].by_key().get("color", ["#1f77b4", "#ff7f0e", "#2ca02c"])

    for ci, chain in enumerate(chains):
        color = palette[ci % len(palette)]
        for anchor in chain:
            rs = int(anchor["ref_start"])
            re = int(anchor.get("ref_end", rs + int(anchor.get("len", 0))))
            qs = int(anchor["qry_start"])
            qe = int(anchor.get("qry_end", qs + int(anchor.get("len", 0))))
            ax.plot([rs, re], [qs, qe], linewidth=1.5, color=color, alpha=0.9)
            line_count += 1

    ax.set_xlim(0, ref_length)
    ax.set_ylim(0, qry_length)
    ax.set_xlabel("Reference (bp)")
    ax.set_ylabel("Query (bp)")
    ax.set_title("Seed chaining (ref vs query)")

    spec = VizSpec(
        kind="seed_chain",
        meta={"ref_length": int(ref_length), "qry_length": int(qry_length), "chains": len(chains)},
        primitives={"line_segments": int(line_count)},
    )
    return finalize(fig, spec, save=save, save_viz_spec=save_viz_spec)
