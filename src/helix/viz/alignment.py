"""Alignment ribbon visualization."""
from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from ._utils import VizSpec, apply_rc, finalize

_CIGAR_RE = re.compile(r"(\d+)([MID=X])")


def _parse_cigar(cigar: str) -> List[Tuple[str, int]]:
    if not cigar:
        return []
    ops = []
    for count, op in _CIGAR_RE.findall(cigar):
        ops.append((op, int(count)))
    return ops


def _path_from_cigar(ref_start: int, read_start: int, cigar: str) -> Tuple[np.ndarray, np.ndarray, List[Dict[str, int]], List[Dict[str, int]]]:
    ref_points: List[int] = [ref_start]
    read_points: List[int] = [read_start]
    insertions: List[Dict[str, int]] = []
    deletions: List[Dict[str, int]] = []

    ref_pos = ref_start
    read_pos = read_start
    for op, count in _parse_cigar(cigar):
        if count <= 0:
            continue
        if op in ("M", "=", "X"):
            for _ in range(count):
                ref_pos += 1
                read_pos += 1
                ref_points.append(ref_pos)
                read_points.append(read_pos)
        elif op == "I":
            start_read = read_pos
            for _ in range(count):
                read_pos += 1
                ref_points.append(ref_pos)
                read_points.append(read_pos)
            insertions.append({"ref": ref_pos, "read_start": start_read, "read_end": read_pos})
        elif op == "D":
            start_ref = ref_pos
            for _ in range(count):
                ref_pos += 1
                ref_points.append(ref_pos)
                read_points.append(read_pos)
            deletions.append({"read": read_pos, "ref_start": start_ref, "ref_end": ref_pos})
    return np.array(ref_points, dtype=float), np.array(read_points, dtype=float), insertions, deletions


def plot_alignment_ribbon(
    *,
    ref_length: int,
    qry_length: int,
    alignment: Dict[str, Any],
    metadata: Optional[Dict[str, Any]] = None,
    title: str | None = None,
    save: Optional[str] = None,
    save_viz_spec: Optional[str] = None,
    extra_meta: Optional[Dict[str, Any]] = None,
):
    """
    Plot a local alignment ribbon using banded Smith-Waterman output.

    Parameters
    ----------
    ref_length : int
        Length of the reference sequence.
    qry_length : int
        Length of the query/read sequence.
    alignment : dict
        Alignment dictionary with ref_start/ref_end/read_start/read_end/score/matches/cigar.
    """

    apply_rc()
    ref_arr, read_arr, insertions, deletions = _path_from_cigar(
        int(alignment.get("ref_start", 0)),
        int(alignment.get("read_start", 0)),
        alignment.get("cigar", ""),
    )
    cigar_ops = _parse_cigar(alignment.get("cigar", ""))
    fig, ax = plt.subplots(figsize=(5.0, 4.0))
    if ref_arr.size and read_arr.size:
        ax.fill_between(ref_arr, read_arr - 0.4, read_arr + 0.4, color="#c6dbef", alpha=0.5, linewidth=0)
        ax.plot(ref_arr, read_arr, color="#084594", linewidth=2.0)
    for ins in insertions:
        ax.vlines(ins["ref"], ins["read_start"], ins["read_end"], color="#cb181d", linewidth=1.3, alpha=0.9)
    for dele in deletions:
        ax.hlines(dele["read"], dele["ref_start"], dele["ref_end"], color="#31a354", linewidth=1.3, alpha=0.9)

    score = alignment.get("score")
    ax.set_xlim(0, max(1, ref_length))
    ax.set_ylim(0, max(1, qry_length))
    ax.set_xlabel("Reference (bp)")
    ax.set_ylabel("Query (bp)")
    ax.set_title(title or f"Alignment ribbon (score={score})")

    meta = {
        "ref_length": int(ref_length),
        "qry_length": int(qry_length),
        "score": float(score) if score is not None else None,
    }
    if metadata:
        meta["metadata"] = metadata

    if extra_meta:
        meta.update(extra_meta)
    spec = VizSpec(
        kind="alignment_ribbon",
        meta=meta,
        primitives={
            "path_points": int(ref_arr.size),
            "insertions": len(insertions),
            "deletions": len(deletions),
            "matches": int(alignment.get("matches", 0)),
            "segments": len(cigar_ops),
        },
    )
    return finalize(fig, spec, save=save, save_viz_spec=save_viz_spec)
