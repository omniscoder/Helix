"""2.5D workflow geometry builders for CRISPR/prime outcomes (v2)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Sequence, Tuple

import math
import numpy as np

from helix import bioinformatics


@dataclass(frozen=True)
class WorkflowLanes:
    """
    Semantic y-positions (world units) for outcome lanes.

    These are *world* coords; the viewer is responsible for mapping them to pixels
    via an orthographic projection. lane_height is a base "world" height which
    is typically overridden to a constant pixel height in the shader.
    """

    y_cut: float = 0.0
    y_del: float = -1.0
    y_ins: float = -2.0
    y_sub: float = -3.0
    y_no_cut: float = 1.0
    y_intended: float = 2.0
    lane_height: float = 0.8


# ---------------------------------------------------------------------------
# Probability → weight: stable, log-scaled mapping in [0, 1]
# ---------------------------------------------------------------------------

def _prob_to_weight(prob: float, floor: float = 1e-5) -> float:
    """
    Map probability in (0, 1] to [0, 1], compressing dynamic range via log10.

    Invariants:
      - prob <= 0        → weight ~ 0
      - prob == floor    → weight = 0
      - prob == 1        → weight = 1
    """
    p = max(float(prob), floor)
    span = math.log10(1.0 / floor)  # = -log10(floor) > 0
    w = (math.log10(p) - math.log10(floor)) / span
    return max(0.0, min(1.0, w))


# ---------------------------------------------------------------------------
# Outcome classification: "visual kinds" that drive lanes & colors
# ---------------------------------------------------------------------------

def _visual_kind(raw_kind: str, outcome: Dict[str, Any], payload: Dict[str, Any]) -> str:
    """
    Collapse raw diff.kind + optional metadata into a small visual vocabulary.

    Returns one of:
      - "deletion"
      - "insertion"
      - "substitution"
      - "no_cut"
      - "intended"
      - "other"
    """
    token = (raw_kind or "").lower()

    # Try to detect "intended" outcome via label / tags / payload hints
    label = (outcome.get("label") or "").lower()
    tags: Sequence[str] = outcome.get("tags") or ()
    tags_lower = {t.lower() for t in tags}

    intended_label = (payload.get("intended_label") or "").lower()
    intended_labels = {s.lower() for s in (payload.get("intended_labels") or [])}

    is_intended = (
        "intended" in label
        or "hdr" in label
        or "scarless" in label
        or "intended" in tags_lower
        or "hdr" in tags_lower
        or label == intended_label
        or (label and label in intended_labels)
    )
    if is_intended:
        return "intended"

    tok = token
    if "del" in tok:
        return "deletion"
    if "ins" in tok:
        return "insertion"
    if "sub" in tok or "snv" in tok or "mismatch" in tok:
        return "substitution"
    if "no_cut" in tok or "nocut" in tok or "none" in tok:
        return "no_cut"

    return "other"


def _kind_to_lane(kind: str, lanes: WorkflowLanes) -> float:
    """Map a *visual* kind string to a y-lane."""
    token = (kind or "no_cut").lower()
    if token == "deletion":
        return lanes.y_del
    if token == "insertion":
        return lanes.y_ins
    if token == "substitution":
        return lanes.y_sub
    if token == "intended":
        return lanes.y_intended
    if token == "no_cut":
        return lanes.y_no_cut
    # Default semantic lane: treat unknown as "substitution-like"
    return lanes.y_sub


def _kind_code(kind: str) -> float:
    """
    Map visual kind into a small float code for the fragment shader.

    These codes must remain stable or your colors will flip.
    """
    token = (kind or "other").lower()
    if token == "deletion":
        return 0.0
    if token == "insertion":
        return 1.0
    if token == "substitution":
        return 2.0
    if token == "no_cut":
        return 3.0
    if token == "intended":
        return 4.0
    return 5.0  # "other"


# ---------------------------------------------------------------------------
# Heat lane: per-base outcome mass
# ---------------------------------------------------------------------------

def _heat_from_outcomes(
    seq_len: int,
    guide: Dict[str, Any],
    outcomes: Iterable[Dict[str, Any]],
) -> np.ndarray:
    """
    Aggregate outcome probability mass across bases.

    heat[i] = sum(probabilities for outcomes whose [start, end) window covers base i),
    normalized so max(heat) == 1 or all zeros if no signal.

    This is intentionally oblivious to kind (del/ins/sub) – it’s "how much editing
    pressure is focused here?".
    """
    heat = np.zeros((seq_len,), dtype="f4")
    for outcome in outcomes:
        prob = float(outcome.get("probability", 0.0))
        if prob <= 0.0:
            continue
        diff = outcome.get("diff") or {}
        kind = (diff.get("kind") or "no_cut").lower()
        if "no_cut" in kind:
            continue
        start = int(diff.get("start", guide.get("start", 0)))
        end = int(diff.get("end", guide.get("end", start + 1)))
        start = max(0, min(seq_len - 1, start))
        end = max(start + 1, min(seq_len, end))
        heat[start:end] += prob
    max_val = float(heat.max())
    if max_val > 0.0:
        heat /= max_val
    return heat


# ---------------------------------------------------------------------------
# Bars (instanced quads) for the top-K outcome mass
# ---------------------------------------------------------------------------

def _bars_from_outcomes(
    seq_len: int,
    guide: Dict[str, Any],
    outcomes: Sequence[Dict[str, Any]],
    payload: Dict[str, Any],
    *,
    lanes: WorkflowLanes,
    keep_mass: float = 0.995,
    min_prob: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """
    Compress outcome distribution into a manageable number of instanced quads.

    Each bar encodes a contiguous [start, end) interval at the lane for its kind,
    plus a log-scaled probability weight in [0, 1].

    Returns:
      inst_events   : (K, 7)  -> [xc, half_width, y_center, lane_height, weight, kind_code, z]
      inst_labels   : (K, 3)  -> [x, y, z] label anchors
      captured_mass : sum(probabilities for rendered bars)
      total_mass    : sum(probabilities for all outcomes)
    """
    # Sort from most likely to least; this makes "keep_mass" do something meaningful
    ordered = sorted(outcomes, key=lambda oc: float(oc.get("probability", 0.0)), reverse=True)

    total_mass = 0.0
    for oc in ordered:
        total_mass += max(0.0, float(oc.get("probability", 0.0)))

    captured = 0.0
    bars: List[List[float]] = []
    label_anchors: List[List[float]] = []

    for outcome in ordered:
        prob = float(outcome.get("probability", 0.0))
        if prob < min_prob:
            continue

        diff = outcome.get("diff") or {}
        raw_kind = str(diff.get("kind", "no_cut"))
        visual_kind = _visual_kind(raw_kind, outcome, payload)

        if visual_kind == "no_cut":
            start = int(guide.get("start", 0))
            end = max(start + 1, int(guide.get("end", start + 1)))
        else:
            start = int(diff.get("start", guide.get("start", 0)))
            end = int(diff.get("end", start + 1))

        start = max(0, min(seq_len - 1, start))
        end = max(start + 1, min(seq_len, end))

        xc = 0.5 * (start + end)
        hw = 0.5 * max(1e-3, end - start)  # half-width in "bases" as world units
        lane_y = _kind_to_lane(visual_kind, lanes)
        weight = _prob_to_weight(prob)
        kind_code = _kind_code(visual_kind)

        # z is a small positive lift that scales with probability for 2.5D ordering
        z = 0.05 + 0.35 * weight

        bars.append([xc, hw, lane_y, lanes.lane_height, weight, kind_code, z])
        label_anchors.append([xc, lane_y + lanes.lane_height * 0.6, z])

        captured += prob
        # Stop once we've accounted for "keep_mass" of total probability.
        # Note: if total_mass < keep_mass due to filtering, we still cut when captured >= total_mass.
        if total_mass > 0.0 and captured >= min(keep_mass * total_mass, total_mass):
            break

    inst = np.array(bars, dtype="f4") if bars else np.zeros((0, 7), dtype="f4")
    anchors = np.array(label_anchors, dtype="f4") if label_anchors else np.zeros((0, 3), dtype="f4")
    return inst, anchors, captured, total_mass


# ---------------------------------------------------------------------------
# Public builders
# ---------------------------------------------------------------------------

def build_workflow2p5d_crispr(
    sequence: str,
    guide: Dict[str, Any],
    payload: Dict[str, Any],
    *,
    lanes: WorkflowLanes | None = None,
    keep_mass: float = 0.995,
    min_prob: float = 1e-6,
) -> Dict[str, np.ndarray | float]:
    """
    Build 2.5D workflow buffers for a CRISPR simulation payload.

    Input:
      sequence : genomic window (string, FASTA headers already stripped)
      guide    : dict with at least "start", "end" (0-based indices into sequence)
      payload  : simulation payload with "outcomes": [...], each outcome having:
                 - "probability": float
                 - "diff": {"kind": ..., "start": int, "end": int, ...}
                 - optional "label", "tags", etc. for intended/HDR classification

    Output keys:
      "rail_xy"        : (2, 3)      -> main rail polyline endpoints
      "tick_xy"        : (T, 3)      -> tick segments every 10 bases (can be empty)
      "heat_values"    : (L,)        -> per-base scalar heat in [0, 1]
      "inst_events"    : (K, 7)      -> bars [xc, half_width, y, lane_h, weight, kind_code, z]
      "inst_labels"    : (K, 3)      -> label anchors [x, y, z]
      "x_extent"       : float       -> domain extent in x (bases)
      "y_extent"       : float       -> total |y| extent for viewport/frustum
      "captured_mass"  : float       -> sum of probabilities in inst_events
      "total_mass"     : float       -> sum of probabilities of all outcomes
    """
    normalized = bioinformatics.normalize_sequence(sequence)
    seq_len = len(normalized)
    if seq_len == 0:
        raise ValueError("Sequence must not be empty.")
    lanes = lanes or WorkflowLanes()

    outcomes = payload.get("outcomes") or []
    heat = _heat_from_outcomes(seq_len, guide, outcomes)
    guide_start = int(guide.get("start", 0) or 0)
    guide_end = int(guide.get("end", guide_start + 1) or guide_start + 1)
    strand = str(guide.get("strand", "+") or "+")
    if strand == "+":
        cut_idx = guide_end - 3
    else:
        cut_idx = guide_start + 3
    cut_idx = max(0, min(seq_len - 1, cut_idx))

    inst_events, inst_labels, captured, total = _bars_from_outcomes(
        seq_len,
        guide,
        outcomes,
        payload,
        lanes=lanes,
        keep_mass=keep_mass,
        min_prob=min_prob,
    )

    # Main rail: x in [0, seq_len-1] at y = y_cut
    rail_xy = np.array(
        [
            [0.0, lanes.y_cut, 0.0],
            [float(max(1, seq_len - 1)), lanes.y_cut, 0.0],
        ],
        dtype="f4",
    )

    # Simple ticks at every 10bp; viewer can choose to ignore or fade at zoom levels
    tick_vertices: List[List[float]] = []
    for idx in range(0, seq_len, 10):
        tick_vertices.append([float(idx), lanes.y_cut - 0.3, 0.0])
        tick_vertices.append([float(idx), lanes.y_cut + 0.3, 0.0])
    tick_xy = np.array(tick_vertices, dtype="f4") if tick_vertices else np.zeros((0, 3), dtype="f4")

    return {
        "rail_xy": rail_xy,
        "tick_xy": tick_xy,
        "heat_values": heat.astype("f4"),
        "inst_events": inst_events,
        "inst_labels": inst_labels,
        "x_extent": float(max(1, seq_len - 1)),
        "y_extent": 4.0,
        "captured_mass": float(captured),
        "total_mass": float(total),
        "cut_index": float(cut_idx),
        "heat_y": float(lanes.y_cut - 1.0),
        "heat_height": float(lanes.lane_height * 0.6),
        "slice_width": 0.15,
    }


def build_workflow2p5d_prime(
    sequence: str,
    payload: Dict[str, Any],
    *,
    lanes: WorkflowLanes | None = None,
) -> Dict[str, np.ndarray | float]:
    """
    Placeholder prime workflow builder.

    Currently reuses the CRISPR mapping with a dummy guide that spans the window.
    In the next iteration, we can:
      - put PBS/RTT on dedicated lanes
      - mark intended prime edit vs byproducts
      - split lanes by nick vs repair vs byproduct class
    """
    guide = {"start": 0, "end": len(sequence), "strand": "+"}
    return build_workflow2p5d_crispr(sequence, guide, payload, lanes=lanes)
