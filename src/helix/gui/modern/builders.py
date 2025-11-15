"""Helpers that translate simulation payloads into ModernGL viz specs."""

from __future__ import annotations

import math
from typing import Any, Mapping, Sequence, Tuple, List

import numpy as np

from helix import bioinformatics
from helix.crispr.model import DigitalGenome
from helix.prime.model import PegRNA, PrimeEditor
from helix.prime.simulator import locate_prime_target_site

from .spec import EditEvent, EditEventType, EditVisualizationSpec


def _clamp_index(value: int, length: int) -> int:
    if length <= 0:
        return 0
    return max(0, min(length - 1, value))


def _safe_range(start: int, end: int, length: int) -> tuple[int, int]:
    start = max(0, min(length, start))
    end = max(start, min(length, end))
    return start, end


def _best_outcome(payload: Mapping[str, Any]) -> Mapping[str, Any] | None:
    outcomes = payload.get("outcomes")
    if not isinstance(outcomes, Sequence) or not outcomes:
        return None
    best = max(outcomes, key=lambda entry: entry.get("probability", 0.0))
    if isinstance(best, Mapping):
        return best
    return None


def build_crispr_viz_spec(
    sequence: str,
    guide: Mapping[str, Any],
    sim_payload: Mapping[str, Any],
) -> EditVisualizationSpec | None:
    """Map a crispr.sim payload to an EditVisualizationSpec."""

    normalized = bioinformatics.normalize_sequence(sequence)
    if not normalized:
        return None
    seq_len = len(normalized)

    guide_start = int(guide.get("start", 0) or 0)
    guide_end = int(guide.get("end", guide_start) or guide_start)
    guide_start, guide_end = _safe_range(guide_start, guide_end, seq_len)

    pam_site = guide.get("pam_site") or {}
    strand = guide.get("strand", "+")
    pam_index = pam_site.get("start")
    if pam_index is None:
        pam_index = guide_end if strand == "+" else guide_start
    pam_index = _clamp_index(int(pam_index), seq_len)
    cut_idx = _clamp_index(int(guide_end - 3 if strand == "+" else guide_start + 3), seq_len)

    strand = str(guide.get("strand", "+"))
    if strand == "+":
        cut_idx = max(guide_start, guide_end - 3)
    else:
        cut_idx = min(guide_end, guide_start + 3)
    cut_idx = _clamp_index(int(cut_idx), seq_len)

    guide_center = _clamp_index(guide_start + (guide_end - guide_start) // 2, seq_len)
    best = _best_outcome(sim_payload) or {}

    events = [
        EditEvent(t=0.2, type=EditEventType.RECOGNITION, index=guide_center),
        EditEvent(t=0.35, type=EditEventType.NICK_PRIMARY, index=cut_idx),
        EditEvent(t=0.45, type=EditEventType.CUT, index=cut_idx, metadata={"label": best.get("label")}),
        EditEvent(
            t=0.9,
            type=EditEventType.REPAIR_COMPLETE,
            index=cut_idx,
            metadata={"label": best.get("label"), "probability": best.get("probability")},
        ),
    ]

    guide_label = guide.get("id") or guide.get("name") or "guide"
    edit_type = f"CRISPR – {guide_label}"
    if best.get("label"):
        edit_type += f" ({best['label']})"

    rail_meta = build_rail_band_metadata(
        normalized,
        (guide_start, guide_end),
        pam_index,
        cut_idx,
        sim_payload,
    )

    metadata = {
        "source": "helix.crispr.sim",
        "guide_id": guide_label,
        "draws": sim_payload.get("draws"),
        "best_outcome": best,
        "pam_site": pam_site,
        "rail_bands": rail_meta,
    }

    return EditVisualizationSpec(
        sequence=normalized,
        pam_index=pam_index,
        guide_range=(guide_start, guide_end),
        edit_type=edit_type,
        events=events,
        metadata=metadata,
    )


def _classify_diff(label: str | None, diff: Mapping[str, Any] | None) -> str:
    token = ""
    if diff:
        token = str(diff.get("edit") or "").lower()
    if not token and label:
        token = str(label).lower()
    if not diff or token == "no_cut" or label == "no_cut":
        return "no_cut"
    if token.startswith("del"):
        return "deletion"
    if token.startswith("ins"):
        return "insertion"
    if token.startswith("scarless") or token.startswith("hdr") or "sub" in token:
        return "substitution"
    return token or "other"


def build_rail_band_metadata(
    sequence: str,
    guide_range: tuple[int, int],
    pam_index: int,
    cut_index: int,
    payload: Mapping[str, Any],
) -> dict[str, Any]:
    seq_len = max(1, len(sequence))
    bands: list[dict[str, Any]] = []
    outcomes = payload.get("outcomes") or []
    for outcome in outcomes:
        probability = float(outcome.get("probability") or 0.0)
        diff = outcome.get("diff")
        label = outcome.get("label") or "outcome"
        start = guide_range[0]
        end = guide_range[1]
        if isinstance(diff, Mapping):
            if diff.get("start") is not None:
                start = int(diff["start"])
            if diff.get("end") is not None:
                end = int(diff["end"])
        start, end = _safe_range(start, end, seq_len)
        if start == end:
            end = min(seq_len, start + 1)
        kind = _classify_diff(label, diff if isinstance(diff, Mapping) else None)
        magnitude = abs(end - start)
        bands.append(
            {
                "label": label,
                "kind": kind,
                "start": start,
                "end": end,
                "probability": probability,
                "magnitude": magnitude,
            }
        )
    bands.sort(key=lambda entry: entry.get("probability", 0.0), reverse=True)
    for idx, band in enumerate(bands):
        band["row"] = idx
    return {
        "sequence": sequence,
        "pam_index": pam_index,
        "guide_range": guide_range,
        "cut_index": cut_index,
        "bands": bands,
        "draws": payload.get("draws"),
    }


def build_prime_viz_spec(
    genome: DigitalGenome,
    peg: PegRNA,
    editor: PrimeEditor,
    sim_payload: Mapping[str, Any],
) -> EditVisualizationSpec | None:
    """Map a prime.edit_sim payload to an EditVisualizationSpec."""

    site = locate_prime_target_site(genome, peg)
    if site is None:
        return None
    chrom_seq = bioinformatics.normalize_sequence(genome.sequences.get(site.chrom, ""))
    if not chrom_seq:
        return None
    flank = max(len(peg.pbs or ""), len(peg.rtt or ""), 10)
    window_start = max(0, site.start - flank)
    window_end = min(len(chrom_seq), site.end + flank)
    if window_end <= window_start:
        return None
    window_seq = chrom_seq[window_start:window_end]
    seq_len = len(window_seq)

    guide_start = site.start - window_start
    guide_end = site.end - window_start
    guide_start, guide_end = _safe_range(guide_start, guide_end, seq_len)
    pam_anchor = guide_end if site.strand >= 0 else guide_start
    pam_index = _clamp_index(pam_anchor, seq_len)

    nick_index = (site.start - window_start) + editor.nick_to_edit_offset
    nick_index = _clamp_index(nick_index, seq_len)
    rtt_seq = bioinformatics.normalize_sequence(peg.rtt or "")
    rtt_len = len(rtt_seq)
    repair_index = _clamp_index(nick_index + max(0, rtt_len - 1), seq_len)

    guide_center = _clamp_index(guide_start + (guide_end - guide_start) // 2, seq_len)
    best = _best_outcome(sim_payload) or {}

    events = [
        EditEvent(t=0.2, type=EditEventType.RECOGNITION, index=guide_center),
        EditEvent(t=0.35, type=EditEventType.NICK_PRIMARY, index=nick_index),
        EditEvent(t=0.55, type=EditEventType.PRIME_INIT, index=nick_index),
        EditEvent(t=0.75, type=EditEventType.RT_SYNTHESIS, start=nick_index, length=rtt_len or None),
        EditEvent(t=1.0, type=EditEventType.FLAP_RESOLUTION, index=repair_index),
        EditEvent(
            t=1.3,
            type=EditEventType.REPAIR_COMPLETE,
            index=repair_index,
            metadata={"description": best.get("description"), "stage": best.get("stage")},
        ),
    ]

    edit_label = peg.name or peg.metadata.get("label") if peg.metadata else None
    if not edit_label:
        edit_label = f"peg-{peg.spacer[:5]}…"

    metadata = {
        "source": "helix.prime.sim",
        "peg": peg.name or peg.metadata,
        "editor": editor.name,
        "chrom": site.chrom,
        "site": {"start": site.start, "end": site.end, "strand": site.strand},
        "best_outcome": best,
    }
    prime_meta = _build_prime_scaffold_metadata(
        sequence=window_seq,
        guide_range=(guide_start, guide_end),
        nick_index=nick_index,
        pbs_length=len(peg.pbs or ""),
        rtt_length=len(rtt_seq),
        payload=sim_payload,
    )
    if prime_meta["outcomes"]:
        metadata["prime_scaffold"] = prime_meta
    metadata["rail_bands"] = build_rail_band_metadata(
        window_seq,
        (guide_start, guide_end),
        pam_index,
        nick_index,
        sim_payload,
    )

    return EditVisualizationSpec(
        sequence=window_seq,
        pam_index=pam_index,
        guide_range=(guide_start, guide_end),
        edit_type=f"Prime – {edit_label}",
        events=events,
        metadata=metadata,
    )


def _build_prime_scaffold_metadata(
    sequence: str,
    guide_range: tuple[int, int],
    nick_index: int,
    pbs_length: int,
    rtt_length: int,
    payload: Mapping[str, Any],
) -> dict[str, Any]:
    seq_len = max(1, len(sequence))
    guide_start, guide_end = guide_range
    pbs_start = max(0, guide_start - max(0, pbs_length))
    pbs_end = guide_start
    rtt_start = guide_end
    rtt_end = min(seq_len, guide_end + max(0, rtt_length))
    outcomes_meta: list[dict[str, Any]] = []
    outcomes = payload.get("outcomes") or []
    for outcome in outcomes:
        diff = outcome.get("diff")
        if isinstance(diff, Mapping):
            start = int(diff.get("start", guide_start))
            end = int(diff.get("end", start))
        else:
            start = guide_start
            end = guide_end
        start, end = _safe_range(start, end, seq_len)
        if start == end:
            end = min(seq_len, start + 1)
        probability = float(outcome.get("probability") or 0.0)
        kind = _classify_diff(outcome.get("label"), diff if isinstance(diff, Mapping) else None)
        outcomes_meta.append(
            {
                "label": outcome.get("label"),
                "kind": kind,
                "start": start,
                "end": end,
                "probability": probability,
            }
        )
    outcomes_meta.sort(key=lambda entry: entry.get("probability", 0.0), reverse=True)
    return {
        "sequence": sequence,
        "guide_range": guide_range,
        "nick_index": _clamp_index(nick_index, seq_len),
        "pbs_span": (pbs_start, pbs_end),
        "rtt_span": (rtt_start, rtt_end),
        "outcomes": outcomes_meta,
    }


# ---------------------------------------------------------------------------
# 3D sci-fi builders
# ---------------------------------------------------------------------------

_SCIFI_COLORS = {
    "bg": (0.02, 0.04, 0.08, 1.0),
    "rail": (0.28, 0.88, 0.86, 1.0),
    "guide": (0.98, 0.78, 0.26, 1.0),
    "deletion": (0.98, 0.65, 0.22, 1.0),
    "insertion": (0.92, 0.33, 0.82, 1.0),
    "substitution": (0.26, 0.86, 0.98, 1.0),
    "no_cut": (0.20, 0.30, 0.40, 0.5),
    "other": (0.78, 0.66, 0.94, 1.0),
}

BP_PER_TURN = 10.5
RISE_NM = 0.34
RADIUS_NM = 1.0
WOBBLE_NM = 0.08
THETA_STEP = 2.0 * math.pi / BP_PER_TURN
PROB_FLOOR = 1e-5
KIND_CODES = {
    "deletion": 0.0,
    "insertion": 1.0,
    "substitution": 2.0,
    "no_cut": 3.0,
    "other": 4.0,
}

# How many arcs / ribbons we draw in 3D before it gets noisy
MAX_CRISPR_RIBBONS = 32
MIN_RIBBON_PROB = 1e-3

MAX_ORBITAL_ARCS = 48
MIN_ORBITAL_PROB = 5e-4


def _label_hash01(label: str | None) -> float:
    """Deterministic hash(label) -> [0, 1) for tiny per-outcome jitters."""
    if not label:
        return 0.0
    h = 0
    for ch in str(label):
        h = (h * 131 + ord(ch)) & 0xFFFFFFFF
    return (h % 10_000) / 10_000.0


def build_crispr_rail_3d_spec(
    sequence: str,
    guide: Mapping[str, Any],
    payload: Mapping[str, Any],
) -> EditVisualizationSpec | None:
    normalized = bioinformatics.normalize_sequence(sequence)
    if not normalized:
        return None
    seq_len = len(normalized)
    guide_start = int(guide.get("start", 0) or 0)
    guide_end = int(guide.get("end", guide_start) or guide_start)
    guide_start, guide_end = _safe_range(guide_start, guide_end, seq_len)
    pam_site = guide.get("pam_site") or {}
    strand = guide.get("strand", "+")
    pam_index = pam_site.get("start")
    if pam_index is None:
        pam_index = guide_end if strand == "+" else guide_start
    pam_index = _clamp_index(int(pam_index), seq_len)
    cut_idx = _clamp_index(int(guide_end - 3 if strand == "+" else guide_start + 3), seq_len)
    cut_idx = _clamp_index(int(guide_end - 3 if guide.get("strand", "+") == "+" else guide_start + 3), seq_len)

    strand_a, strand_b, tangents_a, tangents_b = _dna_double_helix_vertices(seq_len)
    rail_vertices = np.vstack([strand_a, strand_b]).astype("f4")
    rail_ids = [0] * len(strand_a) + [1] * len(strand_b)
    norm = [i / max(1, seq_len - 1) for i in range(seq_len)]
    param_coords = norm + norm

    ribbon_curves: List[np.ndarray] = []
    weight_curves: List[np.ndarray] = []
    kind_curves: List[np.ndarray] = []
    repair_offsets: List[int] = [0]

    outcomes = list(payload.get("outcomes") or [])
    outcomes.sort(key=lambda o: float(o.get("probability") or 0.0), reverse=True)
    outcomes = outcomes[:MAX_CRISPR_RIBBONS]

    for outcome in outcomes:
        prob = max(0.0, float(outcome.get("probability") or 0.0))
        if prob < MIN_RIBBON_PROB:
            continue

        diff = outcome.get("diff") if isinstance(outcome.get("diff"), Mapping) else None
        label = outcome.get("label")
        start_idx, end_idx = _diff_span(diff, guide_start, guide_end, seq_len)
        if start_idx == end_idx:
            end_idx = min(seq_len - 1, start_idx + 1)

        strand_start = strand_a if (start_idx % 2 == 0) else strand_b
        strand_end = strand_b if (end_idx % 2 == 0) else strand_a
        tan_start = tangents_a if strand_start is strand_a else tangents_b
        tan_end = tangents_b if strand_end is strand_b else tangents_a

        kind = _classify_diff(label, diff)
        kind_code = KIND_CODES.get(kind, KIND_CODES["other"])

        jitter = _label_hash01(label)
        if kind == "deletion":
            pan_bias = (0.0 + 0.2 * (jitter - 0.5), 0.0, 1.0)
        elif kind == "insertion":
            pan_bias = (0.2, 0.8 + 0.3 * (jitter - 0.5), 0.4)
        elif kind == "substitution":
            pan_bias = (0.8 + 0.3 * (jitter - 0.5), 0.2, 0.6)
        elif kind == "no_cut":
            pan_bias = (0.0, 0.0, -1.0)
        else:
            pan_bias = (0.3, 0.3, 0.8)

        curve, weights, kinds = _build_ribbon_curve(
            strand_start,
            strand_end,
            tan_start,
            tan_end,
            start_idx,
            end_idx,
            prob,
            kind_code,
            pan_bias=pan_bias,
        )
        if not len(curve):
            continue
        ribbon_curves.append(curve)
        weight_curves.append(weights)
        kind_curves.append(kinds)
        repair_offsets.append(repair_offsets[-1] + len(curve))

    arc_vertices = np.vstack(ribbon_curves).astype("f4") if ribbon_curves else np.zeros((0, 3), dtype="f4")
    arc_weights = np.concatenate(weight_curves).astype("f4") if weight_curves else np.zeros((0,), dtype="f4")
    arc_kinds = np.concatenate(kind_curves).astype("f4") if kind_curves else np.zeros((0,), dtype="f4")

    geometry = {
        "rail_vertices": rail_vertices,
        "rail_ids": rail_ids,
        "param_coords": param_coords,
        "arc_vertices": arc_vertices,
        "arc_offsets": np.array(repair_offsets, dtype="i4"),
        "arc_weights": arc_weights,
        "arc_kinds": arc_kinds,
    }
    z_max = strand_a[-1][2] if strand_a else 4.0
    metadata = {
        "source": "helix.crispr.sim",
        "layout": "rail_3d",
        "custom_geometry": geometry,
        "camera": {
            "pos": (3.5, -5.0, z_max * 0.5),
            "target": (0.0, 0.0, z_max * 0.5),
            "up": (0.0, 0.0, 1.0),
            "fov_deg": 45.0,
        },
        "colors": _SCIFI_COLORS,
        "rail_bands": build_rail_band_metadata(
            normalized,
            (guide_start, guide_end),
            pam_index,
            cut_idx,
            payload,
        ),
    }
    events = [
        EditEvent(t=0.15, type=EditEventType.RECOGNITION, index=guide_start),
        EditEvent(t=0.5, type=EditEventType.CUT, index=cut_idx),
        EditEvent(t=0.9, type=EditEventType.REPAIR_COMPLETE, index=cut_idx),
    ]
    return EditVisualizationSpec(
        sequence=normalized,
        pam_index=pam_index,
        guide_range=(guide_start, guide_end),
        edit_type=f"CRISPR – {guide.get('id') or guide.get('name') or 'guide'} (rail_3d)",
        events=events,
        metadata=metadata,
    )


def build_crispr_orbitals_3d_spec(
    sequence: str,
    guide: Mapping[str, Any],
    payload: Mapping[str, Any],
) -> EditVisualizationSpec | None:
    normalized = bioinformatics.normalize_sequence(sequence)
    if not normalized:
        return None
    seq_len = len(normalized)
    guide_start = int(guide.get("start", 0) or 0)
    guide_end = int(guide.get("end", guide_start) or guide_start)
    guide_start, guide_end = _safe_range(guide_start, guide_end, seq_len)
    pam_site = guide.get("pam_site") or {}
    strand = guide.get("strand", "+")
    pam_index = pam_site.get("start")
    if pam_index is None:
        pam_index = guide_end if strand == "+" else guide_start
    pam_index = _clamp_index(int(pam_index), seq_len)
    cut_idx = _clamp_index(int(guide_end - 3 if strand == "+" else guide_start + 3), seq_len)

    strand_a, strand_b, _, _ = _dna_double_helix_vertices(seq_len, radius=1.05, pitch_per_base=0.2)
    rail_vertices = np.vstack([strand_a, strand_b]).astype("f4")
    rail_ids = [0] * len(strand_a) + [1] * len(strand_b)
    norm = [i / max(1, seq_len - 1) for i in range(seq_len)]
    param_coords = norm + norm

    ring_curves: List[np.ndarray] = []
    weight_curves: List[np.ndarray] = []
    kind_curves: List[np.ndarray] = []
    offsets: List[int] = [0]
    angle_cursor = 0.0

    outcomes = list(payload.get("outcomes") or [])
    outcomes.sort(key=lambda o: float(o.get("probability") or 0.0), reverse=True)
    outcomes = outcomes[:MAX_ORBITAL_ARCS]

    for outcome in outcomes:
        prob = max(0.0, float(outcome.get("probability") or 0.0))
        if prob < MIN_ORBITAL_PROB:
            continue

        diff = outcome.get("diff") if isinstance(outcome.get("diff"), Mapping) else None
        kind = _classify_diff(outcome.get("label"), diff)
        kind_code = KIND_CODES.get(kind, KIND_CODES["other"])
        start_idx, end_idx = _diff_span(diff, guide_start, guide_end, seq_len)
        ins_len = len(str(diff.get("ins") or "")) if diff else 0
        delta_bp = abs((end_idx - start_idx) - ins_len)

        radius = 0.6 + 0.35 * math.log10(delta_bp + 1.5)
        arc_span = max(0.25, 2.0 * math.pi * prob)
        steps = max(32, int(220 * prob))
        weight = _prob_weight(prob)

        base_z_scale = 0.22
        if kind == "insertion":
            z_scale = base_z_scale * 1.6
        elif kind == "deletion":
            z_scale = base_z_scale * 0.8
        elif kind == "substitution":
            z_scale = base_z_scale * 1.2
        else:
            z_scale = base_z_scale

        jitter = _label_hash01(outcome.get("label"))
        local_angle_offset = (jitter - 0.5) * 0.8

        pts = []
        for i in range(steps):
            t = angle_cursor + local_angle_offset + (arc_span * i / max(1, steps - 1))
            x = radius * math.cos(t)
            y = radius * math.sin(t)
            z = z_scale * math.sin(3.0 * t + (0.3 if kind == "insertion" else 0.0))
            pts.append([x, y, z])
        curve = np.array(pts, dtype="f4")
        if not len(curve):
            continue
        ring_curves.append(curve)
        weight_curves.append(np.full((steps,), weight, dtype="f4"))
        kind_curves.append(np.full((steps,), kind_code, dtype="f4"))
        offsets.append(offsets[-1] + len(curve))
        angle_cursor += arc_span + 0.4

    arc_vertices = np.vstack(ring_curves).astype("f4") if ring_curves else np.zeros((0, 3), dtype="f4")
    arc_weights = np.concatenate(weight_curves).astype("f4") if weight_curves else np.zeros((0,), dtype="f4")
    arc_kinds = np.concatenate(kind_curves).astype("f4") if kind_curves else np.zeros((0,), dtype="f4")

    geometry = {
        "rail_vertices": rail_vertices,
        "rail_ids": rail_ids,
        "param_coords": param_coords,
        "arc_vertices": arc_vertices,
        "arc_offsets": np.array(offsets, dtype="i4"),
        "arc_weights": arc_weights,
        "arc_kinds": arc_kinds,
    }
    metadata = {
        "source": "helix.crispr.sim",
        "layout": "orbitals_3d",
        "custom_geometry": geometry,
        "camera": {
            "pos": (0.0, -3.8, 1.6),
            "target": (0.0, 0.0, 0.2),
            "up": (0.0, 0.0, 1.0),
            "fov_deg": 45.0,
        },
        "colors": _SCIFI_COLORS,
        "rail_bands": build_rail_band_metadata(
            normalized,
            (guide_start, guide_end),
            pam_index,
            cut_idx,
            payload,
        ),
    }
    events = [
        EditEvent(t=0.2, type=EditEventType.RECOGNITION, index=guide_start),
        EditEvent(t=0.7, type=EditEventType.REPAIR_COMPLETE, index=guide_end),
    ]
    return EditVisualizationSpec(
        sequence=normalized,
        pam_index=pam_index,
        guide_range=(guide_start, guide_end),
        edit_type=f"CRISPR – {guide.get('id') or guide.get('name') or 'guide'} (orbitals)",
        events=events,
        metadata=metadata,
    )


def build_prime_scaffold_3d_spec(
    genome: DigitalGenome,
    peg: PegRNA,
    editor: PrimeEditor,
    payload: Mapping[str, Any],
) -> EditVisualizationSpec | None:
    site = locate_prime_target_site(genome, peg)
    if site is None:
        return None
    chrom_seq = bioinformatics.normalize_sequence(genome.sequences.get(site.chrom, ""))
    if not chrom_seq:
        return None
    flank = max(len(peg.pbs or ""), len(peg.rtt or ""), 12)
    window_start = max(0, site.start - flank)
    window_end = min(len(chrom_seq), site.end + flank)
    if window_end <= window_start:
        return None
    window_seq = chrom_seq[window_start:window_end]
    seq_len = len(window_seq)

    guide_start = site.start - window_start
    guide_end = site.end - window_start
    guide_start, guide_end = _safe_range(guide_start, guide_end, seq_len)
    pam_index = guide_end
    nick_index = _clamp_index((site.start - window_start) + editor.nick_to_edit_offset, seq_len)

    rail_top, rail_bottom, tangents_top, tangents_bottom = _parallel_rails(seq_len)
    rail_vertices = np.vstack([rail_top, rail_bottom]).astype("f4")
    rail_ids = [0] * len(rail_top) + [1] * len(rail_bottom)
    norm = [i / max(1, seq_len - 1) for i in range(seq_len)]
    param_coords = norm + norm

    ribbon_curves: List[np.ndarray] = []
    weight_curves: List[np.ndarray] = []
    kind_curves: List[np.ndarray] = []
    repair_offsets: List[int] = [0]
    outcomes = payload.get("outcomes") or []
    for outcome in outcomes[:20]:
        prob = max(0.0, float(outcome.get("probability") or 0.0))
        diff = outcome.get("diff") if isinstance(outcome.get("diff"), Mapping) else None
        label = outcome.get("label")
        start_idx, end_idx = _diff_span(diff, guide_start, guide_end, seq_len)
        kind = _classify_diff(outcome.get("label"), diff)
        kind_code = KIND_CODES.get(kind, KIND_CODES["other"])
        curve, weights, kinds = _build_ribbon_curve(
            rail_top,
            rail_bottom,
            tangents_top,
            tangents_bottom,
            start_idx,
            end_idx,
            prob,
            kind_code,
            pan_bias=(0.0, 0.5, 1.0),
        )
        ribbon_curves.append(curve)
        weight_curves.append(weights)
        kind_curves.append(kinds)
        repair_offsets.append(repair_offsets[-1] + len(curve))

    arc_vertices = np.vstack(ribbon_curves).astype("f4") if ribbon_curves else np.zeros((0, 3), dtype="f4")
    arc_weights = np.concatenate(weight_curves).astype("f4") if weight_curves else np.zeros((0,), dtype="f4")
    arc_kinds = np.concatenate(kind_curves).astype("f4") if kind_curves else np.zeros((0,), dtype="f4")

    geometry = {
        "rail_vertices": rail_vertices,
        "rail_ids": rail_ids,
        "param_coords": param_coords,
        "arc_vertices": arc_vertices,
        "arc_offsets": np.array(repair_offsets, dtype="i4"),
        "arc_weights": arc_weights,
        "arc_kinds": arc_kinds,
    }
    metadata = {
        "source": "helix.prime.sim",
        "layout": "prime_scaffold_3d",
        "custom_geometry": geometry,
        "camera": {
            "pos": (0.0, -6.0, 2.5),
            "target": (rail_top[-1][0] * 0.5, 0.0, 0.0),
            "up": (0.0, 0.0, 1.0),
            "fov_deg": 48.0,
        },
        "colors": _SCIFI_COLORS,
        "rail_bands": build_rail_band_metadata(
            window_seq,
            (guide_start, guide_end),
            pam_index,
            nick_index,
            payload,
        ),
    }
    events = [
        EditEvent(t=0.2, type=EditEventType.RECOGNITION, index=guide_start),
        EditEvent(t=0.55, type=EditEventType.PRIME_INIT, index=nick_index),
        EditEvent(t=1.1, type=EditEventType.REPAIR_COMPLETE, index=guide_end),
    ]
    return EditVisualizationSpec(
        sequence=window_seq,
        pam_index=pam_index,
        guide_range=(guide_start, guide_end),
        edit_type=f"Prime – {peg.name or 'pegRNA'} (scaffold_3d)",
        events=events,
        metadata=metadata,
    )


# Helper utilities ----------------------------------------------------------

def _diff_span(diff: Mapping[str, Any] | None, guide_start: int, guide_end: int, seq_len: int) -> tuple[int, int]:
    if not diff:
        return guide_start, guide_end
    start = int(diff.get("start", guide_start))
    end = int(diff.get("end", guide_end))
    return _safe_range(start, end, seq_len)


def _dna_double_helix_vertices(
    n_bp: int,
    radius: float = RADIUS_NM,
    bases_per_turn: float = BP_PER_TURN,
    pitch_per_base: float = RISE_NM,
    wobble_amp: float = WOBBLE_NM,
) -> tuple[
    list[tuple[float, float, float]],
    list[tuple[float, float, float]],
    list[tuple[float, float, float]],
    list[tuple[float, float, float]],
]:
    strand_a: list[tuple[float, float, float]] = []
    strand_b: list[tuple[float, float, float]] = []
    step = 2.0 * math.pi / bases_per_turn
    count = max(1, n_bp)
    for idx in range(count):
        theta = idx * step + _groove_phase(idx)
        wobble = wobble_amp * math.sin(3.0 * theta)
        x = (radius + wobble) * math.cos(theta)
        y = (radius + wobble) * math.sin(theta)
        z = idx * pitch_per_base
        strand_a.append((x, y, z))
        theta_b = theta + math.pi
        wobble_b = wobble_amp * math.sin(3.0 * theta_b)
        x_b = (radius + wobble_b) * math.cos(theta_b)
        y_b = (radius + wobble_b) * math.sin(theta_b)
        strand_b.append((x_b, y_b, z))
    tangents_a = _compute_tangents(strand_a)
    tangents_b = _compute_tangents(strand_b)
    return strand_a, strand_b, tangents_a, tangents_b


def _parallel_rails(
    length: int,
    spacing: float = 1.4,
    step: float = 0.2,
) -> tuple[
    list[tuple[float, float, float]],
    list[tuple[float, float, float]],
    list[tuple[float, float, float]],
    list[tuple[float, float, float]],
]:
    top: list[tuple[float, float, float]] = []
    bottom: list[tuple[float, float, float]] = []
    count = max(1, length)
    for idx in range(count):
        x = idx * step
        top.append((x, spacing, 0.4))
        bottom.append((x, -spacing, -0.5))
    tangents_top = _compute_tangents(top)
    tangents_bottom = _compute_tangents(bottom)
    return top, bottom, tangents_top, tangents_bottom


def _pos_on_strand(verts: Sequence[tuple[float, float, float]], idx: int) -> tuple[float, float, float]:
    if not verts:
        return (0.0, 0.0, 0.0)
    clamped = max(0, min(idx, len(verts) - 1))
    return verts[clamped]


def _bezier3(p0, p1, p2, p3, steps: int) -> list[tuple[float, float, float]]:
    path: list[tuple[float, float, float]] = []
    for i in range(max(2, steps)):
        t = i / max(1, steps - 1)
        u = 1.0 - t
        b0 = u * u * u
        b1 = 3 * u * u * t
        b2 = 3 * u * t * t
        b3 = t * t * t
        x = b0 * p0[0] + b1 * p1[0] + b2 * p2[0] + b3 * p3[0]
        y = b0 * p0[1] + b1 * p1[1] + b2 * p2[1] + b3 * p3[1]
        z = b0 * p0[2] + b1 * p1[2] + b2 * p2[2] + b3 * p3[2]
        path.append((x, y, z))
    return path


def _build_ribbon_curve(
    strand_start: Sequence[tuple[float, float, float]],
    strand_end: Sequence[tuple[float, float, float]],
    tangents_start: Sequence[tuple[float, float, float]],
    tangents_end: Sequence[tuple[float, float, float]],
    start_idx: int,
    end_idx: int,
    probability: float,
    kind_code: float,
    pan_bias: tuple[float, float, float] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not strand_start or not strand_end:
        return np.zeros((0, 3), dtype="f4"), np.zeros((0,), dtype="f4"), np.zeros((0,), dtype="f4")
    start_idx = max(0, min(start_idx, len(strand_start) - 1))
    end_idx = max(0, min(end_idx, len(strand_end) - 1))
    start_point = np.array(strand_start[start_idx], dtype="f4")
    end_point = np.array(strand_end[end_idx], dtype="f4")
    start_tan = np.array(tangents_start[start_idx], dtype="f4") if tangents_start else np.array([0.0, 0.0, 1.0], dtype="f4")
    end_tan = np.array(tangents_end[end_idx], dtype="f4") if tangents_end else np.array([0.0, 0.0, 1.0], dtype="f4")
    weight = _prob_weight(probability)
    span_nm = abs(end_idx - start_idx) * RISE_NM
    alpha = max(0.25, span_nm * (0.25 + 0.35 * weight))
    beta = 0.5 + 3.5 * weight
    bias = pan_bias or (0.0, 0.0, 1.0)
    normal_start = _normal_from_tangent(start_tan, bias)
    normal_end = _normal_from_tangent(end_tan, bias)
    control1 = start_point + _vec_norm(start_tan) * alpha + normal_start * beta
    control2 = end_point - _vec_norm(end_tan) * alpha + normal_end * beta
    steps = int(min(128, max(32, 28 + int(0.5 * abs(end_idx - start_idx)) + int(48 * weight))))
    curve = _bezier3(start_point, control1, control2, end_point, steps)
    weights = np.full((steps,), weight, dtype="f4")
    kinds = np.full((steps,), kind_code, dtype="f4")
    return curve, weights, kinds


def _groove_phase(idx: int) -> float:
    phase = (idx % BP_PER_TURN) / BP_PER_TURN
    stretch = 1.1 if phase < 0.6 else 0.9
    return (phase * 2.0 * math.pi * stretch) - (phase * 2.0 * math.pi)


def _compute_tangents(points: Sequence[tuple[float, float, float]]) -> list[tuple[float, float, float]]:
    tangents: list[tuple[float, float, float]] = []
    if not points:
        return tangents
    for idx in range(len(points)):
        prev_pt = points[idx - 1] if idx > 0 else points[idx]
        next_pt = points[idx + 1] if idx + 1 < len(points) else points[idx]
        tangent = _vec_sub(next_pt, prev_pt)
        tangents.append(_vec_norm(tangent))
    return tangents


def _prob_weight(prob: float) -> float:
    p = max(PROB_FLOOR, min(1.0, float(prob)))
    span = math.log10(1.0 / PROB_FLOOR)
    value = (math.log10(p) - math.log10(PROB_FLOOR)) / span
    return max(0.0, min(1.0, value))


def _normal_from_tangent(tangent: tuple[float, float, float], bias: tuple[float, float, float]) -> tuple[float, float, float]:
    t = _vec_norm(tangent)
    up = _vec_norm(bias)
    cross = _vec_cross(t, up)
    if _vec_length(cross) < 1e-4:
        cross = _vec_cross(t, (1.0, 0.0, 0.0))
    if _vec_length(cross) < 1e-4:
        cross = (0.0, 1.0, 0.0)
    return _vec_norm(cross)


def _vec_sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _vec_norm(v: tuple[float, float, float] | np.ndarray) -> np.ndarray:
    arr = np.asarray(v, dtype="f4")
    length = np.linalg.norm(arr)
    if length < 1e-9:
        return np.zeros(3, dtype="f4")
    return arr / length


def _vec_length(v: tuple[float, float, float] | np.ndarray) -> float:
    arr = np.asarray(v, dtype="f4")
    return float(np.linalg.norm(arr))


def _vec_cross(a: tuple[float, float, float] | np.ndarray, b: tuple[float, float, float] | np.ndarray) -> np.ndarray:
    return np.cross(np.asarray(a, dtype="f4"), np.asarray(b, dtype="f4")).astype("f4")
