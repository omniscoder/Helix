"""CRISPR guide scoring + off-target enumeration."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable, List, Mapping

from .. import bioinformatics
from ..string import edit as string_edit
from .pam import PAM, get_pam, match_pam, reverse_complement_pattern


def _prepare_pam(pam: Dict[str, str] | PAM | str) -> Dict[str, str]:
    if isinstance(pam, str):
        return get_pam(pam)
    if isinstance(pam, PAM):
        return pam.as_dict()
    return pam


def _normalize_guide_sequence(genome: str, guide: Mapping[str, object]) -> str:
    seq = str(guide.get("sequence") or "").strip().upper()
    if seq:
        return seq
    start = int(guide["start"])
    end = int(guide["end"])
    strand = guide.get("strand", "+")
    window = genome[start:end]
    if strand == "+":
        return window
    return bioinformatics.reverse_complement(window)


def _mismatch_profile(guide_seq: str, target_seq: str) -> List[Dict[str, object]]:
    profile: List[Dict[str, object]] = []
    for idx, (g_base, t_base) in enumerate(zip(guide_seq, target_seq)):
        if g_base != t_base:
            profile.append({"position": idx, "guide": g_base, "target": t_base})
    return profile


def enumerate_off_targets(
    genome_seq: str,
    guide: Mapping[str, object],
    pam: Dict[str, str] | PAM | str,
    *,
    max_mm: int = 3,
    max_gap: int = 0,
) -> List[Dict[str, object]]:
    """Enumerate off-target candidates using Myers search (<=64 nt)."""

    if max_mm < 0:
        raise ValueError("max_mm must be >= 0")
    if max_gap != 0:
        raise ValueError("max_gap currently must be 0 (gapless search).")

    genome = bioinformatics.normalize_sequence(genome_seq)
    if not genome:
        return []
    pam_def = _prepare_pam(pam)
    pattern = pam_def.get("pattern")
    if not pattern:
        raise ValueError("PAM definitions require a 'pattern' field.")
    pam_len = len(pattern)

    guide_seq = _normalize_guide_sequence(genome, guide).upper()
    guide_len = len(guide_seq)
    if guide_len == 0:
        return []
    if guide_len > 64:
        raise ValueError("Off-target search currently supports guide lengths <= 64.")

    hits: List[Dict[str, object]] = []
    genome_rc = bioinformatics.reverse_complement(genome)
    rc_pattern = reverse_complement_pattern(pattern)
    rc_pam_def = {"pattern": rc_pattern}

    def _record_hit(start: int, strand: str) -> None:
        end = start + guide_len
        if start < 0 or end > len(genome):
            return
        target_slice = genome[start:end]
        aligned_target = target_slice if strand == "+" else bioinformatics.reverse_complement(target_slice)
        mismatch_profile = _mismatch_profile(guide_seq, aligned_target)
        if len(mismatch_profile) > max_mm:
            return
        if strand == "+":
            pam_pos = end
            pam_ok = pam_pos + pam_len <= len(genome) and match_pam(genome, pam_def, pam_pos)
        else:
            pam_pos = start - pam_len
            pam_ok = pam_pos >= 0 and match_pam(genome, rc_pam_def, pam_pos)
        hits.append(
            {
                "guide_id": guide.get("id"),
                "strand": strand,
                "start": start,
                "end": end,
                "pam_ok": pam_ok,
                "mismatches": mismatch_profile,
                "distance": len(mismatch_profile),
            }
        )

    plus_hits = string_edit.myers_search(guide_seq, genome, max_mm)
    for hit in plus_hits:
        _record_hit(hit["start"], "+")

    minus_hits = string_edit.myers_search(guide_seq, genome_rc, max_mm)
    genome_len = len(genome)
    for hit in minus_hits:
        start_rc = hit["start"]
        if start_rc < 0:
            continue
        start = genome_len - (start_rc + guide_len)
        _record_hit(start, "-")

    hits.sort(key=lambda h: (h["guide_id"], h["start"], h["strand"]))
    return hits


def score_on_target(guide: Mapping[str, object], params: Mapping[str, float]) -> float:
    """Simple GC-content based score (placeholder)."""

    gc = float(guide.get("gc_content", 0.0))
    target = float(params.get("gc_optimum", 0.5))
    weight = float(params.get("gc_weight", 1.0))
    score = max(0.0, 1.0 - abs(gc - target) * weight)
    return round(score, 6)


def score_off_targets(hits: Iterable[Dict[str, object]], weights: Mapping[str, object]) -> List[Dict[str, object]]:
    """Apply multiplicative penalties based on mismatch identity/position + PAM status."""

    mismatch_penalties: Mapping[str, float] = weights.get("mismatch_penalties", {})  # type: ignore[arg-type]
    position_penalties: List[float] = list(weights.get("position_penalties", []))  # type: ignore[arg-type]
    pam_penalty = float(weights.get("pam_penalty", 0.75))
    floor = float(weights.get("floor", 0.0))

    scored: List[Dict[str, object]] = []
    for hit in hits:
        score = 1.0
        mismatches: List[Mapping[str, object]] = hit.get("mismatches", [])
        for mismatch in mismatches:
            ref = str(mismatch.get("target", "N"))
            alt = str(mismatch.get("guide", "N"))
            key = f"{ref}>{alt}"
            score *= float(mismatch_penalties.get(key, weights.get("default_mismatch", 0.8)))
            pos = int(mismatch.get("position", 0))
            if position_penalties:
                idx = min(max(pos, 0), len(position_penalties) - 1)
                score *= float(position_penalties[idx])
        if not hit.get("pam_ok", True):
            score *= pam_penalty
        hit_copy = dict(hit)
        hit_copy["score"] = round(max(score, floor), 6)
        scored.append(hit_copy)
    return scored


DEFAULT_WEIGHTS = {
    "on_target": {"gc_optimum": 0.5, "gc_weight": 1.0},
    "off_target": {
        "default_mismatch": 0.8,
        "mismatch_penalties": {
            "A>G": 0.85,
            "T>C": 0.85,
            "C>T": 0.9,
            "G>A": 0.9,
        },
        "pam_penalty": 0.5,
        "floor": 0.0,
    },
}


def load_weights(path: str | Path | None) -> Dict[str, object]:
    if path is None:
        return json.loads(json.dumps(DEFAULT_WEIGHTS))
    weight_path = Path(path)
    data = json.loads(weight_path.read_text(encoding="utf-8"))
    merged = json.loads(json.dumps(DEFAULT_WEIGHTS))
    for section, values in data.items():
        if isinstance(values, dict):
            merged.setdefault(section, {}).update(values)  # type: ignore[call-arg]
        else:
            merged[section] = values
    return merged
