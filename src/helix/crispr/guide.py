"""Guide discovery and annotation (placeholder scaffolding)."""
from __future__ import annotations

from __future__ import annotations

from typing import Dict, List, Tuple

from .. import bioinformatics
from .pam import PAM, get_pam, match_pam, reverse_complement_pattern

_STRAND_SET = {"+", "-", "both"}


def _normalize_window(window: Tuple[int, int] | None, seq_len: int) -> Tuple[int, int]:
    if window is None:
        return 0, seq_len
    if len(window) != 2:
        raise ValueError("window must be a tuple of (start, end).")
    start, end = map(int, window)
    if start < 0 or end < 0:
        raise ValueError("window boundaries must be >= 0.")
    if start >= end:
        raise ValueError("window start must be < end.")
    return max(0, start), min(seq_len, end)


def _gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for base in seq if base in {"G", "C"})
    return gc / len(seq)


def _prepare_pam(pam: Dict[str, str] | PAM | str) -> Dict[str, str]:
    if isinstance(pam, str):
        return get_pam(pam)
    if isinstance(pam, PAM):
        return pam.as_dict()
    return pam


def _record_guide(
    guides: List[Dict[str, object]],
    *,
    guide_seq: str,
    guide_start: int,
    guide_end: int,
    strand: str,
    pam_start: int,
    pam_end: int,
) -> None:
    guide_id = f"g{len(guides) + 1}"
    guides.append(
        {
            "id": guide_id,
            "start": guide_start,
            "end": guide_end,
            "strand": strand,
            "pam_site": {"start": pam_start, "end": pam_end},
            "gc_content": round(_gc_fraction(guide_seq), 4),
            "sequence": guide_seq,
        }
    )


def find_guides(
    seq: str,
    pam: Dict[str, str] | PAM | str,
    guide_len: int,
    *,
    strand: str = "both",
    window: Tuple[int, int] | None = None,
) -> List[Dict[str, object]]:
    """Return candidate guides with their loci, PAM site, strand, and GC fraction."""

    pam_def = _prepare_pam(pam)
    pattern = pam_def.get("pattern")
    if not pattern:
        raise ValueError("PAM definitions require a 'pattern' field.")
    orientation = pam_def.get("orientation", "3prime").lower()
    if orientation != "3prime":
        raise ValueError("Currently only 3prime PAM orientations are supported.")
    if guide_len <= 0:
        raise ValueError("guide_len must be > 0.")
    strand = strand.lower()
    strand_map = {"+": "+", "-": "-", "both": "both"}
    if strand not in strand_map:
        raise ValueError("strand must be one of '+', '-', or 'both'.")

    normalized_seq = bioinformatics.normalize_sequence(seq)
    seq_len = len(normalized_seq)
    if seq_len == 0:
        return []
    window_start, window_end = _normalize_window(window, seq_len)
    if window_end - window_start < guide_len + len(pattern):
        return []

    guides: List[Dict[str, object]] = []
    pam_len = len(pattern)
    rc_pattern = reverse_complement_pattern(pattern)

    search_plus = strand in {"+", "both"}
    search_minus = strand in {"-", "both"}

    if search_plus:
        start_pos = max(window_start + guide_len, 0)
        stop_pos = max(start_pos, window_end - pam_len + 1)
        for pam_pos in range(start_pos, stop_pos):
            if not match_pam(normalized_seq, pam_def, pam_pos):
                continue
            guide_start = pam_pos - guide_len
            guide_end = pam_pos
            if guide_start < window_start:
                continue
            guide_seq = normalized_seq[guide_start:guide_end]
            _record_guide(
                guides,
                guide_seq=guide_seq,
                guide_start=guide_start,
                guide_end=guide_end,
                strand="+",
                pam_start=pam_pos,
                pam_end=pam_pos + pam_len,
            )

    if search_minus:
        rc_pam = {"pattern": rc_pattern}
        start_pos = max(window_start, 0)
        stop_pos = max(start_pos, window_end - (pam_len + guide_len) + 1)
        for pam_pos in range(start_pos, stop_pos):
            if not match_pam(normalized_seq, rc_pam, pam_pos):
                continue
            guide_start = pam_pos + pam_len
            guide_end = guide_start + guide_len
            if guide_end > window_end:
                continue
            pam_start = pam_pos
            pam_end = pam_pos + pam_len
            raw_seq = normalized_seq[guide_start:guide_end]
            guide_seq = bioinformatics.reverse_complement(raw_seq)
            _record_guide(
                guides,
                guide_seq=guide_seq,
                guide_start=guide_start,
                guide_end=guide_end,
                strand="-",
                pam_start=pam_start,
                pam_end=pam_end,
            )

    guides.sort(key=lambda g: (g["start"], g["strand"]))
    return guides
