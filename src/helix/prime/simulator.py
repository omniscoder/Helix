"""
Prime editing sequence-level simulation helpers.

All operations here work on digital sequences only and are not
wet-lab protocols.
"""
from __future__ import annotations

import hashlib
import json
import logging
import math
import random
import uuid
from bisect import bisect_left
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple

from helix import bioinformatics
from helix.crispr.model import DigitalGenome, TargetSite

from .model import PegRNA, PrimeEditor
from .physics import PrimePhysics
from .priors import resolve_prime_priors

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class PrimeTargetRequest:
    target_id: str
    genome: DigitalGenome
    peg: PegRNA
    editor: PrimeEditor
    search_window: int = 200


@dataclass(frozen=True)
class PrimePrediction:
    target_id: str
    site: TargetSite | None
    predicted_efficiency: float


def _normalize_sequence(value: str, *, allow_ambiguous: bool = True) -> str:
    return bioinformatics.normalize_sequence(value, allow_ambiguous=allow_ambiguous)


def _find_exact_match(sequence: str, query: str) -> List[int]:
    positions: List[int] = []
    start = 0
    while True:
        idx = sequence.find(query, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions


def _best_approximate_match(sequence: str, query: str, limit: int) -> Optional[Tuple[int, int]]:
    best_idx: Optional[int] = None
    best_mismatches: Optional[int] = None
    max_start = max(0, min(len(sequence) - len(query), limit))
    for idx in range(max_start + 1):
        window = sequence[idx : idx + len(query)]
        if len(window) != len(query):
            continue
        mismatches = sum(a != b for a, b in zip(window, query))
        if best_mismatches is None or mismatches < best_mismatches:
            best_idx = idx
            best_mismatches = mismatches
            if mismatches == 0:
                break
    if best_idx is None or best_mismatches is None:
        return None
    return best_idx, best_mismatches


def locate_prime_target_site(
    genome: DigitalGenome,
    peg: PegRNA,
    *,
    search_window: int = 200,
) -> Optional[TargetSite]:
    """
    Identify a primary digital target site for a pegRNA.

    Returns None if no plausible site is found.
    """

    spacer = _normalize_sequence(peg.spacer)
    if not spacer:
        return None
    rc_spacer = bioinformatics.reverse_complement(spacer)

    for chrom, raw_seq in genome.sequences.items():
        sequence = _normalize_sequence(raw_seq)
        if not sequence:
            continue
        for idx in _find_exact_match(sequence, spacer):
            site = TargetSite(
                chrom=chrom,
                start=idx,
                end=idx + len(spacer),
                strand=1,
                sequence=spacer,
                on_target_score=1.0,
            )
            return site
        for idx in _find_exact_match(sequence, rc_spacer):
            site = TargetSite(
                chrom=chrom,
                start=idx,
                end=idx + len(spacer),
                strand=-1,
                sequence=spacer,
                on_target_score=1.0,
            )
            return site

    # Fallback: approximate match within the provided window.
    limit = max(search_window, 0)
    best_site: Optional[TargetSite] = None
    best_mismatches: Optional[int] = None
    for chrom, raw_seq in genome.sequences.items():
        sequence = _normalize_sequence(raw_seq)
        if not sequence:
            continue
        approx = _best_approximate_match(sequence, spacer, limit)
        if approx:
            idx, mismatches = approx
            if best_mismatches is None or mismatches < best_mismatches:
                best_mismatches = mismatches
                best_site = TargetSite(
                    chrom=chrom,
                    start=idx,
                    end=idx + len(spacer),
                    strand=1,
                    sequence=sequence[idx : idx + len(spacer)],
                    on_target_score=max(0.0, 1.0 - mismatches / len(spacer)),
                )
        approx_rc = _best_approximate_match(sequence, rc_spacer, limit)
        if approx_rc:
            idx, mismatches = approx_rc
            if best_mismatches is None or mismatches < best_mismatches:
                seq_slice = sequence[idx : idx + len(spacer)]
                best_mismatches = mismatches
                best_site = TargetSite(
                    chrom=chrom,
                    start=idx,
                    end=idx + len(spacer),
                    strand=-1,
                    sequence=bioinformatics.reverse_complement(seq_slice),
                    on_target_score=max(0.0, 1.0 - mismatches / len(spacer)),
                )
    return best_site


def _extract_window(sequence: str, start: int, end: int) -> str:
    return sequence[start:end]


def apply_rtt_edit(seq: str, site_start: int, site_end: int, rtt: str) -> str:
    """Return a copy of `seq` with [site_start:site_end] replaced by rtt."""
    return seq[:site_start] + rtt + seq[site_end:]


def _apply_rtt(reference: str, offset: int, template: str) -> str:
    if not reference:
        return template
    offset = max(0, min(len(reference), offset))
    window = list(reference)
    for idx, base in enumerate(template):
        pos = offset + idx
        if pos >= len(window):
            break
        window[pos] = base
    return "".join(window)


def _sequence_sha(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def _normalize_prime_priors(priors: Mapping[str, Mapping[str, float]] | None) -> Dict[str, Dict[str, float]]:
    if priors is None:
        priors = resolve_prime_priors("default_indel")
    return {label: dict(config) for label, config in priors.items()}


def _derive_label_seed(seed: int | None, label: str) -> int:
    base = f"{seed}:{label}".encode("utf-8")
    return int.from_bytes(hashlib.blake2b(base, digest_size=8).digest(), "big")


def _peg_field(peg: PegRNA | Mapping[str, object], name: str, default=None):
    if isinstance(peg, Mapping):
        if name in peg:
            return peg[name]
        metadata = peg.get("metadata") if isinstance(peg.get("metadata"), Mapping) else None
        if metadata and name in metadata:
            return metadata[name]
        return default
    if hasattr(peg, name):
        return getattr(peg, name)
    if hasattr(peg, "metadata") and isinstance(peg.metadata, dict):
        return peg.metadata.get(name, default)
    return default


def _nick_position(peg: PegRNA | Mapping[str, object], site_len: int) -> int:
    idx = _peg_field(peg, "nick_index")
    if isinstance(idx, (int, float)):
        nick_idx = int(idx)
    else:
        spacer = str(_peg_field(peg, "spacer", "") or "")
        nick_idx = len(spacer)
    return max(0, min(site_len, nick_idx))


def _sample_prime_diff(
    label: str,
    rng: random.Random,
    nick_pos: int,
    site_len: int,
    config: Mapping[str, float],
    peg: PegRNA | Mapping[str, object],
) -> Dict[str, object] | None:
    rtt_len = len(str(_peg_field(peg, "rtt", "") or ""))
    pbs_len = len(str(_peg_field(peg, "pbs", "") or ""))

    if label in {"left_flap", "scarless_hdr"}:
        length = max(1, int(config.get("length", rtt_len or 1)))
        start = max(0, nick_pos - length)
        end = min(site_len, start + length)
        return {"start": start, "end": end, "edit": "left_flap"}
    if label == "right_flap":
        length = max(1, int(config.get("length", max(1, rtt_len // 2))))
        start = nick_pos
        end = min(site_len, start + length)
        return {"start": start, "end": end, "edit": "right_flap"}
    if label == "indel_loss":
        size = rng.randint(1, max(2, pbs_len or 4))
        start = max(0, nick_pos - size // 2)
        end = min(site_len, start + size)
        return {"start": start, "end": end, "edit": f"del{size}"}
    if label in {"reanneal", "no_edit"}:
        return None
    if label == "right_partial":
        length = max(1, rtt_len // 3 or 1)
        return {"start": nick_pos, "end": min(site_len, nick_pos + length), "edit": "right_partial"}
    return {"start": nick_pos, "end": nick_pos, "edit": label}


def simulate_prime_edit(
    site_seq: str,
    peg: PegRNA | Mapping[str, object],
    priors: Mapping[str, Mapping[str, float]] | None = None,
    *,
    draws: int = 1000,
    seed: int | None = None,
    emit_sequence: bool = False,
    pam_mask: Sequence[float] | None = None,
) -> Dict[str, object]:
    """Sample a multinomial outcome distribution for PRIME editing."""

    if draws <= 0:
        raise ValueError("draws must be > 0.")
    normalized_site = bioinformatics.normalize_sequence(site_seq)
    if not normalized_site:
        raise ValueError("site_seq must contain DNA bases.")
    site_len = len(normalized_site)

    normalized_priors = _normalize_prime_priors(priors)
    if not normalized_priors:
        raise ValueError("prime priors is empty after normalization.")

    nick_pos = _nick_position(peg, site_len)

    pam_factor = 1.0
    if pam_mask is not None:
        if len(pam_mask) != site_len:
            raise ValueError(f"pam_mask length {len(pam_mask)} does not match site length {site_len}.")
        pam_factor = max(0.0, min(1.0, float(pam_mask[nick_pos])))
        if not math.isfinite(pam_factor):
            raise ValueError("pam_mask values must be finite numbers.")

    labels: List[str] = []
    weights: List[float] = []
    for label, cfg in normalized_priors.items():
        weight = float(cfg.get("weight", 0.0))
        if not (math.isfinite(weight) and weight >= 0.0):
            raise ValueError(f"Prime prior '{label}' has invalid weight {weight}.")
        labels.append(label)
        weights.append(weight * pam_factor)

    total_weight = sum(weights)
    if total_weight <= 0.0:
        raise ValueError(
            f"Effective prime weights at nick_pos={nick_pos} are zero; "
            f"PAM profile forbids editing at this site."
        )

    cumulative: List[float] = []
    acc = 0.0
    for weight in weights:
        acc += weight
        cumulative.append(acc)

    rng_counts = random.Random(seed)
    diff_rngs = {label: random.Random(_derive_label_seed(seed, f"prime:{label}")) for label in labels}

    counts = {label: 0 for label in labels}
    diffs: Dict[str, Dict[str, object] | None] = {label: None for label in labels}

    for draw_idx in range(draws):
        r = rng_counts.random() * total_weight
        idx = bisect_left(cumulative, r)
        if idx >= len(labels):
            idx = len(labels) - 1
        label = labels[idx]
        counts[label] += 1
        if diffs[label] is None:
            diffs[label] = _sample_prime_diff(label, diff_rngs[label], nick_pos, site_len, normalized_priors[label], peg)
        if LOGGER.isEnabledFor(logging.DEBUG) and (draw_idx + 1) % 10000 == 0:
            LOGGER.debug("simulate_prime_edit progress: %s/%s draws", draw_idx + 1, draws)

    site_snapshot = {
        "chrom": str(_peg_field(peg, "chrom", "chr")),
        "start": 0,
        "end": site_len,
        "strand": int(_peg_field(peg, "strand", 1) or 1),
        "sequence_sha256": _sequence_sha(normalized_site),
    }
    if emit_sequence:
        site_snapshot["sequence"] = normalized_site

    outcomes = []
    for idx, label in enumerate(labels):
        count = counts[label]
        probability = count / draws if draws else 0.0
        entry = {
            "label": label,
            "count": count,
            "probability": round(probability, 4),
            "diff": diffs[label],
            "site": dict(site_snapshot),
            "edited_sequence": normalized_site,
            "logit_score": weights[idx],
        }
        outcomes.append(
            entry
        )
    outcomes.sort(key=lambda entry: (-entry["count"], entry["label"]))

    site_block: Dict[str, object] = {"length": site_len, "sequence_sha256": _sequence_sha(normalized_site)}
    if emit_sequence:
        site_block["sequence"] = normalized_site

    spacer = str(_peg_field(peg, "spacer", "") or "")
    pbs = str(_peg_field(peg, "pbs", "") or "")
    rtt = str(_peg_field(peg, "rtt", "") or "")
    peg_block: Dict[str, object] = {
        "name": _peg_field(peg, "name"),
        "spacer": spacer,
        "pbs": pbs,
        "rtt": rtt,
        "pbs_length": len(pbs),
        "rtt_length": len(rtt),
    }
    metadata = _peg_field(peg, "metadata")
    if isinstance(metadata, dict) and metadata:
        peg_block["metadata"] = metadata

    run_id = uuid.uuid4().hex
    payload = {
        "schema": {"kind": "prime.edit_sim", "spec_version": "1.0"},
        "meta": {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "seed": seed,
            "run_id": run_id,
            "pam_mask_scale": pam_factor,
        },
        "site": site_block,
        "peg": peg_block,
        "priors": normalized_priors,
        "draws": draws,
        "outcomes": outcomes,
    }
    LOGGER.debug(
        "simulate_prime_edit completed: peg=%s top_outcome=%s run_id=%s",
        peg_block.get("name"),
        outcomes[0]["label"] if outcomes else None,
        run_id,
    )
    return json.loads(json.dumps(payload))
def predict_prime_outcomes_for_targets(requests: Sequence[PrimeTargetRequest]) -> List[PrimePrediction]:
    predictions: List[PrimePrediction] = []
    for request in requests:
        site = locate_prime_target_site(request.genome, request.peg, search_window=request.search_window)
        physics = PrimePhysics(editor=request.editor, peg=request.peg)
        if site is not None:
            mismatches = physics.mismatch_count(site.sequence)
        else:
            mismatches = physics.mismatch_count(request.peg.spacer)
        efficiency = physics.binding_efficiency(mismatches)
        predictions.append(
            PrimePrediction(
                target_id=request.target_id,
                site=site,
                predicted_efficiency=efficiency,
            )
        )
    return predictions
