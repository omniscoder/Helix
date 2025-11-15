"""CRISPR cut/repair simulation."""
from __future__ import annotations

import hashlib
import json
import logging
import math
import random
import uuid
from collections.abc import Mapping, Sequence
from bisect import bisect_left
from datetime import datetime, timezone
from typing import Dict, Tuple

from .. import bioinformatics

LOGGER = logging.getLogger(__name__)

DEFAULT_PRIORS = {
    "no_cut": {"weight": 0.55},
    "small_insertion": {"weight": 0.18, "min": 1, "max": 2},
    "small_deletion": {"weight": 0.17, "min": 1, "max": 5},
    "large_deletion": {"weight": 0.10, "min": 6, "max": 20},
}

CRISPR_PRIOR_PROFILES: Dict[str, Dict[str, Dict[str, float]]] = {
    "default_indel": DEFAULT_PRIORS,
    "scarless": {
        "no_cut": {"weight": 0.7},
        "small_insertion": {"weight": 0.15, "min": 1, "max": 1},
        "small_deletion": {"weight": 0.1, "min": 1, "max": 3},
        "microhomology_del": {"weight": 0.05, "min": 3, "max": 10},
    },
    "indel_heavy": {
        "small_insertion": {"weight": 0.35, "min": 1, "max": 4},
        "small_deletion": {"weight": 0.35, "min": 1, "max": 8},
        "large_deletion": {"weight": 0.25, "min": 8, "max": 40},
        "no_cut": {"weight": 0.05},
    },
    "hdr_only": {
        "no_cut": {"weight": 0.1},
        "scarless_hdr": {"weight": 0.9},
    },
}


def _sequence_sha(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def _normalize_priors(priors: Mapping[str, Mapping[str, float]] | None) -> Dict[str, Dict[str, float]]:
    priors = priors or DEFAULT_PRIORS
    return {label: dict(config) for label, config in priors.items()}


def _derive_label_seed(seed: int | None, label: str) -> int:
    base = f"{seed}:{label}".encode("utf-8")
    return int.from_bytes(hashlib.blake2b(base, digest_size=8).digest(), "big")


def _cut_position(guide: Mapping[str, object]) -> int:
    start = int(guide.get("start", 0))
    end = int(guide.get("end", start))
    strand = guide.get("strand", "+")
    if strand == "+":
        return max(start, end - 3)
    return min(end, start + 3)


def _sample_diff(label: str, rng: random.Random, cut_pos: int, site_len: int, config: Mapping[str, float]) -> Dict[str, object] | None:
    if label == "no_cut":
        return None
    if label == "small_insertion":
        length = rng.randint(int(config.get("min", 1)), int(config.get("max", 2)))
        return {"start": cut_pos, "end": cut_pos, "edit": f"ins{length}"}
    if label in {"small_deletion", "large_deletion"}:
        length = rng.randint(int(config.get("min", 1)), int(config.get("max", 10)))
        start = max(0, cut_pos - length // 2)
        end = min(site_len, start + length)
        return {"start": start, "end": end, "edit": f"del{length}"}
    return {"start": cut_pos, "end": cut_pos, "edit": label}


def resolve_crispr_priors(profile: str | None) -> Dict[str, Dict[str, float]]:
    """Map a profile name to a normalized priors configuration."""

    key = (profile or "default_indel").strip() or "default_indel"
    try:
        template = CRISPR_PRIOR_PROFILES[key]
    except KeyError as exc:
        raise ValueError(f"Unknown CRISPR priors profile: {key}") from exc
    return {label: dict(config) for label, config in template.items()}


def simulate_cut_repair(
    site_seq: str,
    guide: Mapping[str, object],
    priors: Mapping[str, Mapping[str, float]] | None = None,
    *,
    draws: int = 1000,
    seed: int | None = None,
    emit_sequence: bool = False,
    pam_mask: Sequence[float] | None = None,
) -> Dict[str, object]:
    """Sample a simple multinomial outcome distribution for CRISPR cut/repair."""

    if draws <= 0:
        raise ValueError("draws must be > 0.")

    run_id = uuid.uuid4().hex
    LOGGER.debug(
        "simulate_cut_repair start: guide=%s draws=%s seed=%s run_id=%s",
        guide.get("id"),
        draws,
        seed,
        run_id,
    )
    normalized_site = bioinformatics.normalize_sequence(site_seq)
    site_len = len(normalized_site)
    normalized_priors = _normalize_priors(priors)
    if not normalized_priors:
        raise ValueError("priors must contain at least one entry.")

    cut_pos = _cut_position(guide)
    if not (0 <= cut_pos <= site_len):
        raise ValueError(f"Computed cut position {cut_pos} invalid for sequence length {site_len}.")

    pam_factor = 1.0
    if pam_mask is not None:
        if len(pam_mask) != site_len:
            raise ValueError(f"pam_mask length {len(pam_mask)} does not match site length {site_len}.")
        pam_factor = max(0.0, min(1.0, float(pam_mask[cut_pos])))
        if not math.isfinite(pam_factor):
            raise ValueError("pam_mask values must be finite.")

    labels: list[str] = []
    weights: list[float] = []
    for label, config in normalized_priors.items():
        if "weight" not in config:
            raise ValueError(f"Prior '{label}' missing 'weight'.")
        weight = float(config["weight"])
        if not (math.isfinite(weight) and weight >= 0.0):
            raise ValueError(f"Prior '{label}' has invalid weight {weight}.")
        labels.append(label)
        weights.append(weight)

    total_weight = sum(weights)
    if total_weight <= 0.0:
        raise ValueError("Sum of prior weights must be > 0.")

    efficacies = [weight * pam_factor for weight in weights]
    total_efficacy = sum(efficacies)
    if total_efficacy <= 0.0:
        raise ValueError(
            f"Effective weights at cut_pos={cut_pos} are zero; PAM profile forbids cutting at this site."
        )

    cumulative: list[float] = []
    acc = 0.0
    for value in efficacies:
        acc += value
        cumulative.append(acc)

    rng_counts = random.Random(seed)
    diff_rngs = {label: random.Random(_derive_label_seed(seed, label)) for label in labels}

    counts = {label: 0 for label in labels}
    diffs: Dict[str, Dict[str, object] | None] = {label: None for label in labels}
    for draw_idx in range(draws):
        r = rng_counts.random() * total_efficacy
        idx = bisect_left(cumulative, r)
        if idx >= len(labels):
            idx = len(labels) - 1
        label = labels[idx]
        counts[label] += 1
        if diffs[label] is None:
            diffs[label] = _sample_diff(label, diff_rngs[label], cut_pos, site_len, normalized_priors[label])

        if LOGGER.isEnabledFor(logging.DEBUG) and (draw_idx + 1) % 1000 == 0:
            LOGGER.debug(
                "simulate_cut_repair progress: %s/%s draws (run_id=%s)",
                draw_idx + 1,
                draws,
                run_id,
            )

    outcomes = []
    for label in labels:
        count = counts[label]
        probability = count / draws if draws else 0.0
        outcomes.append(
            {
                "label": label,
                "count": count,
                "probability": round(probability, 4),
                "diff": diffs[label],
            }
        )
    outcomes.sort(key=lambda entry: (-entry["count"], entry["label"]))

    site_block = {"length": site_len, "sequence_sha256": _sequence_sha(normalized_site)}
    if emit_sequence:
        site_block["sequence"] = normalized_site

    guide_block = {
        "id": guide.get("id"),
        "start": guide.get("start"),
        "end": guide.get("end"),
        "strand": guide.get("strand"),
        "gc_content": guide.get("gc_content"),
    }
    if emit_sequence:
        guide_block["sequence"] = guide.get("sequence")

    payload = {
        "schema": {"kind": "crispr.sim", "spec_version": "1.0"},
        "meta": {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "seed": seed,
            "run_id": run_id,
            "weights_sum": total_efficacy,
            "pam_mask_scale": pam_factor,
        },
        "site": site_block,
        "guide": guide_block,
        "priors": normalized_priors,
        "draws": draws,
        "outcomes": outcomes,
    }
    LOGGER.debug(
        "simulate_cut_repair completed: guide=%s top_outcome=%s run_id=%s",
        guide.get("id"),
        outcomes[0]["label"] if outcomes else None,
        run_id,
    )
    return json.loads(json.dumps(payload))
