"""CRISPR cut/repair simulation."""
from __future__ import annotations

import hashlib
import json
import random
from datetime import datetime, timezone
from typing import Dict, Mapping, Tuple

from .. import bioinformatics

DEFAULT_PRIORS = {
    "no_cut": {"weight": 0.55},
    "small_insertion": {"weight": 0.18, "min": 1, "max": 2},
    "small_deletion": {"weight": 0.17, "min": 1, "max": 5},
    "large_deletion": {"weight": 0.10, "min": 6, "max": 20},
}


def _sequence_sha(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def _normalize_priors(priors: Mapping[str, Mapping[str, float]] | None) -> Dict[str, Dict[str, float]]:
    priors = priors or DEFAULT_PRIORS
    normalized: Dict[str, Dict[str, float]] = {}
    total = 0.0
    for label, config in priors.items():
        weight = float(config.get("weight", 0.0))
        if weight < 0:
            raise ValueError("Prior weights must be >= 0.")
        total += weight
        normalized[label] = dict(config)
    if total <= 0:
        raise ValueError("Provide priors with positive total weight.")
    for config in normalized.values():
        config["weight"] = config["weight"] / total
    return normalized


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


def simulate_cut_repair(
    site_seq: str,
    guide: Mapping[str, object],
    priors: Mapping[str, Mapping[str, float]] | None = None,
    *,
    draws: int = 1000,
    seed: int | None = None,
    emit_sequence: bool = False,
) -> Dict[str, object]:
    """Sample a simple multinomial outcome distribution for CRISPR cut/repair."""

    if draws <= 0:
        raise ValueError("draws must be > 0.")
    normalized_site = bioinformatics.normalize_sequence(site_seq)
    site_len = len(normalized_site)
    normalized_priors = _normalize_priors(priors)
    rng = random.Random(seed)
    cut_pos = _cut_position(guide)

    labels, weights = zip(*[(label, config["weight"]) for label, config in normalized_priors.items()])
    cumulative = []
    total = 0.0
    for weight in weights:
        total += weight
        cumulative.append(total)

    counts = {label: 0 for label in labels}
    diffs: Dict[str, Dict[str, object] | None] = {label: None for label in labels}
    for _ in range(draws):
        r = rng.random()
        for label, threshold in zip(labels, cumulative):
            if r <= threshold:
                counts[label] += 1
                if diffs[label] is None:
                    diffs[label] = _sample_diff(label, rng, cut_pos, site_len, normalized_priors[label])
                break

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
    outcomes.sort(key=lambda entry: entry["count"], reverse=True)

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
        "meta": {"timestamp": datetime.now(timezone.utc).isoformat(), "seed": seed},
        "site": site_block,
        "guide": guide_block,
        "priors": normalized_priors,
        "draws": draws,
        "outcomes": outcomes,
    }
    return json.loads(json.dumps(payload))
