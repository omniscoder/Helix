"""Prime editing prior profile helpers."""
from __future__ import annotations

from typing import Dict

PRIME_PRIOR_PROFILES: Dict[str, Dict[str, Dict[str, float]]] = {
    "default_indel": {
        "left_flap": {"weight": 0.4},
        "right_flap": {"weight": 0.25},
        "reanneal": {"weight": 0.15},
        "indel_loss": {"weight": 0.15},
        "no_edit": {"weight": 0.05},
    },
    "scarless": {
        "left_flap": {"weight": 0.55},
        "right_flap": {"weight": 0.2},
        "reanneal": {"weight": 0.15},
        "no_edit": {"weight": 0.1},
    },
    "indel_heavy": {
        "indel_loss": {"weight": 0.45},
        "left_flap": {"weight": 0.2},
        "right_flap": {"weight": 0.15},
        "reanneal": {"weight": 0.1},
        "no_edit": {"weight": 0.1},
    },
    "hdr_only": {
        "left_flap": {"weight": 0.6},
        "right_flap": {"weight": 0.3},
        "reanneal": {"weight": 0.1},
    },
}


def resolve_prime_priors(profile: str | None) -> Dict[str, Dict[str, float]]:
    """Return a normalized mapping for the requested profile."""

    key = (profile or "default_indel").strip() or "default_indel"
    try:
        template = PRIME_PRIOR_PROFILES[key]
    except KeyError as exc:
        raise ValueError(f"Unknown prime priors profile: {key}") from exc
    return {label: dict(config) for label, config in template.items()}
