"""Registry of preconfigured multiscale slices."""

from __future__ import annotations

from pathlib import Path

SLICE_REGISTRY = {
    "egfr_grb2": Path("models/egfr_grb2/config.yaml"),
}


def resolve_slice(name: str) -> Path:
    """Return the config path for a named slice."""

    if name not in SLICE_REGISTRY:
        available = ", ".join(sorted(SLICE_REGISTRY))
        raise KeyError(f"Unknown slice '{name}'. Available: {available or 'none'}")
    return SLICE_REGISTRY[name]
