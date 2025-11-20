"""Helix runtime configuration helpers."""

from __future__ import annotations

import os
import logging
from functools import lru_cache

_BACKEND_ENV = "HELIX_CRISPR_BACKEND"
_ALLOW_FALLBACK_ENV = "HELIX_CRISPR_ALLOW_FALLBACK"
_PRIME_BACKEND_ENV = "HELIX_PRIME_BACKEND"
_PRIME_ALLOW_FALLBACK_ENV = "HELIX_PRIME_ALLOW_FALLBACK"

LOGGER = logging.getLogger(__name__)


def _env_backend(env_name: str) -> str | None:
    value = os.getenv(env_name)
    if not value:
        return None
    return value.strip().lower() or None


def _env_bool(name: str, default: bool = False) -> bool:
    raw = os.getenv(name)
    if raw is None:
        return default
    raw = raw.strip().lower()
    if raw in {"", "0", "false", "no"}:
        return False
    if raw in {"1", "true", "yes"}:
        return True
    return default


@lru_cache(maxsize=1)
def native_backend_available() -> bool:
    try:
        from helix_engine import native

        return native.is_available()
    except Exception:
        return False


def resolve_crispr_backend(preferred: str | None, use_gpu: bool) -> tuple[str, bool]:
    """Resolve the backend requested by CLI/env/config."""

    backend = preferred or _env_backend(_BACKEND_ENV)
    allow_fallback = _env_bool(_ALLOW_FALLBACK_ENV)
    if backend is None:
        if use_gpu:
            backend = "gpu"
        else:
            backend = "native-cpu" if native_backend_available() else "cpu-reference"
    backend = backend.lower()
    LOGGER.debug(
        "resolve_crispr_backend backend=%s allow_fallback=%s use_gpu=%s env=%s native_available=%s",
        backend,
        allow_fallback,
        use_gpu,
        preferred or _env_backend(_BACKEND_ENV),
        native_backend_available(),
    )
    return backend, allow_fallback


def resolve_prime_backend(preferred: str | None, use_gpu: bool) -> tuple[str, bool]:
    backend = preferred or _env_backend(_PRIME_BACKEND_ENV)
    allow_fallback = _env_bool(_PRIME_ALLOW_FALLBACK_ENV)
    if backend is None:
        return resolve_crispr_backend(None, use_gpu)
    backend = backend.lower()
    LOGGER.debug(
        "resolve_prime_backend backend=%s allow_fallback=%s use_gpu=%s env=%s",
        backend,
        allow_fallback,
        use_gpu,
        preferred or _env_backend(_PRIME_BACKEND_ENV),
    )
    return backend, allow_fallback


__all__ = [
    "resolve_crispr_backend",
    "resolve_prime_backend",
    "native_backend_available",
    "_BACKEND_ENV",
    "_ALLOW_FALLBACK_ENV",
    "_PRIME_BACKEND_ENV",
    "_PRIME_ALLOW_FALLBACK_ENV",
]
