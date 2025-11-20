"""Python shim for the native Helix physics engine.

The actual implementation lives in :mod:`helix_engine._native`, a pybind11
extension built from the C++ sources under ``helix_engine/``. Until that
extension is compiled, this module falls back to the pure-Python scoring
routine so callers can continue to run tests without native artifacts.
"""
from __future__ import annotations

from typing import Any

import numpy as np

try:  # pragma: no cover - depends on local toolchain
    from . import _native  # type: ignore[attr-defined]
except Exception:  # pragma: no cover - optional extension
    _native = None


def is_available() -> bool:
    """Return True if the compiled native module is importable."""

    return _native is not None


def compute_on_target_score_encoded(
    guide_encoded: np.ndarray,
    window_encoded: np.ndarray,
    weights: np.ndarray,
    *,
    bulge_penalty: float,
    gc_low: float,
    gc_high: float,
    pam_penalty: float,
    min_score: float,
) -> float:
    """Call into the native engine (or Python fallback) with encoded arrays."""

    guide_arr = np.ascontiguousarray(guide_encoded, dtype=np.uint8)
    window_arr = np.ascontiguousarray(window_encoded, dtype=np.uint8)
    weights_arr = np.ascontiguousarray(weights, dtype=np.float32)
    if _native is None:
        from helix.crispr.physics import compute_on_target_score_encoded as _py_compute

        return float(
            _py_compute(
                guide_arr,
                window_arr,
                weights_arr,
                bulge_penalty=bulge_penalty,
                gc_range=(gc_low, gc_high),
                pam_penalty=pam_penalty,
                min_score=min_score,
            )
        )
    return float(
        _native.compute_on_target_score_encoded(
            guide_arr,
            window_arr,
            weights_arr,
            bulge_penalty,
            gc_low,
            gc_high,
            pam_penalty,
            min_score,
        )
    )


def score_pairs_encoded(
    guides_encoded: np.ndarray,
    windows_encoded: np.ndarray,
    weights: np.ndarray,
    *,
    bulge_penalty: float,
    gc_low: float,
    gc_high: float,
    pam_penalties: np.ndarray | None,
    min_score: float,
) -> np.ndarray:
    guides_arr = np.ascontiguousarray(guides_encoded, dtype=np.uint8)
    windows_arr = np.ascontiguousarray(windows_encoded, dtype=np.uint8)
    weights_arr = np.ascontiguousarray(weights, dtype=np.float32)
    pam_arr = None if pam_penalties is None else np.ascontiguousarray(pam_penalties, dtype=np.float32)
    if _native is None:
        from helix.crispr.physics import compute_on_target_score_encoded as _py_compute

        n_guides, n_windows = guides_arr.shape[0], windows_arr.shape[0]
        scores = np.zeros((n_guides, n_windows), dtype=np.float32)
        for gi in range(n_guides):
            for wi in range(n_windows):
                pam_penalty = 0.0
                if pam_arr is not None:
                    pam_penalty = float(pam_arr[gi, wi])
                scores[gi, wi] = _py_compute(
                    guides_arr[gi],
                    windows_arr[wi],
                    weights_arr,
                    bulge_penalty=bulge_penalty,
                    gc_range=(gc_low, gc_high),
                    pam_penalty=pam_penalty,
                    min_score=min_score,
                )
        return scores

    result = _native.score_pairs_encoded(
        guides_arr,
        windows_arr,
        weights_arr,
        bulge_penalty,
        gc_low,
        gc_high,
        pam_arr if pam_arr is not None else None,
        min_score,
    )
    return np.asarray(result, dtype=np.float32)


def cuda_available() -> bool:
    if _native is None or not hasattr(_native, "cuda_available"):
        return False
    try:
        return bool(_native.cuda_available())
    except Exception:
        return False


def score_pairs_encoded_cuda(
    guides_encoded: np.ndarray,
    windows_encoded: np.ndarray,
    weights: np.ndarray,
    *,
    bulge_penalty: float,
    gc_low: float,
    gc_high: float,
    pam_penalties: np.ndarray,
    min_score: float,
) -> np.ndarray:
    if _native is None or not hasattr(_native, "score_pairs_encoded_cuda"):
        raise RuntimeError("Helix CUDA backend is unavailable; build helix_engine with CUDA support.")
    guides_arr = np.ascontiguousarray(guides_encoded, dtype=np.uint8)
    windows_arr = np.ascontiguousarray(windows_encoded, dtype=np.uint8)
    weights_arr = np.ascontiguousarray(weights, dtype=np.float32)
    pam_arr = np.ascontiguousarray(pam_penalties, dtype=np.float32)
    result = _native.score_pairs_encoded_cuda(
        guides_arr,
        windows_arr,
        weights_arr,
        bulge_penalty,
        gc_low,
        gc_high,
        pam_arr,
        min_score,
    )
    return np.asarray(result, dtype=np.float32)


__all__ = [
    "compute_on_target_score_encoded",
    "score_pairs_encoded",
    "is_available",
    "cuda_available",
    "score_pairs_encoded_cuda",
]
