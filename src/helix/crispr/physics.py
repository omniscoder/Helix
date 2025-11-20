"""
CRISPR physics interfaces and backend implementations.

Everything in this module is considered the "engine" side of the CRISPR stack:
it consumes normalized sequences, scores candidate sites, and returns results
to the higher-level simulator facade. Keep the public API limited to
CRISPRPhysicsBase + create_crispr_physics so we can swap the underlying math
for Numba/C++/Rust implementations without changing callers.
"""
from __future__ import annotations

from __future__ import annotations

from abc import ABC, abstractmethod
from collections import OrderedDict
from dataclasses import dataclass
import logging
from typing import Dict, Iterable, List, Mapping, Sequence, Protocol, runtime_checkable

import numpy as np

from .. import bioinformatics
from ..config import native_backend_available, resolve_crispr_backend
from ..engine.encoding import ASCII_TO_BASE, encode_sequence_to_uint8
from .kmm import k_mismatch_positions
from .model import CasSystem, GuideRNA, TargetSite

CRISPR_SCORING_VERSION = "1.0.0"

LOGGER = logging.getLogger(__name__)

_LAST_CRISPR_BACKEND: str | None = None

_IUPAC: dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}

def _encode_sequence_to_uint8(seq: str) -> np.ndarray:
    """Encode an uppercase DNA sequence into tiny uint8 codes (A/C/G/T/N)."""
    return encode_sequence_to_uint8(seq)


def _gc_penalty_encoded(encoded: np.ndarray, lo: float, hi: float) -> float:
    if encoded.size == 0:
        return 0.0
    gc_count = np.count_nonzero((encoded == 1) | (encoded == 2))
    gc_frac = float(gc_count) / float(encoded.size)
    if lo <= gc_frac <= hi:
        return 0.0
    delta = min(abs(gc_frac - lo), abs(gc_frac - hi))
    return delta * 2.0


def _mismatch_cost_encoded(
    guide_encoded: np.ndarray,
    window_encoded: np.ndarray,
    weights: np.ndarray,
    bulge_penalty: float,
) -> float:
    overlap = min(len(window_encoded), len(guide_encoded))
    cost = 0.0
    if overlap:
        mismatches = window_encoded[:overlap] != guide_encoded[:overlap]
        if mismatches.any():
            weight_slice = weights[:overlap]
            cost += float(np.dot(mismatches.astype(np.float32), weight_slice))
    if len(window_encoded) > len(guide_encoded):
        extra = len(window_encoded) - len(guide_encoded)
        cost += 2.0 * extra * bulge_penalty
    elif len(guide_encoded) > len(window_encoded):
        extra = len(guide_encoded) - len(window_encoded)
        cost += extra * bulge_penalty
    return cost


def compute_on_target_score_encoded(
    guide_encoded: np.ndarray,
    window_encoded: np.ndarray,
    weights: np.ndarray,
    *,
    bulge_penalty: float,
    gc_range: tuple[float, float],
    pam_penalty: float,
    min_score: float,
) -> float:
    mismatch_cost = _mismatch_cost_encoded(guide_encoded, window_encoded, weights, bulge_penalty)
    gc_penalty = _gc_penalty_encoded(window_encoded, *gc_range)
    return _compute_score(mismatch_cost, gc_penalty, pam_penalty, len(guide_encoded), min_score)


def score_pairs_encoded(
    guides_encoded: np.ndarray,
    windows_encoded: np.ndarray,
    physics: CrisprPhysics,
    *,
    pam_penalties: np.ndarray | None = None,
) -> np.ndarray:
    """Score batches for a single backend (G guides vs N windows)."""

    guides = np.ascontiguousarray(guides_encoded, dtype=np.uint8)
    windows = np.ascontiguousarray(windows_encoded, dtype=np.uint8)
    if guides.ndim != 2 or windows.ndim != 2:
        raise ValueError("guides_encoded and windows_encoded must be 2D arrays.")
    expected_len = getattr(physics, "_guide_len", None)
    if expected_len is not None and guides.shape[1] != expected_len:
        raise ValueError("Guide encoding length does not match physics instance.")
    if expected_len is not None and windows.shape[1] != expected_len:
        raise ValueError("Window encoding length does not match physics instance.")
    if pam_penalties is not None:
        pam_arr = np.asarray(pam_penalties, dtype=np.float32)
        if pam_arr.shape != (guides.shape[0], windows.shape[0]):
            raise ValueError("pam_penalties shape must match (guides, windows).")
    else:
        pam_arr = np.zeros((guides.shape[0], windows.shape[0]), dtype=np.float32)
    batch_method = getattr(physics, "score_pairs_encoded_batch", None)
    if batch_method is not None:
        return batch_method(guides, windows, pam_arr)
    results = np.zeros((guides.shape[0], windows.shape[0]), dtype=np.float32)
    for row in range(guides.shape[0]):
        for col in range(windows.shape[0]):
            results[row, col] = physics.on_target_score_encoded(
                windows[col],
                pam_penalty=float(pam_arr[row, col]),
            )
    return results


def score_pairs_encoded_multi(
    physics_list: Sequence[CrisprPhysics],
    guides_encoded: np.ndarray,
    windows_encoded: np.ndarray,
    *,
    pam_penalties: np.ndarray | None = None,
) -> np.ndarray:
    """Helper that evaluates multiple guides/backends independently."""

    if not physics_list:
        raise ValueError("physics_list must be non-empty")
    guides = np.ascontiguousarray(guides_encoded, dtype=np.uint8)
    if guides.shape[0] != len(physics_list):
        raise ValueError("guides rows must match physics_list length")
    if pam_penalties is not None:
        pam_arr = np.asarray(pam_penalties, dtype=np.float32)
        if pam_arr.shape != (guides.shape[0], windows_encoded.shape[0]):
            raise ValueError("pam_penalties shape must match (guides, windows)")
    else:
        pam_arr = np.zeros((guides.shape[0], windows_encoded.shape[0]), dtype=np.float32)
    rows = []
    for idx, physics in enumerate(physics_list):
        row = score_pairs_encoded(
            guides[idx : idx + 1],
            windows_encoded,
            physics,
            pam_penalties=pam_arr[idx : idx + 1],
        )
        rows.append(row)
    return np.concatenate(rows, axis=0)


@dataclass
class CRISPRPhysicsResult:
    site: TargetSite
    mismatch_cost: float
    pam_ok: bool
    score: float


@runtime_checkable
class CrisprPhysics(Protocol):
    def score_sites(
        self,
        genome_sequences: Mapping[str, str],
        *,
        max_sites: int | None = None,
    ) -> List[CRISPRPhysicsResult]:
        ...

    def on_target_score_encoded(
        self,
        window_encoded: np.ndarray,
        *,
        pam_penalty: float = 0.0,
    ) -> float:
        ...


def _canonical_backend_name(name: str) -> str:
    lowered = name.lower()
    if lowered in {"cpu", "cpu-reference"}:
        return "cpu-reference"
    if lowered in {"gpu", "gpu-cuda"}:
        return "gpu"
    if lowered == "native-cpu":
        return "native-cpu"
    return lowered


def _record_backend(backend_name: str, physics: CrisprPhysics) -> CrisprPhysics:
    global _LAST_CRISPR_BACKEND
    _LAST_CRISPR_BACKEND = backend_name
    setattr(physics, "backend_name", backend_name)
    return physics


def get_last_crispr_backend() -> str | None:
    """Return the backend used by the most recent create_crispr_physics call."""

    return _LAST_CRISPR_BACKEND


class CRISPRPhysicsBase(CrisprPhysics, ABC):
    """
    Base class for CRISPR physics implementations.
    """

    def __init__(self, cas: CasSystem, guide: GuideRNA):
        self.cas = cas
        self.guide = GuideRNA(
            sequence=bioinformatics.normalize_sequence(guide.sequence),
            pam=guide.pam,
            name=guide.name,
            metadata=dict(guide.metadata),
        )
        if not self.guide.sequence:
            raise ValueError("Guide sequence must be non-empty.")
        self.backend_name = "unknown"

    @abstractmethod
    def score_sites(
        self,
        genome_sequences: Mapping[str, str],
        *,
        max_sites: int | None = None,
    ) -> List[CRISPRPhysicsResult]:
        """
        Return candidate sites ranked by physics score.
        """

    def on_target_score_encoded(
        self,
        window_encoded: np.ndarray,
        *,
        pam_penalty: float = 0.0,
    ) -> float:
        """
        Score a single target window that already satisfies PAM constraints.
        """
        raise NotImplementedError


def _seed_weight_vector(length: int, seed_length: int) -> Sequence[float]:
    weights = [1.0] * length
    seed_len = max(0, min(seed_length, length))
    seed_start = length - seed_len
    for idx in range(seed_start, length):
        weights[idx] = 1.5
    return weights


def _pam_ok(pam_pattern: str, pam_seq: str) -> tuple[bool, float]:
    if not pam_pattern:
        return True, 0.0
    if len(pam_seq) != len(pam_pattern):
        return False, 1.0
    mismatches = 0
    for pattern_char, base in zip(pam_pattern.upper(), pam_seq.upper()):
        allowed = _IUPAC.get(pattern_char, pattern_char)
        if base not in allowed:
            mismatches += 1
            if mismatches > 1:
                return False, 1.0
    return True, float(mismatches)


def _gc_penalty(seq: str, lo: float, hi: float) -> float:
    if not seq:
        return 0.0
    gc = bioinformatics.gc_content(seq)
    if lo <= gc <= hi:
        return 0.0
    delta = min(abs(gc - lo), abs(gc - hi))
    return delta * 2.0


def _base_candidate_positions(
    sequence: str,
    guide_seq: str,
    max_mismatches: int | None,
) -> Iterable[int]:
    limit = len(sequence) - len(guide_seq) + 1
    if limit <= 0:
        return []
    if max_mismatches is None:
        return range(0, limit)
    return k_mismatch_positions(sequence, guide_seq, max_mismatches)


def _pam_filter_positions(
    sequence: str,
    positions: Iterable[int],
    guide_len: int,
    pam_pattern: str,
) -> List[int]:
    pam_len = len(pam_pattern)
    if pam_len == 0:
        return list(positions)
    seq_upper = sequence.upper()
    filtered: List[int] = []
    for start in positions:
        pam_start = start + guide_len
        if pam_start + pam_len > len(seq_upper):
            continue
        pam_seq = seq_upper[pam_start : pam_start + pam_len]
        ok, _ = _pam_ok(pam_pattern, pam_seq)
        if ok:
            filtered.append(start)
    return filtered


def _compute_score(
    mismatch_cost: float,
    gc_penalty: float,
    pam_penalty: float,
    guide_len: int,
    min_score: float,
) -> float:
    total_cost = mismatch_cost + gc_penalty + pam_penalty
    span = max(1.0, float(guide_len))
    score = max(min_score, 1.0 - (total_cost / span))
    return min(1.0, score)


class CRISPRPhysicsCPU(CRISPRPhysicsBase):
    """
    Reference CPU implementation used for correctness and fallback.
    """

    def __init__(self, cas: CasSystem, guide: GuideRNA):
        super().__init__(cas, guide)
        pam_pattern = guide.pam or (cas.pam_rules[0].pattern if cas.pam_rules else "")
        self.pam_pattern = pam_pattern.upper()
        seed_len = min(len(self.guide.sequence), 12)
        weights = tuple(_seed_weight_vector(len(self.guide.sequence), seed_len))
        self._weights = weights
        self._weight_array = np.asarray(weights, dtype=np.float32)
        self._guide_encoded = _encode_sequence_to_uint8(self.guide.sequence)
        self._guide_len = len(self._guide_encoded)
        self._window_cache: OrderedDict[str, np.ndarray] = OrderedDict()
        self._window_cache_limit = 64
        self.seed_length = seed_len
        self.bulge_penalty = max(0.5, cas.weight_mismatch_penalty * 2.0)
        self.pam_weak_penalty = max(0.5, cas.weight_pam_penalty or 1.0)
        self.pam_fail_penalty = max(self.pam_weak_penalty * 2.0, 2.0)
        self.gc_opt_range = (0.4, 0.8)
        self.min_score = 1e-6

    def _encode_window_cached(self, seq: str) -> np.ndarray:
        cached = self._window_cache.get(seq)
        if cached is not None:
            self._window_cache.move_to_end(seq)
            return cached
        encoded = _encode_sequence_to_uint8(seq)
        self._window_cache[seq] = encoded
        self._window_cache.move_to_end(seq)
        if len(self._window_cache) > self._window_cache_limit:
            self._window_cache.popitem(last=False)
        return encoded

    def score_sites(
        self,
        genome_sequences: Mapping[str, str],
        *,
        max_sites: int | None = None,
    ) -> List[CRISPRPhysicsResult]:
        results: List[CRISPRPhysicsResult] = []
        for chrom, seq in genome_sequences.items():
            normalized = bioinformatics.normalize_sequence(seq)
            results.extend(self._score_strand(chrom, normalized, strand=1))
            results.extend(self._score_strand(chrom, normalized, strand=-1))
        results.sort(key=lambda r: r.score, reverse=True)
        if max_sites is not None:
            results = results[:max_sites]
        return results

    def on_target_score_encoded(
        self,
        window_encoded: np.ndarray,
        *,
        pam_penalty: float = 0.0,
    ) -> float:
        return compute_on_target_score_encoded(
            self._guide_encoded,
            window_encoded,
            self._weight_array,
            bulge_penalty=self.bulge_penalty,
            gc_range=self.gc_opt_range,
            pam_penalty=pam_penalty,
            min_score=self.min_score,
        )

    def _score_strand(self, chrom: str, seq: str, strand: int) -> List[CRISPRPhysicsResult]:
        guide_len = len(self.guide.sequence)
        pam_pattern = self.pam_pattern
        max_mm = self.cas.max_mismatches
        if strand == -1:
            seq_to_scan = bioinformatics.reverse_complement(seq)
        else:
            seq_to_scan = seq
        base_positions = _base_candidate_positions(seq_to_scan, self.guide.sequence, max_mm)
        positions = _pam_filter_positions(seq_to_scan, base_positions, guide_len, pam_pattern)
        results: List[CRISPRPhysicsResult] = []
        for start in positions:
            pam_seq = ""
            if pam_pattern:
                pam_start = start + guide_len
                pam_seq = seq_to_scan[pam_start : pam_start + len(pam_pattern)]
            pam_ok, pam_penalty_raw = _pam_ok(pam_pattern, pam_seq)
            if not pam_ok:
                continue
            target_seq = seq_to_scan[start : start + guide_len]
            target_encoded = self._encode_window_cached(target_seq)
            score = self.on_target_score_encoded(
                target_encoded,
                pam_penalty=pam_penalty_raw * self.pam_weak_penalty,
            )
            if score <= 0:
                continue
            if strand == -1:
                orig_start = len(seq) - (start + guide_len)
                orig_end = orig_start + guide_len
                sequence = bioinformatics.reverse_complement(target_seq)
            else:
                orig_start = start
                orig_end = start + guide_len
                sequence = target_seq
            site = TargetSite(
                chrom=chrom,
                start=orig_start,
                end=orig_end,
                strand=strand,
                sequence=sequence,
                on_target_score=score,
            )
            results.append(
                CRISPRPhysicsResult(
                    site=site,
                    mismatch_cost=_mismatch_cost_encoded(
                        self._guide_encoded,
                        target_encoded,
                        self._weight_array,
                        self.bulge_penalty,
                    ),
                    pam_ok=True,
                    score=score,
                )
            )
        return results


class NativeCrisprPhysics(CRISPRPhysicsCPU):
    """Bridge to the pybind11-backed native scoring engine."""

    def __init__(self, cas: CasSystem, guide: GuideRNA):
        from helix_engine import native as native_engine

        if not native_engine.is_available():
            raise RuntimeError(
                "helix_engine._native is unavailable. Build the native Helix engine "
                "to use backend='native-cpu'."
            )
        self._native_engine = native_engine
        super().__init__(cas, guide)

    def on_target_score_encoded(
        self,
        window_encoded: np.ndarray,
        *,
        pam_penalty: float = 0.0,
    ) -> float:
        window_arr = np.ascontiguousarray(window_encoded, dtype=np.uint8)
        return self._native_engine.compute_on_target_score_encoded(
            self._guide_encoded,
            window_arr,
            self._weight_array,
            bulge_penalty=self.bulge_penalty,
            gc_low=self.gc_opt_range[0],
            gc_high=self.gc_opt_range[1],
            pam_penalty=pam_penalty,
            min_score=self.min_score,
        )

    def score_pairs_encoded_batch(
        self,
        guides_encoded: np.ndarray,
        windows_encoded: np.ndarray,
        pam_penalties: np.ndarray,
    ) -> np.ndarray:
        guides_arr = np.ascontiguousarray(guides_encoded, dtype=np.uint8)
        windows_arr = np.ascontiguousarray(windows_encoded, dtype=np.uint8)
        pam_arr = np.ascontiguousarray(pam_penalties, dtype=np.float32)
        if guides_arr.shape[1] != self._guide_len:
            raise ValueError("Guide encoding length does not match native physics guide length.")
        if windows_arr.shape[1] != self._guide_len:
            raise ValueError("Window encoding length does not match native physics guide length.")
        return self._native_engine.score_pairs_encoded(
            guides_arr,
            windows_arr,
            self._weight_array,
            bulge_penalty=self.bulge_penalty,
            gc_low=self.gc_opt_range[0],
            gc_high=self.gc_opt_range[1],
            pam_penalties=pam_arr,
            min_score=self.min_score,
        )


class CudaCrisprPhysics(CRISPRPhysicsCPU):
    """GPU-backed CRISPR physics using the native CUDA kernels."""

    def __init__(self, cas: CasSystem, guide: GuideRNA):
        from helix_engine import native as native_engine

        if not getattr(native_engine, "cuda_available", lambda: False)():
            raise RuntimeError("Helix CUDA backend is unavailable.")
        self._native_engine = native_engine
        super().__init__(cas, guide)

    def score_pairs_encoded_batch(
        self,
        guides_encoded: np.ndarray,
        windows_encoded: np.ndarray,
        pam_penalties: np.ndarray,
    ) -> np.ndarray:
        guides_arr = np.ascontiguousarray(guides_encoded, dtype=np.uint8)
        windows_arr = np.ascontiguousarray(windows_encoded, dtype=np.uint8)
        pam_arr = np.ascontiguousarray(pam_penalties, dtype=np.float32)
        if guides_arr.shape[1] != self._guide_len:
            raise ValueError("Guide encoding length does not match CUDA physics guide length.")
        if windows_arr.shape[1] != self._guide_len:
            raise ValueError("Window encoding length does not match CUDA physics guide length.")
        return self._native_engine.score_pairs_encoded_cuda(
            guides_arr,
            windows_arr,
            self._weight_array,
            bulge_penalty=self.bulge_penalty,
            gc_low=self.gc_opt_range[0],
            gc_high=self.gc_opt_range[1],
            pam_penalties=pam_arr,
            min_score=self.min_score,
        )


def create_crispr_physics(
    cas: CasSystem,
    guide: GuideRNA,
    *,
    use_gpu: bool = False,
    backend: str | None = None,
) -> CrisprPhysics:
    backend, allow_fallback = resolve_crispr_backend(backend, use_gpu)
    requested_backend = _canonical_backend_name(backend)
    choice = backend
    while True:
        if choice in {"gpu", "gpu-cuda"}:
            try:
                physics = CudaCrisprPhysics(cas, guide)
            except RuntimeError as exc:
                if not allow_fallback:
                    raise
                LOGGER.warning(
                    "Falling back from backend='%s' to cpu-reference: %s",
                    requested_backend,
                    exc,
                )
                choice = "cpu-reference"
                continue
            return _record_backend("gpu", physics)
        if choice == "native-cpu":
            if not native_backend_available():
                if allow_fallback:
                    LOGGER.warning(
                        "Falling back from backend='%s' to cpu-reference: native engine unavailable.",
                        requested_backend,
                    )
                    choice = "cpu-reference"
                    continue
                raise RuntimeError(
                    "Native CRISPR engine is unavailable. Build helix_engine._native or choose "
                    "backend='cpu-reference'."
                )
            return _record_backend("native-cpu", NativeCrisprPhysics(cas, guide))
        if choice in {"cpu-reference", "cpu"}:
            return _record_backend("cpu-reference", CRISPRPhysicsCPU(cas, guide))
        raise ValueError(f"Unknown CRISPR physics backend: {choice}")
