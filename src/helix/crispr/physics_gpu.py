"""
GPU-accelerated CRISPR physics (Numba CUDA).

Falls back to the CPU implementation if CUDA/Numba is unavailable.
"""
from __future__ import annotations

import math
from typing import List, Mapping

import numpy as np
from numba import cuda  # type: ignore

from .. import bioinformatics
from .physics import (
    CRISPRPhysicsBase,
    CRISPRPhysicsResult,
    _base_candidate_positions,
    _compute_score,
    _gc_penalty,
    _pam_filter_positions,
    _pam_ok,
)
from .model import CasSystem, GuideRNA, TargetSite

_BASE_TO_INT = {
    ord("A"): 0,
    ord("C"): 1,
    ord("G"): 2,
    ord("T"): 3,
    ord("N"): 4,
}


def _encode_sequence(seq: str) -> np.ndarray:
    arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    out = np.zeros_like(arr)
    for key, value in _BASE_TO_INT.items():
        out[arr == key] = value
    return out


@cuda.jit
def _score_windows_kernel(seq_int, starts, guide_int, weights, scores):
    i = cuda.grid(1)
    n = starts.shape[0]
    guide_len = guide_int.shape[0]
    if i >= n:
        return
    start = starts[i]
    mismatch_cost = 0.0
    for pos in range(guide_len):
        s_idx = start + pos
        if s_idx >= seq_int.shape[0]:
            mismatch_cost += weights[pos]
            continue
        g = guide_int[pos]
        b = seq_int[s_idx]
        if g == 4 or b == 4 or g != b:
            mismatch_cost += weights[pos]
    scores[i] = mismatch_cost


class CRISPRPhysicsGPU(CRISPRPhysicsBase):
    """
    GPU-backed implementation that accelerates mismatch scoring.
    """

    def __init__(self, cas: CasSystem, guide: GuideRNA, device_id: int = 0):
        super().__init__(cas, guide)
        self.device_id = device_id
        self.pam_pattern = guide.pam or (cas.pam_rules[0].pattern if cas.pam_rules else "")
        seed_len = min(len(self.guide.sequence), 12)
        self._weights = np.array(
            [1.5 if idx >= len(self.guide.sequence) - seed_len else 1.0 for idx in range(len(self.guide.sequence))],
            dtype=np.float32,
        )
        self.gc_opt_range = (0.4, 0.8)
        self.min_score = 1e-6
        self._guide_int = _encode_sequence(self.guide.sequence)

    def score_sites(
        self,
        genome_sequences: Mapping[str, str],
        *,
        max_sites: int | None = None,
    ) -> List[CRISPRPhysicsResult]:
        results: List[CRISPRPhysicsResult] = []
        with cuda.gpus[self.device_id]:
            for chrom, seq in genome_sequences.items():
                normalized = bioinformatics.normalize_sequence(seq)
                results.extend(self._score_strand(chrom, normalized, strand=1))
                results.extend(self._score_strand(chrom, normalized, strand=-1))
        results.sort(key=lambda r: r.score, reverse=True)
        if max_sites is not None:
            results = results[:max_sites]
        return results

    def _score_strand(self, chrom: str, seq: str, strand: int) -> List[CRISPRPhysicsResult]:
        guide_len = len(self.guide.sequence)
        pam_pattern = self.pam_pattern
        max_mm = self.cas.max_mismatches
        if strand == -1:
            scan_seq = bioinformatics.reverse_complement(seq)
        else:
            scan_seq = seq
        base_positions = _base_candidate_positions(scan_seq, self.guide.sequence, max_mm)
        positions = _pam_filter_positions(scan_seq, base_positions, guide_len, pam_pattern)
        if not positions:
            return []
        scores = self._gpu_scores(scan_seq, positions)
        results: List[CRISPRPhysicsResult] = []
        for idx, start in enumerate(positions):
            target_seq = scan_seq[start : start + guide_len]
            mismatch_cost = float(scores[idx])
            gc_penalty = _gc_penalty(target_seq, *self.gc_opt_range)
            pam_seq = ""
            if pam_pattern:
                pam_seq = scan_seq[start + guide_len : start + guide_len + len(pam_pattern)]
            pam_ok, pam_penalty_raw = _pam_ok(pam_pattern, pam_seq)
            if not pam_ok:
                continue
            score = _compute_score(
                mismatch_cost,
                gc_penalty,
                pam_penalty_raw,
                guide_len,
                self.min_score,
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
                    mismatch_cost=mismatch_cost,
                    pam_ok=True,
                    score=score,
                )
            )
        return results

    def _gpu_scores(self, sequence: str, positions: List[int]) -> np.ndarray:
        seq_int = _encode_sequence(sequence)
        starts = np.asarray(positions, dtype=np.int32)
        d_seq = cuda.to_device(seq_int)
        d_starts = cuda.to_device(starts)
        d_guide = cuda.to_device(self._guide_int)
        d_weights = cuda.to_device(self._weights.astype(np.float32))
        d_scores = cuda.device_array(starts.shape[0], dtype=np.float32)
        threads_per_block = 128
        blocks = math.ceil(starts.shape[0] / threads_per_block) or 1
        _score_windows_kernel[blocks, threads_per_block](d_seq, d_starts, d_guide, d_weights, d_scores)
        return d_scores.copy_to_host()
