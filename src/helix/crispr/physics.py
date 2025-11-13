"""
CRISPR physics interfaces and CPU reference implementation.

These helpers remain purely computational: they operate on digital sequences,
compute mismatch/PAM penalties, and feed candidate sites back to higher-level
simulators. Swapping the backend (CPU vs. GPU) only requires changing the
physics factory.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Sequence

from .. import bioinformatics
from .kmm import k_mismatch_positions
from .model import CasSystem, GuideRNA, TargetSite

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


@dataclass
class CRISPRPhysicsResult:
    site: TargetSite
    mismatch_cost: float
    pam_ok: bool
    score: float


class CRISPRPhysicsBase(ABC):
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
        self._weights = tuple(_seed_weight_vector(len(self.guide.sequence), seed_len))
        self.seed_length = seed_len
        self.bulge_penalty = max(0.5, cas.weight_mismatch_penalty * 2.0)
        self.pam_weak_penalty = max(0.5, cas.weight_pam_penalty or 1.0)
        self.pam_fail_penalty = max(self.pam_weak_penalty * 2.0, 2.0)
        self.gc_opt_range = (0.4, 0.8)
        self.min_score = 1e-6

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
            mismatch_cost = self._mismatch_cost(target_seq)
            gc_penalty = _gc_penalty(target_seq, *self.gc_opt_range)
            score = _compute_score(
                mismatch_cost,
                gc_penalty,
                pam_penalty_raw * self.pam_weak_penalty,
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

    def _mismatch_cost(self, target_seq: str) -> float:
        cost = 0.0
        for idx, base in enumerate(target_seq):
            if idx >= len(self.guide.sequence):
                cost += self.bulge_penalty
                continue
            if base != self.guide.sequence[idx]:
                weight = self._weights[idx] if idx < len(self._weights) else 1.0
                cost += weight
        extra = abs(len(target_seq) - len(self.guide.sequence))
        if extra:
            cost += extra * self.bulge_penalty
        return cost


def create_crispr_physics(
    cas: CasSystem,
    guide: GuideRNA,
    *,
    use_gpu: bool = False,
) -> CRISPRPhysicsBase:
    if use_gpu:
        try:
            from .physics_gpu import CRISPRPhysicsGPU  # pragma: no cover - optional
        except Exception as exc:  # pragma: no cover - GPU optional
            raise RuntimeError("GPU backend requested, but CUDA/Numba is unavailable.") from exc
        return CRISPRPhysicsGPU(cas, guide)
    return CRISPRPhysicsCPU(cas, guide)
