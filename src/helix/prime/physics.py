"""
Prime editing physics helpers.

This module keeps the pegRNA binding + RTT extension heuristics separate from
the simulators so we can iterate on physics-inspired behaviour in one place.
All numbers here are purely computational knobs – they are not calibrated to
any specific dataset yet.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import zip_longest
from typing import Dict, List, Optional, Sequence, Tuple

from .. import bioinformatics
from .model import PegRNA, PrimeEditor

PRIME_SCORING_VERSION = "1.0.0"


@dataclass(frozen=True)
class PrimeBranchDistribution:
    """Compact container for branch logits coming out of the prime physics model."""

    efficiency: float
    left: float
    right: float
    reanneal: float
    indel: float
    no_edit: float

    def as_dict(self) -> Dict[str, float]:
        return {
            "efficiency": self.efficiency,
            "left": self.left,
            "right": self.right,
            "reanneal": self.reanneal,
            "indel": self.indel,
            "no_edit": self.no_edit,
        }


@dataclass
class PrimePhysics:
    """
    Parameterized heuristics for prime editing simulations.

    The class exposes interpretable knobs:
      - anchor_strength: spacer/PBS binding quality
      - extension_strength: RTT copying quality
      - flap distribution: probability mass over left/right/reanneal/indel/no-edit
    """

    editor: PrimeEditor
    peg: PegRNA
    preferred_pbs_range: Tuple[int, int] = (10, 16)
    preferred_rtt_range: Tuple[int, int] = (10, 34)
    anchor_mismatch_penalty: float = 0.2
    anchor_length_penalty: float = 0.02
    extension_length_penalty: float = 0.01
    gc_penalty_scale: float = 0.5
    min_score: float = 1e-6

    _pbs_len: int = 0
    _rtt_len: int = 0
    _spacer: str = ""

    def __post_init__(self) -> None:
        self._spacer = bioinformatics.normalize_sequence(self.peg.spacer)
        self._pbs_len = len(bioinformatics.normalize_sequence(self.peg.pbs))
        self._rtt_len = len(bioinformatics.normalize_sequence(self.peg.rtt))

    @classmethod
    def from_config(cls, editor: PrimeEditor, peg: PegRNA) -> "PrimePhysics":
        return cls(editor=editor, peg=peg)

    def _range_penalty(self, value: int, preferred: Tuple[int, int], scale: float) -> float:
        lo, hi = preferred
        if lo <= value <= hi:
            return 0.0
        delta = min(abs(value - lo), abs(value - hi))
        return delta * scale

    def anchor_strength(self, mismatch_count: int) -> float:
        penalty = mismatch_count * self.anchor_mismatch_penalty
        penalty += self._range_penalty(self._pbs_len, self.preferred_pbs_range, self.anchor_length_penalty)
        return max(self.min_score, 1.0 - penalty)

    def extension_strength(self, local_gc: Optional[float] = None) -> float:
        penalty = self._range_penalty(self._rtt_len, self.preferred_rtt_range, self.extension_length_penalty)
        if local_gc is not None:
            delta = max(0.0, abs(local_gc - 0.45) - 0.15)
            penalty += delta * self.gc_penalty_scale
        return max(self.min_score, 1.0 - penalty)

    def binding_efficiency(self, mismatch_count: int, *, local_gc: Optional[float] = None) -> float:
        anchor = self.anchor_strength(mismatch_count)
        extension = self.extension_strength(local_gc)
        eff = anchor * extension * max(self.editor.efficiency_scale, self.min_score)
        return max(self.min_score, min(1.0, eff))

    def log_efficiency(self, mismatch_count: int, *, local_gc: Optional[float] = None) -> float:
        eff = self.binding_efficiency(mismatch_count, local_gc=local_gc)
        return math.log(max(eff, self.min_score))

    def branch_distribution(
        self,
        mismatch_count: int,
        *,
        local_gc: Optional[float] = None,
    ) -> PrimeBranchDistribution:
        eff = self.binding_efficiency(mismatch_count, local_gc=local_gc)
        flap_bias = min(max(self.editor.flap_balance, 0.0), 1.0)
        left = eff * flap_bias
        right = eff - left
        reanneal = max(self.min_score, (1.0 - eff) * max(self.editor.reanneal_bias, 0.0))
        indel = max(self.min_score, self.editor.indel_bias + mismatch_count * 0.05)
        remaining = max(self.min_score, 1.0 - (left + right + reanneal))
        no_edit = max(self.min_score, remaining)
        return PrimeBranchDistribution(
            efficiency=eff,
            left=left,
            right=right,
            reanneal=reanneal,
            indel=indel,
            no_edit=no_edit,
        )

    def mismatch_count(self, sequence: str) -> int:
        """Helper for scoring routines to compare the spacer against an observed target sequence."""

        seq = bioinformatics.normalize_sequence(sequence)
        n = min(len(seq), len(self._spacer))
        mismatches = sum(1 for idx in range(n) if seq[idx] != self._spacer[idx])
        mismatches += abs(len(seq) - len(self._spacer))
        return mismatches


_NN_DG: Dict[Tuple[str, str], float] = {
    ("A", "T"): -1.0,
    ("T", "A"): -1.0,
    ("C", "G"): -1.5,
    ("G", "C"): -1.5,
    ("G", "T"): -0.5,
    ("T", "G"): -0.5,
    ("A", "A"): 0.5,
    ("C", "C"): 0.5,
    ("G", "G"): 0.5,
    ("T", "T"): 0.5,
}

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _complement(seq: str) -> str:
    return "".join(_COMPLEMENT.get(ch.upper(), ch.upper()) for ch in seq)


def _pair_delta_g(base_a: str, base_b: str) -> float:
    if not base_a or not base_b:
        return 0.0
    a = base_a.upper()
    b = base_b.upper()
    if b == _COMPLEMENT.get(a):
        return _NN_DG.get((a, b), -1.0)
    return _NN_DG.get((a, b), 0.75)


def compute_pbs_dG(pbs_seq: str, target_slice: str) -> float:
    seq = bioinformatics.normalize_sequence(pbs_seq)
    target = _complement(bioinformatics.normalize_sequence(target_slice))
    total = 0.0
    for a, b in zip(seq, target):
        total += _pair_delta_g(a, b)
    total += abs(len(seq) - len(target)) * 0.5
    return total


def compute_rt_path_dG(rt_template: str, target_slice: str) -> List[float]:
    seq = bioinformatics.normalize_sequence(rt_template)
    target = bioinformatics.normalize_sequence(target_slice)
    cumulative: List[float] = []
    total = 0.0
    for a, b in zip_longest(seq, target, fillvalue="N"):
        total += _pair_delta_g(a, b)
        cumulative.append(total)
    return cumulative


def compute_flap_ddG(edited_seq: str, wt_seq: str) -> float:
    edited = bioinformatics.normalize_sequence(edited_seq)
    wt = bioinformatics.normalize_sequence(wt_seq)
    matches = sum(1 for a, b in zip(edited, wt) if a == b)
    mismatches = sum(1 for a, b in zip(edited, wt) if a != b)
    mismatches += abs(len(edited) - len(wt))
    return mismatches * 0.5 - matches * 0.2


def compute_microhomology(edited: str, wt: str) -> int:
    edited_norm = bioinformatics.normalize_sequence(edited)
    wt_norm = bioinformatics.normalize_sequence(wt)
    limit = min(len(edited_norm), len(wt_norm))
    best = 0
    for size in range(limit, 0, -1):
        if edited_norm[-size:] == wt_norm[:size]:
            best = size
            break
    return best


def estimate_mmr_flag(mismatch_positions: Sequence[int], tolerance: int) -> bool:
    return len(mismatch_positions) > max(tolerance, 0)


@dataclass
class PrimePhysicsScore:
    pbs_dG: float
    rt_cum_dG: List[float]
    flap_ddG: float
    microhomology: int
    mmr_flag: bool
    nick_distance: int
    P_RT: float
    P_flap: float
    E_pred: float


def _sigmoid(value: float) -> float:
    if value >= 0:
        z = math.exp(-value)
        return 1.0 / (1.0 + z)
    z = math.exp(value)
    return z / (1.0 + z)


def _mismatch_positions(seq_a: str, seq_b: str) -> List[int]:
    a_norm = bioinformatics.normalize_sequence(seq_a)
    b_norm = bioinformatics.normalize_sequence(seq_b)
    return [idx for idx, (aa, bb) in enumerate(zip(a_norm, b_norm)) if aa != bb]


def score_prime_design(
    peg: PegRNA,
    target_sequence: str,
    editor: PrimeEditor,
) -> PrimePhysicsScore:
    target = bioinformatics.normalize_sequence(target_sequence)
    spacer = bioinformatics.normalize_sequence(peg.spacer or "")
    pbs_seq = bioinformatics.normalize_sequence(peg.pbs or "")
    rtt_seq = bioinformatics.normalize_sequence(peg.rtt or "")
    if not target:
        return PrimePhysicsScore(
            pbs_dG=0.0,
            rt_cum_dG=[],
            flap_ddG=0.0,
            microhomology=0,
            mmr_flag=False,
            nick_distance=0,
            P_RT=0.5,
            P_flap=0.5,
            E_pred=0.25,
        )
    spacer_idx = target.find(spacer) if spacer else -1
    if spacer_idx == -1:
        spacer_idx = len(target) // 2
    pbs_target_start = max(spacer_idx - len(pbs_seq), 0)
    pbs_target = target[pbs_target_start : pbs_target_start + len(pbs_seq)]
    rtt_target_start = spacer_idx
    rtt_target = target[rtt_target_start : rtt_target_start + len(rtt_seq)]
    pbs_dg = compute_pbs_dG(pbs_seq, pbs_target) if pbs_seq else 0.0
    rt_cum = compute_rt_path_dG(rtt_seq, rtt_target) if rtt_seq else []
    flap_ddg = compute_flap_ddG(rtt_seq, rtt_target) if rtt_seq else 0.0
    microhomology = compute_microhomology(rtt_seq, rtt_target) if rtt_seq else 0
    mismatch_positions = _mismatch_positions(rtt_seq, rtt_target)
    tolerance = editor.mismatch_tolerance if editor and editor.mismatch_tolerance is not None else 3
    mmr_flag = estimate_mmr_flag(mismatch_positions, tolerance)
    nick_distance = abs((editor.nick_to_edit_offset or 0) - len(rtt_seq))
    rt_total = rt_cum[-1] if rt_cum else 0.0
    p_rt = _sigmoid((-pbs_dg - rt_total) * 0.3)
    p_flap = _sigmoid((microhomology - abs(flap_ddg)) * 0.4)
    e_pred = p_rt * p_flap * max(editor.efficiency_scale, 1e-3)
    e_pred = max(0.0, min(1.0, e_pred))
    return PrimePhysicsScore(
        pbs_dG=pbs_dg,
        rt_cum_dG=rt_cum,
        flap_ddG=flap_ddg,
        microhomology=microhomology,
        mmr_flag=mmr_flag,
        nick_distance=nick_distance,
        P_RT=p_rt,
        P_flap=p_flap,
        E_pred=e_pred,
    )
