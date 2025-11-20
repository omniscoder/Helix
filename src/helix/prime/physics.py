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
from typing import Dict, Optional, Tuple

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
