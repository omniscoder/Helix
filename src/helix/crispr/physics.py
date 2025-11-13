"""
Physics-inspired scoring helpers for CRISPR simulations.

These utilities stay entirely within the digital genome world; they do not
describe laboratory protocols. The goal is to separate the "how do we score
candidate sites?" logic from the simulators so that more detailed physics
models (PAM penalties, seed weighting, bulges, GC bias, etc.) can be swapped
in without rewriting the scanners.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence, Tuple

from .. import bioinformatics
from .model import CasSystem, GuideRNA

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
class CRISPRPhysics:
    """
    Tunable scoring surface for candidate CRISPR cut sites.

    Parameters
    ----------
    guide_sequence:
        Spacer sequence (5'→3'). Automatically normalized to uppercase DNA.
    pam_pattern:
        IUPAC string describing the PAM requirement. Empty string disables PAM
        checks (useful for near-PAM-less editors).
    seed_length:
        Number of PAM-proximal bases that receive heavier mismatch penalties.
    mismatch_weights:
        Optional per-position weights. If omitted, the class generates a seed-
        weighted profile automatically.
    bulge_penalty:
        Cost added per inserted/deleted base when the alignment length differs
        from the guide length.
    pam_weak_penalty / pam_fail_penalty:
        Penalties for single-symbol PAM mismatches vs outright failures.
    gc_opt_range:
        Acceptable GC range before GC penalties kick in (inclusive bounds).
    """

    guide_sequence: str
    pam_pattern: str = "NGG"
    seed_length: int = 10
    mismatch_weights: Sequence[float] | None = None
    bulge_penalty: float = 2.0
    pam_weak_penalty: float = 1.5
    pam_fail_penalty: float = 4.0
    gc_opt_range: Tuple[float, float] = (0.4, 0.8)
    min_score: float = 1e-6

    _weights: Tuple[float, ...] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self.guide_sequence = bioinformatics.normalize_sequence(self.guide_sequence)
        self.pam_pattern = (self.pam_pattern or "").upper()
        if not self.guide_sequence:
            raise ValueError("CRISPRPhysics requires a non-empty guide sequence.")
        if self.mismatch_weights:
            weights = tuple(float(w) for w in self.mismatch_weights)
            if len(weights) < len(self.guide_sequence):
                pad = (1.0,) * (len(self.guide_sequence) - len(weights))
                self._weights = weights + pad
            else:
                self._weights = weights[: len(self.guide_sequence)]
        else:
            self._weights = self._seed_weight_vector(len(self.guide_sequence))

    @classmethod
    def from_system(cls, cas: CasSystem, guide: GuideRNA) -> "CRISPRPhysics":
        """
        Convenience constructor that derives reasonable defaults from a Cas system.
        """

        pam_pattern = guide.pam or (cas.pam_rules[0].pattern if cas.pam_rules else "")
        seed_len = min(len(bioinformatics.normalize_sequence(guide.sequence)), 12)
        bulge_penalty = max(0.5, cas.weight_mismatch_penalty * 2.0)
        pam_weak = max(0.5, cas.weight_pam_penalty or 1.0)
        pam_fail = max(pam_weak * 2.0, 2.0)
        return cls(
            guide_sequence=guide.sequence,
            pam_pattern=pam_pattern,
            seed_length=seed_len,
            mismatch_weights=None,
            bulge_penalty=bulge_penalty,
            pam_weak_penalty=pam_weak,
            pam_fail_penalty=pam_fail,
        )

    def _seed_weight_vector(self, length: int) -> Tuple[float, ...]:
        weights = [1.0] * length
        seed_len = max(0, min(self.seed_length, length))
        seed_start = length - seed_len
        for idx in range(seed_start, length):
            weights[idx] = 1.5
        return tuple(weights)

    def _pam_penalty(self, pam_seq: str) -> tuple[bool, float]:
        if not self.pam_pattern:
            return True, 0.0
        if len(pam_seq) != len(self.pam_pattern):
            return False, self.pam_fail_penalty
        mismatches = 0
        for pattern_char, base in zip(self.pam_pattern, pam_seq.upper()):
            allowed = _IUPAC.get(pattern_char, pattern_char)
            if base not in allowed:
                mismatches += 1
                if mismatches > 1:
                    return False, self.pam_fail_penalty
        if mismatches == 0:
            return True, 0.0
        return True, self.pam_weak_penalty * mismatches

    def _mismatch_cost(self, target_seq: str) -> float:
        n = min(len(target_seq), len(self.guide_sequence))
        cost = 0.0
        for idx in range(n):
            if target_seq[idx] != self.guide_sequence[idx]:
                weight = self._weights[idx] if idx < len(self._weights) else 1.0
                cost += weight
        extra = abs(len(target_seq) - len(self.guide_sequence))
        if extra:
            cost += extra * self.bulge_penalty
        return cost

    def _gc_penalty(self, seq: str) -> float:
        if not seq:
            return 0.0
        gc = bioinformatics.gc_content(seq)
        lo, hi = self.gc_opt_range
        if lo <= gc <= hi:
            return 0.0
        delta = min(abs(gc - lo), abs(gc - hi))
        return delta * 2.0

    def score_site(self, target_seq: str, pam_seq: str, strand: int = 1) -> float:
        """
        Return a normalized score (0..1) for a guide/PAM pairing on the given strand.

        The strand argument is currently informational but kept for future
        orientation-specific effects.
        """

        del strand  # strand-specific parameters can be incorporated later.
        normalized_target = bioinformatics.normalize_sequence(target_seq)
        pam_ok, pam_cost = self._pam_penalty(pam_seq.upper())
        if not pam_ok:
            return 0.0
        mismatch_cost = self._mismatch_cost(normalized_target)
        gc_cost = self._gc_penalty(normalized_target)
        total_cost = mismatch_cost + gc_cost + pam_cost
        span = max(1.0, float(len(self.guide_sequence)))
        score = max(self.min_score, 1.0 - (total_cost / span))
        return min(1.0, score)
