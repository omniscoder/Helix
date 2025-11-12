"""
Prime editing simulation models for Helix.

These types support in-silico modeling of edit outcomes only.
They are not instructions for constructing or using any real-world system.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional

from helix.crispr.model import CasSystem, GuideRNA, TargetSite


@dataclass
class PegRNA:
    """
    Digital representation of a prime editing guide RNA for simulation.

    All sequence fields are abstract strings used for modeling.
    """

    spacer: str
    pbs: str  # primer binding site
    rtt: str  # reverse transcription template
    name: Optional[str] = None
    metadata: Dict[str, str] = field(default_factory=dict)


@dataclass
class PrimeEditor:
    """
    Parameterization of a prime editing system for simulation.

    This wraps a CasSystem plus additional numerical knobs that
    influence simulated outcome distributions.
    """

    name: str
    cas: CasSystem
    nick_to_edit_offset: int = 0
    """
    Abstract offset between nick and expected edit window in the
    digital sequence model.
    """

    # Generic parameters for simulation; interpretations are up to the model.
    efficiency_scale: float = 1.0
    indel_bias: float = 0.0
    mismatch_tolerance: int = 3


@dataclass
class PrimeEditOutcome:
    """
    A possible edit outcome in the digital genome model.

    Each outcome is defined by:
      - a target site
      - a description of the edited sequence
      - a probability weight (before normalization)
    """

    site: TargetSite
    edited_sequence: str
    logit_score: float  # unnormalized log-probability or score
    description: str = ""
