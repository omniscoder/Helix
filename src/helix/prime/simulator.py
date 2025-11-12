"""
Prime editing sequence-level simulation for Helix.

All operations here work on digital sequences only and are not
wet-lab protocols.
"""
from __future__ import annotations

from typing import List, Optional

from helix.crispr.model import DigitalGenome, TargetSite

from .model import PegRNA, PrimeEditOutcome, PrimeEditor


def locate_prime_target_site(
    genome: DigitalGenome,
    peg: PegRNA,
    *,
    search_window: int = 200,
) -> Optional[TargetSite]:
    """
    Identify a primary digital target site for a pegRNA.

    This function is responsible for aligning the spacer sequence to
    the digital genome and returning a site suitable for prime editing
    simulation.

    Returns None if no plausible site is found.
    """

    # TODO: implement alignment/search against DigitalGenome.
    raise NotImplementedError("locate_prime_target_site is not yet implemented.")


def simulate_prime_edit(
    genome: DigitalGenome,
    editor: PrimeEditor,
    peg: PegRNA,
    *,
    max_outcomes: int = 16,
) -> List[PrimeEditOutcome]:
    """
    Simulate prime editing outcomes in a digital genome.

    Conceptual responsibilities:
      - resolve the target site (or sites)
      - apply an abstract edit model using peg.rtt and editor parameters
      - generate a set of possible edited sequences
      - assign scores / logits that downstream code can normalize

    This function does not alter the input genome; it only returns
    hypothetical outcomes for further analysis.
    """

    # TODO: call locate_prime_target_site and construct outcome set.
    raise NotImplementedError("simulate_prime_edit is not yet implemented.")
