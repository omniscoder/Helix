"""Digital genome primitives for Helix."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

from helix.edit.events import EditEvent

from .diff import apply_diffs


@dataclass(frozen=True)
class DigitalGenome:
    """
    Immutable base genome representation.

    Stores chromosome -> sequence mappings (uppercase DNA strings).
    """

    sequences: Dict[str, str]

    def view(self) -> "DigitalGenomeView":
        """Return a zero-diff view rooted at this genome."""
        return DigitalGenomeView(base=self, diffs=[])


@dataclass(frozen=True)
class DigitalGenomeView:
    """
    Immutable view over a DigitalGenome, defined by a sequence of EditEvents.
    """

    base: DigitalGenome
    diffs: List[EditEvent] = field(default_factory=list)

    def apply(self, event: EditEvent) -> "DigitalGenomeView":
        """Return a new view with the edit appended."""
        return DigitalGenomeView(base=self.base, diffs=[*self.diffs, event])

    def materialize_chrom(self, chrom: str) -> str:
        """Apply relevant diffs to a chromosome and return the resulting sequence."""
        seq = self.base.sequences[chrom]
        relevant = [event for event in self.diffs if event.chrom == chrom]
        if not relevant:
            return seq
        return apply_diffs(seq, relevant)

    def materialize_all(self) -> Dict[str, str]:
        """Return a dict of all chromosome sequences for this view."""
        return {chrom: self.materialize_chrom(chrom) for chrom in self.base.sequences}
