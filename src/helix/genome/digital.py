"""Digital genome primitives for Helix."""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List

from helix.edit.events import EditEvent

from .diff import apply_diffs


@dataclass(frozen=True)
class DigitalGenome:
    """
    Immutable base genome representation.

    Stores chromosome -> sequence mappings (uppercase DNA strings).
    """

    sequences: Dict[str, str]

    @classmethod
    def from_fasta(cls, path: str | Path) -> "DigitalGenome":
        """
        Load a FASTA (single or multi) into a DigitalGenome.
        """

        fasta_path = Path(path)
        if not fasta_path.exists():
            raise FileNotFoundError(f"Genome FASTA '{fasta_path}' not found.")
        sequences: Dict[str, str] = {}
        current_name: str | None = None
        buffer: List[str] = []
        with fasta_path.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_name is not None:
                        sequences[current_name] = "".join(buffer).upper()
                    header = line[1:].strip()
                    current_name = header.split()[0] if header else f"sequence_{len(sequences) + 1}"
                    buffer = []
                else:
                    buffer.append(line)
            if current_name is not None:
                sequences[current_name] = "".join(buffer).upper()
        if not sequences:
            raise ValueError(f"No sequences were found in {fasta_path}.")
        return cls(sequences=sequences)

    def slice_region(self, region: str) -> "DigitalGenome":
        """
        Slice a region of the form 'chrom:start-end' (1-based, inclusive).
        Returns a new DigitalGenome with a single chromosome entry.
        """

        region = region.strip()
        match = re.match(r"([^:]+):\s*(\d+)\s*-\s*(\d+)", region)
        if not match:
            raise ValueError(f"Invalid region specifier '{region}'. Expected 'chrom:start-end'.")
        chrom, start_str, end_str = match.groups()
        if chrom not in self.sequences:
            raise KeyError(f"Chromosome '{chrom}' not found in genome.")
        start_1 = int(start_str)
        end_1 = int(end_str)
        if start_1 <= 0 or end_1 < start_1:
            raise ValueError(f"Invalid region bounds in '{region}'.")
        sequence = self.sequences[chrom]
        start0 = max(0, start_1 - 1)
        end0 = min(len(sequence), end_1)
        sliced = sequence[start0:end0]
        if not sliced:
            raise ValueError(f"Region '{region}' produced an empty slice.")
        return DigitalGenome(sequences={chrom: sliced})

    def subset(self, chromosomes: Iterable[str]) -> "DigitalGenome":
        """Return a new DigitalGenome containing only the selected chromosome names."""

        selected: Dict[str, str] = {}
        for chrom in chromosomes:
            if chrom not in self.sequences:
                raise KeyError(f"Chromosome '{chrom}' not found in genome.")
            selected[chrom] = self.sequences[chrom]
        if not selected:
            raise ValueError("Subset operation produced an empty genome.")
        return DigitalGenome(sequences=selected)

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
