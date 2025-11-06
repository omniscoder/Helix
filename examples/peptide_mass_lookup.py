"""Compute simple peptide masses using Helix's amino acid table.

Run with:
    python examples/peptide_mass_lookup.py PEPTIDE
"""
import sys
from typing import Iterable

from amino_acids import mass_index_table


def peptide_mass(sequence: str) -> int:
    """Return the total mass for a peptide sequence."""
    try:
        masses: Iterable[int] = (mass_index_table[aa] for aa in sequence.upper())
        return sum(masses)
    except KeyError as exc:
        unknown = exc.args[0]
        raise ValueError(f"Unknown residue '{unknown}'") from exc


def main(argv: list[str]) -> None:
    if len(argv) < 2:
        print("Usage: python examples/peptide_mass_lookup.py PEPTIDE")
        print("Example: python examples/peptide_mass_lookup.py GAVLIM")
        sys.exit(1)

    peptide = argv[1]
    mass = peptide_mass(peptide)
    print(f"Peptide: {peptide}")
    print(f"Total mass: {mass} Da")


if __name__ == "__main__":
    main(sys.argv)
