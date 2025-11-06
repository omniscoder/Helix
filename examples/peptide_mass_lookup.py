"""Compute simple peptide masses using Helix's amino acid table."""
from __future__ import annotations

import argparse
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Add up monoisotopic masses for a peptide.")
    parser.add_argument(
        "peptides",
        nargs="+",
        help="One or more peptide sequences (e.g. GAVLIM).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    for peptide in args.peptides:
        try:
            mass = peptide_mass(peptide)
        except ValueError as err:
            print(f"{peptide}\tERROR: {err}")
            continue
        print(f"{peptide}\t{mass} Da")


if __name__ == "__main__":
    main()
