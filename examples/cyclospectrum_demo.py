"""Play with peptide spectra and leaderboard sequencing."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List

from helix.cyclospectrum import (
    cyclic_spectrum,
    leaderboard_cyclopeptide_sequencing,
    linear_spectrum,
    score_peptide,
    theoretical_spectrum,
)


def _parse_masses(text: str | None) -> List[int]:
    if not text:
        return []
    tokens = text.replace(",", " ").split()
    if not tokens:
        return []
    try:
        return [int(token) for token in tokens]
    except ValueError as exc:
        raise SystemExit(f"Invalid spectrum token: {exc}") from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Helix cyclo-spectrum helpers from the CLI.")
    parser.add_argument("--peptide", help="Optional peptide sequence to analyse.")
    parser.add_argument(
        "--linear",
        action="store_true",
        help="Use the linear spectrum when printing --peptide (default: cyclic).",
    )
    parser.add_argument(
        "--spectrum",
        help="Comma/space-separated list of experimental masses (e.g. '0,113,128').",
    )
    parser.add_argument(
        "--spectrum-file",
        type=Path,
        help="Path to a text file containing whitespace-separated masses.",
    )
    parser.add_argument(
        "--leaderboard",
        type=int,
        default=5,
        help="Size of the leaderboard when sequencing (default: 5).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    spectrum = _parse_masses(args.spectrum)
    if args.spectrum_file:
        file_masses = _parse_masses(args.spectrum_file.read_text(encoding="utf-8"))
        if spectrum:
            spectrum.extend(file_masses)
        else:
            spectrum = file_masses

    if args.peptide:
        spec = theoretical_spectrum(args.peptide, cyclic=not args.linear)
        mode = "cyclic" if not args.linear else "linear"
        print(f"{mode.title()} spectrum for {args.peptide}:")
        print(" ".join(str(m) for m in spec))
        if spectrum:
            print(f"Score vs provided spectrum: {score_peptide(args.peptide, spectrum, cyclic=not args.linear)}")

    if spectrum:
        print(f"\nRunning leaderboard sequencing (top {args.leaderboard})...")
        hits = leaderboard_cyclopeptide_sequencing(spectrum, leaderboard_size=args.leaderboard)
        if not hits:
            print("No candidates found. Try increasing --leaderboard or adjusting the spectrum.")
        else:
            for peptide, score in hits:
                print(f"{peptide}\tscore={score}")
    elif not args.peptide:
        raise SystemExit("Provide at least --peptide or --spectrum/--spectrum-file.")


if __name__ == "__main__":
    main()
