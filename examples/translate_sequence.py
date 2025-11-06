"""Translate a DNA/RNA sequence using the Helix codon helper."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

from codon import translate_rna


def load_sequence(sequence: Optional[str], path: Optional[Path]) -> str:
    if sequence and path:
        raise ValueError("Provide either a sequence or --input file, not both.")
    if path:
        contents = path.read_text(encoding="utf-8")
        return "".join(contents.split())
    if sequence:
        return sequence
    raise ValueError("No sequence provided. Use positional argument or --input.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Translate a DNA/RNA sequence.")
    parser.add_argument(
        "sequence",
        nargs="?",
        help="Inline sequence to translate. If omitted, supply --input.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        help="Path to a text/FASTA file containing a single sequence.",
    )
    parser.add_argument(
        "--no-stop",
        action="store_true",
        help="Keep translating after stop codons (default breaks at the first stop).",
    )
    parser.add_argument(
        "--stop-symbol",
        default="*",
        help="Character used when emitting stop codons (default: '*').",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    raw_sequence = load_sequence(args.sequence, args.input)
    stop_at_stop = not args.no_stop
    protein = translate_rna(raw_sequence, stop_symbol=args.stop_symbol, stop_at_stop=stop_at_stop)
    print(protein)


if __name__ == "__main__":
    main()
