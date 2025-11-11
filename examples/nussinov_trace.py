"""Visualize the Helix Nussinov folding result for a sequence."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

from nussinov_algorithm import nussinov

DEFAULT_SEQUENCE = "GGGAAACCC"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute Nussinov dot-bracket notation for a sequence.")
    parser.add_argument(
        "sequence",
        nargs="?",
        help="Inline RNA/DNA string. Defaults to a short demo hairpin if omitted.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        help="Optional path to a FASTA/plaintext file containing a single sequence.",
    )
    parser.add_argument(
        "--min-loop",
        type=int,
        default=3,
        help="Minimum unpaired loop length (default: 3).",
    )
    parser.add_argument(
        "--no-wobble",
        action="store_true",
        help="Disable G-U wobble pairs.",
    )
    parser.add_argument(
        "--dot-output",
        type=Path,
        help="Optional path to write just the dot-bracket string.",
    )
    return parser.parse_args()


def load_sequence(sequence: Optional[str], path: Optional[Path]) -> str:
    if sequence and path:
        raise ValueError("Provide either a positional sequence or --input, not both.")
    if path:
        return path.read_text(encoding="utf-8")
    return sequence or DEFAULT_SEQUENCE


def main() -> None:
    args = parse_args()
    raw = load_sequence(args.sequence, args.input)
    result = nussinov(
        raw,
        min_loop_length=args.min_loop,
        allow_wobble_pairs=not args.no_wobble,
    )
    print(f"Sequence ({len(result.sequence)} nt): {result.sequence}")
    print(f"Dot-bracket ({result.score()} pairs):\n{result.structure}")
    if args.dot_output:
        args.dot_output.write_text(result.structure, encoding="utf-8")
        print(f"Dot-bracket saved to {args.dot_output}")
    if result.pairs:
        print("\nBase pairs (0-indexed):")
        for i, j in result.pairs:
            print(f"{i:>3} - {j:<3} ({result.sequence[i]}-{result.sequence[j]})")


if __name__ == "__main__":
    main()
