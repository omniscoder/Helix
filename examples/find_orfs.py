"""Scan for open reading frames using Helix codon helpers."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Optional

import bioinformatics
from codon import (
    detect_frameshifts,
    frameshifts_to_csv,
    find_orfs,
    orfs_to_csv,
    orfs_to_fasta,
    Orf,
    Frameshift,
)


def load_sequence(sequence: Optional[str], path: Optional[Path]) -> str:
    if sequence and path:
        raise ValueError("Provide either a sequence or --input file, not both.")
    if path:
        contents = path.read_text(encoding="utf-8")
        return bioinformatics.normalize_sequence(contents)
    if sequence:
        return bioinformatics.normalize_sequence(sequence)
    raise ValueError("No sequence provided. Use positional argument or --input.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Identify ORFs (forward and reverse) in a sequence.")
    parser.add_argument(
        "sequence",
        nargs="?",
        help="Inline DNA/RNA string. If omitted, supply --input.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        help="Path to a text/FASTA file containing a single sequence.",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=30,
        help="Minimum ORF length in nucleotides (default: 30).",
    )
    parser.add_argument(
        "--include-partial",
        action="store_true",
        help="Include ORFs that run off the sequence end without a stop codon.",
    )
    parser.add_argument(
        "--detect-frameshifts",
        action="store_true",
        help="Report candidate frameshifts by chaining adjacent ORFs.",
    )
    parser.add_argument(
        "--gap-tolerance",
        type=int,
        default=3,
        help="Maximum gap between ORFs when detecting frameshifts (default: 3).",
    )
    parser.add_argument(
        "--orf-fasta",
        type=Path,
        help="Optional path to write ORF peptides as FASTA.",
    )
    parser.add_argument(
        "--orf-csv",
        type=Path,
        help="Optional path to write ORF metadata as CSV.",
    )
    parser.add_argument(
        "--frameshift-csv",
        type=Path,
        help="Optional path to write frameshift candidates as CSV.",
    )
    return parser.parse_args()


def format_orf(start: int, end: int, frame: int, strand: str, peptide: str) -> str:
    length_nt = end - start
    length_aa = len(peptide)
    return (
        f"start={start:>4} end={end:>4} strand={strand} frame={frame} "
        f"length_nt={length_nt:>4} length_aa={length_aa:>3} peptide={peptide}"
    )


def main() -> None:
    args = parse_args()
    sequence = load_sequence(args.sequence, args.input)
    orfs = find_orfs(sequence, min_length=args.min_length, include_partial=args.include_partial)
    if not orfs:
        print("No ORFs found.")
    else:
        for orf in orfs:
            print(format_orf(orf.start, orf.end, orf.frame, orf.strand, orf.peptide))
        if args.orf_fasta:
            fasta_text = orfs_to_fasta(orfs)
            args.orf_fasta.write_text(fasta_text, encoding="utf-8")
            print(f"Wrote ORF peptides to {args.orf_fasta}")
        if args.orf_csv:
            rows = orfs_to_csv(orfs)
            with args.orf_csv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            print(f"Wrote ORF table to {args.orf_csv}")

    if args.detect_frameshifts:
        events = detect_frameshifts(
            sequence,
            min_orf_length=args.min_length,
            gap_tolerance=args.gap_tolerance,
        )
        if not events:
            print("No candidate frameshifts.")
        else:
            print("\nFrameshift candidates:")
            for event in events:
                peptides = "/".join(event.peptides)
                print(
                    f"start={event.start:>4} end={event.end:>4} strand={event.strand} "
                    f"frames={event.frames} shift={event.shift} gap={event.gap} peptides={peptides}"
                )
            if args.frameshift_csv:
                rows = frameshifts_to_csv(events)
                with args.frameshift_csv.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
                    writer.writeheader()
                    writer.writerows(rows)
                print(f"Wrote frameshift table to {args.frameshift_csv}")


if __name__ == "__main__":
    main()
