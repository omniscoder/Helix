"""Convert examples/prime_demo_genome.json into a FASTA file for CLI demos."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Write a FASTA from a JSON genome config.")
    parser.add_argument("--input", type=Path, default=Path("examples/prime_demo_genome.json"))
    parser.add_argument("--out", type=Path, required=True, help="Output FASTA path")
    args = parser.parse_args()

    data = json.loads(args.input.read_text())
    with args.out.open("w", encoding="utf-8") as handle:
        for chrom in data.get("chromosomes", []):
            name = chrom.get("name", "chr")
            seq = chrom.get("sequence", "")
            handle.write(f">{name}\n{seq}\n")
    print(f"Wrote FASTA to {args.out}")


if __name__ == "__main__":
    main()
