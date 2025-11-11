"""Lightweight IO helpers (FASTA, etc.)."""
from __future__ import annotations

from pathlib import Path
from typing import List, Tuple


def read_fasta(path: str | Path) -> List[Tuple[str, str]]:
    """Return a list of (header, sequence) tuples from a FASTA file."""
    header = None
    chunks: List[str] = []
    records: List[Tuple[str, str]] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            records.append((header, "".join(chunks)))
    return records


def read_plain_sequence(path: str | Path) -> str:
    """Return a contiguous uppercase DNA/RNA string from a FASTA or raw text file."""
    records = read_fasta(path)
    if records:
        return "".join(seq for _, seq in records).upper()
    return Path(path).read_text(encoding="utf-8").strip().upper()
