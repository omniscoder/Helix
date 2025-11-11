"""Lightweight RNA translation helpers for Helix."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterator, List, Tuple


CODON_TABLE: Dict[str, str] = {
    # Phenylalanine / Leucine
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    # Leucine
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    # Isoleucine / Methionine
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    # Valine
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    # Serine
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "AGU": "S",
    "AGC": "S",
    # Proline
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    # Threonine
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    # Alanine
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    # Tyrosine / Stop
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    # Histidine / Glutamine
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    # Asparagine / Lysine
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    # Aspartate / Glutamate
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    # Cysteine / Stop / Tryptophan
    "UGU": "C",
    "UGC": "C",
    "UGA": "*",
    "UGG": "W",
    # Arginine
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGA": "R",
    "AGG": "R",
    # Glycine
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def _clean_rna(rna: str) -> str:
    """Normalize to uppercase RNA alphabet (replace thymine with uracil)."""
    cleaned: list[str] = []
    for base in rna.upper():
        if base in {"\n", "\r", "\t", " "}:
            continue
        if base == "-":
            continue
        if base == "T":
            base = "U"
        if base not in {"A", "U", "G", "C"}:
            raise ValueError(f"Unexpected base '{base}' in RNA sequence.")
        cleaned.append(base)
    return "".join(cleaned)


def translate_rna(
    rna: str,
    *,
    stop_symbol: str = "*",
    stop_at_stop: bool = True,
) -> str:
    """Translate an RNA string into amino acids.

    Args:
        rna: String containing RNA or DNA characters.
        stop_symbol: Character to emit when a stop codon is encountered.
        stop_at_stop: If True, translation stops before appending a stop symbol.

    Returns:
        Protein sequence as a string of one-letter amino-acid codes.
    """
    cleaned = _clean_rna(rna)
    if len(cleaned) % 3 != 0:
        raise ValueError("RNA length must be divisible by 3 for translation.")

    amino_acids: list[str] = []
    for i in range(0, len(cleaned), 3):
        codon = cleaned[i : i + 3]
        aa = CODON_TABLE.get(codon)
        if aa is None:
            raise KeyError(f"Unknown codon '{codon}'.")

        if aa == "*":
            if stop_at_stop:
                break
            amino_acids.append(stop_symbol)
        else:
            amino_acids.append(aa)

    return "".join(amino_acids)


def translate_sequence_file(path: str) -> str:
    """Utility to translate a file containing a single RNA/DNA sequence."""
    with open(path, encoding="utf-8") as handle:
        contents = "".join(handle.read().split())
    return translate_rna(contents)


START_CODON = "AUG"
STOP_CODONS = {"UAA", "UAG", "UGA"}
BASE_COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G"}


@dataclass
class Orf:
    start: int
    end: int
    frame: int
    strand: str
    peptide: str

    def length_nt(self) -> int:
        return self.end - self.start

    def length_aa(self) -> int:
        return len(self.peptide)


@dataclass
class Frameshift:
    start: int
    end: int
    strand: str
    frames: Tuple[int, int]
    shift: int
    gap: int
    peptides: Tuple[str, str]

    def combined_peptide(self) -> str:
        return "".join(self.peptides)


def _iter_codons(seq: str, frame: int = 0) -> Iterator[Tuple[int, str]]:
    for i in range(frame, len(seq) - 2, 3):
        yield i, seq[i : i + 3]


def find_orfs(sequence: str, *, min_length: int = 30, include_partial: bool = False) -> List[Orf]:
    """Find open reading frames (ORFs) across forward and reverse frames.

    Args:
        sequence: DNA/RNA sequence (case-insensitive).
        min_length: Minimum nucleotide length for an ORF (default 30 nt).
        include_partial: If True, include ORFs that lack a terminating stop codon.
    """
    cleaned = _clean_rna(sequence)
    orfs: List[Orf] = []
    for frame in range(3):
        start_index: int | None = None
        for i, codon in _iter_codons(cleaned, frame):
            if start_index is None:
                if codon == START_CODON:
                    start_index = i
                continue

            if codon in STOP_CODONS:
                end = i + 3
                if end - start_index >= min_length:
                    peptide = translate_rna(cleaned[start_index:end])
                    orfs.append(Orf(start=start_index, end=end, frame=frame, strand="+", peptide=peptide))
                start_index = None

        if include_partial and start_index is not None:
            segment = cleaned[start_index:]
            usable_len = len(segment) - (len(segment) % 3)
            if usable_len >= 3:
                peptide = translate_rna(segment[:usable_len], stop_at_stop=False)
                if len(peptide) * 3 >= min_length:
                    orfs.append(
                        Orf(start=start_index, end=start_index + usable_len, frame=frame, strand="+", peptide=peptide)
                    )

    rc = "".join(BASE_COMPLEMENT[base] for base in reversed(cleaned))
    offset = len(cleaned)
    for frame in range(3):
        start_index: int | None = None
        for i, codon in _iter_codons(rc, frame):
            if start_index is None:
                if codon == START_CODON:
                    start_index = i
                continue

            if codon in STOP_CODONS:
                end = i + 3
                if end - start_index >= min_length:
                    peptide = translate_rna(rc[start_index:end])
                    start = offset - end
                    stop = offset - start_index
                    orfs.append(Orf(start=start, end=stop, frame=frame, strand="-", peptide=peptide))
                start_index = None

        if include_partial and start_index is not None:
            segment = rc[start_index:]
            usable_len = len(segment) - (len(segment) % 3)
            if usable_len >= 3:
                peptide = translate_rna(segment[:usable_len], stop_at_stop=False)
                if len(peptide) * 3 >= min_length:
                    start = offset - (start_index + usable_len)
                    stop = offset - start_index
                    orfs.append(Orf(start=start, end=stop, frame=frame, strand="-", peptide=peptide))

    return sorted(orfs, key=lambda orf: (orf.start, orf.end, orf.strand, orf.frame))


def orfs_to_fasta(orfs: List[Orf], *, line_width: int = 60) -> str:
    def wrap(seq: str) -> str:
        return "\n".join(seq[i:i + line_width] for i in range(0, len(seq), line_width))

    chunks: List[str] = []
    for idx, orf in enumerate(orfs, start=1):
        header = f">orf_{idx}_start{orf.start}_end{orf.end}_strand{orf.strand}_frame{orf.frame}"
        chunks.append(f"{header}\n{wrap(orf.peptide)}")
    return "\n".join(chunks)


def orfs_to_csv(orfs: List[Orf]) -> List[dict]:
    rows = []
    for idx, orf in enumerate(orfs, start=1):
        rows.append(
            {
                "id": f"orf_{idx}",
                "start": orf.start,
                "end": orf.end,
                "strand": orf.strand,
                "frame": orf.frame,
                "length_nt": orf.length_nt(),
                "length_aa": orf.length_aa(),
                "peptide": orf.peptide,
            }
        )
    return rows


def detect_frameshifts(
    sequence: str,
    *,
    min_orf_length: int = 30,
    gap_tolerance: int = 3,
    shift_sizes: Tuple[int, ...] = (1, 2),
) -> List[Frameshift]:
    """Detect candidate frameshifts by chaining adjacent ORFs on the same strand.

    Args:
        sequence: DNA/RNA sequence (case-insensitive).
        min_orf_length: Minimum nucleotide length for ORFs considered in chaining.
        gap_tolerance: Maximum nucleotide gap allowed between ORFs (default 3).
        shift_sizes: Allowed frame offsets (mod 3) between chained ORFs.
    """
    if any(size % 3 == 0 for size in shift_sizes):
        raise ValueError("shift_sizes must not include multiples of 3.")

    orfs = find_orfs(sequence, min_length=min_orf_length, include_partial=False)
    events: List[Frameshift] = []

    for i, current in enumerate(orfs):
        for nxt in orfs[i + 1 :]:
            if nxt.strand != current.strand:
                continue
            if nxt.start < current.end:
                continue
            gap = nxt.start - current.end
            if gap > gap_tolerance:
                break
            shift = (nxt.frame - current.frame) % 3
            if shift == 0 or shift not in shift_sizes:
                continue

            events.append(
                Frameshift(
                    start=current.start,
                    end=nxt.end,
                    strand=current.strand,
                    frames=(current.frame, nxt.frame),
                    shift=shift,
                    gap=gap,
                    peptides=(current.peptide, nxt.peptide),
                )
            )
    return events


def frameshifts_to_csv(frameshifts: List[Frameshift]) -> List[dict]:
    rows = []
    for idx, fs in enumerate(frameshifts, start=1):
        rows.append(
            {
                "id": f"frameshift_{idx}",
                "start": fs.start,
                "end": fs.end,
                "strand": fs.strand,
                "frames": f"{fs.frames[0]}/{fs.frames[1]}",
                "shift": fs.shift,
                "gap": fs.gap,
                "peptides": "/".join(fs.peptides),
            }
        )
    return rows

    return orfs
