"""Practical protein sequence helpers for Helix.

Hydropathy scales are resolved dynamically; unavailable tables are logged and skipped.
"""
from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

try:  # pragma: no cover - import guard
    from Bio import SeqIO
    from Bio.SeqUtils import ProtParamData
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    BIOPYTHON_AVAILABLE = True
except ImportError:  # pragma: no cover - fallback when Biopython missing
    SeqIO = None
    ProtParamData = None
    ProteinAnalysis = None
    BIOPYTHON_AVAILABLE = False


VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
LOGGER = logging.getLogger(__name__)
if BIOPYTHON_AVAILABLE:
    _scales: Dict[str, Dict[str, float]] = {}
    for key in ("kd", "hs", "flexibility"):
        attr = getattr(ProtParamData, key, None)
        if attr is not None:
            _scales[key] = attr
        else:
            LOGGER.info("ProtParamData.%s unavailable; hydropathy scale skipped.", key)
    SCALE_LIBRARY = _scales
else:
    SCALE_LIBRARY = {}


def _require_biopython() -> None:
    if not BIOPYTHON_AVAILABLE:
        raise ImportError(
            "Biopython is required for protein utilities. Install it with 'pip install biopython'."
        )


@dataclass(frozen=True)
class ProteinSummary:
    sequence: str
    length: int
    molecular_weight: float
    aromaticity: float
    instability_index: float
    gravy: float
    charge_at_pH7: float
    aa_percentages: Dict[str, float]


@dataclass(frozen=True)
class HydropathyWindow:
    start: int
    end: int
    score: float


def clean_sequence(raw: str) -> str:
    """Uppercase and strip whitespace/FASTA headers."""
    if not raw:
        raise ValueError("Empty protein sequence.")
    sequence: List[str] = []
    for line in raw.splitlines():
        if line.startswith(">"):
            continue
        for char in line.strip().upper():
            if char in {" ", "\t"}:
                continue
            if char not in VALID_AA:
                raise ValueError(f"Invalid residue '{char}'.")
            sequence.append(char)
    if not sequence:
        raise ValueError("No residues found in the provided sequence.")
    return "".join(sequence)


def summarize_sequence(sequence: str) -> ProteinSummary:
    _require_biopython()
    cleaned = clean_sequence(sequence)
    analysis = ProteinAnalysis(cleaned)
    return ProteinSummary(
        sequence=cleaned,
        length=len(cleaned),
        molecular_weight=analysis.molecular_weight(),
        aromaticity=analysis.aromaticity(),
        instability_index=analysis.instability_index(),
        gravy=analysis.gravy(),
        charge_at_pH7=analysis.charge_at_pH(7.0),
        aa_percentages=analysis.get_amino_acids_percent(),
    )


def hydropathy_profile(sequence: str, *, window: int = 9, step: int = 1, scale: str = "kd") -> List[HydropathyWindow]:
    _require_biopython()
    if window <= 0:
        raise ValueError("window must be > 0")
    if step <= 0:
        raise ValueError("step must be > 0")
    cleaned = clean_sequence(sequence)
    if window > len(cleaned):
        raise ValueError("window larger than sequence length")

    scale_dict = SCALE_LIBRARY.get(scale.lower())
    if scale_dict is None:
        raise KeyError(f"Unknown scale '{scale}'. Available: {', '.join(SCALE_LIBRARY)}")

    analysis = ProteinAnalysis(cleaned)
    scores = analysis.protein_scale(scale_dict, window=window)
    windows: List[HydropathyWindow] = []
    for idx in range(0, len(scores), step):
        windows.append(HydropathyWindow(start=idx, end=idx + window, score=scores[idx]))
    return windows


def load_fasta_sequences(path: Path) -> Dict[str, str]:
    _require_biopython()
    records: Dict[str, str] = {}
    for record in SeqIO.parse(str(path), "fasta"):
        records[record.id] = clean_sequence(str(record.seq))
    return records


def _select_sequence(sequence: Optional[str], path: Optional[Path]) -> str:
    if path or sequence:
        _require_biopython()
    if sequence and path:
        raise ValueError("Provide either a sequence or --input, not both.")
    if sequence:
        return clean_sequence(sequence)
    if path:
        sequences = load_fasta_sequences(path)
        if not sequences:
            raise ValueError(f"No FASTA records found in {path}.")
        first_id = next(iter(sequences))
        print(f"Loaded record '{first_id}' from {path}")
        return sequences[first_id]
    raise ValueError("Provide a sequence or --input path.")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize a protein sequence and hydropathy profile.")
    parser.add_argument(
        "sequence",
        nargs="?",
        help="Inline amino-acid string.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        help="Path to a FASTA file. The first record will be analysed.",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=9,
        help="Window size for hydropathy (default: 9).",
    )
    parser.add_argument(
        "--step",
        type=int,
        default=1,
        help="Step size when sampling the hydropathy profile (default: 1).",
    )
    parser.add_argument(
        "--scale",
        default="kd",
        help=f"Hydropathy scale ({', '.join(SCALE_LIBRARY)}) (default: kd).",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=5,
        help="Display the top-N hydrophobic windows (default: 5).",
    )
    return parser.parse_args()


def main() -> None:
    _require_biopython()
    args = _parse_args()
    sequence = _select_sequence(args.sequence, args.input)
    summary = summarize_sequence(sequence)
    print(f"Length: {summary.length}")
    print(f"Molecular weight: {summary.molecular_weight:.2f} Da")
    print(f"Aromaticity: {summary.aromaticity:.3f}")
    print(f"Instability index: {summary.instability_index:.2f}")
    print(f"GRAVY: {summary.gravy:.3f}")
    print(f"Charge @ pH 7.0: {summary.charge_at_pH7:.2f}")

    windows = []
    if summary.length >= args.window:
        windows = hydropathy_profile(
            summary.sequence,
            window=args.window,
            step=args.step,
            scale=args.scale,
        )
        sorted_windows = sorted(windows, key=lambda w: w.score, reverse=True)
        print(f"\nTop {min(args.top, len(sorted_windows))} hydrophobic windows (scale={args.scale}):")
        for window in sorted_windows[: args.top]:
            print(f"{window.start:>4}-{window.end:<4}\tscore={window.score:.3f}")
    else:
        print("Hydropathy profile skipped (sequence shorter than the requested window).")


if __name__ == "__main__":
    main()
