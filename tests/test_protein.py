from pathlib import Path
import sys

import pytest

pytest.importorskip("Bio")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from helix.protein import clean_sequence, hydropathy_profile, load_fasta_sequences, summarize_sequence


def test_clean_sequence_handles_fasta_and_whitespace():
    raw = ">sample\nac d e\nFGH"
    assert clean_sequence(raw) == "ACDEFGH"


def test_summarize_sequence_reports_basic_metrics():
    seq = "ACDEFGHIKLMNPQRSTVWY"
    summary = summarize_sequence(seq)
    assert summary.length == len(seq)
    assert summary.sequence == seq
    assert summary.molecular_weight > 2000
    assert summary.aa_percentages["A"] == pytest.approx(1 / len(seq))


def test_hydropathy_profile_sampling():
    seq = "ACDEFGHIKLMNPQRSTVWY"
    windows = hydropathy_profile(seq, window=5, step=2)
    assert len(windows) == 8
    assert windows[0].start == 0
    assert windows[0].end == 5


def test_load_fasta_sequences(tmp_path: Path):
    fasta = tmp_path / "demo.fasta"
    fasta.write_text(">p1\nACDE\n>p2\nFGHI", encoding="utf-8")
    records = load_fasta_sequences(fasta)
    assert records["p1"] == "ACDE"
    assert records["p2"] == "FGHI"
