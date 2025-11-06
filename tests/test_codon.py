from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from codon import (
    Frameshift,
    Orf,
    detect_frameshifts,
    find_orfs,
    frameshifts_to_csv,
    orfs_to_csv,
    orfs_to_fasta,
    translate_rna,
    START_CODON,
)


@pytest.mark.parametrize(
    "rna,expected",
    [
        ("AUGGCCUUU", "MAF"),
        ("aug-gcu uaa", "MA"),
        ("AUGUGAUAA", "M"),
    ],
)
def test_translate_rna_basic(rna, expected):
    assert translate_rna(rna) == expected


def test_translate_rna_stop_symbol():
    assert translate_rna("AUGUGAUAA", stop_at_stop=False) == "M**"


def test_translate_rna_invalid_length():
    with pytest.raises(ValueError):
        translate_rna("AUGG")


def test_translate_rna_invalid_base():
    with pytest.raises(ValueError):
        translate_rna("AUGXCC")


def test_find_orfs_basic():
    seq = "AUGGCCUUUUAA"
    orfs = find_orfs(seq, min_length=9)
    assert len(orfs) == 1
    orf = orfs[0]
    assert orf.start == 0
    assert orf.end == len(seq)
    assert orf.peptide == "MAF"
    assert orf.frame == 0
    assert orf.strand == "+"


def test_find_orfs_multiple_frames():
    seq = "AUGAAAUAA" + "C" + "AUGAAAUAA"
    orfs = find_orfs(seq, min_length=9)
    frames = {(orf.frame, orf.strand) for orf in orfs}
    assert frames == {(0, "+"), (1, "+")}


def test_find_orfs_partial():
    seq = "CCAUGGCCUUU"
    orfs = find_orfs(seq, min_length=9, include_partial=True)
    assert len(orfs) == 1
    assert orfs[0].peptide == "MAF"
    assert orfs[0].end == len(seq.replace("T", "U"))
    assert orfs[0].strand == "+"


def test_find_orfs_reverse_strand():
    seq = "CTATTCCATTCATTA"
    orfs = find_orfs(seq, min_length=6)
    assert any(orf.strand == "-" for orf in orfs)
    rev = [orf for orf in orfs if orf.strand == "-"][0]
    assert rev.peptide.startswith("M")


def test_detect_frameshifts():
    seq = "AUGAAAUAAGCAUGGGCUAA"
    events = detect_frameshifts(seq, min_orf_length=9, gap_tolerance=2)
    assert events
    event = events[0]
    assert event.shift in {1, 2}
    assert event.strand == "+"


def test_orfs_export_helpers():
    seq = "AUGGCCUUUUAA"
    orfs = find_orfs(seq, min_length=9)
    fasta = orfs_to_fasta(orfs, line_width=10)
    assert fasta.startswith(">orf_1")
    rows = orfs_to_csv(orfs)
    assert rows[0]["length_aa"] == len(orfs[0].peptide)


def test_frameshift_export_helpers():
    seq = "AUGAAAUAAGCAUGGGCUAA"
    events = detect_frameshifts(seq, min_orf_length=9, gap_tolerance=2)
    rows = frameshifts_to_csv(events)
    assert rows[0]["id"].startswith("frameshift")
