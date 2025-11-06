from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from triage import compute_triage_report


def test_compute_triage_report_basic():
    seq = "AUGGCCUUUUAA"
    report = compute_triage_report(seq, k=3, max_diff=0, min_orf_length=9)
    assert report.sequence == seq
    assert len(report.skew) == len(seq) + 1
    assert report.clusters
    assert any(cluster.canonical == "UUU" for cluster in report.clusters)
    assert report.orfs  # ORF should be detected for the sequence
