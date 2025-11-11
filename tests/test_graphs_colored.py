from helix.graphs import build_colored_dbg, pseudoalign


def test_colored_dbg_presence_and_pseudoalign():
    reads_by_sample = {
        "s1": ["AAGGTT"],
        "s2": ["CCGGTT"],
    }
    colored = build_colored_dbg(reads_by_sample, k=3)
    assert "AA" in colored.presence
    assert "s1" in colored.presence["AA"]
    hits = pseudoalign("AAGGTT", colored)
    assert "s1" in hits
    assert "s2" not in hits
