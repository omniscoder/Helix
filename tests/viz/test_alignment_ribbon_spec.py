from helix.viz.alignment import plot_alignment_ribbon


def test_alignment_ribbon_spec(tmp_path):
    alignment = {
        "ref_start": 10,
        "read_start": 5,
        "cigar": "5M2I3M1D2M",
        "matches": 10,
        "score": 42,
    }
    _, spec = plot_alignment_ribbon(
        ref_length=100,
        qry_length=80,
        alignment=alignment,
        save=str(tmp_path / "aln.png"),
        save_viz_spec=str(tmp_path / "aln.viz.json"),
    )
    assert spec["kind"] == "alignment_ribbon"
    assert spec["meta"]["ref_length"] == 100
    assert spec["primitives"]["insertions"] == 1
    assert spec["primitives"]["deletions"] == 1
    assert spec["primitives"]["path_points"] > 0
