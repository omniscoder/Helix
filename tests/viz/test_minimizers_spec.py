from helix.viz.minimizers import plot_minimizer_density


def test_minimizer_density_viz_spec_stable(tmp_path):
    seq_len = 1000
    mins = [(10, "AAA", 1), (50, "TTT", 2), (200, "GGA", 3), (205, "CCC", 4), 900]
    fig, spec = plot_minimizer_density(
        sequence_length=seq_len,
        minimizers=mins,
        bin_count=20,
        save=str(tmp_path / "min.png"),
        save_viz_spec=str(tmp_path / "min.viz.json"),
    )
    assert spec["kind"] == "minimizer_density"
    assert spec["meta"]["sequence_length"] == seq_len
    assert spec["meta"]["bin_count"] == 20
    assert spec["primitives"]["total_minimizers"] == len(mins)
    assert spec["primitives"]["points"] == 20
    assert 1 <= spec["primitives"]["density_max"] <= 2
