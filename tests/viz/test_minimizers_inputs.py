from helix.viz.minimizers import plot_minimizer_density


def test_minimizer_mixed_input_forms(tmp_path):
    minimizers = [
        5,
        (10, "AAA"),
        (15, "CCC", 99),
        {"pos": "20"},
        {"position": 25},
        {"ref_start": 30},
        {"ignored": 40},
    ]
    _, spec = plot_minimizer_density(
        sequence_length=100,
        minimizers=minimizers,
        bin_count=10,
        save=str(tmp_path / "mixed.png"),
        save_viz_spec=str(tmp_path / "mixed.viz.json"),
    )
    assert spec["primitives"]["total_minimizers"] == 6


def test_minimizer_out_of_range_and_empty(tmp_path):
    minimizers = [-5, 150]
    _, spec = plot_minimizer_density(
        sequence_length=100,
        minimizers=minimizers,
        bin_count=10,
        save=str(tmp_path / "oor.png"),
        save_viz_spec=str(tmp_path / "oor.viz.json"),
    )
    assert spec["primitives"]["total_minimizers"] == 0
    assert spec["primitives"]["density_max"] == 0


def test_minimizer_bin_count_clamped(tmp_path):
    minimizers = [0, 1, 2]
    _, spec = plot_minimizer_density(
        sequence_length=5,
        minimizers=minimizers,
        bin_count=200,
        save=str(tmp_path / "clamp.png"),
        save_viz_spec=str(tmp_path / "clamp.viz.json"),
    )
    assert spec["meta"]["bin_count"] == 5
    assert spec["primitives"]["points"] == 5
