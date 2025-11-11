from helix.viz.motif import plot_motif_logo


def test_motif_logo_spec(tmp_path):
    pwm = [
        {"A": 0.7, "C": 0.1, "G": 0.1, "T": 0.1},
        {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
        {"A": 0.05, "C": 0.05, "G": 0.8, "T": 0.1},
    ]
    _, spec = plot_motif_logo(
        pwm=pwm,
        title="Test motif",
        save=str(tmp_path / "motif.png"),
        save_viz_spec=str(tmp_path / "motif.viz.json"),
    )
    assert spec["kind"] == "motif_logo"
    assert spec["meta"]["columns"] == 3
    assert spec["primitives"]["max_information"] <= 2.0
    assert spec["primitives"]["mean_information"] > 0.0
