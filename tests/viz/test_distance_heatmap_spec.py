from helix.viz.distance import plot_distance_heatmap


def test_distance_heatmap_spec(tmp_path):
    matrix = [
        [0.0, 0.1, 0.2],
        [0.1, 0.0, 0.15],
        [0.2, 0.15, 0.0],
    ]
    labels = ["A", "B", "C"]
    _, spec = plot_distance_heatmap(
        matrix=matrix,
        labels=labels,
        method="minhash",
        save=str(tmp_path / "dist.png"),
        save_viz_spec=str(tmp_path / "dist.viz.json"),
    )
    assert spec["kind"] == "distance_heatmap"
    assert spec["meta"]["n"] == 3
    assert spec["primitives"]["min"] == 0.0
    assert spec["primitives"]["max"] == 0.2
