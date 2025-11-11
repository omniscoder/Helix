from helix.viz.rna import plot_rna_dotplot


def test_rna_dotplot_viz_spec_stable(tmp_path):
    posterior = [
        [0.0, 0.7, 0.0, 0.0],
        [0.7, 0.0, 0.4, 0.0],
        [0.0, 0.4, 0.0, 0.2],
        [0.0, 0.0, 0.2, 0.0],
    ]
    fig, spec = plot_rna_dotplot(
        posterior=posterior,
        save=str(tmp_path / "rna.png"),
        save_viz_spec=str(tmp_path / "rna.viz.json"),
    )
    assert spec["kind"] == "rna_dotplot"
    assert spec["meta"]["n"] == 4
    # Only upper-triangle entries are visualized, so we count 3 unique pairs.
    assert spec["primitives"]["nonzero_pairs"] == 3
    quantiles = spec["primitives"]["quantiles"]
    assert len(quantiles) == 6
    assert 0.0 <= quantiles[0] <= quantiles[-1] <= 1.0
