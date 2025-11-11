import pytest

from helix.viz.rna import plot_rna_dotplot


def test_rna_dotplot_rejects_non_square():
    posterior = [
        [0.0, 0.5],
        [0.5],
    ]
    with pytest.raises(AssertionError, match="posterior must be square"):
        plot_rna_dotplot(posterior=posterior)


def test_rna_dotplot_all_zero_matrix(tmp_path):
    posterior = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ]
    _, spec = plot_rna_dotplot(
        posterior=posterior,
        save=str(tmp_path / "zero.png"),
        save_viz_spec=str(tmp_path / "zero.viz.json"),
    )
    assert spec["primitives"]["nonzero_pairs"] == 0
    assert all(q == 0.0 for q in spec["primitives"]["quantiles"])
