from helix.viz.seed_chain import plot_seed_chain


def test_seed_chain_viz_spec_stable(tmp_path):
    chains = [
        [
            {"ref_start": 100, "ref_end": 140, "qry_start": 80, "qry_end": 120},
            {"ref_start": 150, "ref_end": 170, "qry_start": 130, "qry_end": 150},
        ],
        [
            {"ref_start": 400, "ref_end": 450, "qry_start": 390, "qry_end": 440},
        ],
    ]
    fig, spec = plot_seed_chain(
        ref_length=1000,
        qry_length=1000,
        chains=chains,
        save=str(tmp_path / "chains.png"),
        save_viz_spec=str(tmp_path / "chains.viz.json"),
    )
    assert spec["kind"] == "seed_chain"
    assert spec["meta"]["chains"] == 2
    assert spec["primitives"]["line_segments"] == 3
