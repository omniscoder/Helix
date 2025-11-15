"""Smoke tests for the sci-fi geometry builders."""

from __future__ import annotations

from helix.crispr.pam import build_crispr_pam_mask
from helix.crispr.simulate import resolve_crispr_priors, simulate_cut_repair
from helix.gui.modern.builders import (
    build_crispr_rail_3d_spec,
    build_crispr_orbitals_3d_spec,
)


def _sample_crispr_payload(sequence: str) -> tuple[dict, dict]:
    guide = {"id": "g1", "start": 4, "end": 20, "strand": "+", "gc_content": 0.5}
    priors = resolve_crispr_priors("default_indel")
    pam = build_crispr_pam_mask(sequence, guide, "SpCas9_NGG", 1.0)
    payload = simulate_cut_repair(
        site_seq=sequence,
        guide=guide,
        priors=priors,
        draws=25,
        seed=123,
        emit_sequence=True,
        pam_mask=pam,
    )
    return guide, payload


def test_crispr_rail_3d_outputs_geometry_payloads() -> None:
    sequence = "ACGT" * 10
    guide, payload = _sample_crispr_payload(sequence)
    spec = build_crispr_rail_3d_spec(sequence, guide, payload)
    assert spec is not None
    meta = spec.metadata
    assert meta is not None and "custom_geometry" in meta
    geom = meta["custom_geometry"]
    assert geom["rail_vertices"].shape[1] == 3
    assert geom["rail_vertices"].shape[0] == len(sequence) * 2
    weights = geom["arc_weights"]
    kinds = geom["arc_kinds"]
    assert weights.shape == kinds.shape == (geom["arc_vertices"].shape[0],)
    assert "camera" in meta and "pos" in meta["camera"]


def test_crispr_orbitals_3d_has_ring_data() -> None:
    sequence = "ACGT" * 8
    guide, payload = _sample_crispr_payload(sequence)
    spec = build_crispr_orbitals_3d_spec(sequence, guide, payload)
    assert spec is not None
    geom = spec.metadata["custom_geometry"]
    assert geom["arc_vertices"].size > 0
    assert geom["arc_offsets"].ndim == 1
