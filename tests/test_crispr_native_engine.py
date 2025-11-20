"""Parity tests for the native CRISPR physics backend.

The entire module is skipped unless the helix_engine._native extension can be
imported. That keeps CI green on pure-Python environments while still giving us
strong guarantees whenever the native backend is available locally.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("helix_engine._native")

from helix_engine import native as native_engine

from helix.crispr.model import CasSystem, CasSystemType, GuideRNA, PAMRule, DigitalGenome
from helix.crispr.physics import (
    _encode_sequence_to_uint8,
    compute_on_target_score_encoded as py_compute,
    create_crispr_physics,
)
from helix.crispr.simulator import find_candidate_sites


GUIDE_SEQ = "ACCCAGGAAACCCGGGTTTT"
PAM_PATTERN = "NGG"


def _demo_cas() -> CasSystem:
    return CasSystem(
        name="demo-cas9",
        system_type=CasSystemType.CAS9,
        pam_rules=[PAMRule(pattern=PAM_PATTERN)],
        cut_offset=3,
    )


def _demo_genome() -> DigitalGenome:
    plus = "TTT" + GUIDE_SEQ + "AGGTTT"
    rc = "TTT" + GUIDE_SEQ[::-1] + "TTT"
    return DigitalGenome({"chr_plus": plus, "chr_rc": rc})


def test_native_scalar_matches_reference() -> None:
    cas = _demo_cas()
    guide = GuideRNA(sequence=GUIDE_SEQ)
    physics = create_crispr_physics(cas, guide, backend="cpu-reference")
    window_seq = GUIDE_SEQ
    window_encoded = _encode_sequence_to_uint8(window_seq)
    score_py = py_compute(
        physics._guide_encoded,  # type: ignore[attr-defined]
        window_encoded,
        physics._weight_array,  # type: ignore[attr-defined]
        bulge_penalty=physics.bulge_penalty,  # type: ignore[attr-defined]
        gc_range=physics.gc_opt_range,  # type: ignore[attr-defined]
        pam_penalty=0.0,
        min_score=physics.min_score,  # type: ignore[attr-defined]
    )
    score_native = native_engine.compute_on_target_score_encoded(
        physics._guide_encoded,  # type: ignore[attr-defined]
        window_encoded,
        physics._weight_array,  # type: ignore[attr-defined]
        bulge_penalty=physics.bulge_penalty,  # type: ignore[attr-defined]
        gc_low=physics.gc_opt_range[0],  # type: ignore[attr-defined]
        gc_high=physics.gc_opt_range[1],  # type: ignore[attr-defined]
        pam_penalty=0.0,
        min_score=physics.min_score,  # type: ignore[attr-defined]
    )
    assert pytest.approx(score_py, rel=0, abs=1e-9) == score_native


def test_native_batch_matches_reference() -> None:
    cas = _demo_cas()
    guide = GuideRNA(sequence=GUIDE_SEQ)
    physics = create_crispr_physics(cas, guide, backend="cpu-reference")
    guides = np.expand_dims(physics._guide_encoded, axis=0)  # type: ignore[attr-defined]
    windows = np.stack(
        [
            _encode_sequence_to_uint8(GUIDE_SEQ),
            _encode_sequence_to_uint8("ACCCAGGAAACCCGGGTTTA"),
            _encode_sequence_to_uint8("CCCCAGGAAACCCGGGTTTT"),
        ]
    )
    pam_penalties = np.zeros((guides.shape[0], windows.shape[0]), dtype=np.float32)
    native_scores = native_engine.score_pairs_encoded(
        guides,
        windows,
        physics._weight_array,  # type: ignore[attr-defined]
        bulge_penalty=physics.bulge_penalty,  # type: ignore[attr-defined]
        gc_low=physics.gc_opt_range[0],  # type: ignore[attr-defined]
        gc_high=physics.gc_opt_range[1],  # type: ignore[attr-defined]
        pam_penalties=pam_penalties,
        min_score=physics.min_score,  # type: ignore[attr-defined]
    )
    ref_scores = np.zeros_like(native_scores)
    for row in range(guides.shape[0]):
        for col in range(windows.shape[0]):
            ref_scores[row, col] = py_compute(
                guides[row],
                windows[col],
                physics._weight_array,  # type: ignore[attr-defined]
                bulge_penalty=physics.bulge_penalty,  # type: ignore[attr-defined]
                gc_range=physics.gc_opt_range,  # type: ignore[attr-defined]
                pam_penalty=0.0,
                min_score=physics.min_score,  # type: ignore[attr-defined]
            )
    assert np.allclose(native_scores, ref_scores)


def test_find_candidate_sites_matches_reference() -> None:
    cas = _demo_cas()
    guide = GuideRNA(sequence=GUIDE_SEQ)
    genome = _demo_genome()
    physics_native = create_crispr_physics(cas, guide, backend="native-cpu")
    physics_ref = create_crispr_physics(cas, guide, backend="cpu-reference")
    sites_native = find_candidate_sites(genome, cas, guide, physics=physics_native, max_sites=4)
    sites_ref = find_candidate_sites(genome, cas, guide, physics=physics_ref, max_sites=4)
    assert len(sites_native) == len(sites_ref)
    for native_site, ref_site in zip(sites_native, sites_ref):
        assert native_site.sequence == ref_site.sequence
        assert pytest.approx(ref_site.on_target_score, rel=0, abs=1e-9) == native_site.on_target_score

