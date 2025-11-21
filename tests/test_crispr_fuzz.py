"""Fuzz and edge-case tests for CRISPR physics backends."""

from __future__ import annotations

import numpy as np
import pytest

from helix.engine.encoding import encode_sequence_to_uint8
from helix.crispr.model import CasSystem, CasSystemType, GuideRNA, PAMRule
from helix.crispr.physics import create_crispr_physics, score_pairs_encoded

CAS = CasSystem(
    name="fuzz-cas9",
    system_type=CasSystemType.CAS9,
    pam_rules=[PAMRule(pattern="NGG")],
    cut_offset=3,
    max_mismatches=3,
    weight_mismatch_penalty=1.0,
    weight_pam_penalty=2.0,
)
GUIDE = GuideRNA(sequence="ACCCAGGAAACCCGGGTTTT")
GUIDE_ENCODED = encode_sequence_to_uint8(GUIDE.sequence)

EDGE_CASE_WINDOWS = [
    "A" * len(GUIDE.sequence),
    "G" * len(GUIDE.sequence),
    "AC" * (len(GUIDE.sequence) // 2),
    "GT" * (len(GUIDE.sequence) // 2),
    "AG" * (len(GUIDE.sequence) // 2),
]


def _make_physics(backend: str):
    try:
        physics = create_crispr_physics(CAS, GUIDE, backend=backend)
    except RuntimeError as exc:  # pragma: no cover - availability
        pytest.skip(str(exc))
    actual_backend = getattr(physics, "backend_name", backend)
    if backend == "gpu" and actual_backend != "gpu":
        pytest.skip("CUDA backend unavailable")
    return physics


@pytest.mark.slow
@pytest.mark.parametrize("backend", ["native-cpu", "gpu"])
def test_crispr_random_fuzz_large_batch(backend: str) -> None:
    test_physics = _make_physics(backend)
    ref_physics = create_crispr_physics(CAS, GUIDE, backend="cpu-reference")

    rng = np.random.default_rng(42)
    guide_count = 4
    window_count = 8192
    guide_array = np.repeat(GUIDE_ENCODED[None, :], guide_count, axis=0)
    windows = rng.integers(0, 4, size=(window_count, GUIDE_ENCODED.size), dtype=np.uint8)
    pam = np.zeros((guide_count, window_count), dtype=np.float32)

    ref_scores = score_pairs_encoded(guide_array, windows, ref_physics, pam_penalties=pam)
    test_scores = score_pairs_encoded(guide_array, windows, test_physics, pam_penalties=pam)
    np.testing.assert_allclose(test_scores, ref_scores, rtol=1e-6, atol=1e-7)


@pytest.mark.parametrize("backend", ["cpu-reference", "native-cpu", "gpu"])
def test_crispr_edge_cases_no_nan(backend: str) -> None:
    physics = _make_physics(backend)
    ref_physics = create_crispr_physics(CAS, GUIDE, backend="cpu-reference")
    edge_windows = np.stack([encode_sequence_to_uint8(seq) for seq in EDGE_CASE_WINDOWS], dtype=np.uint8)
    guide_array = GUIDE_ENCODED[None, :]
    pam = np.zeros((1, edge_windows.shape[0]), dtype=np.float32)

    scores = score_pairs_encoded(guide_array, edge_windows, physics, pam_penalties=pam)
    ref_scores = score_pairs_encoded(guide_array, edge_windows, ref_physics, pam_penalties=pam)
    assert np.all(np.isfinite(scores)), f"Non-finite scores for backend={backend}"
    np.testing.assert_allclose(scores, ref_scores, rtol=1e-7, atol=1e-8)
