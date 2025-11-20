"""Golden fixtures for CRISPR physics backends."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from helix.crispr.model import CasSystem, CasSystemType, GuideRNA, PAMRule
from helix.crispr.physics import (
    CRISPRPhysicsCPU,
    compute_on_target_score_encoded,
    create_crispr_physics,
)

FIXTURE_DIR = Path(__file__).resolve().parent / "data" / "engine_crispr"
FIXTURES = json.loads((FIXTURE_DIR / "fixtures.json").read_text())
GUIDES_ENC = np.load(FIXTURE_DIR / "guides_encoded.npy")
WINDOWS_ENC = np.load(FIXTURE_DIR / "windows_encoded.npy")
SCORES_REFERENCE = np.load(FIXTURE_DIR / "scores_cpu_reference.npy")

CAS = CasSystem(
    name="demo",
    system_type=CasSystemType.CAS9,
    pam_rules=[PAMRule(pattern="NGG")],
    cut_offset=3,
    max_mismatches=3,
    weight_mismatch_penalty=1.0,
    weight_pam_penalty=2.0,
)


def _reference_scores() -> np.ndarray:
    guide = GuideRNA(sequence=FIXTURES["guides"][0])
    physics = CRISPRPhysicsCPU(CAS, guide)
    scores = np.zeros_like(SCORES_REFERENCE)
    for g_idx, guide_vec in enumerate(GUIDES_ENC):
        for w_idx, window_vec in enumerate(WINDOWS_ENC):
            scores[g_idx, w_idx] = compute_on_target_score_encoded(
                guide_vec,
                window_vec,
                physics._weight_array,  # type: ignore[attr-defined]
                bulge_penalty=physics.bulge_penalty,  # type: ignore[attr-defined]
                gc_range=physics.gc_opt_range,  # type: ignore[attr-defined]
                pam_penalty=0.0,
                min_score=physics.min_score,  # type: ignore[attr-defined]
            )
    return scores


def test_cpu_reference_matches_golden() -> None:
    scores = _reference_scores()
    assert np.allclose(scores, SCORES_REFERENCE)


def test_native_matches_golden_if_available() -> None:
    try:
        physics = create_crispr_physics(CAS, GuideRNA(sequence=FIXTURES["guides"][0]), backend="native-cpu")
    except RuntimeError:
        pytest.skip("Native backend not available")
    scores = physics.score_pairs_encoded_batch(GUIDES_ENC, WINDOWS_ENC, np.zeros_like(SCORES_REFERENCE))  # type: ignore[attr-defined]
    assert np.allclose(scores, SCORES_REFERENCE)
