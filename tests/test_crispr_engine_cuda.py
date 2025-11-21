"""CUDA backend parity tests."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

from helix.crispr.model import CasSystem, CasSystemType, GuideRNA, PAMRule
from helix.crispr.physics import create_crispr_physics, score_pairs_encoded, _encode_sequence_to_uint8

from tests.test_helix_cli import run_cli

try:  # optional CUDA availability
    from helix_engine import native as native_engine
except Exception:  # pragma: no cover - module missing
    native_engine = None


def _has_cuda() -> bool:
    return bool(native_engine and hasattr(native_engine, "cuda_available") and native_engine.cuda_available())


pytestmark = pytest.mark.skipif(not _has_cuda(), reason="CUDA backend unavailable")

FIXTURE_DIR = Path(__file__).resolve().parent / "data" / "engine_crispr"
GUIDES_ENC = np.load(FIXTURE_DIR / "guides_encoded.npy")
WINDOWS_ENC = np.load(FIXTURE_DIR / "windows_encoded.npy")
SCORES_REFERENCE = np.load(FIXTURE_DIR / "scores_cpu_reference.npy")
GUIDE_SEQ = json.loads((FIXTURE_DIR / "fixtures.json").read_text())["guides"][0]

CAS = CasSystem(
    name="demo",
    system_type=CasSystemType.CAS9,
    pam_rules=[PAMRule(pattern="NGG")],
    cut_offset=3,
    max_mismatches=3,
    weight_mismatch_penalty=1.0,
    weight_pam_penalty=2.0,
)


def test_cuda_backend_matches_golden(monkeypatch):
    monkeypatch.setenv("HELIX_CRISPR_ALLOW_FALLBACK", "0")
    physics = create_crispr_physics(CAS, GuideRNA(sequence=GUIDE_SEQ), backend="gpu")
    scores = physics.score_pairs_encoded_batch(GUIDES_ENC, WINDOWS_ENC, np.zeros_like(SCORES_REFERENCE))  # type: ignore[attr-defined]
    assert np.allclose(scores, SCORES_REFERENCE)


def test_cuda_backend_matches_random_microbatch(monkeypatch):
    monkeypatch.setenv("HELIX_CRISPR_ALLOW_FALLBACK", "0")
    guide = GuideRNA(sequence=GUIDE_SEQ)
    gpu_physics = create_crispr_physics(CAS, guide, backend="gpu")
    cpu_physics = create_crispr_physics(CAS, guide, backend="cpu-reference")
    guide_encoded = _encode_sequence_to_uint8(GUIDE_SEQ)
    guides = np.repeat(guide_encoded[None, :], 3, axis=0)
    rng = np.random.default_rng(1337)
    windows = rng.integers(0, 4, size=(5, len(GUIDE_SEQ)), dtype=np.uint8)
    scores_gpu = score_pairs_encoded(guides, windows, gpu_physics)
    scores_cpu = score_pairs_encoded(guides, windows, cpu_physics)
    assert np.allclose(scores_gpu, scores_cpu, atol=5e-7, rtol=1e-5)


def test_cli_crispr_simulate_reports_gpu(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("HELIX_CRISPR_BACKEND", "gpu")
    monkeypatch.setenv("HELIX_CRISPR_ALLOW_FALLBACK", "0")
    genome_path = tmp_path / "genome.fna"
    genome_path.write_text(">chr\nACCCAGGAAACCCGGGTTTT\n", encoding="utf-8")
    out_path = tmp_path / "cuts.json"
    run_cli(
        "crispr",
        "genome-sim",
        "--genome",
        str(genome_path),
        "--guide-sequence",
        "ACCCAGGAAACCCGGGTTTT",
        "--json",
        str(out_path),
    )
    payload = json.loads(out_path.read_text())
    assert payload["meta"].get("crispr_engine_backend") == "gpu"
