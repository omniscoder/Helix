"""CRISPR backend selection error handling tests."""

from __future__ import annotations

import pytest

from helix.cli import _attach_crispr_engine_meta
from helix.crispr.model import CasSystem, CasSystemType, GuideRNA, PAMRule
from helix.crispr.physics import create_crispr_physics, get_last_crispr_backend


CAS = CasSystem(
    name="cas9",
    system_type=CasSystemType.CAS9,
    pam_rules=[PAMRule(pattern="NGG")],
    cut_offset=3,
    max_mismatches=3,
    weight_mismatch_penalty=1.0,
    weight_pam_penalty=2.0,
)
GUIDE = GuideRNA(sequence="ACCCAGGAAACCCGGGTTTT")


def test_native_cpu_backend_errors_without_native(monkeypatch):
    monkeypatch.setattr("helix_engine.native.is_available", lambda: False)
    monkeypatch.setattr("helix.crispr.physics.native_backend_available", lambda: False)
    monkeypatch.delenv("HELIX_CRISPR_ALLOW_FALLBACK", raising=False)
    with pytest.raises(RuntimeError, match="Native CRISPR engine is unavailable"):
        create_crispr_physics(CAS, GUIDE, backend="native-cpu")


def test_gpu_backend_errors_without_cuda_when_no_fallback(monkeypatch):
    monkeypatch.setattr("helix_engine.native.cuda_available", lambda: False)
    monkeypatch.delenv("HELIX_CRISPR_ALLOW_FALLBACK", raising=False)
    with pytest.raises(RuntimeError, match="Helix CUDA backend is unavailable"):
        create_crispr_physics(CAS, GUIDE, backend="gpu")


def test_gpu_backend_fallback_updates_metadata(monkeypatch, caplog):
    monkeypatch.setenv("HELIX_CRISPR_BACKEND", "gpu")
    monkeypatch.setenv("HELIX_CRISPR_ALLOW_FALLBACK", "1")
    monkeypatch.setattr("helix_engine.native.cuda_available", lambda: False)
    caplog.set_level("WARNING", logger="helix.crispr.physics")
    physics = create_crispr_physics(CAS, GUIDE, use_gpu=True)
    assert physics.backend_name == "cpu-reference"
    assert "Falling back from backend='gpu'" in caplog.text
    assert get_last_crispr_backend() == "cpu-reference"
    meta: dict[str, str] = {}
    _attach_crispr_engine_meta(meta, use_gpu=True)
    assert meta["crispr_engine_backend"] == "cpu-reference"
