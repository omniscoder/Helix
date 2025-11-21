"""Shared helpers for engine benchmark reporting."""

from __future__ import annotations

import platform
from dataclasses import dataclass
from time import perf_counter
from typing import Any, Mapping, Sequence

import numpy as np

from .. import bioinformatics
from ..crispr.model import CasSystem, CasSystemType, GuideRNA, PAMRule, DigitalGenome as CrisprDigitalGenome
from ..crispr.physics import create_crispr_physics, score_pairs_encoded_multi
from ..config import native_backend_available, resolve_prime_backend
from ..prime.model import PegRNA, PrimeEditor
from ..prime.simulator import PrimeTargetRequest, predict_prime_outcomes_for_targets
from .. import __version__ as HELIX_VERSION
from ..prime.physics import PRIME_SCORING_VERSION
from ..crispr.physics import CRISPR_SCORING_VERSION

_BENCH_BASES = np.array(list("ACGT"))


def _random_sequences(count: int, length: int, rng: np.random.Generator) -> list[str]:
    encoded = rng.integers(0, 4, size=(count, length))
    return ["".join(_BENCH_BASES[row]) for row in encoded]


def _encode_sequence_batch(sequences: Sequence[str]) -> np.ndarray:
    from ..engine.encoding import encode_sequence_to_uint8

    return np.stack([encode_sequence_to_uint8(seq) for seq in sequences])


def _query_gpu_name() -> str | None:
    import subprocess

    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            check=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    line = result.stdout.strip().splitlines()
    if not line:
        return None
    return line[0].strip() or None


def collect_benchmark_env() -> dict[str, Any]:
    env: dict[str, Any] = {
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "cuda_available": False,
        "gpu_name": None,
        "native_backend_available": native_backend_available(),
        "backends_built": ["cpu-reference"],
    }
    try:
        from helix_engine import native as native_engine  # type: ignore
    except Exception:
        native_engine = None
    if native_engine is not None:
        native_ok = bool(getattr(native_engine, "is_available", lambda: False)())
        if native_ok:
            env["backends_built"].append("native-cpu")
        cuda_ok = native_ok and bool(getattr(native_engine, "cuda_available", lambda: False)())
        env["cuda_available"] = cuda_ok
        if cuda_ok:
            env["backends_built"].append("gpu")
            env["gpu_name"] = _query_gpu_name()
    env["backends_built"] = sorted(set(env["backends_built"]))
    return env


def run_crispr_benchmarks(
    backends: Sequence[str],
    shapes: Sequence[tuple[int, int, int]],
    seed: int,
) -> list[dict[str, Any]]:
    rng = np.random.default_rng(seed)
    cas = CasSystem(
        name="bench-cas9",
        system_type=CasSystemType.CAS9,
        pam_rules=[PAMRule(pattern="NGG")],
        cut_offset=3,
    )
    results: list[dict[str, Any]] = []
    for g, n, l in shapes:
        guide_strs = _random_sequences(g, l, rng)
        guides_encoded = _encode_sequence_batch(guide_strs)
        window_strs = _random_sequences(n, l, rng)
        windows_encoded = _encode_sequence_batch(window_strs)
        pam_penalties = np.zeros((g, n), dtype=np.float32)
        for backend in backends:
            entry: dict[str, Any] = {
                "backend_requested": backend,
                "backend_used": None,
                "g": g,
                "n": n,
                "l": l,
                "shape": f"{g}x{n}x{l}",
            }
            try:
                physics_list = [create_crispr_physics(cas, GuideRNA(sequence=seq), backend=backend) for seq in guide_strs]
            except RuntimeError as exc:
                entry["error"] = str(exc)
                results.append(entry)
                continue
            start = perf_counter()
            score_pairs_encoded_multi(physics_list, guides_encoded, windows_encoded, pam_penalties=pam_penalties)
            elapsed = perf_counter() - start
            mpairs = (g * n) / max(elapsed, 1e-9) / 1e6
            entry["backend_used"] = getattr(physics_list[0], "backend_name", backend)
            entry["mpairs_per_s"] = float(mpairs)
            entry["elapsed_seconds"] = float(elapsed)
            results.append(entry)
    return results


def run_prime_benchmarks(
    workloads: Sequence[tuple[int, int, int]],
    seed: int,
) -> list[dict[str, Any]]:
    rng = np.random.default_rng(seed)
    cas = CasSystem(name="primer", system_type=CasSystemType.CAS9, pam_rules=[PAMRule(pattern="NGG")], cut_offset=3)
    editor = PrimeEditor(name="PE-bench", cas=cas, nick_to_edit_offset=1, efficiency_scale=0.8, indel_bias=0.1)
    backend, _ = resolve_prime_backend(None, use_gpu=False)
    results: list[dict[str, Any]] = []
    for targets, genome_len, spacer_len in workloads:
        genome_seq = "".join(rng.choice(list("ACGT")) for _ in range(genome_len))
        genome = CrisprDigitalGenome({"chr_bench": genome_seq})
        requests: list[PrimeTargetRequest] = []
        max_start = max(0, genome_len - spacer_len)
        for idx in range(targets):
            start = rng.integers(0, max_start + 1)
            spacer = genome_seq[start : start + spacer_len]
            peg = PegRNA(spacer=spacer, pbs="GAAAC", rtt="TTTTAA", name=f"peg_{idx}")
            requests.append(
                PrimeTargetRequest(
                    target_id=f"target_{idx}",
                    genome=genome,
                    peg=peg,
                    editor=editor,
                )
            )
        start_time = perf_counter()
        # Inline predict to avoid extra serialization; use simple loop to avoid large allocations.
        predict_prime_outcomes_for_targets(requests)
        elapsed = perf_counter() - start_time
        preds_per_s = targets / max(elapsed, 1e-9)
        results.append(
            {
                "backend_requested": backend,
                "backend_used": backend,
                "workload": f"{targets}x{genome_len}x{spacer_len}",
                "targets": targets,
                "genome_len": genome_len,
                "spacer_len": spacer_len,
                "predictions": targets,
                "predictions_per_s": float(preds_per_s),
                "elapsed_seconds": float(elapsed),
            }
        )
    return results


@dataclass
class BenchmarkConfig:
    backends: Sequence[str]
    crispr_shapes: Sequence[tuple[int, int, int]]
    prime_workloads: Sequence[tuple[int, int, int]]
    seed: int = 0


def run_benchmark(config: BenchmarkConfig) -> dict[str, Any]:
    crispr_results = run_crispr_benchmarks(config.backends, config.crispr_shapes, config.seed)
    prime_results = run_prime_benchmarks(config.prime_workloads, config.seed)
    return {
        "helix_version": HELIX_VERSION,
        "scoring_versions": {
            "crispr": CRISPR_SCORING_VERSION,
            "prime": PRIME_SCORING_VERSION,
        },
        "env": collect_benchmark_env(),
        "seed": config.seed,
        "config": {
            "backends": list(config.backends),
            "crispr_shapes": [f"{g}x{n}x{l}" for g, n, l in config.crispr_shapes],
            "prime_workloads": [f"{t}x{g}x{s}" for t, g, s in config.prime_workloads],
        },
        "benchmarks": {
            "crispr": crispr_results,
            "prime": prime_results,
        },
    }
