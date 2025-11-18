from __future__ import annotations

import json
import uuid
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from time import perf_counter
from typing import Any, Dict, Mapping, Protocol

from .metrics_utils import run_metrics_to_summary
from .ogn.client import MockOGNClient, OGNJobResult, OGNJobSpec
from .studio.run_metrics import compute_run_metrics
from .crispr.simulate import resolve_crispr_priors, simulate_cut_repair
from .crispr.pam import build_crispr_pam_mask, build_prime_pam_mask
from .prime.priors import resolve_prime_priors
from .prime.simulator import simulate_prime_edit
from .studio.pcr_engine import PCRConfig, simulate_pcr


@dataclass(frozen=True)
class RunRequest:
    run_kind: str
    session_id: str
    run_id: str
    snapshot_payload: Mapping[str, Any]
    parameters: Mapping[str, Any]
    labels: Dict[str, str] = field(default_factory=dict)


@dataclass
class EngineResult:
    snapshot: Mapping[str, Any]
    metrics: Mapping[str, Any] | None = None
    status: str = "succeeded"
    info: Dict[str, Any] | None = None


class Engine(Protocol):
    def execute(self, request: RunRequest) -> EngineResult:  # pragma: no cover - interface
        ...


class LocalEngine:
    def __init__(self, *, version: str = "local-0.1.0") -> None:
        self.version = version

    def execute(self, request: RunRequest) -> EngineResult:
        state = _extract_state(request.snapshot_payload)
        run_kind = request.run_kind.upper()
        if run_kind == "CRISPR":
            snapshot = self._run_crispr(state, request)
        elif run_kind == "PRIME":
            snapshot = self._run_prime(state, request)
        elif run_kind == "PCR":
            snapshot = self._run_pcr(state, request)
        else:
            raise RuntimeError(f"Run kind '{request.run_kind}' is not supported by LocalEngine.")
        metrics = compute_run_metrics(snapshot)
        return EngineResult(
            snapshot=snapshot,
            metrics=run_metrics_to_summary(metrics),
            status="succeeded",
            info={"engine_version": self.version},
        )

    def _run_crispr(self, state: Mapping[str, Any], request: RunRequest) -> Dict[str, Any]:
        config = dict(state.get("config") or {})
        run_config = _merge_run_config(config.get("run_config"), request.parameters)
        guide = config.get("guide")
        if not isinstance(guide, Mapping):
            raise RuntimeError("CRISPR snapshot is missing guide information.")
        sequence = _extract_sequence(state)
        draws = int(run_config.get("draws", 1000))
        seed = _maybe_int(run_config.get("seed"))
        pam_profile = run_config.get("pam_profile") or "SpCas9_NGG"
        pam_softness = float(run_config.get("pam_softness", 0.0))
        priors_profile = run_config.get("priors_profile") or "default_indel"
        priors = resolve_crispr_priors(priors_profile)
        pam_mask = build_crispr_pam_mask(sequence, guide, pam_profile, pam_softness)
        start = perf_counter()
        payload = simulate_cut_repair(
            site_seq=sequence,
            guide=guide,
            priors=priors,
            draws=draws,
            seed=seed,
            emit_sequence=True,
            pam_mask=pam_mask,
        )
        sim_ms = _elapsed_ms(start)
        return _finalize_snapshot(state, request, run_config, payload, perf={"sim_ms": sim_ms})

    def _run_prime(self, state: Mapping[str, Any], request: RunRequest) -> Dict[str, Any]:
        config = dict(state.get("config") or {})
        run_config = _merge_run_config(config.get("run_config"), request.parameters)
        peg = state.get("peg") or config.get("peg")
        if not isinstance(peg, Mapping):
            raise RuntimeError("Prime snapshot is missing pegRNA information.")
        sequence = _extract_sequence(state)
        draws = int(run_config.get("draws", 1000))
        seed = _maybe_int(run_config.get("seed"))
        pam_profile = run_config.get("pam_profile") or "SpCas9_NGG"
        pam_softness = float(run_config.get("pam_softness", 0.0))
        priors_profile = run_config.get("priors_profile") or "default_indel"
        priors = resolve_prime_priors(priors_profile)
        pam_mask = build_prime_pam_mask(sequence, peg, pam_profile, pam_softness)
        start = perf_counter()
        payload = simulate_prime_edit(
            site_seq=sequence,
            peg=peg,
            priors=priors,
            draws=draws,
            seed=seed,
            emit_sequence=True,
            pam_mask=pam_mask,
        )
        sim_ms = _elapsed_ms(start)
        return _finalize_snapshot(state, request, run_config, payload, perf={"sim_ms": sim_ms})

    def _run_pcr(self, state: Mapping[str, Any], request: RunRequest) -> Dict[str, Any]:
        config = dict(state.get("pcr_config") or {})
        run_cfg = _merge_run_config(config, request.parameters)
        try:
            pcr_cfg = PCRConfig(**run_cfg)
        except TypeError as exc:
            raise RuntimeError(f"Invalid PCR configuration: {exc}")
        start = perf_counter()
        result = simulate_pcr(pcr_cfg)
        sim_ms = _elapsed_ms(start)
        state_copy = json.loads(json.dumps(state))
        state_copy["run_kind"] = request.run_kind
        state_copy["run_id"] = request.run_id
        state_copy["pcr_config"] = run_cfg
        state_copy["pcr_result"] = asdict(result)
        state_copy.setdefault("config", {})["run_config"] = run_cfg
        state_copy["config"].setdefault("perf", {})["sim_ms"] = sim_ms
        return {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "state": state_copy,
        }


class OGNEngine:
    def __init__(
        self,
        *,
        client: MockOGNClient | None = None,
        queue: str = "gpu-a100",
        priority: str = "normal",
        labels: Dict[str, str] | None = None,
    ) -> None:
        self._client = client or MockOGNClient()
        self._queue = queue
        self._priority = priority
        self._labels = dict(labels or {})

    def _build_spec(self, request: RunRequest) -> OGNJobSpec:
        job_id = f"ogn-{request.session_id}-{request.run_id}-{uuid.uuid4().hex[:6]}"
        labels = {"run_kind": request.run_kind, "queue": self._queue}
        labels.update(self._labels)
        labels.update(request.labels)
        return OGNJobSpec(
            job_id=job_id,
            run_kind=request.run_kind,
            session_id=request.session_id,
            snapshot_payload=request.snapshot_payload,
            parameters=request.parameters,
            priority=self._priority,
            labels=labels,
        )

    def execute(self, request: RunRequest) -> EngineResult:
        spec = self._build_spec(request)
        self._client.submit(spec)
        status, result = self._client.wait_for_completion(spec.job_id)
        info = {
            "ogn_job_id": spec.job_id,
            "ogn_status": status.status,
            "ogn_engine_version": status.engine_version,
            "ogn_gpu_model": status.gpu_model,
            "ogn_queue": self._queue,
        }
        if status.status != "succeeded":
            info["ogn_error_code"] = status.error_code
            info["ogn_error_message"] = status.error_message
            raise RuntimeError(
                f"OGN job {spec.job_id} failed: {status.error_message or status.error_code or 'unknown error'}"
            )
        metrics = result.metrics if isinstance(result, OGNJobResult) else None
        return EngineResult(snapshot=result.snapshot, metrics=metrics or None, status=status.status, info=info)


def _extract_state(snapshot_payload: Mapping[str, Any]) -> Mapping[str, Any]:
    state = snapshot_payload.get("state") if isinstance(snapshot_payload, Mapping) else None
    if isinstance(state, Mapping):
        return state
    if isinstance(snapshot_payload, Mapping):
        return snapshot_payload
    raise RuntimeError("Snapshot payload is missing 'state'.")


def _extract_sequence(state: Mapping[str, Any]) -> str:
    genome = state.get("genome")
    if not isinstance(genome, Mapping) or not genome:
        raise RuntimeError("Snapshot is missing genome sequences.")
    for sequence in genome.values():
        if isinstance(sequence, str) and sequence:
            return sequence
    raise RuntimeError("Genome sequences are invalid or empty.")


def _merge_run_config(base: Mapping[str, Any] | None, overrides: Mapping[str, Any]) -> Dict[str, Any]:
    config = dict(base or {})
    override_block = overrides.get("run_config") if isinstance(overrides, Mapping) else None
    if isinstance(override_block, Mapping):
        for key, value in override_block.items():
            config[key] = value
    else:
        for key in ("draws", "seed", "pam_profile", "pam_softness", "priors_profile"):
            if key in overrides:
                config[key] = overrides[key]
    return config


def _finalize_snapshot(
    base_state: Mapping[str, Any],
    request: RunRequest,
    run_config: Mapping[str, Any],
    sim_payload: Mapping[str, Any],
    *,
    perf: Mapping[str, Any],
) -> Dict[str, Any]:
    state_copy = json.loads(json.dumps(base_state))
    state_copy["run_kind"] = request.run_kind
    state_copy["run_id"] = request.run_id
    state_copy["outcomes"] = list(sim_payload.get("outcomes", []))
    config = dict(state_copy.get("config") or {})
    config["run_config"] = dict(run_config)
    config["sim_payload"] = sim_payload
    config.setdefault("perf", {}).update(perf)
    state_copy["config"] = config
    return {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "state": state_copy,
    }


def _elapsed_ms(start: float) -> float:
    return max(0.0, (perf_counter() - start) * 1000.0)


def _maybe_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None
