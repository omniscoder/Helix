from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Mapping, Tuple


def _iso_now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass(frozen=True)
class OGNJobSpec:
    job_id: str
    run_kind: str
    session_id: str
    snapshot_payload: Mapping[str, Any]
    parameters: Mapping[str, Any]
    priority: str = "normal"
    created_at: str = field(default_factory=_iso_now)
    labels: Dict[str, str] = field(default_factory=dict)


@dataclass
class OGNJobStatus:
    job_id: str
    status: str
    engine_version: str | None = None
    gpu_model: str | None = None
    started_at: str | None = None
    finished_at: str | None = None
    progress: float | None = None
    error_code: str | None = None
    error_message: str | None = None


@dataclass
class OGNJobResult:
    job_id: str
    snapshot: Mapping[str, Any]
    metrics: Mapping[str, Any]
    artifacts: Dict[str, str]


class MockOGNClient:
    def __init__(self, root: Path | None = None) -> None:
        self.root = root or Path("mock_ogn_runs")
        self.root.mkdir(parents=True, exist_ok=True)
        self._statuses: Dict[str, OGNJobStatus] = {}
        self._results: Dict[str, OGNJobResult] = {}

    def submit(self, spec: OGNJobSpec) -> OGNJobStatus:
        status = OGNJobStatus(job_id=spec.job_id, status="queued")
        self._statuses[spec.job_id] = status
        self._execute(spec)
        return self._statuses[spec.job_id]

    def _execute(self, spec: OGNJobSpec) -> None:
        status = self._statuses[spec.job_id]
        status.status = "running"
        status.started_at = _iso_now()
        job_dir = self.root / spec.job_id
        job_dir.mkdir(parents=True, exist_ok=True)
        snapshot_path = job_dir / "snapshot.json"
        snapshot_text = json.dumps(spec.snapshot_payload, indent=2) + "\n"
        snapshot_path.write_text(snapshot_text, encoding="utf-8")
        status.finished_at = _iso_now()
        status.status = "succeeded"
        status.engine_version = "mock-ogn-0"
        status.gpu_model = "mock-gpu"
        result = OGNJobResult(
            job_id=spec.job_id,
            snapshot=json.loads(json.dumps(spec.snapshot_payload)),
            metrics={},
            artifacts={"snapshot": str(snapshot_path)},
        )
        self._results[spec.job_id] = result
        (job_dir / "status.json").write_text(json.dumps(status.__dict__, default=str), encoding="utf-8")

    def get_status(self, job_id: str) -> OGNJobStatus:
        return self._statuses[job_id]

    def wait_for_completion(self, job_id: str) -> Tuple[OGNJobStatus, OGNJobResult]:
        status = self._statuses[job_id]
        result = self._results[job_id]
        return status, result
