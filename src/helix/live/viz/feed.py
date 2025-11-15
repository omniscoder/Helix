"""Streaming/replay helpers for the helix live viz client."""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

from ..realtime import RealtimeClient


@dataclass
class LiveFrame:
    time: Optional[float]
    snapshot: Dict[str, Dict[str, Any]]
    runtime: Dict[str, Any]
    node_meta: Dict[str, Dict[str, Any]]
    schema: str
    payload: Dict[str, Any]


class _SnapshotState:
    def __init__(self) -> None:
        self.snapshot: Dict[str, Dict[str, Any]] = {}
        self.node_meta: Dict[str, Dict[str, Any]] = {}

    def apply(self, payload: Mapping[str, Any]) -> LiveFrame:
        if payload.get("node_meta"):
            self.node_meta = {k: dict(v) for k, v in payload["node_meta"].items()}
        if payload.get("snapshot"):
            self.snapshot = {node: dict(values) for node, values in payload["snapshot"].items()}
        else:
            self.snapshot = _apply_delta(self.snapshot, payload.get("delta") or {})
        runtime = payload.get("runtime") or {}
        schema = payload.get("schema") or payload.get("kind") or ""
        return LiveFrame(
            time=payload.get("time"),
            snapshot=self.snapshot,
            runtime=dict(runtime),
            node_meta=self.node_meta,
            schema=schema,
            payload=dict(payload),
        )


def _apply_delta(snapshot: Mapping[str, Dict[str, Any]], delta: Mapping[str, Mapping[str, Any]]):
    current = {node: dict(values) for node, values in snapshot.items()}
    for node, values in delta.get("added", {}).items():
        current[node] = dict(values)
    for node, values in delta.get("updated", {}).items():
        current[node] = dict(values)
    for node in delta.get("removed", {}):
        current.pop(node, None)
    return current


class BaseFeed:
    def poll(self) -> Optional[LiveFrame]:  # pragma: no cover - interface
        raise NotImplementedError

    def send_command(self, command: Mapping[str, Any]) -> None:  # pragma: no cover - interface
        raise NotImplementedError

    def close(self) -> None:  # pragma: no cover - interface
        raise NotImplementedError


class RealtimeFeed(BaseFeed):
    """Fronts a RealtimeClient and emits normalized frames."""

    def __init__(self, endpoint: str):
        self._client = RealtimeClient(endpoint)
        self._state = _SnapshotState()

    def poll(self) -> Optional[LiveFrame]:
        payload = self._client.poll()
        if not payload:
            return None
        return self._state.apply(payload)

    def send_command(self, command: Mapping[str, Any]) -> None:
        self._client.send_command(dict(command))

    def close(self) -> None:
        self._client.close()


class BundleFeed(BaseFeed):
    """Replays bundles/deltas.json payloads for offline debugging."""

    def __init__(self, bundle_path: Path, *, loop: bool = True, interval: float = 0.05):
        deltas_path = bundle_path / "deltas.json"
        if not deltas_path.exists():
            raise FileNotFoundError(f"{deltas_path} was not found.")
        data = json.loads(deltas_path.read_text(encoding="utf-8"))
        self._frames: List[Mapping[str, Any]] = list(data)
        self._loop = loop
        self._interval = interval
        self._state = _SnapshotState()
        self._next_idx = 0
        self._next_deadline = time.time()
        self._paused = False
        self._step_pending = False

    def poll(self) -> Optional[LiveFrame]:
        if not self._frames:
            return None
        if self._paused and not self._step_pending:
            return None
        now = time.time()
        if not self._paused and now < self._next_deadline:
            return None
        payload = self._frames[self._next_idx]
        self._advance_index()
        if self._paused and self._step_pending:
            self._step_pending = False
        self._next_deadline = now + self._interval
        return self._state.apply(payload)

    def send_command(self, command: Mapping[str, Any]) -> None:
        # Offline replays ignore commands to keep semantics clear.
        return None

    def close(self) -> None:
        return None

    def pause(self) -> None:
        self._paused = True

    def resume(self) -> None:
        self._paused = False
        self._next_deadline = time.time()

    def toggle(self) -> None:
        if self._paused:
            self.resume()
        else:
            self.pause()

    def step_once(self) -> None:
        self._paused = True
        self._step_pending = True

    @property
    def paused(self) -> bool:
        return self._paused

    def _advance_index(self) -> None:
        self._next_idx += 1
        if self._loop and self._next_idx >= len(self._frames):
            self._next_idx = 0
        elif self._next_idx >= len(self._frames):
            self._next_idx = len(self._frames) - 1
