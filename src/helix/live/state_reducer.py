"""Delta extraction utilities."""

from __future__ import annotations

from typing import Dict, Mapping, MutableMapping, Optional


def diff_snapshots(
    current: Mapping[str, Mapping[str, float]], previous: Optional[Mapping[str, Mapping[str, float]]]
) -> Dict[str, Dict[str, float]]:
    """Compute delta between snapshots."""

    delta = {"added": {}, "removed": {}, "updated": {}}
    previous = previous or {}
    for key, value in current.items():
        if key not in previous:
            delta["added"][key] = value
        elif value != previous[key]:
            delta["updated"][key] = value
    for key in previous:
        if key not in current:
            delta["removed"][key] = previous[key]
    return delta


class StateReducer:
    """Maintains last snapshot and emits diffs at UI cadence."""

    DEFAULT_SCHEMA = "helix.live.delta.v1"

    def __init__(
        self,
        hz: float = 60.0,
        realtime_queue=None,
        *,
        node_meta: Optional[Mapping[str, Mapping[str, object]]] = None,
        schema: str = DEFAULT_SCHEMA,
    ):
        self.hz = hz
        self._last: Optional[Mapping[str, Mapping[str, float]]] = None
        self.history: list = []
        self._realtime_queue = realtime_queue
        self._emitted_full = False
        self._node_meta: Mapping[str, Mapping[str, object]] = node_meta or {}
        self._schema = schema

    def push(self, snapshot: Mapping[str, Mapping[str, float]], meta: Optional[Mapping[str, object]] = None) -> None:
        delta = diff_snapshots(snapshot, self._last)
        self._last = snapshot
        payload = {
            "time": meta.get("time") if meta else None,
            "delta": delta,
            "runtime": meta.get("runtime") if meta else None,
            "schema": self._schema,
        }
        include_snapshot = not self._emitted_full
        if meta and "include_snapshot" in meta:
            include_snapshot = bool(meta["include_snapshot"])
        if include_snapshot:
            payload["snapshot"] = snapshot
            self._emitted_full = True
            if self._node_meta:
                payload["node_meta"] = self._node_meta
        self.emit(payload)

    def emit(self, payload: Mapping[str, object]) -> None:
        self.history.append(payload)
        if self._realtime_queue:
            self._realtime_queue.put(dict(payload))
