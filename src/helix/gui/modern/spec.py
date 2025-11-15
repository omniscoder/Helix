"""Data classes describing CRISPR/prime editing visualization specs."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


class EditEventType(str, Enum):
    """Canonical visualization events mapped from the CRISPR simulator."""

    RECOGNITION = "recognition"
    NICK_PRIMARY = "nick_primary"
    NICK_SECONDARY = "nick_secondary"
    CUT = "cut"
    PRIME_INIT = "prime_init"
    RT_SYNTHESIS = "rt_synthesis"
    FLAP_RESOLUTION = "flap_resolution"
    REPAIR_COMPLETE = "repair_complete"
    CUSTOM = "custom"

    @classmethod
    def from_token(cls, value: str | "EditEventType") -> "EditEventType":
        if isinstance(value, cls):
            return value
        token = str(value).lower().strip()
        for member in cls:
            if member.value == token:
                return member
        return cls.CUSTOM


@dataclass(slots=True)
class EditEvent:
    """Single timeline event used to drive helix shading and overlays."""

    t: float
    type: EditEventType
    index: int | None = None
    start: int | None = None
    length: int | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_payload(cls, payload: Mapping[str, Any]) -> "EditEvent":
        try:
            t = float(payload.get("t") or payload.get("time"))
        except Exception as exc:  # pragma: no cover - defensive
            raise ValueError(f"Timeline event missing 't': {payload}") from exc
        return cls(
            t=t,
            type=EditEventType.from_token(payload.get("type", "custom")),
            index=payload.get("index"),
            start=payload.get("start") or payload.get("begin"),
            length=payload.get("length") or payload.get("size"),
            metadata=dict(payload.get("metadata") or {}),
        )


@dataclass(slots=True)
class EditVisualizationSpec:
    """User-friendly bundle describing one visualized edit sequence."""

    sequence: str
    pam_index: int
    guide_range: tuple[int, int]
    edit_type: str
    events: list[EditEvent]
    description: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.sequence = self.sequence.upper()
        g0, g1 = self.guide_range
        if g0 < 0 or g1 < g0 or g1 > len(self.sequence):
            raise ValueError(
                f"guide_range {self.guide_range} invalid for sequence of length {len(self.sequence)}"
            )
        if self.pam_index < 0 or self.pam_index >= len(self.sequence):
            raise ValueError(
                f"pam_index {self.pam_index} invalid for sequence of length {len(self.sequence)}"
            )
        self.events.sort(key=lambda e: e.t)

    @property
    def duration(self) -> float:
        return self.events[-1].t if self.events else 0.0

    @classmethod
    def from_payload(cls, payload: Mapping[str, Any]) -> "EditVisualizationSpec":
        sequence = str(payload.get("sequence", "")).strip()
        if not sequence:
            raise ValueError("Visualization spec missing 'sequence'")
        pam_token = payload.get("pam_index")
        if pam_token is None:
            pam_token = payload.get("pam")
        if pam_token is None:
            raise ValueError("Visualization spec missing 'pam_index'")
        pam_index = int(pam_token)
        guide_range_seq = payload.get("guide_range") or payload.get("guide") or [0, 0]
        if isinstance(guide_range_seq, Sequence):
            try:
                g0 = int(guide_range_seq[0])
                g1 = int(guide_range_seq[1])
            except Exception as exc:  # pragma: no cover - defensive
                raise ValueError("guide_range must be a 2-item sequence") from exc
            guide_range = (min(g0, g1), max(g0, g1))
        else:
            raise ValueError("guide_range must be provided as a 2-element sequence")
        edit_type = str(payload.get("edit_type", "unknown"))
        description = payload.get("description")
        metadata = dict(payload.get("metadata") or {})
        events_payload = payload.get("edit_events") or payload.get("events") or []
        if not isinstance(events_payload, Iterable):
            raise ValueError("edit_events must be a list-like structure")
        events = [EditEvent.from_payload(entry) for entry in events_payload]
        return cls(
            sequence=sequence,
            pam_index=pam_index,
            guide_range=guide_range,
            edit_type=edit_type,
            events=events,
            description=description,
            metadata=metadata,
        )

    def to_payload(self) -> dict[str, Any]:
        return {
            "sequence": self.sequence,
            "pam_index": self.pam_index,
            "guide_range": list(self.guide_range),
            "edit_type": self.edit_type,
            "description": self.description,
            "metadata": dict(self.metadata),
            "edit_events": [
                {
                    "t": event.t,
                    "type": event.type.value,
                    "index": event.index,
                    "start": event.start,
                    "length": event.length,
                    "metadata": dict(event.metadata),
                }
                for event in self.events
            ],
        }


def load_viz_spec(source: str | Path | Mapping[str, Any] | EditVisualizationSpec) -> EditVisualizationSpec:
    """Load an :class:`EditVisualizationSpec` from JSON text, a path, or a mapping."""

    if isinstance(source, EditVisualizationSpec):
        return source

    payload: Mapping[str, Any]
    if isinstance(source, Mapping):
        payload = dict(source)
    else:
        candidate = Path(str(source)).expanduser()
        text: str
        if candidate.exists():
            text = candidate.read_text(encoding="utf-8")
        else:
            text = str(source)
        try:
            payload = json.loads(text)
        except json.JSONDecodeError as exc:  # pragma: no cover - defensive
            raise ValueError(f"Failed to parse visualization spec: {exc}") from exc
    return EditVisualizationSpec.from_payload(payload)


def load_viz_specs(
    source: str | Path | Mapping[str, Any] | Sequence[Mapping[str, Any]]
) -> list[EditVisualizationSpec]:
    """Vectorized loader that understands bundles of visualization specs."""

    if isinstance(source, Sequence) and not isinstance(source, (str, bytes, bytearray)):
        if not source:
            return []
        first = source[0]
        if isinstance(first, Mapping):
            return [EditVisualizationSpec.from_payload(payload) for payload in source]  # type: ignore[arg-type]

    if isinstance(source, Path) or (isinstance(source, str) and Path(source).expanduser().exists()):
        path = Path(source).expanduser()
        text = path.read_text(encoding="utf-8")
    elif isinstance(source, str):
        text = source
    elif isinstance(source, Mapping):
        return [EditVisualizationSpec.from_payload(source)]
    else:
        raise TypeError(f"Unsupported viz spec source: {type(source)}")

    data = json.loads(text)
    if isinstance(data, list):
        payloads = data
    elif isinstance(data, Mapping):
        if "edits" in data and isinstance(data["edits"], list):
            payloads = data["edits"]  # type: ignore[assignment]
        else:
            payloads = [data]
    else:
        raise ValueError("Visualization JSON must be an object or list")
    return [EditVisualizationSpec.from_payload(payload) for payload in payloads]
