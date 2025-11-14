"""Helpers for building runtime nodes."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Mapping

from ..core.graph import GraphNode, Port


def port(name: str, direction: str, unit: str = "arb", positive: bool = True) -> Port:
    return Port(name=name, direction=direction, unit=unit, positive=positive)


@dataclass
class BaseLiveNode(GraphNode):
    """Convenience base with safe state helpers."""

    state: Dict[str, Any] = field(default_factory=dict)

    def get_state(self, key: str, default: float = 0.0) -> float:
        return float(self.state.get(key, default))

    def set_state(self, key: str, value: float) -> None:
        self.state[key] = float(value)

    def step(self, t: float, dt: float, inputs: Mapping[str, float]) -> Dict[str, float]:  # type: ignore[override]
        return super().step(t, dt, inputs)
