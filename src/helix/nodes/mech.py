"""Tissue mechanics node."""

from __future__ import annotations

from typing import Mapping

from .base import BaseLiveNode, port


class MechNode(BaseLiveNode):
    def __init__(self, name: str, stiffness: float = 1.0) -> None:
        ports = {
            "stress_in": port("stress_in", "in", unit="Pa"),
            "stress_out": port("stress_out", "out", unit="Pa"),
        }
        super().__init__(name=name, kind="mech", ports=ports, state={"stress": 0.0})
        self.stiffness = stiffness
        self.metadata["hash"] = f"mech:{name}:{stiffness}"

    def step(self, t: float, dt: float, inputs: Mapping[str, float]):  # type: ignore[override]
        stress = inputs.get("stress_in", 0.0) * self.stiffness
        self.set_state("stress", stress)
        return {"stress_out": stress}
