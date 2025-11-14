"""Coupler nodes lift/restrict data between islands."""

from __future__ import annotations

from typing import Mapping

from .base import BaseLiveNode, port


class CouplerNode(BaseLiveNode):
    def __init__(self, name: str) -> None:
        ports = {
            "input": port("input", "in"),
            "output": port("output", "out"),
        }
        super().__init__(name=name, kind="coupler", ports=ports)
        self.metadata["hash"] = f"coupler:{name}"

    def step(self, t: float, dt: float, inputs: Mapping[str, float]):  # type: ignore[override]
        return {"output": inputs.get("input", 0.0)}
