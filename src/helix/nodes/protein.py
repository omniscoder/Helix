"""Protein nodes convert MD/ΔΔG into rate multipliers."""

from __future__ import annotations

import math
from typing import Mapping

from .base import BaseLiveNode, port


class ProteinNode(BaseLiveNode):
    def __init__(self, name: str, base_multiplier: float = 1.0) -> None:
        ports = {
            "energy": port("energy", "in", unit="kcal/mol"),
            "rate_multiplier": port("rate_multiplier", "out"),
        }
        super().__init__(name=name, kind="protein", ports=ports)
        self.metadata["hash"] = f"protein:{name}:{base_multiplier}"
        self.base_multiplier = base_multiplier

    def step(self, t: float, dt: float, inputs: Mapping[str, float]):  # type: ignore[override]
        energy = inputs.get("energy", 0.0)
        multiplier = self.base_multiplier * math.exp(-energy / 10.0)
        self.set_state("multiplier", multiplier)
        return {"rate_multiplier": multiplier}
