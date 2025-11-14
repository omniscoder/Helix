"""Gene regulatory network node."""

from __future__ import annotations

import math
from typing import Mapping

from .base import BaseLiveNode, port


class GRNNode(BaseLiveNode):
    def __init__(self, name: str, hill: float = 2.0, threshold: float = 0.5) -> None:
        ports = {
            "activation": port("activation", "in"),
            "expression": port("expression", "out"),
        }
        super().__init__(name=name, kind="grn", ports=ports, state={"expression": 0.0})
        self.hill = hill
        self.threshold = threshold
        self.metadata["hash"] = f"grn:{name}:{hill}:{threshold}"

    def step(self, t: float, dt: float, inputs: Mapping[str, float]):  # type: ignore[override]
        activation = inputs.get("activation", 0.0)
        numerator = activation ** self.hill
        denom = numerator + self.threshold ** self.hill
        expression = numerator / denom if denom else 0.0
        self.set_state("expression", expression)
        return {"expression": expression}
