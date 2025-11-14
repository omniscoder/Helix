"""Observer nodes emit metrics/losses."""

from __future__ import annotations

from typing import Mapping

from .base import BaseLiveNode, port


class ObserverNode(BaseLiveNode):
    def __init__(self, name: str) -> None:
        ports = {
            "signal": port("signal", "in"),
            "metric": port("metric", "out"),
        }
        super().__init__(name=name, kind="observer", ports=ports)
        self.metadata["hash"] = f"observer:{name}"

    def step(self, t: float, dt: float, inputs: Mapping[str, float]):  # type: ignore[override]
        metric = inputs.get("signal", 0.0)
        self.set_state("metric", metric)
        return {"metric": metric}
