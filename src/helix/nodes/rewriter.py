"""Rewriter nodes mutate graph parameters mid-run."""

from __future__ import annotations

from typing import Mapping

from .base import BaseLiveNode, port


class RewriterNode(BaseLiveNode):
    def __init__(self, name: str) -> None:
        ports = {
            "mutation": port("mutation", "in"),
            "ack": port("ack", "out"),
        }
        super().__init__(name=name, kind="rewriter", ports=ports)
        self.metadata["hash"] = f"rewriter:{name}"

    def step(self, t: float, dt: float, inputs: Mapping[str, float]):  # type: ignore[override]
        mutation = inputs.get("mutation", 0.0)
        self.set_state("last_mutation", mutation)
        return {"ack": mutation}
