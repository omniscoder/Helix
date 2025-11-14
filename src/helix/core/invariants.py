"""Invariant checks for Helix graphs."""

from __future__ import annotations

from typing import Mapping

from .graph import GraphIR, Port


class InvariantViolation(RuntimeError):
    """Base error for invariant violations."""


class UnitMismatch(InvariantViolation):
    """Raised when edge ports use incompatible units."""


class NegativeSignal(InvariantViolation):
    """Raised when a positive-only port emits a negative value."""


def validate_edge_units(graph: GraphIR) -> None:
    """Ensure every connection maintains consistent units."""

    for edge in graph.edges:
        src_node, src_port = edge.source
        dst_node, dst_port = edge.target
        src = graph.node(src_node).ports[src_port]
        dst = graph.node(dst_node).ports[dst_port]
        if src.unit != dst.unit:
            raise UnitMismatch(
                f"Unit mismatch {src_node}.{src_port} ({src.unit}) -> {dst_node}.{dst_port} ({dst.unit})"
            )


def validate_snapshot(graph: GraphIR, snapshot: Mapping[str, Mapping[str, float]]) -> None:
    """Check positivity invariants for a runtime snapshot."""

    for node in graph.iter_nodes():
        if node.name not in snapshot:
            continue
        values = snapshot[node.name]
        for port in node.ports.values():
            if port.direction != "out" or not port.positive:
                continue
            value = values.get(port.name)
            if value is None:
                continue
            if value < 0:
                raise NegativeSignal(f"{node.name}.{port.name} emitted negative value {value}")


def assert_invariants(graph: GraphIR, snapshot: Mapping[str, Mapping[str, float]]) -> None:
    """Run all invariant checks."""

    validate_edge_units(graph)
    validate_snapshot(graph, snapshot)
