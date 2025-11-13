"""Edit rule (\"physics\") registry for the DAG engine."""
from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Dict, Iterable, Protocol, Tuple

from .dag import EditNode
from .events import EditEvent

if TYPE_CHECKING:  # pragma: no cover - typing only
    from .simulate import SimulationContext


class EditRule(Protocol):
    """Minimal interface for an edit rule plug-in."""

    name: str

    def propose(
        self,
        node: EditNode,
        context: "SimulationContext",
    ) -> Iterable[Tuple[EditEvent, float, Dict[str, Any]]]:
        """Yield (event, log_prob_delta, metadata) tuples."""


_RULES: Dict[str, EditRule] = {}


def register_rule(rule: EditRule) -> None:
    """Register an EditRule implementation."""
    _RULES[rule.name] = rule


def get_rule(name: str) -> EditRule:
    """Return a registered rule by name."""
    if name not in _RULES:
        raise KeyError(f"Unknown EditRule '{name}'.")
    return _RULES[name]


def edit_rule(name: str) -> Callable[[Callable[..., Iterable[Tuple[EditEvent, float, Dict[str, Any]]]]], Callable[..., Iterable[Tuple[EditEvent, float, Dict[str, Any]]]]]:
    """
    Decorator to register a simple function as an EditRule.

    The decorated function receives (node, context) and returns an iterable of proposals.
    """

    def decorator(fn: Callable[[EditNode, "SimulationContext"], Iterable[Tuple[EditEvent, float, Dict[str, Any]]]]):
        class _FunctionRule:
            def __init__(self, func: Callable[[EditNode, "SimulationContext"], Iterable[Tuple[EditEvent, float, Dict[str, Any]]]]):
                self.name = name
                self._func = func

            def propose(
                self,
                node: EditNode,
                context: "SimulationContext",
            ) -> Iterable[Tuple[EditEvent, float, Dict[str, Any]]]:
                return self._func(node, context)

        register_rule(_FunctionRule(fn))
        return fn

    return decorator

# Backwards compatibility alias
rule = edit_rule
