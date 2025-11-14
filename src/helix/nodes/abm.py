"""Agent-based model node."""

from __future__ import annotations

from typing import List, Mapping, Sequence, Tuple

from .base import BaseLiveNode, port


class ABMNode(BaseLiveNode):
    """Toy ABM node that exposes per-agent positions and pERK-like scores."""

    def __init__(self, name: str, agents: int = 100, *, layout: Sequence[int] | None = None) -> None:
        ports = {
            "field_signal": port("field_signal", "in"),
            "field_map": port("field_map", "in", positive=False),
            "agent_delta": port("agent_delta", "out"),
            "agents": port("agents", "out"),
            "positions": port("positions", "out", positive=False),
            "perk_levels": port("perk_levels", "out", positive=False),
        }
        super().__init__(name=name, kind="abm", ports=ports, state={"agents": float(agents)})
        self.metadata["hash"] = f"abm:{name}:{agents}"
        self._layout = self._normalize_layout(layout)
        self._positions: List[Tuple[float, float]] = self._seed_positions(max(agents, self._layout[0] * self._layout[1]))
        self._perk_levels: List[float] = [0.1 for _ in self._positions]

    def _normalize_layout(self, layout: Sequence[int] | None) -> Tuple[int, int]:
        if not layout:
            return (32, 32)
        if len(layout) == 1:
            return (int(layout[0]), int(layout[0]))
        return (max(1, int(layout[0])), max(1, int(layout[1])))

    def _seed_positions(self, count: int) -> List[Tuple[float, float]]:
        cols, rows = self._layout
        total = max(1, cols * rows)
        positions: List[Tuple[float, float]] = []
        for idx in range(count):
            x = (idx % cols) / max(cols - 1, 1)
            y = (idx // cols) / max(rows - 1, 1)
            positions.append((x, y))
            if len(positions) >= total:
                break
        return positions

    def step(self, t: float, dt: float, inputs: Mapping[str, float]):  # type: ignore[override]
        agents = self.get_state("agents")
        global_signal = float(inputs.get("field_signal", 0.0))
        field_map = inputs.get("field_map")
        local_drive = self._sample_field(field_map)
        perk_levels = []
        for idx, perk in enumerate(self._perk_levels):
            drive = local_drive[idx] if idx < len(local_drive) else global_signal
            target = 0.5 * global_signal + 0.5 * drive
            perk = perk + (target - perk) * 0.5 * dt
            perk = max(0.0, min(2.0, perk))
            perk_levels.append(perk)
        self._perk_levels = perk_levels
        growth_rate = 0.05 + 0.1 * global_signal
        delta = growth_rate * dt
        agents = max(0.0, agents + delta * agents)
        self.set_state("agents", agents)
        return {
            "agent_delta": delta,
            "agents": agents,
            "positions": [list(pos) for pos in self._positions],
            "perk_levels": perk_levels,
        }

    def _sample_field(self, field_map) -> List[float]:
        if not isinstance(field_map, dict):
            return [0.0 for _ in self._positions]
        data = field_map.get("data")
        shape = field_map.get("shape") or [len(data or []), 1]
        if not isinstance(data, list):
            return [0.0 for _ in self._positions]
        if not isinstance(shape, list) or len(shape) < 2:
            width = len(data)
            height = 1
        else:
            width = max(1, int(shape[0]))
            height = max(1, int(shape[1]))
        samples: List[float] = []
        for x, y in self._positions:
            ix = min(width - 1, max(0, int(x * (width - 1))))
            iy = min(height - 1, max(0, int(y * (height - 1))))
            idx = iy * width + ix
            if 0 <= idx < len(data):
                samples.append(float(data[idx]))
            else:
                samples.append(0.0)
        return samples
