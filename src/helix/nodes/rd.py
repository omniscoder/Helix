"""Reaction-diffusion field node."""

from __future__ import annotations

from typing import List, Mapping, Sequence, Tuple

from .base import BaseLiveNode, port


class FieldNode(BaseLiveNode):
    """Simple reaction-diffusion tile that also emits a 2D heatmap."""

    def __init__(
        self,
        name: str,
        diffusion: float = 0.1,
        *,
        shape: Sequence[int] | None = None,
        baseline: float = 0.2,
        stripe_value: float = 1.0,
        stripe_span: float = 0.2,
    ) -> None:
        ports = {
            "signal": port("signal", "in"),
            "control": port("control", "in"),
            "tile": port("tile", "out"),
            "map": port("map", "out", positive=False),
        }
        super().__init__(name=name, kind="field", ports=ports, state={"tile": float(baseline)})
        self.diffusion = float(diffusion)
        self.baseline = float(baseline)
        self.stripe_value = float(stripe_value)
        self.stripe_span = max(0.01, min(0.5, float(stripe_span)))
        width, height = self._normalize_shape(shape)
        self.shape: Tuple[int, int] = (width, height)
        self._grid: List[float] = [self.baseline] * (width * height)
        self._control_value = stripe_value
        self.metadata["hash"] = f"field:{name}:{diffusion}:{width}x{height}"
        self._seed_gradient()

    def _normalize_shape(self, shape: Sequence[int] | None) -> Tuple[int, int]:
        if not shape:
            return (32, 32)
        if len(shape) == 2:
            w, h = int(shape[0]), int(shape[1])
        else:
            w, h = int(shape[0]), int(shape[0])
        return max(w, 8), max(h, 8)

    def _seed_gradient(self) -> None:
        width, height = self.shape
        stripe_cols = max(1, int(width * self.stripe_span))
        for y in range(height):
            for x in range(width):
                idx = y * width + x
                if x < stripe_cols:
                    self._grid[idx] = self.stripe_value
                else:
                    self._grid[idx] = self.baseline

    def step(self, t: float, dt: float, inputs: Mapping[str, float]):  # type: ignore[override]
        signal = float(inputs.get("signal", 0.0))
        control = float(inputs.get("control", self._control_value))
        self._control_value = control
        self._apply_control(control)
        self._diffuse(dt)
        self._consume(signal, dt)
        tile = sum(self._grid) / len(self._grid)
        self.set_state("tile", tile)
        return {"tile": tile, "map": {"shape": list(self.shape), "data": list(self._grid)}}

    def _apply_control(self, control: float) -> None:
        width, height = self.shape
        stripe_cols = max(1, int(width * self.stripe_span))
        for y in range(height):
            for x in range(stripe_cols):
                idx = y * width + x
                self._grid[idx] = max(control, 0.0)

    def _diffuse(self, dt: float) -> None:
        width, height = self.shape
        source = self._grid
        target = source.copy()
        diff = self.diffusion * dt
        for y in range(height):
            for x in range(width):
                idx = y * width + x
                neighbors = [source[idx]]
                if x > 0:
                    neighbors.append(source[idx - 1])
                if x < width - 1:
                    neighbors.append(source[idx + 1])
                if y > 0:
                    neighbors.append(source[idx - width])
                if y < height - 1:
                    neighbors.append(source[idx + width])
                avg = sum(neighbors) / len(neighbors)
                target[idx] = max(0.0, source[idx] + (avg - source[idx]) * diff)
        self._grid = target

    def _consume(self, signal: float, dt: float) -> None:
        if signal == 0.0:
            return
        width, height = self.shape
        delta = signal * 0.01 * dt
        for idx in range(width * height):
            self._grid[idx] = max(0.0, self._grid[idx] - delta)
