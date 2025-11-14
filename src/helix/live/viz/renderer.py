"""Moderngl helpers for the helix live visualizer."""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import moderngl
import numpy as np

from .feed import LiveFrame


class HeatmapRenderer:
    """Draws a fullscreen quad textured with the field heatmap."""

    def __init__(self, ctx: moderngl.Context):
        self.ctx = ctx
        self.program = ctx.program(
            vertex_shader="""
                #version 330
                in vec2 in_pos;
                in vec2 in_uv;
                out vec2 v_uv;
                void main() {
                    v_uv = in_uv;
                    gl_Position = vec4(in_pos, 0.0, 1.0);
                }
            """,
            fragment_shader="""
                #version 330
                uniform sampler2D data_tex;
                in vec2 v_uv;
                out vec4 fragColor;

                vec3 colormap(float value) {
                    float t = clamp(value, 0.0, 1.0);
                    return vec3(
                        smoothstep(0.0, 1.0, t),
                        smoothstep(0.2, 1.0, t),
                        1.0 - smoothstep(0.0, 1.0, t)
                    );
                }

                void main() {
                    float v = texture(data_tex, v_uv).r;
                    vec3 color = colormap(v);
                    fragColor = vec4(color, 1.0);
                }
            """,
        )
        vertices = np.array(
            [
                -1.0,
                -1.0,
                0.0,
                0.0,
                1.0,
                -1.0,
                1.0,
                0.0,
                -1.0,
                1.0,
                0.0,
                1.0,
                -1.0,
                -1.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
            ],
            dtype="f4",
        )
        self.vao = ctx.simple_vertex_array(self.program, ctx.buffer(vertices.tobytes()), "in_pos", "in_uv")
        self.texture = ctx.texture((1, 1), 1, data=np.zeros((1, 1), dtype="f4").tobytes())
        self.texture.filter = (moderngl.LINEAR, moderngl.LINEAR)

    def update(self, shape: Sequence[int], data: Sequence[float]) -> None:
        width = max(1, int(shape[0]))
        height = max(1, int(shape[1]))
        array = np.array(data, dtype="f4")
        if array.size != width * height:
            array = np.resize(array, width * height)
        array = array.reshape((height, width))
        self.texture = self.ctx.texture((width, height), 1, data=array.tobytes())
        self.texture.filter = (moderngl.LINEAR, moderngl.LINEAR)

    def render(self) -> None:
        self.texture.use(location=0)
        self.program["data_tex"] = 0
        self.vao.render()


class AgentRenderer:
    """Renders agents as instanced points colorized by pERK levels."""

    def __init__(self, ctx: moderngl.Context):
        self.ctx = ctx
        self.program = ctx.program(
            vertex_shader="""
                #version 330
                in vec2 in_pos;
                in vec3 in_color;
                out vec3 v_color;
                void main() {
                    v_color = in_color;
                    gl_Position = vec4(in_pos, 0.0, 1.0);
                    gl_PointSize = 8.0;
                }
            """,
            fragment_shader="""
                #version 330
                in vec3 v_color;
                out vec4 fragColor;
                void main() {
                    vec2 coord = gl_PointCoord - vec2(0.5);
                    if (length(coord) > 0.5) {
                        discard;
                    }
                    fragColor = vec4(v_color, 1.0);
                }
            """,
        )
        self.buffer = ctx.buffer(reserve=0)
        self.vao = ctx.vertex_array(
            self.program,
            [
                (self.buffer, "2f 3f", "in_pos", "in_color"),
            ],
        )
        self.agent_count = 0

    def update(self, positions: Sequence[Sequence[float]], colors: Sequence[float]) -> None:
        if not positions:
            self.agent_count = 0
            return
        positions_arr = np.array(positions, dtype="f4")
        perk_arr = np.array(colors, dtype="f4")
        if perk_arr.shape[0] != positions_arr.shape[0]:
            perk_arr = np.resize(perk_arr, positions_arr.shape[0])
        clip = np.empty_like(positions_arr)
        clip[:, 0] = positions_arr[:, 0] * 2.0 - 1.0
        clip[:, 1] = 1.0 - positions_arr[:, 1] * 2.0
        colors_arr = _colormap(perk_arr)
        packed = np.hstack([clip, colors_arr]).astype("f4")
        data = packed.tobytes()
        self.buffer.orphan(len(data))
        self.buffer.write(data)
        self.agent_count = positions_arr.shape[0]

    def render(self) -> None:
        if self.agent_count <= 0:
            return
        self.ctx.enable(moderngl.BLEND)
        self.ctx.blend_func = moderngl.SRC_ALPHA, moderngl.ONE_MINUS_SRC_ALPHA
        self.vao.render(mode=moderngl.POINTS, vertices=self.agent_count)
        self.ctx.disable(moderngl.BLEND)


def _colormap(values: np.ndarray) -> np.ndarray:
    normalized = np.clip(values, 0.0, 2.0) / 2.0
    r = normalized
    g = np.sqrt(normalized)
    b = 1.0 - normalized
    return np.stack([r, g, b], axis=1)


class SceneRenderer:
    """High-level scene aggregator fed by LiveFrame snapshots."""

    def __init__(self, ctx: moderngl.Context):
        self.ctx = ctx
        self.field = HeatmapRenderer(ctx)
        self.agents = AgentRenderer(ctx)
        self.metrics: Dict[str, float] = {}
        self._field_node: Optional[str] = None
        self._agent_node: Optional[str] = None
        self._perk_node: Optional[str] = None
        self._cell_node: Optional[str] = None

    def apply_frame(self, frame: LiveFrame) -> None:
        snapshot = frame.snapshot
        if frame.node_meta and not self._field_node:
            self._detect_nodes(frame.node_meta)
        if self._field_node and self._field_node in snapshot:
            field_payload = snapshot[self._field_node].get("map")
            if isinstance(field_payload, dict):
                shape = field_payload.get("shape") or (32, 32)
                data = field_payload.get("data") or []
                if data:
                    self.field.update(shape, data)
        if self._agent_node and self._agent_node in snapshot:
            agent_payload = snapshot[self._agent_node]
            positions = agent_payload.get("positions") or []
            perks = agent_payload.get("perk_levels") or []
            self.agents.update(positions, perks)

        self.metrics["time"] = float(frame.time or self.metrics.get("time", 0.0))
        if self._perk_node:
            value = snapshot.get(self._perk_node, {}).get("metric")
            if isinstance(value, (int, float)):
                self.metrics["pERK"] = float(value)
        if self._cell_node:
            cells = snapshot.get(self._cell_node, {}).get("agents")
            if isinstance(cells, (int, float)):
                self.metrics["cell_count"] = float(cells)

    def draw(self) -> None:
        self.field.render()
        self.agents.render()

    def _detect_nodes(self, node_meta: Dict[str, Dict[str, object]]) -> None:
        for name, meta in node_meta.items():
            kind = str(meta.get("kind", "")).lower()
            lowered = name.lower()
            if kind == "field" and not self._field_node:
                self._field_node = name
            elif kind == "abm" and not self._agent_node:
                self._agent_node = name
            elif kind == "observer":
                if ("perk" in lowered or "erk" in lowered) and not self._perk_node:
                    self._perk_node = name
                elif ("cell" in lowered or "agent" in lowered) and not self._cell_node:
                    self._cell_node = name
