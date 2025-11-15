"""Moderngl helpers for the helix live visualizer."""

from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import base64
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
                uniform float u_min;
                uniform float u_max;
                in vec2 v_uv;
                out vec4 fragColor;

                vec3 colormap(float value) {
                    float t = clamp((value - u_min) / (u_max - u_min), 0.0, 1.0);
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
        self.texture: Optional[moderngl.Texture] = None
        self.size = (0, 0)
        self.value_range = (0.0, 1.0)

    def update(self, shape: Sequence[int], data: Sequence[float]) -> None:
        width = max(1, int(shape[0]))
        height = max(1, int(shape[1]))
        array = np.array(data, dtype="float32")
        if array.size != width * height:
            array = np.resize(array, width * height)
        array = array.reshape((height, width))
        self._ensure_texture(width, height)
        self.texture.write(array.tobytes())

    def update_from_tiles(self, spec: Mapping[str, Any]) -> None:
        nx = int(spec.get("nx") or self.size[0] or 1)
        ny = int(spec.get("ny") or self.size[1] or 1)
        self._ensure_texture(nx, ny)
        for tile in spec.get("tiles") or []:
            data_b64 = tile.get("data")
            if not data_b64:
                continue
            raw = base64.b64decode(data_b64)
            tile_nx = int(tile.get("nx") or nx)
            tile_ny = int(tile.get("ny") or ny)
            arr = np.frombuffer(raw, dtype="<f4")
            if arr.size != tile_nx * tile_ny:
                arr = np.resize(arr, tile_nx * tile_ny)
            arr = arr.reshape((tile_ny, tile_nx))
            x0 = int(tile.get("x0", 0))
            y0 = int(tile.get("y0", 0))
            self.texture.write(arr.tobytes(), viewport=(x0, y0, tile_nx, tile_ny))
        min_val = float(spec.get("min", self.value_range[0]))
        max_val = float(spec.get("max", self.value_range[1]))
        if max_val <= min_val:
            max_val = min_val + 1e-6
        self.value_range = (min_val, max_val)

    def render(self) -> None:
        if not self.texture:
            return
        self.texture.use(location=0)
        self.program["data_tex"] = 0
        self.program["u_min"].value = self.value_range[0]
        self.program["u_max"].value = self.value_range[1]
        self.vao.render()

    def _ensure_texture(self, width: int, height: int) -> None:
        if self.texture and self.size == (width, height):
            return
        self.texture = self.ctx.texture(
            (width, height),
            1,
            data=np.zeros((height, width), dtype="float32").tobytes(),
            dtype="f4",
        )
        self.texture.filter = (moderngl.LINEAR, moderngl.LINEAR)
        self.size = (width, height)


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
        self.buffer = ctx.buffer(reserve=4)
        self._build_vao()
        self.agent_count = 0

    def _build_vao(self) -> None:
        self.vao = self.ctx.vertex_array(
            self.program,
            [
                (self.buffer, "2f 3f", "in_pos", "in_color"),
            ],
        )

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
        required = packed.nbytes
        if required > self.buffer.size:
            self.buffer.release()
            self.buffer = self.ctx.buffer(reserve=required)
            self._build_vao()
        self.buffer.write(packed.tobytes())
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
        self._live_agents: Dict[int, Dict[str, float]] = {}
        self._field_size = (1.0, 1.0)

    def apply_frame(self, frame: LiveFrame) -> None:
        if frame.payload.get("kind") == "live_delta":
            self._apply_live_delta(frame.payload)
            return
        self._apply_snapshot_frame(frame)

    def _apply_snapshot_frame(self, frame: LiveFrame) -> None:
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
        if frame.time is not None:
            self.metrics["time"] = float(frame.time)

    def _apply_live_delta(self, payload: Mapping[str, Any]) -> None:
        metrics = payload.get("metrics") or {}
        for key, value in metrics.items():
            if isinstance(value, (int, float)):
                self.metrics[key] = float(value)
        if "t" in payload:
            self.metrics["time"] = float(payload["t"])
        fields = payload.get("fields") or {}
        field_spec = None
        for field_entry in fields.values():
            field_spec = field_entry
            break
        if field_spec:
            nx = int(field_spec.get("nx") or self._field_size[0] or 1)
            ny = int(field_spec.get("ny") or self._field_size[1] or 1)
            self._field_size = (float(nx), float(ny))
            self.field.update_from_tiles(field_spec)
        agents_spec = payload.get("agents")
        if agents_spec:
            self._update_live_agents(agents_spec, self._field_size)
        elif not field_spec:
            self.agents.update([], [])

    def _update_live_agents(self, spec: Mapping[str, Any], bounds: Tuple[float, float]) -> None:
        width = max(bounds[0], 1.0)
        height = max(bounds[1], 1.0)
        if spec.get("full"):
            self._live_agents.clear()
        for entry in spec.get("add", []):
            if "id" not in entry:
                continue
            agent_id = int(entry["id"])
            self._live_agents[agent_id] = {
                "x": float(entry.get("x", 0.0)),
                "y": float(entry.get("y", 0.0)),
                "perk": float(entry.get("perk", 0.0)),
            }
        for entry in spec.get("update", []):
            agent_id = entry.get("id")
            if agent_id is None:
                continue
            agent = self._live_agents.setdefault(int(agent_id), {})
            if "x" in entry:
                agent["x"] = float(entry["x"])
            if "y" in entry:
                agent["y"] = float(entry["y"])
            if "perk" in entry:
                agent["perk"] = float(entry["perk"])
        for removed in spec.get("remove", []):
            self._live_agents.pop(int(removed), None)

        if not self._live_agents:
            self.agents.update([], [])
            return
        positions = []
        perks = []
        for agent in self._live_agents.values():
            norm_x = float(agent.get("x", 0.0)) / width
            norm_y = float(agent.get("y", 0.0)) / height
            positions.append([norm_x, norm_y])
            perks.append(float(agent.get("perk", 0.0)))
        self.agents.update(positions, perks)

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
