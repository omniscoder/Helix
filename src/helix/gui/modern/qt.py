"""PyQt shell that wraps the ModernGL helix core."""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Sequence, Any, Mapping

import moderngl
import numpy as np
from OpenGL import GL
from PySide6.QtCore import Qt, QTimer, Signal, QPointF, QThread
from PySide6.QtGui import QSurfaceFormat, QWindow
from PySide6.QtOpenGL import QOpenGLWindow
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QPushButton,
    QTextBrowser,
    QSlider,
    QVBoxLayout,
    QWidget,
)

from ..theme import apply_helix_theme
from ...config import native_backend_available, resolve_crispr_backend
from ...engine.benchmark import BenchmarkConfig, run_benchmark
from ...crispr.physics import CRISPR_SCORING_VERSION
from ...prime.physics import PRIME_SCORING_VERSION
from .engine import AnimationPhase, HelixVizEngine
from .spec import EditVisualizationSpec, load_viz_spec, load_viz_specs


class BenchmarkWorker(QThread):
    completed = Signal(dict)
    failed = Signal(str)

    def __init__(self, config: BenchmarkConfig, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._config = config

    def run(self) -> None:
        try:
            payload = run_benchmark(self._config)
        except Exception as exc:  # pragma: no cover - UI helper
            self.failed.emit(str(exc))
            return
        self.completed.emit(payload)


class HelixModernWidget(QOpenGLWindow):
    """QOpenGLWidget that pushes geometry buffers into ModernGL programs."""

    timeUpdated = Signal(float)
    phaseChanged = Signal(str)
    frameDrawn = Signal()

    def __init__(self, parent: QWindow | None = None) -> None:
        super().__init__(parent=parent)
        fmt = QSurfaceFormat()
        fmt.setRenderableType(QSurfaceFormat.OpenGL)
        fmt.setVersion(3, 3)
        fmt.setProfile(QSurfaceFormat.CoreProfile)
        fmt.setDepthBufferSize(24)
        fmt.setStencilBufferSize(8)
        fmt.setSwapBehavior(QSurfaceFormat.DoubleBuffer)
        fmt.setSamples(4)
        fmt.setSwapInterval(1)
        self.setFormat(fmt)
        self.engine: HelixVizEngine | None = None
        self.ctx: moderngl.Context | None = None
        self._rail_program: moderngl.Program | None = None
        self._arc_program: moderngl.Program | None = None
        self._rail_vbo: moderngl.Buffer | None = None
        self._vao: moderngl.VertexArray | None = None
        self._arc_vbo: moderngl.Buffer | None = None
        self._arc_vao: moderngl.VertexArray | None = None
        self._arc_ranges: list[tuple[int, int]] = []
        self._arc_firsts: np.ndarray = np.zeros((0,), dtype="i4")
        self._arc_counts: np.ndarray = np.zeros((0,), dtype="i4")
        self._bars_prog: moderngl.Program | None = None
        self._bars_unit_vbo: moderngl.Buffer | None = None
        self._bars_inst0_vbo: moderngl.Buffer | None = None
        self._bars_inst1_vbo: moderngl.Buffer | None = None
        self._bars_vao: moderngl.VertexArray | None = None
        self._heat_prog: moderngl.Program | None = None
        self._heat_vbo: moderngl.Buffer | None = None
        self._heat_vao: moderngl.VertexArray | None = None
        self._heat_vertex_count = 0
        self._flat_prog: moderngl.Program | None = None
        self._slice_vbo: moderngl.Buffer | None = None
        self._slice_vao: moderngl.VertexArray | None = None
        self._workflow_buffers: dict[str, np.ndarray] | None = None
        self._workflow_instance_count = 0
        self._workflow_rail_vbo: moderngl.Buffer | None = None
        self._workflow_tick_vbo: moderngl.Buffer | None = None
        self._mode = "3d"
        self._playing = True
        self._show_template = True
        self._show_flap = True
        self._timer = QTimer(self)
        self._timer.timeout.connect(self._tick)
        self._timer.start(16)
        self._has_drawn_frame = False
        self._last_fbo = None
        self._custom_camera: dict[str, Any] | None = None
        self._spec: EditVisualizationSpec | None = None

        # Camera interaction state
        self._cam_target = np.array([0.0, 0.0, 0.0], dtype="f4")
        self._cam_distance = 6.0
        self._cam_yaw = 0.0
        self._cam_pitch = 0.3
        self._cam_sensitivity = 0.005
        self._zoom_sensitivity = 0.2
        self._pan_sensitivity = 0.002
        self._last_mouse_pos: QPointF | None = None
        self._drag_mode: str | None = None
        if hasattr(self, "setMouseTracking"):
            self.setMouseTracking(True)  # type: ignore[attr-defined]
        if hasattr(self, "setFocusPolicy"):
            self.setFocusPolicy(Qt.StrongFocus)  # type: ignore[attr-defined]

    # Public API -------------------------------------------------------
    def set_spec(self, spec: EditVisualizationSpec) -> None:
        print(
            "[HelixModernWidget] set_spec called:",
            getattr(spec, "edit_type", "<unknown>"),
            "sequence len:",
            len(spec.sequence),
            flush=True,
        )
        self._spec = spec
        self._mode = "3d"
        self._workflow_buffers = None
        self.engine = HelixVizEngine(spec)
        cam_meta = (spec.metadata or {}).get("camera") if spec.metadata else None
        if isinstance(cam_meta, Mapping):
            self._custom_camera = dict(cam_meta)
            self._sync_camera_state_from_meta(cam_meta)
        else:
            self._custom_camera = None
            self._reset_camera_defaults()
        self._has_drawn_frame = False
        print(
            "[HelixModernWidget] engine geometry:",
            self.engine.helix_geometry.rail_vertices.shape,
            "arcs:",
            self.engine.repair_geometry.vertices.shape,
            flush=True,
        )
        if __debug__:
            self._assert_geometry()
        if self.ctx is not None:
            self._upload_buffers()
        self.update()
        if self.engine:
            self.phaseChanged.emit(self.engine.state.phase.value)

    def set_workflow_view(self, buffers: dict[str, np.ndarray] | None) -> None:
        """Switch into/out of the 2.5D workflow mode."""
        self._workflow_buffers = buffers
        self._mode = "workflow_2p5d" if buffers else "3d"
        if self.ctx is not None:
            self._upload_workflow_buffers()
        self.update()

    def set_playing(self, playing: bool) -> None:
        self._playing = playing

    def set_template_visible(self, visible: bool) -> None:
        self._show_template = visible
        self.update()

    def set_flap_visible(self, visible: bool) -> None:
        self._show_flap = visible
        self.update()

    def scrub_to(self, normalized_t: float) -> None:
        if self.engine is None:
            return
        target_time = max(0.0, min(1.0, normalized_t)) * self.engine.loop_duration
        state = self.engine.set_time(target_time)
        self.timeUpdated.emit(state.time)
        self.phaseChanged.emit(state.phase.value)
        self.update()

    # Qt overrides -----------------------------------------------------
    def initializeGL(self) -> None:  # pragma: no cover - requires OpenGL context
        print("[HelixModernWidget] initializeGL", flush=True)
        self.ctx = moderngl.create_context()
        print("[HelixModernWidget] GL version:", self.ctx.version_code, flush=True)
        self.ctx.enable(moderngl.BLEND)
        self.ctx.enable(moderngl.DEPTH_TEST)
        self.ctx.blend_func = (moderngl.SRC_ALPHA, moderngl.ONE)
        GL.glEnable(GL.GL_FRAMEBUFFER_SRGB)
        self._rail_program = self.ctx.program(
            vertex_shader=self._rail_vertex_shader(),
            fragment_shader=self._rail_fragment_shader(),
        )
        self._arc_program = self.ctx.program(
            vertex_shader=self._arc_vertex_shader(),
            fragment_shader=self._arc_fragment_shader(),
        )
        self._bars_prog = self.ctx.program(
            vertex_shader=self._bars_vertex_shader(),
            fragment_shader=self._bars_fragment_shader(),
        )
        self._heat_prog = self.ctx.program(
            vertex_shader=self._heat_vertex_shader(),
            fragment_shader=self._heat_fragment_shader(),
        )
        self._flat_prog = self.ctx.program(
            vertex_shader=self._flat_vertex_shader(),
            fragment_shader=self._flat_fragment_shader(),
        )
        self._upload_buffers()

    def resizeGL(self, width: int, height: int) -> None:  # pragma: no cover - requires OpenGL context
        if self.ctx is not None:
            scale = self.devicePixelRatioF()
            self.ctx.viewport = (0, 0, int(width * scale), int(height * scale))

    def paintGL(self) -> None:  # pragma: no cover - requires OpenGL context
        if self.ctx is None:
            return
        if self._mode == "workflow_2p5d" and self._workflow_buffers:
            self._render_workflow()
            if not self._has_drawn_frame:
                self._has_drawn_frame = True
                self.frameDrawn.emit()
            return
        qt_fbo = self.ctx.detect_framebuffer()
        qt_fbo.use()
        bound = GL.glGetIntegerv(GL.GL_DRAW_FRAMEBUFFER_BINDING)
        if bound != self._last_fbo:
            print(
                "[HelixModernWidget] draw framebuffer:",
                bound,
                "default:",
                self.defaultFramebufferObject(),
                flush=True,
            )
            self._last_fbo = bound
        self.ctx.clear(0.02, 0.02, 0.05, 1.0, 1.0)
        if not self.engine or not self._vao:
            return
        state = self.engine.state
        aspect = self.width() / max(1, self.height())
        mvp = self._mvp_matrix(aspect, orbit_time=state.time)
        self._rail_program["u_mvp"].write(mvp.tobytes())
        self._rail_program["u_time"].value = state.time
        self._rail_program["u_phase"].value = self._phase_value(state.phase)
        self._rail_program["u_template"].value = 1.0 if self._show_template else 0.0
        total_vertices = self.engine.helix_geometry.rail_vertices.shape[0]
        strand_vertices = total_vertices // 2
        self._vao.render(mode=moderngl.LINE_STRIP, vertices=strand_vertices, first=0)
        self._vao.render(mode=moderngl.LINE_STRIP, vertices=strand_vertices, first=strand_vertices)

        if self._arc_vao and self._arc_ranges:
            self._arc_program["u_mvp"].write(mvp.tobytes())
            self._arc_program["u_time"].value = state.time
            self._arc_program["u_flap"].value = 1.0 if self._show_flap else 0.0
            if self._arc_firsts.size and hasattr(self._arc_vao, "render"):
                try:
                    self._arc_vao.render(
                        mode=moderngl.LINE_STRIP,
                        vertices=self._arc_counts,
                        first=self._arc_firsts,
                    )
                except TypeError:
                    for start, count in self._arc_ranges:
                        self._arc_vao.render(mode=moderngl.LINE_STRIP, first=start, vertices=count)
            else:
                for start, count in self._arc_ranges:
                    self._arc_vao.render(mode=moderngl.LINE_STRIP, first=start, vertices=count)
        if not self._has_drawn_frame:
            self._has_drawn_frame = True
            self.frameDrawn.emit()

    # Internal helpers -------------------------------------------------
    def _tick(self) -> None:
        if self.engine and self._playing:
            state = self.engine.advance(1.0 / 60.0)
            self.timeUpdated.emit(state.time)
            self.phaseChanged.emit(state.phase.value)
        self.update()

    def _upload_buffers(self) -> None:
        if self.ctx is None or self.engine is None:
            return
        helix_geom = self.engine.helix_geometry
        repair_geom = self.engine.repair_geometry
        print(
            "[HelixModernWidget] uploading buffers:",
            f"rail_vertices={helix_geom.rail_vertices.shape}",
            f"rail_ids={helix_geom.rail_ids.shape}",
            f"param_coords={helix_geom.param_coords.shape}",
            f"repair_vertices={repair_geom.vertices.shape}",
            f"repair_offsets={repair_geom.offsets.shape}",
            flush=True,
        )
        payload = self.engine.buffer_payload()
        rail_vertices = payload["rail_vertices"]
        rail_ids = payload["rail_ids"][:, None]
        params = payload["param_coords"][:, None]
        combined = np.hstack([rail_vertices, rail_ids, params]).astype("f4")
        if self._rail_vbo is None:
            self._rail_vbo = self.ctx.buffer(combined.tobytes())
        else:
            self._rail_vbo.orphan(combined.nbytes)
            self._rail_vbo.write(combined.tobytes())
        self._vao = self.ctx.simple_vertex_array(
            self._rail_program,
            self._rail_vbo,
            "in_pos",
            "in_rail",
            "in_param",
        )

        arc_vertices = payload["arc_vertices"]
        arc_offsets = payload["arc_offsets"]
        arc_weights = payload.get("arc_weights")
        arc_kinds = payload.get("arc_kinds")
        self._arc_ranges = []
        self._arc_firsts = np.zeros((0,), dtype="i4")
        self._arc_counts = np.zeros((0,), dtype="i4")
        if arc_vertices.size:
            weights = arc_weights
            kinds = arc_kinds
            if weights is None or not weights.size:
                weights = np.ones((arc_vertices.shape[0],), dtype="f4")
            if kinds is None or not kinds.size:
                kinds = np.zeros((arc_vertices.shape[0],), dtype="f4")
            weights = weights.astype("f4").reshape(-1, 1)
            kinds = kinds.astype("f4").reshape(-1, 1)
            interleaved = np.hstack([arc_vertices.astype("f4"), weights, kinds])
            if self._arc_vbo is None:
                self._arc_vbo = self.ctx.buffer(interleaved.tobytes())
            else:
                self._arc_vbo.orphan(interleaved.nbytes)
                self._arc_vbo.write(interleaved.tobytes())
            self._arc_vao = self.ctx.simple_vertex_array(
                self._arc_program,
                self._arc_vbo,
                "in_pos",
                "in_weight",
                "in_kind",
            )
            firsts: list[int] = []
            counts: list[int] = []
            for i in range(len(arc_offsets) - 1):
                start = int(arc_offsets[i])
                end = int(arc_offsets[i + 1])
                if end > start:
                    length = end - start
                    self._arc_ranges.append((start, length))
                    firsts.append(start)
                    counts.append(length)
            if firsts:
                self._arc_firsts = np.array(firsts, dtype="i4")
                self._arc_counts = np.array(counts, dtype="i4")
        else:
            self._arc_vao = None

    def _upload_workflow_buffers(self) -> None:
        if self.ctx is None or not self._workflow_buffers:
            self._bars_vao = None
            self._workflow_instance_count = 0
            return
        buffers = self._workflow_buffers
        inst = buffers.get("inst_events")
        if inst is not None and inst.size:
            data = inst.astype("f4")
            col0 = data[:, :4].astype("f4")
            col1 = data[:, 4:].astype("f4")
            if self._bars_inst0_vbo is None:
                self._bars_inst0_vbo = self.ctx.buffer(col0.tobytes())
            else:
                self._bars_inst0_vbo.orphan(col0.nbytes)
                self._bars_inst0_vbo.write(col0.tobytes())
            if self._bars_inst1_vbo is None:
                self._bars_inst1_vbo = self.ctx.buffer(col1.tobytes())
            else:
                self._bars_inst1_vbo.orphan(col1.nbytes)
                self._bars_inst1_vbo.write(col1.tobytes())
            if self._bars_unit_vbo is None:
                unit = np.array([
                    [-0.5, -0.5],
                    [0.5, -0.5],
                    [0.5, 0.5],
                    [-0.5, -0.5],
                    [0.5, 0.5],
                    [-0.5, 0.5],
                ], dtype="f4")
                self._bars_unit_vbo = self.ctx.buffer(unit.tobytes())
            self._bars_vao = self.ctx.vertex_array(
            self._bars_prog,
                [
                    (self._bars_unit_vbo,  "2f",   "in_unit"),   # per-vertex
                    (self._bars_inst0_vbo, "4f/i", "in_box0"),   # per-instance
                    (self._bars_inst1_vbo, "3f/i", "in_box1"),   # per-instance
                ],
            )

            self._workflow_instance_count = data.shape[0]
        else:
            self._bars_vao = None
            self._workflow_instance_count = 0

        heat = buffers.get("heat_values")
        if heat is not None and getattr(self, "_heat_prog", None) is not None and heat.size >= 2:
            values = np.asarray(heat, dtype="f4").flatten()
            y = float(buffers.get("heat_y", -1.0))
            height = float(buffers.get("heat_height", 0.6))
            verts = np.zeros((values.shape[0] * 2, 4), dtype="f4")
            for i, val in enumerate(values):
                x = float(i)
                verts[2 * i] = (x, y, -0.5, val)
                verts[2 * i + 1] = (x, y + height, -0.5, val)
            if self._heat_vbo is None:
                self._heat_vbo = self.ctx.buffer(verts.tobytes())
            else:
                self._heat_vbo.orphan(verts.nbytes)
                self._heat_vbo.write(verts.tobytes())
            self._heat_vao = self.ctx.vertex_array(
                self._heat_prog,
                [
                    (self._heat_vbo, "3f 1f", "in_pos", "in_value"),
                ],
            )
            self._heat_vertex_count = verts.shape[0]
        else:
            self._heat_vao = None
            self._heat_vertex_count = 0

        cut_index = buffers.get("cut_index")
        if cut_index is not None and self._flat_prog is not None:
            cut = float(cut_index)
            width = float(buffers.get("slice_width", 0.12))
            extent = float(buffers.get("y_extent", 4.0)) + 1.0
            y_min = -extent
            y_max = extent
            data = np.array(
                [
                    [cut - width * 0.5, y_min, -0.4],
                    [cut - width * 0.5, y_max, -0.4],
                    [cut + width * 0.5, y_min, -0.4],
                    [cut + width * 0.5, y_max, -0.4],
                ],
                dtype="f4",
            )
            if self._slice_vbo is None:
                self._slice_vbo = self.ctx.buffer(data.tobytes())
            else:
                self._slice_vbo.orphan(data.nbytes)
                self._slice_vbo.write(data.tobytes())
            self._slice_vao = self.ctx.simple_vertex_array(self._flat_prog, self._slice_vbo, "in_pos")
        else:
            self._slice_vao = None

        rail = buffers.get("rail_xy")
        if rail is not None and rail.size:
            data = rail.astype("f4")
            if self._workflow_rail_vbo is None:
                self._workflow_rail_vbo = self.ctx.buffer(data.tobytes())
            else:
                self._workflow_rail_vbo.orphan(data.nbytes)
                self._workflow_rail_vbo.write(data.tobytes())
        else:
            self._workflow_rail_vbo = None

        ticks = buffers.get("tick_xy")
        if ticks is not None and ticks.size:
            data = ticks.astype("f4")
            if self._workflow_tick_vbo is None:
                self._workflow_tick_vbo = self.ctx.buffer(data.tobytes())
            else:
                self._workflow_tick_vbo.orphan(data.nbytes)
                self._workflow_tick_vbo.write(data.tobytes())
        else:
            self._workflow_tick_vbo = None

    def _mvp_matrix(self, aspect: float, orbit_time: float) -> np.ndarray:
        if self._custom_camera:
            return self._camera_mvp(aspect, self._custom_camera)
        radius = 2.2
        orbit = orbit_time * 0.25
        eye = np.array([
            radius * math.cos(orbit),
            radius * 0.35,
            radius * math.sin(orbit),
        ], dtype="f4")
        target = np.array([0.0, 0.0, 0.0], dtype="f4")
        up = np.array([0.0, 1.0, 0.0], dtype="f4")
        view = HelixModernWidget._look_at(eye, target, up)
        proj = HelixModernWidget._perspective(math.radians(40.0), aspect, 0.1, 20.0)
        mvp = proj @ view
        return mvp.astype("f4")

    def _camera_mvp(self, aspect: float, camera: Mapping[str, Any]) -> np.ndarray:
        center, radius = self._scene_bounds()
        default_pos = center + np.array([radius * 1.5, -radius * 1.2, radius * 0.6], dtype="f4")
        pos = np.array(camera.get("pos", default_pos), dtype="f4")
        target = np.array(camera.get("target", center), dtype="f4")
        up = np.array(camera.get("up", (0.0, 0.0, 1.0)), dtype="f4")
        fov = math.radians(float(camera.get("fov_deg", 45.0)))
        view = HelixModernWidget._look_at(pos, target, up)
        cam_dist = max(1e-3, float(np.linalg.norm(pos - target)))
        near = max(0.05, cam_dist - radius * 2.0)
        far = max(near + 5.0, cam_dist + radius * 2.5)
        proj = HelixModernWidget._perspective(fov, aspect, near, far)
        return (proj @ view).astype("f4")

    def _reset_camera_defaults(self) -> None:
        self._cam_target = np.array([0.0, 0.0, 0.0], dtype="f4")
        self._cam_distance = 6.0
        self._cam_yaw = 0.0
        self._cam_pitch = 0.3
        self._update_camera_meta()

    def reset_camera(self) -> None:
        """Public slot to reset the interactive camera."""
        self._reset_camera_defaults()
        self.update()

    def _sync_camera_state_from_meta(self, camera: Mapping[str, Any]) -> None:
        target = camera.get("target", (0.0, 0.0, 0.0))
        pos = camera.get("pos")
        if pos is None:
            self._reset_camera_defaults()
            return
        tx, ty, tz = target
        dx = pos[0] - tx
        dy = pos[1] - ty
        dz = pos[2] - tz
        distance = max(0.5, math.sqrt(dx * dx + dy * dy + dz * dz))
        yaw = math.atan2(dy, dx)
        pitch = math.asin(max(-0.999, min(0.999, dz / distance)))
        self._cam_target = np.array([tx, ty, tz], dtype="f4")
        self._cam_distance = distance
        self._cam_yaw = yaw
        self._cam_pitch = pitch

    def _update_camera_meta(self) -> None:
        cp = math.cos(self._cam_pitch)
        sp = math.sin(self._cam_pitch)
        cy = math.cos(self._cam_yaw)
        sy = math.sin(self._cam_yaw)
        tx, ty, tz = self._cam_target
        px = tx + self._cam_distance * cp * cy
        py = ty + self._cam_distance * cp * sy
        pz = tz + self._cam_distance * sp
        camera = {
            "pos": (float(px), float(py), float(pz)),
            "target": (float(tx), float(ty), float(tz)),
            "up": (0.0, 0.0, 1.0),
            "fov_deg": 45.0,
        }
        self._custom_camera = camera
        if self._spec is not None:
            meta = self._spec.metadata or {}
            meta["camera"] = camera
            self._spec.metadata = meta

    def _scene_bounds(self) -> tuple[np.ndarray, float]:
        if not self.engine:
            return np.zeros(3, dtype="f4"), 5.0
        vertices = getattr(self.engine.helix_geometry, "rail_vertices", None)
        if vertices is None or not vertices.size:
            return np.zeros(3, dtype="f4"), 5.0
        mn = vertices.min(axis=0)
        mx = vertices.max(axis=0)
        center = (mn + mx) * 0.5
        radius = float(np.linalg.norm(mx - center))
        return center.astype("f4"), max(radius, 1.0)

    def _ortho_mvp(self, left: float, right: float, bottom: float, top: float, near: float = -5.0, far: float = 5.0) -> np.ndarray:
        m = np.identity(4, dtype="f4")
        m[0, 0] = 2.0 / (right - left)
        m[1, 1] = 2.0 / (top - bottom)
        m[2, 2] = -2.0 / (far - near)
        m[0, 3] = -(right + left) / (right - left)
        m[1, 3] = -(top + bottom) / (top - bottom)
        m[2, 3] = -(far + near) / (far - near)
        return m

    def _render_workflow(self) -> None:
        if self.ctx is None or not self._workflow_buffers:
            return
        buffers = self._workflow_buffers
        x_extent = float(buffers.get("x_extent", 1.0))
        y_extent = float(buffers.get("y_extent", 4.0))
        left, right = -1.0, x_extent + 1.0
        bottom, top = -y_extent - 0.5, y_extent + 1.5
        mvp = self._ortho_mvp(left, right, bottom, top)
        px_scale = np.array([
            (right - left) / max(1.0, self.width() * self.devicePixelRatioF()),
            (top - bottom) / max(1.0, self.height() * self.devicePixelRatioF()),
        ], dtype="f4")
        if self._heat_vao and self._heat_vertex_count:
            self._heat_prog["u_mvp"].write(mvp.tobytes())
            self._heat_prog["u_color_low"].value = (0.08, 0.14, 0.28)
            self._heat_prog["u_color_high"].value = (0.98, 0.53, 0.18)
            self._heat_vao.render(mode=moderngl.TRIANGLE_STRIP, vertices=self._heat_vertex_count)
        if self._bars_vao and self._workflow_instance_count:
            self._bars_prog["u_mvp"].write(mvp.tobytes())
            self._bars_prog["u_px2world"].write(px_scale.tobytes())
            self._bars_vao.render(mode=moderngl.TRIANGLES, instances=self._workflow_instance_count)
        if self._slice_vao:
            self._flat_prog["u_mvp"].write(mvp.tobytes())
            self._flat_prog["u_color"].value = (1.0, 0.55, 0.18, 0.28)
            self._slice_vao.render(mode=moderngl.TRIANGLE_STRIP)
        if self._workflow_rail_vbo is not None:
            vao = self.ctx.simple_vertex_array(self._rail_program, self._workflow_rail_vbo, "in_pos")
            self._rail_program["u_mvp"].write(mvp.tobytes())
            self._rail_program["u_time"].value = 0.0
            self._rail_program["u_phase"].value = 0.0
            self._rail_program["u_template"].value = 1.0
            vao.render(mode=moderngl.LINE_STRIP)
        if self._workflow_tick_vbo is not None:
            vao = self.ctx.simple_vertex_array(self._rail_program, self._workflow_tick_vbo, "in_pos")
            vao.render(mode=moderngl.LINES)

    @staticmethod
    def _look_at(eye: np.ndarray, target: np.ndarray, up: np.ndarray) -> np.ndarray:
        f = target - eye
        f /= max(np.linalg.norm(f), 1e-6)
        u = up / max(np.linalg.norm(up), 1e-6)
        s = np.cross(f, u)
        s /= max(np.linalg.norm(s), 1e-6)
        u = np.cross(s, f)
        m = np.identity(4, dtype="f4")
        m[0, :3] = s
        m[1, :3] = u
        m[2, :3] = -f
        translate = np.identity(4, dtype="f4")
        translate[:3, 3] = -eye
        return m @ translate

    @staticmethod
    def _perspective(fovy: float, aspect: float, z_near: float, z_far: float) -> np.ndarray:
        f = 1.0 / math.tan(fovy / 2.0)
        m = np.zeros((4, 4), dtype="f4")
        m[0, 0] = f / aspect
        m[1, 1] = f
        m[2, 2] = (z_far + z_near) / (z_near - z_far)
        m[2, 3] = (2 * z_far * z_near) / (z_near - z_far)
        m[3, 2] = -1.0
        return m

    @staticmethod
    def _phase_value(phase: AnimationPhase) -> float:
        mapping = {
            AnimationPhase.APPROACH: 0.0,
            AnimationPhase.RECOGNITION: 0.2,
            AnimationPhase.CUT: 0.4,
            AnimationPhase.PRIME: 0.6,
            AnimationPhase.REPAIR: 0.8,
            AnimationPhase.RESOLVE: 1.0,
        }
        return mapping.get(phase, 0.0)

    # Interaction events -----------------------------------------------
    def mousePressEvent(self, event) -> None:  # type: ignore[override]
        self._last_mouse_pos = event.position()
        if event.button() == Qt.LeftButton:
            self._drag_mode = "orbit"
        elif event.button() == Qt.MiddleButton or (event.button() == Qt.RightButton and event.modifiers() & Qt.ControlModifier):
            self._drag_mode = "pan"
        else:
            self._drag_mode = None
        super().mousePressEvent(event)

    def mouseReleaseEvent(self, event) -> None:  # type: ignore[override]
        self._drag_mode = None
        self._last_mouse_pos = None
        super().mouseReleaseEvent(event)

    def mouseMoveEvent(self, event) -> None:  # type: ignore[override]
        if self._last_mouse_pos is None or self._drag_mode is None:
            super().mouseMoveEvent(event)
            return
        pos = event.position()
        dx = pos.x() - self._last_mouse_pos.x()
        dy = pos.y() - self._last_mouse_pos.y()
        self._last_mouse_pos = pos
        if self._drag_mode == "orbit":
            self._cam_yaw += dx * self._cam_sensitivity
            self._cam_pitch -= dy * self._cam_sensitivity
            limit = math.radians(89.0)
            self._cam_pitch = max(-limit, min(limit, self._cam_pitch))
        elif self._drag_mode == "pan":
            cy = math.cos(self._cam_yaw)
            sy = math.sin(self._cam_yaw)
            cp = math.cos(self._cam_pitch)
            sp = math.sin(self._cam_pitch)
            forward = np.array([cp * cy, cp * sy, sp], dtype="f4")
            up = np.array([0.0, 0.0, 1.0], dtype="f4")
            right = np.cross(forward, up)
            pan_scale = self._cam_distance * self._pan_sensitivity
            self._cam_target -= dx * pan_scale * right
            self._cam_target += dy * pan_scale * up
        self._update_camera_meta()
        self.update()

    def wheelEvent(self, event) -> None:  # type: ignore[override]
        delta = event.angleDelta().y() / 120.0
        factor = math.exp(-delta * self._zoom_sensitivity)
        self._cam_distance = max(0.5, self._cam_distance * factor)
        self._update_camera_meta()
        self.update()
        super().wheelEvent(event)

    def keyPressEvent(self, event) -> None:  # type: ignore[override]
        if event.key() == Qt.Key_R:
            self._reset_camera_defaults()
            self.update()
            return
        super().keyPressEvent(event)

    def _assert_geometry(self) -> None:
        if not self.engine:
            return
        try:
            rv = self.engine.helix_geometry.rail_vertices
            if rv.size == 0:
                return
            n = rv.shape[0] // 2
            if n == 0:
                return
            strand_a = rv[:n]
            strand_b = rv[n:]
            assert np.all(np.diff(strand_a[:, 2]) >= -1e-6), "Strand A z not monotonic"
            assert np.all(np.diff(strand_b[:, 2]) >= -1e-6), "Strand B z not monotonic"
            assert np.allclose(strand_a[:, 2], strand_b[:, 2], atol=1e-3), "Strand z mismatch"
            assert np.allclose(strand_a[:, 0] + strand_b[:, 0], 0.0, atol=0.1), "Strands not opposite in X"
            assert np.allclose(strand_a[:, 1] + strand_b[:, 1], 0.0, atol=0.1), "Strands not opposite in Y"
        except AssertionError as exc:
            print(f"[HelixModernWidget] Geometry invariant failed: {exc}", flush=True)

    @staticmethod
    def _rail_vertex_shader() -> str:
        return """
        #version 330
        in vec3 in_pos;
        in float in_rail;
        in float in_param;
        uniform mat4 u_mvp;
        uniform float u_time;
        uniform float u_phase;
        out float v_param;
        out float v_rail;
        out float v_phase;
        void main() {
            float ripple = sin(in_param * 6.2831 + u_time * 0.8) * 0.02;
            vec3 pos = in_pos;
            pos.xy *= 1.05 + ripple;
            gl_Position = u_mvp * vec4(pos, 1.0);
            v_param = in_param;
            v_rail = in_rail;
            v_phase = u_phase;
        }
        """

    @staticmethod
    def _rail_fragment_shader() -> str:
        return """
        #version 330
        in float v_param;
        in float v_rail;
        in float v_phase;
        uniform float u_template;
        out vec4 f_color;
        vec3 gradient(float t) {
            vec3 pre = vec3(0.04, 0.65, 0.64);
            vec3 cut = vec3(1.0, 1.0, 1.0);
            vec3 repair = vec3(0.96, 0.62, 0.18);
            return mix(pre, repair, t);
        }
        void main() {
            float glow = smoothstep(0.9, 1.0, v_phase);
            vec3 base = gradient(v_param);
            vec3 rail_offset = mix(vec3(0.0), vec3(0.02, 0.08, 0.1), v_rail);
            vec3 color = base + rail_offset + glow * 0.2;
            float alpha = mix(0.15, 1.0, u_template);
            f_color = vec4(color, alpha);
        }
        """

    @staticmethod
    def _arc_vertex_shader() -> str:
        return """
        #version 330
        in vec3  in_pos;
        in float in_weight;
        in float in_kind;
        uniform mat4 u_mvp;
        uniform float u_time;
        out float v_weight;
        out float v_kind;
        void main() {
            vec3 pos = in_pos;
            pos *= 1.02 + sin(u_time * 1.5) * 0.03;
            gl_Position = u_mvp * vec4(pos, 1.0);
            v_weight = in_weight;
            v_kind = in_kind;
        }
        """

    @staticmethod
    def _arc_fragment_shader() -> str:
        return """
        #version 330
        uniform float u_flap;
        in float v_weight;
        in float v_kind;
        out vec4 f_color;
        vec3 kind_color(float k) {
            if (k < 0.5) return vec3(0.98, 0.65, 0.22);
            if (k < 1.5) return vec3(0.92, 0.33, 0.82);
            if (k < 2.5) return vec3(0.26, 0.86, 0.98);
            if (k < 3.5) return vec3(0.20, 0.40, 0.55);
            return vec3(0.76, 0.65, 0.92);
        }
        void main() {
            float w = clamp(v_weight, 0.0, 1.0);
            float alpha = mix(0.18, 0.95, w) * u_flap;
            vec3 color = kind_color(v_kind);
            f_color = vec4(color, alpha);
        }
        """

    @staticmethod
    def _bars_vertex_shader() -> str:
        return """
        #version 330
        in vec2 in_unit;
        in vec4 in_box0;
        in vec3 in_box1;
        uniform mat4 u_mvp;
        uniform vec2 u_px2world;
        out float v_weight;
        out float v_kind;
        void main() {
            float xc = in_box0.x;
            float hw = max(in_box0.y, 1e-5);
            float yc = in_box0.z;
            float base_h = in_box0.w;
            float w = in_box1.x;
            float kind = in_box1.y;
            float z = in_box1.z;
            float px_h = 14.0 * u_px2world.y;
            float world_h = max(base_h, px_h);
            float x = xc + in_unit.x * (2.0 * hw);
            float y = yc + in_unit.y * world_h;
            gl_Position = u_mvp * vec4(x, y, z, 1.0);
            v_weight = w;
            v_kind = kind;
        }
        """

    @staticmethod
    def _bars_fragment_shader() -> str:
        return """
        #version 330
        in float v_weight;
        in float v_kind;
        out vec4 f_color;
        vec3 kind_color(float kind) {
            if (kind < 0.5) return vec3(0.98, 0.65, 0.22);
            if (kind < 1.5) return vec3(0.92, 0.33, 0.82);
            if (kind < 2.5) return vec3(0.26, 0.86, 0.98);
            if (kind < 3.5) return vec3(0.20, 0.40, 0.55);
            if (kind < 4.5) return vec3(0.35, 0.88, 0.75);
            return vec3(0.65, 0.65, 0.72);
        }
        void main() {
            float w = clamp(v_weight, 0.0, 1.0);
            vec3 color = kind_color(v_kind);
            float alpha = mix(0.2, 0.95, w);
            f_color = vec4(color, alpha);
        }
        """

    @staticmethod
    def _heat_vertex_shader() -> str:
        return """
        #version 330
        in vec3 in_pos;
        in float in_value;
        uniform mat4 u_mvp;
        out float v_value;
        void main() {
            gl_Position = u_mvp * vec4(in_pos, 1.0);
            v_value = in_value;
        }
        """

    @staticmethod
    def _heat_fragment_shader() -> str:
        return """
        #version 330
        in float v_value;
        uniform vec3 u_color_low;
        uniform vec3 u_color_high;
        out vec4 f_color;
        void main() {
            float t = clamp(v_value, 0.0, 1.0);
            vec3 color = mix(u_color_low, u_color_high, t);
            float alpha = mix(0.12, 0.55, t);
            f_color = vec4(color, alpha);
        }
        """

    @staticmethod
    def _flat_vertex_shader() -> str:
        return """
        #version 330
        in vec3 in_pos;
        uniform mat4 u_mvp;
        void main() {
            gl_Position = u_mvp * vec4(in_pos, 1.0);
        }
        """

    @staticmethod
    def _flat_fragment_shader() -> str:
        return """
        #version 330
        uniform vec4 u_color;
        out vec4 f_color;
        void main() {
            f_color = u_color;
        }
        """


class HelixControlPanel(QWidget):
    """Side panel with edit selector, transport, and toggles."""

    specSelected = Signal(int)
    timeScrubbed = Signal(float)
    playToggled = Signal(bool)
    templateToggled = Signal(bool)
    flapToggled = Signal(bool)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._block_slider = False
        self._benchmark_worker: BenchmarkWorker | None = None
        self._benchmark_config = BenchmarkConfig(
            backends=["cpu-reference", "native-cpu", "gpu"],
            crispr_shapes=[(1, 512, 20), (96, 4096, 20), (256, 16384, 20)],
            prime_workloads=[(32, 2000, 20), (64, 4000, 20)],
            seed=0,
        )
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(10)

        self.phase_label = QLabel("Phase: approach", self)
        layout.addWidget(self.phase_label)

        self.spec_combo = QComboBox(self)
        self.spec_combo.currentIndexChanged.connect(self.specSelected.emit)
        layout.addWidget(self.spec_combo)

        self.play_button = QPushButton("Pause", self)
        self.play_button.setCheckable(True)
        self.play_button.setChecked(True)
        self.play_button.toggled.connect(self._on_play_toggled)
        layout.addWidget(self.play_button)

        slider_label = QLabel("Time", self)
        layout.addWidget(slider_label)
        self.time_slider = QSlider(Qt.Horizontal, self)
        self.time_slider.setRange(0, 1000)
        self.time_slider.valueChanged.connect(self._on_slider)
        layout.addWidget(self.time_slider)

        self.template_toggle = QCheckBox("Show template strand", self)
        self.template_toggle.setChecked(True)
        self.template_toggle.toggled.connect(self.templateToggled.emit)
        layout.addWidget(self.template_toggle)

        self.flap_toggle = QCheckBox("Show flap/template path", self)
        self.flap_toggle.setChecked(True)
        self.flap_toggle.toggled.connect(self.flapToggled.emit)
        layout.addWidget(self.flap_toggle)

        self.engine_box = QGroupBox("Engine Health", self)
        engine_layout = QVBoxLayout(self.engine_box)
        self.backend_label = QLabel("CRISPR engine: --", self.engine_box)
        engine_layout.addWidget(self.backend_label)
        self.scoring_label = QLabel("Scoring versions: -- / --", self.engine_box)
        engine_layout.addWidget(self.scoring_label)
        self.benchmark_status = QLabel("Benchmark not run", self.engine_box)
        engine_layout.addWidget(self.benchmark_status)
        self.run_benchmark_button = QPushButton("Run Benchmark", self.engine_box)
        self.run_benchmark_button.clicked.connect(self._on_run_benchmark)
        engine_layout.addWidget(self.run_benchmark_button)
        self.benchmark_output = QTextBrowser(self.engine_box)
        self.benchmark_output.setReadOnly(True)
        self.benchmark_output.setMinimumHeight(120)
        engine_layout.addWidget(self.benchmark_output)
        layout.addWidget(self.engine_box)

        self.physics_box = QGroupBox("Prime Physics Score", self)
        phys_layout = QFormLayout(self.physics_box)
        self.pbs_label = QLabel("--", self.physics_box)
        self.pbs_hint_label = QLabel("", self.physics_box)
        phys_layout.addRow("PBS ΔG", self.pbs_label)
        phys_layout.addRow("PBS hint", self.pbs_hint_label)
        self.rtt_label = QLabel("--", self.physics_box)
        phys_layout.addRow("RT path ΔG", self.rtt_label)
        self.flap_label = QLabel("--", self.physics_box)
        phys_layout.addRow("Flap ΔΔG", self.flap_label)
        self.micro_label = QLabel("--", self.physics_box)
        phys_layout.addRow("Microhomology", self.micro_label)
        self.nick_label = QLabel("--", self.physics_box)
        phys_layout.addRow("Nick distance", self.nick_label)
        self.p_rt_label = QLabel("--", self.physics_box)
        phys_layout.addRow("P_RT", self.p_rt_label)
        self.p_flap_label = QLabel("--", self.physics_box)
        phys_layout.addRow("P_flap", self.p_flap_label)
        self.e_pred_label = QLabel("--", self.physics_box)
        phys_layout.addRow("E_pred", self.e_pred_label)
        layout.addWidget(self.physics_box)
        self.physics_box.hide()

        layout.addStretch(1)

    def set_specs(self, labels: Sequence[str]) -> None:
        block = self.spec_combo.blockSignals(True)
        self.spec_combo.clear()
        self.spec_combo.addItems(labels)
        self.spec_combo.blockSignals(block)

    def set_time_norm(self, value: float) -> None:
        clamped = max(0.0, min(1.0, value))
        self._block_slider = True
        self.time_slider.setValue(int(clamped * 1000))
        self._block_slider = False

    def set_phase(self, phase: str) -> None:
        self.phase_label.setText(f"Phase: {phase}")

    def set_engine_status(self, text: str) -> None:
        self.backend_label.setText(text)

    def set_scoring_versions(self, crispr: str, prime: str) -> None:
        self.scoring_label.setText(f"Scoring versions: CRISPR {crispr} | Prime {prime}")

    def set_prime_physics_score(self, score: Mapping[str, Any] | None) -> None:
        if not score:
            self.physics_box.hide()
            return
        self.physics_box.show()
        pbs_dg = float(score.get("pbs_dG", 0.0))
        self.pbs_label.setText(f"{pbs_dg:.2f} kcal/mol")
        self.pbs_hint_label.setText(self._pbs_hint(pbs_dg))
        rt_cum = score.get("rt_cum_dG") or []
        rt_final = float(rt_cum[-1]) if rt_cum else 0.0
        self.rtt_label.setText(f"{rt_final:.2f} kcal/mol")
        flap_ddg = float(score.get("flap_ddG", 0.0))
        self.flap_label.setText(f"{flap_ddg:.2f}")
        micro = int(score.get("microhomology", 0))
        self.micro_label.setText(f"{micro} nt")
        self.nick_label.setText(str(int(score.get("nick_distance", 0))))
        p_rt = float(score.get("P_RT", 0.0))
        p_flap = float(score.get("P_flap", 0.0))
        self.p_rt_label.setText(f"{p_rt:.3f}")
        self.p_flap_label.setText(f"{p_flap:.3f}")
        e_pred = float(score.get("E_pred", 0.0))
        self.e_pred_label.setText(f"{e_pred:.3f}")
        self.e_pred_label.setStyleSheet(f"color: {self._score_color(e_pred)};")

    def _on_slider(self, value: int) -> None:
        if self._block_slider:
            return
        self.timeScrubbed.emit(value / 1000.0)

    def _on_play_toggled(self, checked: bool) -> None:
        self.play_button.setText("Pause" if checked else "Play")
        self.playToggled.emit(checked)

    def _on_run_benchmark(self) -> None:
        if self._benchmark_worker is not None:
            return
        self.run_benchmark_button.setEnabled(False)
        self.benchmark_status.setText("Running benchmark…")
        self.benchmark_output.clear()
        self._benchmark_worker = BenchmarkWorker(self._benchmark_config, self)
        self._benchmark_worker.completed.connect(self._on_benchmark_complete)
        self._benchmark_worker.failed.connect(self._on_benchmark_failed)
        self._benchmark_worker.finished.connect(self._on_benchmark_finished)
        self._benchmark_worker.start()

    def _on_benchmark_complete(self, payload: Mapping[str, Any]) -> None:
        self.benchmark_status.setText("Benchmark complete")
        self._update_benchmark_output(payload)

    def _on_benchmark_failed(self, message: str) -> None:
        self.benchmark_status.setText(f"Benchmark failed: {message}")

    def _on_benchmark_finished(self) -> None:
        self.run_benchmark_button.setEnabled(True)
        self._benchmark_worker = None

    def _update_benchmark_output(self, payload: Mapping[str, Any]) -> None:
        lines = ["CRISPR:"]
        for entry in payload.get("benchmarks", {}).get("crispr", []):
            if "error" in entry:
                lines.append(
                    f"  {entry.get('backend_requested')}: {entry['error']}"
                )
            else:
                lines.append(
                    f"  {entry.get('backend_used') or entry.get('backend_requested')}: {entry.get('mpairs_per_s', 0):.2f} MPairs/s"
                )
        lines.append("Prime:")
        for entry in payload.get("benchmarks", {}).get("prime", []):
            lines.append(
                f"  {entry.get('backend_used') or entry.get('backend_requested')}: {entry.get('predictions_per_s', 0):.2f} preds/sec"
            )
        self.benchmark_output.setPlainText("\n".join(lines))

    def _pbs_hint(self, dg: float) -> str:
        if dg <= -10:
            return "Strong binding"
        if dg <= -3:
            return "Within target"
        return "Weak binding"

    def _score_color(self, value: float) -> str:
        clamped = max(0.0, min(1.0, value))
        red = int((1.0 - clamped) * 200 + 55)
        green = int(clamped * 200 + 55)
        return f"rgb({red},{green},90)"


class HelixModernWindow(QMainWindow):
    """Main window combining the ModernGL viewport and control panel."""

    def __init__(self, specs: Sequence[EditVisualizationSpec]) -> None:
        super().__init__()
        self.setWindowTitle("Helix – ModernGL CRISPR Viz")
        self.viewer = HelixModernWidget()
        self.viewer_container = QWidget.createWindowContainer(self.viewer, self)
        self.panel = HelixControlPanel(self)
        self.panel.specSelected.connect(self._on_spec_selected)
        self.panel.timeScrubbed.connect(self.viewer.scrub_to)
        self.panel.playToggled.connect(self.viewer.set_playing)
        self.panel.templateToggled.connect(self.viewer.set_template_visible)
        self.panel.flapToggled.connect(self.viewer.set_flap_visible)
        self.viewer.timeUpdated.connect(self._on_time_update)
        self.viewer.phaseChanged.connect(self.panel.set_phase)

        main = QWidget(self)
        layout = QHBoxLayout(main)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.viewer_container, stretch=3)
        layout.addWidget(self.panel, stretch=1)
        self.setCentralWidget(main)

        self.specs = list(specs)
        labels = [spec.metadata.get("label") or spec.metadata.get("name") or spec.edit_type for spec in self.specs]
        self.panel.set_specs(labels)
        backend, allow = resolve_crispr_backend(None, use_gpu=False)
        native_ok = native_backend_available()
        status_parts = [f"CRISPR engine: {backend}"]
        status_parts.append(f"native {'available' if native_ok else 'missing'}")
        status_parts.append(f"fallback {'on' if allow else 'off'}")
        self.panel.set_engine_status(" | ".join(status_parts))
        self.panel.set_scoring_versions(CRISPR_SCORING_VERSION, PRIME_SCORING_VERSION)
        if self.specs:
            self.viewer.set_spec(self.specs[0])
            self._update_prime_physics(self.specs[0])

    def _on_spec_selected(self, index: int) -> None:
        if 0 <= index < len(self.specs):
            spec = self.specs[index]
            self.viewer.set_spec(spec)
            self._update_prime_physics(spec)

    def _on_time_update(self, time_value: float) -> None:
        if not self.viewer.engine:
            return
        norm = time_value / self.viewer.engine.loop_duration
        self.panel.set_time_norm(norm)

    def _update_prime_physics(self, spec: EditVisualizationSpec) -> None:
        metadata = getattr(spec, "metadata", {}) or {}
        score = metadata.get("physics_score")
        self.panel.set_prime_physics_score(score if isinstance(score, Mapping) else None)


def run_modern_viz(spec_source: str | Path | Sequence[EditVisualizationSpec] | None = None) -> None:
    """Launch the PyQt + ModernGL helix visualizer."""

    app = QApplication.instance()
    owns_app = False
    if app is None:
        app = QApplication(sys.argv)
        owns_app = True
    apply_helix_theme(app)

    if spec_source is None:
        specs = [load_viz_spec(
            {
                "sequence": "ACGTACGTACGTACGTACGTACGT",
                "pam_index": 10,
                "guide_range": [4, 24],
                "edit_type": "prime",
                "edit_events": [
                    {"t": 0.2, "type": "recognition", "index": 8},
                    {"t": 0.5, "type": "nick_primary", "index": 10},
                    {"t": 0.8, "type": "rt_synthesis", "start": 11, "length": 6},
                    {"t": 1.2, "type": "flap_resolution", "index": 12},
                    {"t": 1.6, "type": "repair_complete", "index": 12},
                ],
            }
        )]
    elif isinstance(spec_source, Sequence) and spec_source and isinstance(spec_source[0], EditVisualizationSpec):  # type: ignore[index]
        specs = list(spec_source)  # type: ignore[assignment]
    else:
        specs = load_viz_specs(spec_source)  # type: ignore[arg-type]

    window = HelixModernWindow(specs)
    window.resize(1280, 800)
    window.show()

    if owns_app:
        sys.exit(app.exec())
    else:  # pragma: no cover - embedded Qt contexts
        app.exec()


if __name__ == "__main__":  # pragma: no cover - manual launch helper
    source = Path(sys.argv[1]).expanduser() if len(sys.argv) > 1 else None
    run_modern_viz(source)
