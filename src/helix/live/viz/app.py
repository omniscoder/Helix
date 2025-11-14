"""GLFW + moderngl realtime viewer with a custom Helix HUD."""

from __future__ import annotations

import time
import uuid
from collections import deque
from pathlib import Path
from typing import Deque, Dict, Optional

try:  # pragma: no cover - optional dependency guard
    import glfw  # type: ignore
    import moderngl  # type: ignore

    from .channel import VizChannel
    from .hud import HelixHUD, HudEvent
    from .renderer import SceneRenderer
except Exception as exc:  # pragma: no cover - import guard
    RUNTIME_ERROR = exc
else:  # pragma: no cover - import guard
    RUNTIME_ERROR = None


class LiveVizApp:
    """Minimal standalone viewer for helix.live delta/control streams."""

    def __init__(self, *, endpoint: Optional[str], bundle: Optional[Path], target_hz: float = 60.0):
        if RUNTIME_ERROR is not None:
            raise RuntimeError(
                "Helix live viz requires the 'realtime' extra (glfw + moderngl). "
                "Install with `pip install \"veri-helix[realtime]\"`."
            ) from RUNTIME_ERROR
        self.channel = VizChannel(endpoint=endpoint, bundle=bundle)
        self.target_dt = 1.0 / max(target_hz, 1.0)
        self.window = None
        self.ctx: Optional[moderngl.Context] = None
        self.renderer: Optional[SceneRenderer] = None
        self.hud: Optional[HelixHUD] = None
        self.metric_history: Dict[str, Deque[float]] = {
            "pERK": deque(maxlen=512),
            "cell_count": deque(maxlen=512),
        }
        self._last_frame = time.time()
        self._paused = False
        self._space_prev = False

    def run(self) -> None:
        self._init_window()
        try:
            while not glfw.window_should_close(self.window):
                self._tick()
        finally:
            self.channel.close()
            glfw.terminate()

    # Internal helpers -----------------------------------------------------

    def _init_window(self) -> None:
        if not glfw.init():
            raise RuntimeError("Failed to initialize GLFW.")
        glfw.window_hint(glfw.CONTEXT_VERSION_MAJOR, 3)
        glfw.window_hint(glfw.CONTEXT_VERSION_MINOR, 3)
        glfw.window_hint(glfw.OPENGL_PROFILE, glfw.OPENGL_CORE_PROFILE)
        self.window = glfw.create_window(1280, 720, "Helix Live", None, None)
        if not self.window:
            glfw.terminate()
            raise RuntimeError("Failed to create GLFW window.")
        glfw.make_context_current(self.window)
        self.ctx = moderngl.create_context()
        self.renderer = SceneRenderer(self.ctx)
        self.hud = HelixHUD(self.ctx)

    def _tick(self) -> None:
        now = time.time()
        remaining = self.target_dt - (now - self._last_frame)
        if remaining > 0:
            time.sleep(remaining)
        self._last_frame = time.time()
        glfw.poll_events()
        frame = self.channel.poll()
        if frame and self.renderer:
            self.renderer.apply_frame(frame)
            self._update_metric_history()
        if self.ctx and self.renderer and self.hud:
            width, height = glfw.get_framebuffer_size(self.window)
            metrics = dict(self.renderer.metrics)
            self.hud.build_controls(
                width=width,
                height=height,
                metrics=metrics,
                metric_history={k: list(v) for k, v in self.metric_history.items()},
                slice_name=self.channel.slice,
                variant=self.channel.variant,
                variant_options=self.channel.variant_options,
                paused=self._paused,
                interactive=self.channel.interactive,
            )
            self._handle_input(width, height)
            self.ctx.viewport = (0, 0, width, height)
            self.ctx.clear(0.04, 0.05, 0.08, 1.0)
            self.renderer.draw()
            self.hud.draw()
            glfw.swap_buffers(self.window)

    def _handle_input(self, width: int, height: int) -> None:
        if not self.hud:
            return
        cursor_x, cursor_y = glfw.get_cursor_pos(self.window)
        cursor_y = float(height) - float(cursor_y)
        left_down = glfw.get_mouse_button(self.window, glfw.MOUSE_BUTTON_LEFT) == glfw.PRESS
        events = self.hud.handle_mouse(
            float(cursor_x),
            float(cursor_y),
            left_down,
            interactive=self.channel.interactive,
        )
        for event in events:
            self._dispatch_hud_event(event)

        space_pressed = glfw.get_key(self.window, glfw.KEY_SPACE) == glfw.PRESS
        if (
            self.channel.interactive
            and space_pressed
            and not self._space_prev
        ):
            self._toggle_pause()
        self._space_prev = space_pressed

    def _dispatch_hud_event(self, event: HudEvent) -> None:
        if not self.channel.interactive:
            return
        if event.id == "pause" and event.kind == "click":
            self._toggle_pause()
        elif event.id == "grb2_scale" and event.kind == "change" and event.value is not None:
            self._send_control(
                {
                    "type": "set_param",
                    "target": "egfr_rules",
                    "param": "grb2_scale",
                    "value": float(event.value),
                }
            )
            self.renderer and self.renderer.metrics.__setitem__("grb2_scale", float(event.value))
        elif event.id == "variant" and event.kind == "change" and event.value:
            self._send_control(
                {
                    "type": "set_variant",
                    "variant": str(event.value),
                }
            )

    def _toggle_pause(self) -> None:
        control_type = "resume" if self._paused else "pause"
        self._send_control({"type": control_type})
        self._paused = not self._paused

    def _send_control(self, payload: Dict[str, object]) -> None:
        base = {
            "schema_version": 1,
            "kind": "live_control",
            "slice": self.channel.slice,
            "variant": self.channel.variant,
            "run_id": self.channel.run_id or "helix_run",
            "event_id": str(uuid.uuid4()),
            "ts_wall": time.time(),
        }
        base.update(payload)
        self.channel.send_control(base)

    def _update_metric_history(self) -> None:
        if not self.renderer:
            return
        metrics = self.renderer.metrics
        if "pERK" in metrics:
            self.metric_history["pERK"].append(float(metrics["pERK"]))
        if "cell_count" in metrics:
            self.metric_history["cell_count"].append(float(metrics["cell_count"]))
