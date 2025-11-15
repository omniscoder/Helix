"""Custom HUD primitives for the realtime viewer."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import moderngl
import numpy as np

# font8x8_basic subset (U+0020-U+007F). Public domain by Daniel Hepper.
FONT_8x8: Dict[str, Tuple[int, ...]] = {
    " ": (0, 0, 0, 0, 0, 0, 0, 0),
    "!": (24, 60, 60, 24, 24, 0, 24, 0),
    '"': (54, 54, 0, 0, 0, 0, 0, 0),
    "#": (54, 54, 127, 54, 127, 54, 54, 0),
    "$": (12, 62, 3, 30, 48, 31, 12, 0),
    "%": (0, 99, 51, 24, 12, 102, 99, 0),
    "&": (28, 54, 28, 110, 59, 51, 110, 0),
    "'": (6, 6, 3, 0, 0, 0, 0, 0),
    "(": (24, 12, 6, 6, 6, 12, 24, 0),
    ")": (6, 12, 24, 24, 24, 12, 6, 0),
    "*": (0, 102, 60, 255, 60, 102, 0, 0),
    "+": (0, 12, 12, 63, 12, 12, 0, 0),
    ",": (0, 0, 0, 0, 0, 12, 12, 6),
    "-": (0, 0, 0, 63, 0, 0, 0, 0),
    ".": (0, 0, 0, 0, 0, 12, 12, 0),
    "/": (96, 48, 24, 12, 6, 3, 1, 0),
    "0": (62, 99, 115, 123, 111, 103, 62, 0),
    "1": (12, 14, 12, 12, 12, 12, 63, 0),
    "2": (30, 51, 48, 28, 6, 51, 63, 0),
    "3": (30, 51, 48, 28, 48, 51, 30, 0),
    "4": (56, 60, 54, 51, 127, 48, 120, 0),
    "5": (63, 3, 31, 48, 48, 51, 30, 0),
    "6": (28, 6, 3, 31, 51, 51, 30, 0),
    "7": (63, 51, 48, 24, 12, 12, 12, 0),
    "8": (30, 51, 51, 30, 51, 51, 30, 0),
    "9": (30, 51, 51, 62, 48, 24, 14, 0),
    ":": (0, 12, 12, 0, 0, 12, 12, 0),
    ";": (0, 12, 12, 0, 0, 12, 12, 6),
    "<": (24, 12, 6, 3, 6, 12, 24, 0),
    "=": (0, 0, 63, 0, 0, 63, 0, 0),
    ">": (6, 12, 24, 48, 24, 12, 6, 0),
    "?": (30, 51, 48, 24, 12, 0, 12, 0),
    "@": (62, 99, 123, 123, 123, 3, 30, 0),
    "A": (12, 30, 51, 51, 63, 51, 51, 0),
    "B": (63, 102, 102, 62, 102, 102, 63, 0),
    "C": (60, 102, 3, 3, 3, 102, 60, 0),
    "D": (31, 54, 102, 102, 102, 54, 31, 0),
    "E": (127, 70, 22, 30, 22, 70, 127, 0),
    "F": (127, 70, 22, 30, 22, 6, 15, 0),
    "G": (60, 102, 3, 3, 115, 102, 124, 0),
    "H": (51, 51, 51, 63, 51, 51, 51, 0),
    "I": (30, 12, 12, 12, 12, 12, 30, 0),
    "J": (120, 48, 48, 48, 51, 51, 30, 0),
    "K": (103, 102, 54, 30, 54, 102, 103, 0),
    "L": (15, 6, 6, 6, 70, 102, 127, 0),
    "M": (99, 119, 127, 127, 107, 99, 99, 0),
    "N": (99, 103, 111, 123, 115, 99, 99, 0),
    "O": (28, 54, 99, 99, 99, 54, 28, 0),
    "P": (63, 102, 102, 62, 6, 6, 15, 0),
    "Q": (30, 51, 51, 51, 59, 30, 56, 0),
    "R": (63, 102, 102, 62, 54, 102, 103, 0),
    "S": (30, 51, 7, 14, 56, 51, 30, 0),
    "T": (63, 45, 12, 12, 12, 12, 30, 0),
    "U": (51, 51, 51, 51, 51, 51, 63, 0),
    "V": (51, 51, 51, 51, 51, 30, 12, 0),
    "W": (99, 99, 99, 107, 127, 119, 99, 0),
    "X": (99, 99, 54, 28, 28, 54, 99, 0),
    "Y": (51, 51, 51, 30, 12, 12, 30, 0),
    "Z": (127, 99, 49, 24, 76, 102, 127, 0),
    "[": (30, 6, 6, 6, 6, 6, 30, 0),
    "\\": (3, 6, 12, 24, 48, 96, 64, 0),
    "]": (30, 24, 24, 24, 24, 24, 30, 0),
    "^": (8, 28, 54, 99, 0, 0, 0, 0),
    "_": (0, 0, 0, 0, 0, 0, 0, 255),
    "`": (12, 12, 24, 0, 0, 0, 0, 0),
    "a": (0, 0, 30, 48, 62, 51, 110, 0),
    "b": (7, 6, 6, 62, 102, 102, 59, 0),
    "c": (0, 0, 30, 51, 3, 51, 30, 0),
    "d": (56, 48, 48, 62, 51, 51, 110, 0),
    "e": (0, 0, 30, 51, 63, 3, 30, 0),
    "f": (28, 54, 6, 15, 6, 6, 15, 0),
    "g": (0, 0, 110, 51, 51, 62, 48, 31),
    "h": (7, 6, 54, 110, 102, 102, 103, 0),
    "i": (12, 0, 14, 12, 12, 12, 30, 0),
    "j": (48, 0, 48, 48, 48, 51, 51, 30),
    "k": (7, 6, 102, 54, 30, 54, 103, 0),
    "l": (14, 12, 12, 12, 12, 12, 30, 0),
    "m": (0, 0, 51, 127, 127, 107, 99, 0),
    "n": (0, 0, 31, 51, 51, 51, 51, 0),
    "o": (0, 0, 30, 51, 51, 51, 30, 0),
    "p": (0, 0, 59, 102, 102, 62, 6, 15),
    "q": (0, 0, 110, 51, 51, 62, 48, 120),
    "r": (0, 0, 59, 110, 102, 6, 15, 0),
    "s": (0, 0, 62, 3, 30, 48, 31, 0),
    "t": (8, 12, 62, 12, 12, 44, 24, 0),
    "u": (0, 0, 51, 51, 51, 51, 110, 0),
    "v": (0, 0, 51, 51, 51, 30, 12, 0),
    "w": (0, 0, 99, 107, 127, 127, 54, 0),
    "x": (0, 0, 99, 54, 28, 54, 99, 0),
    "y": (0, 0, 51, 51, 51, 62, 48, 31),
    "z": (0, 0, 63, 25, 12, 38, 63, 0),
    "{": (56, 12, 12, 7, 12, 12, 56, 0),
    "|": (24, 24, 24, 0, 24, 24, 24, 0),
    "}": (7, 12, 12, 56, 12, 12, 7, 0),
    "~": (110, 59, 0, 0, 0, 0, 0, 0),
}


@dataclass
class Rect:
    x: float
    y: float
    w: float
    h: float

    def contains(self, px: float, py: float) -> bool:
        return self.x <= px <= self.x + self.w and self.y <= py <= self.y + self.h


@dataclass
class HUDControl:
    id: str
    rect: Rect
    kind: str
    label: str
    value: Optional[float | bool | str] = None
    min: Optional[float] = None
    max: Optional[float] = None
    options: Optional[Sequence[str]] = None
    data: Optional[Sequence[float]] = None


@dataclass
class HudEvent:
    id: str
    kind: str
    value: Optional[float | str | bool] = None


def _clamp(value: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, value))


class BitmapFont:
    """Tiny bitmap font renderer suitable for HUD annotations."""

    def __init__(self, ctx: moderngl.Context):
        self.ctx = ctx
        self._glyph_w = 8
        self._glyph_h = 8
        self._chars = sorted(FONT_8x8.keys())
        self._atlas, self._uv = self._build_atlas()

        self._program = ctx.program(
            vertex_shader="""
                #version 330
                in vec2 in_pos;
                in vec2 in_uv;
                uniform vec2 u_screen;
                out vec2 v_uv;
                void main() {
                    vec2 clip = vec2(
                        (in_pos.x / u_screen.x) * 2.0 - 1.0,
                        (in_pos.y / u_screen.y) * 2.0 - 1.0
                    );
                    gl_Position = vec4(clip, 0.0, 1.0);
                    v_uv = in_uv;
                }
            """,
            fragment_shader="""
                #version 330
                uniform sampler2D u_font;
                uniform vec4 u_color;
                in vec2 v_uv;
                out vec4 fragColor;
                void main() {
                    float alpha = texture(u_font, v_uv).r;
                    fragColor = vec4(u_color.rgb, u_color.a * alpha);
                }
            """,
        )
        self._buffer = ctx.buffer(reserve=256)
        self._build_vao()
        self._texture = ctx.texture(self._atlas.shape[::-1], 1, data=self._atlas.tobytes())
        self._texture.filter = (moderngl.NEAREST, moderngl.NEAREST)
        self._screen = (1.0, 1.0)

    def set_screen(self, width: float, height: float) -> None:
        self._screen = (max(width, 1.0), max(height, 1.0))

    def _build_vao(self) -> None:
        self._vao = self.ctx.vertex_array(
            self._program,
            [(self._buffer, "2f 2f", "in_pos", "in_uv")],
        )

    def draw_text(self, text: str, x: float, y: float, color: Tuple[float, float, float, float], size: float = 14.0) -> None:
        if not text:
            return
        scale = size / self._glyph_h
        cursor_x = x
        cursor_y = y
        vertices: List[float] = []
        for char in text:
            if char == "\n":
                cursor_x = x
                cursor_y += size + 4.0
                continue
            glyph_uv = self._uv.get(char)
            if glyph_uv is None:
                glyph_uv = self._uv.get("?")
                if glyph_uv is None:
                    continue
            w = self._glyph_w * scale
            h = self._glyph_h * scale
            x0 = cursor_x
            y0 = cursor_y
            x1 = x0 + w
            y1 = y0 + h
            u0, v0, u1, v1 = glyph_uv
            vertices.extend(
                [
                    x0,
                    y0,
                    u0,
                    v1,
                    x1,
                    y0,
                    u1,
                    v1,
                    x0,
                    y1,
                    u0,
                    v0,
                    x0,
                    y1,
                    u0,
                    v0,
                    x1,
                    y0,
                    u1,
                    v1,
                    x1,
                    y1,
                    u1,
                    v0,
                ]
            )
            cursor_x += w + scale
        if not vertices:
            return
        data = np.array(vertices, dtype="f4")
        self._ensure_capacity(data.nbytes)
        self._buffer.write(data.tobytes())
        self._program["u_screen"].value = self._screen
        self._program["u_color"].value = color
        self._texture.use(location=0)
        self._program["u_font"] = 0
        self._vao.render()

    def _build_atlas(self) -> Tuple[np.ndarray, Dict[str, Tuple[float, float, float, float]]]:
        width = self._glyph_w * len(self._chars)
        atlas = np.zeros((self._glyph_h, width), dtype=np.uint8)
        lookup: Dict[str, Tuple[float, float, float, float]] = {}
        for idx, char in enumerate(self._chars):
            cols = FONT_8x8[char]
            x_offset = idx * self._glyph_w
            for y, row in enumerate(cols):
                for bit in range(self._glyph_w):
                    if row & (1 << bit):
                        column = x_offset + (self._glyph_w - 1 - bit)
                        atlas[self._glyph_h - 1 - y, column] = 255
            u0 = x_offset / width
            u1 = (x_offset + self._glyph_w) / width
            lookup[char] = (u0, 0.0, u1, 1.0)
        return atlas, lookup

    def _ensure_capacity(self, required: int) -> None:
        if required <= self._buffer.size:
            return
        self._buffer.release()
        self._buffer = self.ctx.buffer(reserve=required)
        self._build_vao()


class HelixHUD:
    """Minimal immediate-mode HUD with bespoke controls."""

    PANEL_COLOR = (0.07, 0.09, 0.12, 0.9)
    BUTTON_COLOR = (0.22, 0.36, 0.56, 0.9)
    BUTTON_ACTIVE_COLOR = (0.12, 0.56, 0.38, 0.95)
    SLIDER_BG = (0.18, 0.2, 0.26, 0.9)
    SLIDER_FG = (0.93, 0.55, 0.2, 0.95)
    TEXT_COLOR = (0.88, 0.9, 0.94, 0.95)
    WARN_COLOR = (0.93, 0.42, 0.3, 0.95)

    def __init__(self, ctx: moderngl.Context):
        self.ctx = ctx
        self.controls: List[HUDControl] = []
        self._values: Dict[str, float | str] = {"grb2_scale": 1.0}
        self._program = ctx.program(
            vertex_shader="""
                #version 330
                in vec2 in_pos;
                uniform vec2 u_screen;
                void main() {
                    vec2 clip = vec2(
                        (in_pos.x / u_screen.x) * 2.0 - 1.0,
                        (in_pos.y / u_screen.y) * 2.0 - 1.0
                    );
                    gl_Position = vec4(clip, 0.0, 1.0);
                }
            """,
            fragment_shader="""
                #version 330
                uniform vec4 u_color;
                out vec4 fragColor;
                void main() {
                    fragColor = u_color;
                }
            """,
        )
        self._buffer = ctx.buffer(reserve=256)
        self._rebuild_vao()
        self._font = BitmapFont(ctx)
        self._hover: Optional[str] = None
        self._mouse_down_prev = False
        self._screen = (1280.0, 720.0)
        self._open_dropdown: Optional[str] = None
        self._dropdown_hover_index: Optional[int] = None

    def build_controls(
        self,
        *,
        width: int,
        height: int,
        metrics: Dict[str, float],
        metric_history: Dict[str, Sequence[float]],
        slice_name: str,
        variant: Optional[str],
        variant_options: Sequence[str],
        paused: bool,
        interactive: bool,
        bundle_mode: bool,
        bundle_playing: bool,
    ) -> None:
        self._screen = (max(width, 1), max(height, 1))
        self._font.set_screen(*self._screen)
        panel = Rect(24, 24, 300, 220)
        controls: List[HUDControl] = [HUDControl("panel", panel, "panel", "")]
        cursor_y = panel.y + panel.h - 50
        control_height = 34
        if interactive:
            button_rect = Rect(panel.x + 16, cursor_y, panel.w - 32, control_height)
            label = "RESUME" if paused else "PAUSE"
            controls.append(HUDControl("pause", button_rect, "button", label, value=paused))
            cursor_y -= control_height + 10
            slider_rect = Rect(panel.x + 16, cursor_y, panel.w - 32, control_height)
            slider_value = float(metrics.get("grb2_scale", self._values.get("grb2_scale", 1.0)))
            self._values["grb2_scale"] = slider_value
            slider_label = f"GRB2 SCALE {slider_value:0.2f}"
            controls.append(
                HUDControl(
                    "grb2_scale",
                    slider_rect,
                    "slider",
                    slider_label,
                    value=slider_value,
                    min=0.0,
                    max=3.0,
                )
            )
            cursor_y -= control_height + 10
            if variant_options:
                dropdown_rect = Rect(panel.x + 16, cursor_y, panel.w - 32, control_height)
                current_variant = variant or variant_options[0]
                self._values["variant"] = current_variant
                dropdown_label = f"VARIANT {current_variant.upper()} ▾"
                controls.append(
                    HUDControl(
                        "variant",
                        dropdown_rect,
                        "dropdown",
                        dropdown_label,
                        value=current_variant,
                        options=list(variant_options),
                    )
                )
                cursor_y -= control_height + 10
        elif bundle_mode:
            play_rect = Rect(panel.x + 16, cursor_y, panel.w - 32, control_height)
            play_label = "PAUSE" if bundle_playing else "PLAY"
            controls.append(
                HUDControl("bundle_toggle", play_rect, "button", play_label, value=bundle_playing)
            )
            cursor_y -= control_height + 10
            step_rect = Rect(panel.x + 16, cursor_y, panel.w - 32, control_height)
            controls.append(HUDControl("bundle_step", step_rect, "button", "STEP", value=False))
            cursor_y -= control_height + 10

        stats_rect = Rect(panel.x + 16, panel.y + 16, panel.w - 32, control_height + 20)
        perk_val = metrics.get("pERK", 0.0)
        cells_val = metrics.get("cell_count", 0.0)
        sim_time = metrics.get("time", 0.0)
        lines = [
            f"SLICE {slice_name.upper()}",
            f"VAR {variant.upper() if variant else '-'}",
            f"PERK {perk_val:0.2f}",
            f"CELL {int(cells_val):d}",
            f"TIME {sim_time:0.2f}s",
        ]
        controls.append(HUDControl("metrics", stats_rect, "label", "\n".join(lines)))

        plot_rect = Rect(panel.x, panel.y + panel.h + 20, panel.w, 120)
        controls.append(
            HUDControl(
                "perk_plot",
                plot_rect,
                "plot",
                "PERK",
                data=list(metric_history.get("pERK", [])),
            )
        )
        self.controls = controls

    def handle_mouse(self, mouse_x: float, mouse_y: float, mouse_down: bool, *, interactive: bool, bundle_mode: bool) -> List[HudEvent]:
        if not (interactive or bundle_mode):
            self._mouse_down_prev = mouse_down
            self._hover = None
            self._open_dropdown = None
            self._dropdown_hover_index = None
            return []
        events: List[HudEvent] = []
        self._hover = None

        if self._open_dropdown is not None:
            control = next((c for c in self.controls if c.id == self._open_dropdown), None)
            if control is None or not control.options:
                self._open_dropdown = None
                self._dropdown_hover_index = None
            else:
                menu_rect, item_height = self._dropdown_menu_rect(control)
                self._dropdown_hover_index = None
                if menu_rect.contains(mouse_x, mouse_y):
                    rel_y = mouse_y - menu_rect.y
                    idx = int(rel_y // max(item_height, 1e-3))
                    if 0 <= idx < len(control.options or []):
                        self._dropdown_hover_index = idx
                if mouse_down and not self._mouse_down_prev:
                    if menu_rect.contains(mouse_x, mouse_y) and self._dropdown_hover_index is not None:
                        selected = control.options[self._dropdown_hover_index]
                        self._values[control.id] = selected
                        events.append(HudEvent(control.id, "change", value=selected))
                    self._open_dropdown = None
                    self._dropdown_hover_index = None
            self._mouse_down_prev = mouse_down
            return events

        for control in self.controls:
            if control.kind not in {"button", "slider", "dropdown"}:
                continue
            if not control.rect.contains(mouse_x, mouse_y):
                continue
            self._hover = control.id
            if control.kind == "button":
                if mouse_down and not self._mouse_down_prev:
                    events.append(HudEvent(control.id, "click"))
            elif control.kind == "slider" and control.min is not None and control.max is not None and interactive:
                if mouse_down:
                    ratio = _clamp((mouse_x - control.rect.x) / max(control.rect.w, 1e-6), 0.0, 1.0)
                    value = control.min + ratio * (control.max - control.min)
                    self._values[control.id] = value
                    events.append(HudEvent(control.id, "change", value=value))
            elif control.kind == "dropdown" and control.options and interactive:
                if mouse_down and not self._mouse_down_prev:
                    if self._open_dropdown == control.id:
                        self._open_dropdown = None
                        self._dropdown_hover_index = None
                    else:
                        self._open_dropdown = control.id
                        self._dropdown_hover_index = None
            break
        self._mouse_down_prev = mouse_down
        return events

    def draw(self) -> None:
        for control in self.controls:
            if control.kind == "panel":
                self._draw_rect(control.rect, self.PANEL_COLOR)
            elif control.kind == "button":
                color = self.BUTTON_ACTIVE_COLOR if control.value else self.BUTTON_COLOR
                if self._hover == control.id:
                    color = (
                        min(color[0] + 0.1, 1.0),
                        min(color[1] + 0.1, 1.0),
                        min(color[2] + 0.1, 1.0),
                        color[3],
                    )
                self._draw_rect(control.rect, color)
                self._draw_label(control)
            elif control.kind == "slider":
                self._draw_rect(control.rect, self.SLIDER_BG)
                fill = float(self._values.get(control.id, control.value or 0.0))
                min_v = control.min or 0.0
                max_v = control.max or 1.0
                ratio = _clamp((fill - min_v) / (max_v - min_v), 0.0, 1.0)
                fill_rect = Rect(control.rect.x, control.rect.y, control.rect.w * ratio, control.rect.h)
                self._draw_rect(fill_rect, self.SLIDER_FG)
                self._draw_label(control)
            elif control.kind == "dropdown":
                base_color = self.BUTTON_COLOR
                if self._hover == control.id:
                    base_color = (
                        min(base_color[0] + 0.1, 1.0),
                        min(base_color[1] + 0.1, 1.0),
                        min(base_color[2] + 0.1, 1.0),
                        base_color[3],
                    )
                self._draw_rect(control.rect, base_color)
                self._draw_label(control)
            elif control.kind == "label":
                self._draw_label(control, align_left=True, color=self.TEXT_COLOR, size=13)
            elif control.kind == "plot":
                self._draw_rect(control.rect, (0.05, 0.07, 0.1, 0.85))
                self._draw_plot(control)
        if self._open_dropdown is not None:
            control = next((c for c in self.controls if c.id == self._open_dropdown), None)
            if control is not None and control.options:
                self._draw_dropdown_menu(control)

    def _draw_rect(self, rect: Rect, color: Tuple[float, float, float, float]) -> None:
        x0, y0 = rect.x, rect.y
        x1, y1 = rect.x + rect.w, rect.y + rect.h
        vertices = np.array(
            [
                x0,
                y0,
                x1,
                y0,
                x0,
                y1,
                x0,
                y1,
                x1,
                y0,
                x1,
                y1,
            ],
            dtype="f4",
        )
        self._ensure_buffer(vertices.nbytes)
        self._buffer.write(vertices.tobytes())
        self._program["u_screen"].value = self._screen
        self._program["u_color"].value = color
        self._vao.render()

    def _draw_label(
        self,
        control: HUDControl,
        *,
        align_left: bool = False,
        color: Tuple[float, float, float, float] = TEXT_COLOR,
        size: float = 15.0,
    ) -> None:
        text = control.label.upper()
        if not text:
            return
        padding = 10
        if align_left:
            x = control.rect.x + 6
        else:
            x = control.rect.x + padding
        y = control.rect.y + 8
        self._font.draw_text(text, x, y, color, size=size)

    def _draw_plot(self, control: HUDControl) -> None:
        values = list(control.data or [])
        if not values:
            return
        v_min = min(values)
        v_max = max(values)
        if abs(v_max - v_min) < 1e-9:
            v_max = v_min + 1.0
        rect = control.rect
        step = rect.w / max(len(values) - 1, 1)
        points = []
        for idx, value in enumerate(values[-256:]):
            ratio = (value - v_min) / (v_max - v_min)
            x = rect.x + idx * step
            y = rect.y + ratio * rect.h
            points.append((x, y))
        if len(points) < 2:
            return
        flat = np.array(points, dtype="f4").ravel()
        self._ensure_buffer(flat.nbytes)
        self._buffer.write(flat.tobytes())
        self._program["u_screen"].value = self._screen
        self._program["u_color"].value = self.SLIDER_FG
        self.ctx.line_width = 2
        self._vao.render(mode=moderngl.LINE_STRIP)

    def _dropdown_menu_rect(self, control: HUDControl) -> tuple[Rect, float]:
        item_height = control.rect.h * 0.8
        spacing = 4.0
        total_height = len(control.options or []) * (item_height + spacing) - spacing
        x = control.rect.x
        y = control.rect.y + control.rect.h + 6.0
        return Rect(x, y, control.rect.w, total_height), item_height

    def _draw_dropdown_menu(self, control: HUDControl) -> None:
        options = control.options or []
        if not options:
            return
        menu_rect, item_height = self._dropdown_menu_rect(control)
        self._draw_rect(Rect(0.0, 0.0, self._screen[0], self._screen[1]), (0.0, 0.0, 0.0, 0.25))
        self._draw_rect(menu_rect, (0.05, 0.07, 0.10, 0.96))
        current = str(self._values.get(control.id, options[0]))
        for idx, opt in enumerate(options):
            y = menu_rect.y + idx * (item_height + 4.0)
            item_rect = Rect(menu_rect.x + 2.0, y, menu_rect.w - 4.0, item_height)
            color = self.BUTTON_ACTIVE_COLOR if opt == current else self.BUTTON_COLOR
            if self._dropdown_hover_index == idx:
                color = (
                    min(color[0] + 0.1, 1.0),
                    min(color[1] + 0.1, 1.0),
                    min(color[2] + 0.1, 1.0),
                    color[3],
                )
            self._draw_rect(item_rect, color)
            fake = HUDControl(
                id=f"{control.id}_opt_{idx}",
                rect=item_rect,
                kind="label",
                label=opt.upper(),
            )
            self._draw_label(fake, align_left=True, color=self.TEXT_COLOR, size=14.0)

    def _ensure_buffer(self, required: int) -> None:
        if required <= self._buffer.size:
            return
        self._buffer.release()
        self._buffer = self.ctx.buffer(reserve=required)
        self._rebuild_vao()

    def _rebuild_vao(self) -> None:
        self._vao = self.ctx.simple_vertex_array(self._program, self._buffer, "in_pos")
