"""2D rail band visualization widget."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from PySide6.QtCore import Qt, QTimer, QRectF
from PySide6.QtGui import QColor, QFont, QPainter, QPen
from PySide6.QtWidgets import QWidget

from helix.gui.modern.spec import EditVisualizationSpec


class RailBandWidget(QWidget):
    """Simple QWidget renderer for CRISPR rail bands."""

    _COLORS = {
        "deletion": QColor("#f97316"),
        "insertion": QColor("#ec4899"),
        "substitution": QColor("#0ea5e9"),
        "no_cut": QColor("#22d3ee"),
        "other": QColor("#a78bfa"),
    }

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumHeight(360)
        self._sequence = ""
        self._bands: list[Mapping[str, Any]] = []
        self._guide_range = (0, 0)
        self._pam_index = 0
        self._cut_index = 0
        self._title = ""
        self._message = "No visualization data."
        self._anim_progress = 1.0
        self._anim_timer = QTimer(self)
        self._anim_timer.setInterval(16)
        self._anim_timer.timeout.connect(self._advance_animation)

    def set_spec(self, spec: EditVisualizationSpec | None) -> None:
        if spec is None or not spec.metadata:
            self._sequence = ""
            self._bands = []
            self._title = ""
            self._message = "No visualization data."
            self._anim_progress = 1.0
            self._anim_timer.stop()
            self.update()
            return
        meta = spec.metadata.get("rail_bands")
        if not isinstance(meta, Mapping):
            self._sequence = ""
            self._bands = []
            self._title = ""
            self._message = "Rail band metadata unavailable."
            self._anim_progress = 1.0
            self._anim_timer.stop()
            self.update()
            return
        self._sequence = str(meta.get("sequence") or spec.sequence or "")
        self._bands = list(meta.get("bands") or [])
        self._guide_range = tuple(meta.get("guide_range") or spec.guide_range)
        self._pam_index = int(meta.get("pam_index", spec.pam_index))
        self._cut_index = int(meta.get("cut_index", self._guide_range[0]))
        self._title = spec.edit_type
        if not self._bands:
            self._message = "No outcomes to display."
            self._anim_progress = 1.0
            self._anim_timer.stop()
        else:
            self._message = ""
            self._start_animation()
        self.update()

    def paintEvent(self, event) -> None:  # type: ignore[override]
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing, True)
        rect = self.rect()
        painter.fillRect(rect, QColor("#050910"))
        if not self._sequence or not self._bands:
            painter.setPen(QColor("#94a3b8"))
            painter.setFont(QFont("Inter", 14))
            painter.drawText(rect, Qt.AlignCenter, self._message or "No data.")
            return

        margin_x = 60.0
        margin_top = 50.0
        seq_len = max(1, len(self._sequence) - 1)
        width = rect.width() - margin_x * 2
        height = rect.height() - margin_top - 60.0

        def x_for(idx: int) -> float:
            idx = max(0, min(len(self._sequence) - 1, idx))
            return margin_x + (idx / max(1, seq_len)) * width

        rail_y = margin_top
        painter.setPen(QPen(QColor("#38bdf8"), 3))
        painter.drawLine(x_for(0), rail_y, x_for(len(self._sequence) - 1), rail_y)
        cut_x = x_for(self._cut_index)

        # Guide region
        g0, g1 = self._guide_range
        painter.setPen(Qt.NoPen)
        guide_color = QColor("#22d3ee")
        guide_color.setAlpha(60)
        painter.setBrush(guide_color)
        painter.drawRect(x_for(g0), rail_y - 12, max(2.0, x_for(g1) - x_for(g0)), 24)

        # PAM marker
        pam_color = QColor("#be123c")
        painter.setBrush(pam_color)
        painter.drawRect(x_for(self._pam_index) - 2, rail_y - 18, 4, 36)

        # Cut line
        painter.setPen(QPen(QColor("#f97316"), 2))
        painter.drawLine(cut_x, rail_y - 30, cut_x, rail_y + height * 0.5)
        self._draw_axis_annotations(painter, QRectF(margin_x, margin_top - 20, width, 20))
        self._draw_probability_scale(painter, QRectF(margin_x, rail_y + height + 10, width, 40))
        painter.setFont(QFont("Inter", 10, QFont.Medium))
        painter.setPen(QColor("#fbbf24"))
        painter.drawText(QRectF(cut_x + 6, rail_y - 48, 140, 18), Qt.AlignLeft | Qt.AlignVCenter, "Cut site (0 bp)")

        self._draw_reference_strip(painter, x_for, rail_y, len(self._sequence))

        painter.setFont(QFont("Inter", 16, QFont.Bold))
        painter.setPen(QColor("#e2e8f0"))
        painter.drawText(rect.adjusted(16, 10, -16, 0), Qt.AlignLeft | Qt.AlignTop, self._title)

        painter.setFont(QFont("Inter", 12))
        band_height = 20
        gap = 18
        for band in self._bands:
            start = int(band.get("start", 0))
            end = int(band.get("end", start + 1))
            label = str(band.get("label") or "")
            prob = float(band.get("probability") or 0.0)
            kind = str(band.get("kind") or "other")
            row = int(band.get("row", 0))
            y = rail_y + 40 + row * (band_height + gap)
            base_start = x_for(start)
            base_end = x_for(end)
            anim_start = cut_x + (base_start - cut_x) * self._anim_progress
            anim_end = cut_x + (base_end - cut_x) * self._anim_progress
            x_left = min(anim_start, anim_end)
            bar_width = max(2.0, abs(anim_end - anim_start))
            band_rect = QRectF(x_left, y, bar_width, band_height)

            color = QColor(self._COLORS.get(kind, self._COLORS["other"]))
            alpha = int(90 + min(prob, 1.0) * 160)
            color.setAlpha(max(60, min(alpha, 255)))

            glow = QColor(color)
            glow.setAlpha(int(glow.alpha() * 0.45))
            painter.setPen(Qt.NoPen)
            painter.setBrush(glow)
            painter.drawRoundedRect(band_rect.adjusted(-2.0, -2.0, 2.0, 2.0), 6, 6)

            painter.setBrush(color)
            painter.drawRoundedRect(band_rect, 4, 4)
            self._draw_band_texture(painter, band_rect, kind, prob)

            painter.setPen(QPen(QColor("#0f172a")))
            painter.drawRoundedRect(band_rect, 4, 4)
            painter.setPen(QColor("#cbd5f5"))
            painter.drawText(
                int(band_rect.left()),
                int(band_rect.bottom()) + 14,
                f"{label} ({prob * 100:.1f}%)",
            )

    def _draw_reference_strip(self, painter: QPainter, x_for, rail_y: float, seq_len: int) -> None:
        painter.save()
        strip_y = rail_y + 26.0
        base_pen = QPen(QColor("#0f172a"), 1)
        painter.setPen(base_pen)
        painter.drawLine(x_for(0), strip_y, x_for(seq_len - 1), strip_y)
        label_pen = QPen(QColor("#475569"), 1)
        font = QFont("Inter", 9)
        painter.setFont(font)
        for idx in range(0, seq_len, 5):
            x = x_for(idx)
            major = idx % 10 == 0
            tick_len = 8 if major else 4
            painter.setPen(base_pen if major else QPen(QColor("#0b1120"), 1))
            painter.drawLine(x, strip_y, x, strip_y + tick_len)
            if major:
                painter.setPen(label_pen)
                offset = idx - self._cut_index
                label = f"{offset:+d}"
                painter.drawText(int(x - 10), int(strip_y + tick_len + 12), label)
        painter.restore()

    def _draw_axis_annotations(self, painter: QPainter, rect: QRectF) -> None:
        painter.save()
        painter.setPen(QColor("#94a3b8"))
        painter.setFont(QFont("Inter", 10))
        painter.drawText(
            rect.adjusted(0, 0, 0, 0),
            Qt.AlignLeft | Qt.AlignVCenter,
            "Position (bp relative to cut)",
        )
        painter.restore()

    def _draw_probability_scale(self, painter: QPainter, rect: QRectF) -> None:
        painter.save()
        painter.setPen(QColor("#94a3b8"))
        painter.setFont(QFont("Inter", 10))
        painter.drawText(rect, Qt.AlignLeft | Qt.AlignTop, "Bar length ∝ probability")
        line_y = rect.top() + 18
        scale_width = min(180.0, rect.width())
        painter.setPen(QPen(QColor("#475569"), 1))
        painter.drawLine(rect.left(), line_y, rect.left() + scale_width, line_y)
        ticks = [0.0, 0.5, 1.0]
        for frac in ticks:
            x = rect.left() + scale_width * frac
            painter.drawLine(x, line_y - 4, x, line_y + 4)
            label = f"{int(frac * 100)}%"
            painter.drawText(QRectF(x - 12, line_y + 6, 40, 14), Qt.AlignLeft | Qt.AlignTop, label)
        painter.restore()

    def _draw_band_texture(self, painter: QPainter, rect: QRectF, kind: str, prob: float) -> None:
        painter.save()
        painter.setClipRect(rect)
        stroke = QColor("#0f172a")
        stroke.setAlpha(int(80 + prob * 100))
        painter.setPen(QPen(stroke, 1))
        if kind == "deletion":
            step = 7
            span = int(rect.width()) + int(rect.height()) + step
            for offset in range(-int(rect.height()), span, step):
                x1 = rect.left() + offset
                painter.drawLine(x1, rect.bottom(), x1 + rect.height(), rect.top())
        elif kind == "insertion":
            step = 6
            x = rect.left()
            while x <= rect.right():
                painter.drawLine(x, rect.top(), x, rect.bottom())
                x += step
        elif kind == "substitution":
            spacing = 6
            y = rect.top()
            while y <= rect.bottom():
                x = rect.left()
                while x <= rect.right():
                    painter.drawPoint(x + 1, y + 1)
                    x += spacing
                y += spacing
        elif kind == "no_cut":
            painter.setPen(QPen(stroke, 0.8))
            painter.drawLine(rect.left(), rect.center().y(), rect.right(), rect.center().y())
        else:
            step = 8
            y = rect.top()
            while y <= rect.bottom():
                painter.drawLine(rect.left(), y, rect.right(), y)
                y += step
        painter.restore()

    def _start_animation(self) -> None:
        self._anim_progress = 0.0
        if not self._anim_timer.isActive():
            self._anim_timer.start()

    def _advance_animation(self) -> None:
        self._anim_progress = min(1.0, self._anim_progress + 0.08)
        if self._anim_progress >= 1.0:
            self._anim_timer.stop()
        self.update()
