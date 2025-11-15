"""Application-wide Qt theme helpers."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QPalette
from PySide6.QtWidgets import QApplication

_BASE_BG = "#0f172a"
_ELEVATED_BG = "#1e293b"
_TEXT = "#e2e8f0"
_SUBDUED_TEXT = "#94a3b8"
_ACCENT = "#22d3ee"
_ACCENT_DARK = "#0ea5a4"


def _build_palette() -> QPalette:
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor(_BASE_BG))
    palette.setColor(QPalette.WindowText, QColor(_TEXT))
    palette.setColor(QPalette.Base, QColor(_ELEVATED_BG))
    palette.setColor(QPalette.AlternateBase, QColor("#152238"))
    palette.setColor(QPalette.ToolTipBase, QColor(_ELEVATED_BG))
    palette.setColor(QPalette.ToolTipText, QColor(_TEXT))
    palette.setColor(QPalette.Text, QColor(_TEXT))
    palette.setColor(QPalette.Button, QColor(_ELEVATED_BG))
    palette.setColor(QPalette.ButtonText, QColor(_TEXT))
    palette.setColor(QPalette.Highlight, QColor(_ACCENT))
    palette.setColor(QPalette.HighlightedText, QColor("#0c111d"))
    palette.setColor(QPalette.PlaceholderText, QColor(_SUBDUED_TEXT))
    palette.setColor(QPalette.Disabled, QPalette.Text, QColor("#4b576d"))
    palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor("#4b576d"))
    palette.setColor(QPalette.Disabled, QPalette.WindowText, QColor("#4b576d"))
    return palette


_GLOBAL_STYLESHEET = f"""
* {{
    font-family: "Inter", "Segoe UI", "Helvetica Neue", Arial, sans-serif;
    font-size: 12px;
}}
QWidget {{
    color: {_TEXT};
    background-color: {_BASE_BG};
}}
QLabel#phaseLabel {{
    font-weight: 600;
}}
QPushButton {{
    background-color: {_ELEVATED_BG};
    border: 1px solid #243447;
    border-radius: 4px;
    padding: 6px 12px;
}}
QPushButton:hover {{
    border-color: {_ACCENT};
}}
QPushButton:pressed {{
    background-color: {_ACCENT_DARK};
    border-color: {_ACCENT};
}}
QPushButton:disabled {{
    background-color: #1b2533;
    color: {_SUBDUED_TEXT};
    border-color: #1b2533;
}}
QLineEdit, QPlainTextEdit, QListWidget, QSpinBox, QDoubleSpinBox {{
    background-color: #{_BASE_BG.strip('#')}ee;
    border: 1px solid #1f2933;
    border-radius: 4px;
    padding: 4px;
}}
QTabWidget::pane {{
    border: 1px solid #1f2a37;
    background: {_BASE_BG};
}}
QTabBar::tab {{
    background: {_ELEVATED_BG};
    border: 1px solid #1f2a37;
    border-bottom: none;
    padding: 6px 12px;
}}
QTabBar::tab:selected {{
    background: {_ACCENT_DARK};
}}
QScrollBar:vertical {{
    background: {_BASE_BG};
    width: 12px;
    margin: 2px 0 2px 0;
}}
QScrollBar::handle:vertical {{
    background: {_ACCENT_DARK};
    min-height: 20px;
    border-radius: 6px;
}}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
    height: 0;
}}
"""


def apply_helix_theme(app: QApplication) -> None:
    """Apply the Helix Qt palette/stylesheet if not already enabled."""

    if app.property("helix_theme_applied"):
        return
    app.setPalette(_build_palette())
    app.setStyleSheet(_GLOBAL_STYLESHEET)
    app.setProperty("helix_theme_applied", True)
