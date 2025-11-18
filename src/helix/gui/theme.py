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
    font-size: 13px;
    letter-spacing: 0.01em;
}}
QWidget {{
    color: {_TEXT};
    background-color: {_BASE_BG};
}}
QLabel#phaseLabel {{
    font-weight: 600;
}}
QLabel#StudioHero {{
    font-size: 16px;
    font-weight: 600;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    color: {_ACCENT};
}}
QWidget#ViewerPlaceholder {{
    border: 1px dashed #233043;
    border-radius: 8px;
    background-color: #0b1221;
}}
QPushButton {{
    background-color: {_ELEVATED_BG};
    border: 1px solid #243447;
    border-radius: 6px;
    padding: 8px 14px;
    font-weight: 500;
}}
QPushButton:hover {{
    border-color: {_ACCENT};
    color: {_ACCENT};
}}
QPushButton:pressed {{
    background-color: {_ACCENT_DARK};
    border-color: {_ACCENT};
    color: #05131b;
}}
QPushButton:disabled {{
    background-color: #1b2533;
    color: {_SUBDUED_TEXT};
    border-color: #1b2533;
}}
QLineEdit, QPlainTextEdit, QListWidget, QSpinBox, QDoubleSpinBox {{
    background-color: rgba(15, 23, 42, 0.85);
    border: 1px solid #253247;
    border-radius: 6px;
    padding: 6px 10px;
}}
QLineEdit:focus, QPlainTextEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus {{
    border-color: {_ACCENT};
}}
QPlainTextEdit#LogPanel {{
    font-family: "JetBrains Mono", "Fira Mono", monospace;
    background-color: #0c1324;
    border-radius: 8px;
}}
QTabWidget::pane {{
    border: 1px solid #1f2a37;
    border-radius: 8px;
    background: {_BASE_BG};
}}
QTabBar::tab {{
    background: {_ELEVATED_BG};
    border: 1px solid #1f2a37;
    border-bottom: none;
    min-width: 120px;
    padding: 8px 16px;
    margin-right: 2px;
    text-transform: uppercase;
    letter-spacing: 0.08em;
}}
QTabBar::tab:selected {{
    background: {_ACCENT_DARK};
    color: #041017;
    font-weight: 600;
}}
QSplitter::handle {{
    background-color: #192133;
    border: 1px solid #101624;
    margin: 0 4px;
}}
QSplitter::handle:horizontal {{
    width: 6px;
}}
QSplitter::handle:vertical {{
    height: 6px;
}}
QTreeWidget, QTableWidget {{
    background-color: rgba(11, 17, 29, 0.9);
    border: 1px solid #212b3d;
    border-radius: 6px;
}}
QHeaderView::section {{
    background-color: {_ELEVATED_BG};
    padding: 6px;
    border: none;
    font-weight: 500;
    text-transform: uppercase;
    letter-spacing: 0.05em;
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
QToolTip {{
    background-color: {_ELEVATED_BG};
    color: {_TEXT};
    border: 1px solid {_ACCENT_DARK};
    padding: 6px;
    border-radius: 6px;
}}
"""


def apply_helix_theme(app: QApplication) -> None:
    """Apply the Helix Qt palette/stylesheet if not already enabled."""

    if app.property("helix_theme_applied"):
        return
    app.setPalette(_build_palette())
    app.setStyleSheet(_GLOBAL_STYLESHEET)
    app.setProperty("helix_theme_applied", True)
