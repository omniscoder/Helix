"""Qt GUI entrypoint for Helix realtime simulations."""
from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import QUrl
from PySide6.QtWidgets import QApplication, QMainWindow
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWebEngineWidgets import QWebEngineView

from .bridge import SimulationBridge


def run_gui() -> None:
    """Launch the PySide6 GUI with the realtime simulator."""
    app = QApplication.instance()
    owns_app = False
    if app is None:
        app = QApplication(sys.argv)
        owns_app = True

    window = QMainWindow()
    window.setWindowTitle("Helix – Genome Editing Digital Twin")

    view = QWebEngineView(window)
    window.setCentralWidget(view)

    channel = QWebChannel()
    bridge = SimulationBridge(parent=view)
    channel.registerObject("HelixBridge", bridge)
    view.page().setWebChannel(channel)

    html_path = Path(__file__).resolve().parent / "resources" / "qt_realtime.html"
    view.load(QUrl.fromLocalFile(str(html_path)))

    window.resize(1200, 800)
    window.show()

    if owns_app:
        sys.exit(app.exec())
    else:  # pragma: no cover - embedding scenarios
        app.exec()
