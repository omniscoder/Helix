from __future__ import annotations

import sys

from PySide6.QtWidgets import QApplication

from helix.gui.simbuilder.app import configure_opengl_surface_format
from helix.gui.theme import apply_helix_theme
from .main_window import HelixStudioMainWindow
from .session import SessionModel


def main() -> None:
    app = QApplication(sys.argv)
    configure_opengl_surface_format()
    apply_helix_theme(app)
    session = SessionModel()
    window = HelixStudioMainWindow(session)
    window.show()
    sys.exit(app.exec())
