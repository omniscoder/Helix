from __future__ import annotations

import json
from pathlib import Path

from typing import Mapping

from PySide6.QtCore import Qt
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QDockWidget,
    QFileDialog,
    QMainWindow,
    QMessageBox,
    QSplitter,
    QTabWidget,
    QTextEdit,
    QWidget,
)

from .gl_panel import GLViewport
from .session import SessionModel
from .session_inspector import SessionInspector
from .workflow_panel import Workflow2p5DPanel
from .run_history import RunHistoryWidget
from .run_compare import RunCompareWidget
from helix.gui.simbuilder.app import create_simbuilder_widget


class HelixStudioMainWindow(QMainWindow):
    def __init__(self, session: SessionModel | None = None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.session = session or SessionModel(self)
        self.setWindowTitle("Helix Studio")
        self._create_menus()

        # --- central: GL viewport (and later 2.5D / rails tabs) ---
        central_split = QSplitter(Qt.Vertical, self)

        tab = QTabWidget()
        self.gl_view = GLViewport(self.session)
        tab.addTab(self.gl_view, "Realtime / DAG")
        self.workflow_panel = Workflow2p5DPanel(self.session)
        tab.addTab(self.workflow_panel, "2.5D Workflow")

        central_split.addWidget(tab)

        # bottom log area
        self.log_widget = QTextEdit()
        self.log_widget.setReadOnly(True)
        central_split.addWidget(self.log_widget)

        central_split.setStretchFactor(0, 3)
        central_split.setStretchFactor(1, 1)
        self.setCentralWidget(central_split)

        self.session.logMessage.connect(self.log_widget.append)
        self.session.statusChanged.connect(self.statusBar().showMessage)
        self.session.errorOccurred.connect(self._show_error)

        # --- left dock: SimBuilder as a panel ---
        sim_dock = QDockWidget("Sim Builder", self)
        sim_dock.setObjectName("SimBuilderDock")
        self.sim_panel = create_simbuilder_widget(self.session, enable_gl_viewer=False)
        sim_dock.setWidget(self.sim_panel)
        self.addDockWidget(Qt.LeftDockWidgetArea, sim_dock)

        inspector_dock = QDockWidget("Session Inspector", self)
        inspector_dock.setObjectName("SessionInspectorDock")
        inspector_dock.setWidget(SessionInspector(self.session))
        self.addDockWidget(Qt.RightDockWidgetArea, inspector_dock)

        history_dock = QDockWidget("Run History", self)
        history_dock.setObjectName("RunHistoryDock")
        self.history_widget = RunHistoryWidget(self.session, on_compare=self._on_compare_runs)
        history_dock.setWidget(self.history_widget)
        self.addDockWidget(Qt.RightDockWidgetArea, history_dock)

        compare_dock = QDockWidget("Compare Runs", self)
        compare_dock.setObjectName("RunCompareDock")
        self.compare_widget = RunCompareWidget(self)
        compare_dock.setWidget(self.compare_widget)
        self.addDockWidget(Qt.RightDockWidgetArea, compare_dock)

        self._init_layout()
        self.statusBar().showMessage("Ready")
        self._setup_shortcuts()

    def _init_layout(self) -> None:
        self.resize(1800, 1000)

    def _create_menus(self) -> None:
        menu = self.menuBar().addMenu("&File")
        save_action = menu.addAction("Save Session…")
        save_action.triggered.connect(self._save_session)
        save_action.setShortcut(QKeySequence("Ctrl+S"))
        save_as_action = menu.addAction("Save Session As…")
        save_as_action.triggered.connect(self._save_session)
        save_as_action.setShortcut(QKeySequence("Ctrl+Shift+S"))
        load_action = menu.addAction("Load Session…")
        load_action.triggered.connect(self._load_session)
        load_action.setShortcut(QKeySequence("Ctrl+O"))

        help_menu = self.menuBar().addMenu("&Help")
        shortcuts_action = help_menu.addAction("Keyboard Shortcuts")
        shortcuts_action.triggered.connect(self._show_shortcuts_dialog)

    def _save_session(self) -> None:
        path_str, _ = QFileDialog.getSaveFileName(
            self,
            "Save Helix Session",
            filter="Helix Session (*.helix.json);;JSON (*.json)",
        )
        if not path_str:
            return
        path = Path(path_str)
        payload = self.session.to_payload()
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        self.session.set_status(f"Saved session to {path}")

    def _load_session(self) -> None:
        path_str, _ = QFileDialog.getOpenFileName(
            self,
            "Load Helix Session",
            filter="Helix Session (*.helix.json);;JSON (*.json)",
        )
        if not path_str:
            return
        path = Path(path_str)
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            self._show_error(f"Failed to load session: {exc}")
            return
        try:
            payload["viz_dirty"] = True
            self.session.load_payload(payload)
        except Exception as exc:
            self._show_error(f"Failed to restore session: {exc}")
            return
        self.session.set_status(f"Loaded session from {path}")

    def _show_error(self, message: str) -> None:
        QMessageBox.critical(self, "Helix Studio", message)

    def _setup_shortcuts(self) -> None:
        self._shortcuts: list[QShortcut] = []
        rerun = QShortcut(QKeySequence("Ctrl+R"), self)
        rerun.activated.connect(self._rerun_active_sim)
        self._shortcuts.append(rerun)
        for i in range(1, 10):
            shortcut = QShortcut(QKeySequence(f"Ctrl+{i}"), self)
            shortcut.activated.connect(lambda pos=i - 1: self._select_history_position(pos))
            self._shortcuts.append(shortcut)
        zero_shortcut = QShortcut(QKeySequence("Ctrl+0"), self)
        zero_shortcut.activated.connect(lambda: self._select_history_position(9))
        self._shortcuts.append(zero_shortcut)

    def _rerun_active_sim(self) -> None:
        if hasattr(self.sim_panel, "run_active_tab"):
            self.sim_panel.run_active_tab()

    def _select_history_position(self, position: int) -> None:
        if hasattr(self, "history_widget"):
            self.history_widget.select_run_by_position(position)

    def _on_compare_runs(self, lhs: Mapping[str, object], rhs: Mapping[str, object]) -> None:
        self.compare_widget.compare_snapshots(lhs, rhs)
        self.session.set_status("Comparing selected runs")

    def _show_shortcuts_dialog(self) -> None:
        text = (
            "Global\n"
            "Ctrl+R – Re-run active SimBuilder tab\n"
            "Ctrl+S / Ctrl+Shift+S – Save session / Save as\n"
            "Ctrl+O – Load session\n\n"
            "Run History\n"
            "Ctrl+F – Focus search\n"
            "Ctrl+Shift+I – Toggle intended filter\n"
            "Ctrl+Shift+K – Cycle run kinds\n"
            "Ctrl+Shift+F – Clear filters\n"
            "Ctrl+1..0 – Jump to run rows"
        )
        QMessageBox.information(self, "Keyboard Shortcuts", text)
