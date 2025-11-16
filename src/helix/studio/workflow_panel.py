from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from PySide6.QtWidgets import QLabel, QTextEdit, QVBoxLayout, QWidget

from .session import ExperimentState, SessionModel
from .workflow_stats import summarize_workflow_meta


class Workflow2p5DPanel(QWidget):
    def __init__(self, session: SessionModel, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._session = session

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(8)

        self._status = QLabel("Workflow will appear here after simulation.", self)
        self._status.setWordWrap(True)
        layout.addWidget(self._status)

        self._summary = QTextEdit(self)
        self._summary.setReadOnly(True)
        self._summary.setPlaceholderText("Workflow 2.5D geometry details will appear here.")
        layout.addWidget(self._summary)

        self._session.stateChanged.connect(self._on_state_changed)
        self._on_state_changed(self._session.state)

    def _on_state_changed(self, state: ExperimentState) -> None:
        config = state.config if isinstance(state.config, Mapping) else {}
        workflow_meta = config.get("workflow_view")
        if not isinstance(workflow_meta, Mapping):
            self._status.setText("Workflow will appear here after simulation.")
            self._summary.setPlainText("Run a CRISPR or Prime simulation to populate workflow geometry.")
            return

        lines = summarize_workflow_meta(workflow_meta)

        if not lines:
            self._summary.setPlainText("Workflow metadata present but no geometry stats available yet.")
        else:
            self._summary.setPlainText("\n".join(lines))

        sim_kind = config.get("sim_type", "unknown")
        self._status.setText(f"{sim_kind} workflow geometry ready.")
