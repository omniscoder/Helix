from __future__ import annotations

import os
from collections.abc import Mapping
from typing import Any

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QLabel, QTextEdit, QVBoxLayout, QWidget

from helix.gui.modern.qt import HelixModernWidget

from .gl_panel import _workflow_from_serializable
from .session import ExperimentState, SessionModel, RunKind
from .workflow_stats import summarize_workflow_meta
from .pcr_workflow import summarize_pcr_workflow_meta


class Workflow2p5DPanel(QWidget):
    def __init__(self, session: SessionModel, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._session = session
        self._viewer: HelixModernWidget | None = None
        self._viewer_container: QWidget | None = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(8)

        self._status = QLabel("Workflow will appear here after simulation.", self)
        self._status.setWordWrap(True)
        layout.addWidget(self._status)

        app = QApplication.instance()
        platform = app.platformName().lower() if app else ""
        is_wayland = "wayland" in platform or bool(os.environ.get("WAYLAND_DISPLAY"))
        self._embed_supported = bool(platform) and not is_wayland

        if self._embed_supported:
            self._viewer = HelixModernWidget()
            self._viewer_container = QWidget.createWindowContainer(self._viewer, self)
            self._viewer_container.setMinimumHeight(320)
            layout.addWidget(self._viewer_container, stretch=1)
            self._viewer_container.hide()
        else:
            fallback = QLabel(
                "Workflow preview requires an X11/OpenGL backend. "
                "Launch Helix Studio with QT_QPA_PLATFORM=xcb or view the workflow in SimBuilder.",
                self,
            )
            fallback.setWordWrap(True)
            fallback.setAlignment(Qt.AlignLeft | Qt.AlignTop)
            layout.addWidget(fallback)

        self._summary = QTextEdit(self)
        self._summary.setReadOnly(True)
        self._summary.setPlaceholderText("Workflow 2.5D geometry details will appear here.")
        layout.addWidget(self._summary)

        self._session.stateChanged.connect(self._on_state_changed)
        self._on_state_changed(self._session.state)

    def _on_state_changed(self, state: ExperimentState) -> None:
        config = state.config if isinstance(state.config, Mapping) else {}
        workflow_meta = self._workflow_meta_from_state(state, config)

        if state.run_kind is RunKind.PCR and state.pcr_result is not None:
            lines: list[str] = []
            if isinstance(workflow_meta, Mapping) and workflow_meta.get("kind") == "PCR":
                lines = summarize_pcr_workflow_meta(workflow_meta)
            else:
                lines = [
                    f"Amplicon length: {state.pcr_result.amplicon_length} bp",
                    f"Final mass: {state.pcr_result.final_amplicon_mass_ng:.2f} ng",
                    f"Final mutation rate: {state.pcr_result.final_mutation_rate:.4f}",
                ]
            self._set_viewer_buffers(None)
            self._status.setText("PCR workflow ready.")
            self._summary.setPlainText("\n".join(lines))
            return

        if not isinstance(workflow_meta, Mapping):
            self._set_viewer_buffers(None)
            self._status.setText("Workflow will appear here after simulation.")
            self._summary.setPlainText("Run a simulation to populate workflow geometry.")
            return

        buffers = self._extract_workflow_buffers(workflow_meta)
        self._set_viewer_buffers(buffers)

        lines = summarize_workflow_meta(workflow_meta)

        if not lines:
            self._summary.setPlainText("Workflow metadata present but no geometry stats available yet.")
        else:
            self._summary.setPlainText("\n".join(lines))

        sim_kind = config.get("sim_type", "unknown")
        self._status.setText(f"{sim_kind} workflow geometry ready.")

    def _workflow_meta_from_state(
        self,
        state: ExperimentState,
        config: Mapping[str, Any],
    ) -> Mapping[str, Any] | None:
        workflow_meta = config.get("workflow_view")
        if isinstance(workflow_meta, Mapping):
            return workflow_meta
        spec_payload = state.viz_spec_payload
        if isinstance(spec_payload, Mapping):
            metadata = spec_payload.get("metadata")
            if isinstance(metadata, Mapping):
                workflow_meta = metadata.get("workflow_view")
                if isinstance(workflow_meta, Mapping):
                    return workflow_meta
        return None

    def _extract_workflow_buffers(self, workflow_meta: Mapping[str, Any]) -> dict[str, Any] | None:
        if workflow_meta.get("kind") == "PCR":
            return None
        try:
            buffers = _workflow_from_serializable(workflow_meta)
            return buffers if buffers else None
        except Exception:
            return None

    def _set_viewer_buffers(self, buffers: dict[str, Any] | None) -> None:
        if self._viewer is not None and self._viewer_container is not None:
            if buffers:
                self._viewer_container.show()
                self._viewer.set_workflow_view(buffers)
            else:
                self._viewer.set_workflow_view(None)
                self._viewer_container.hide()
