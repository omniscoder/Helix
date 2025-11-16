from __future__ import annotations

from PySide6.QtWidgets import QApplication

from helix.studio.session import SessionModel
from helix.studio.gl_panel import GLViewport
from helix.studio.workflow_panel import Workflow2p5DPanel


def _get_app() -> QApplication:
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    return app


def test_session_snapshot_drives_panels() -> None:
    app = _get_app()
    session = SessionModel()
    gl = GLViewport(session)
    workflow = Workflow2p5DPanel(session)
    payload = {
        "version": 1,
        "state": {
            "run_kind": "CRISPR",
            "run_id": 1,
            "viz_dirty": True,
            "config": {
                "run_config": {"draws": 10},
                "workflow_view": {
                    "inst_events": [[0.0, 1.0, 0.0, 0.5, 0.3, 0.0, 0.0]],
                    "heat_values": [0.1, 0.2, 0.3],
                },
            },
            "outcomes": [
                {"label": "intended", "probability": 0.6, "diff": {"kind": "none"}},
                {"label": "off", "probability": 0.4, "diff": {"kind": "deletion"}},
            ],
            "viz_spec_payload": {
                "sequence": "ACGT" * 8,
                "pam_index": 5,
                "guide_range": [2, 22],
                "edit_type": "crispr.sim",
                "edit_events": [],
                "metadata": {},
            },
            "genome": {"chr": "ACGT" * 8},
        },
    }
    session.load_payload(payload)
    app.processEvents()
    assert not gl._overlay.isHidden()
    assert "Outcomes" in gl._overlay.text()
    assert workflow._summary.toPlainText()
