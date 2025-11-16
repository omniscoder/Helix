from __future__ import annotations

from PySide6.QtWidgets import QApplication

from helix.studio.session import RunKind, SessionModel


def _get_app() -> QApplication:
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    return app


def test_session_dirty_cycle() -> None:
    _get_app()
    session = SessionModel()
    session.update(run_kind=RunKind.CRISPR, viz_dirty=True)
    assert session.state.viz_dirty is True
    assert session.state.run_kind == RunKind.CRISPR
    session.mark_viz_clean()
    assert session.state.viz_dirty is False


def test_session_run_recorded_and_roundtrip() -> None:
    _get_app()
    session = SessionModel()
    recorded: list[dict] = []
    session.runRecorded.connect(lambda payload: recorded.append(payload))
    session.update(run_kind=RunKind.CRISPR, viz_dirty=True)
    assert recorded, "runRecorded should emit when viz_dirty run updates occur"
    payload = session.to_payload()
    payload["state"]["viz_dirty"] = True
    session.load_payload(payload)
    assert session.state.run_kind == RunKind.CRISPR
    assert session.state.viz_dirty is True
