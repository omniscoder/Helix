from __future__ import annotations

from PySide6.QtWidgets import QApplication

from helix.studio.run_history import RunHistoryWidget
from helix.studio.session import SessionModel


def _get_app() -> QApplication:
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    return app


def _snapshot_payload(run_kind: str, run_id: int, genome_label: str, intended: int, total: int) -> dict:
    outcomes = []
    for _ in range(intended):
        outcomes.append({"label": "intended", "probability": 0.5, "diff": {"kind": "none"}})
    for _ in range(total - intended):
        outcomes.append({"label": "off", "probability": 0.1, "diff": {"kind": "del"}})
    return {
        "timestamp": f"t{run_id}",
        "state": {
            "run_kind": run_kind,
            "run_id": run_id,
            "viz_dirty": False,
            "config": {"run_config": {"draws": 10}},
            "genome_source": "file",
            "genome_uri": genome_label,
            "genome": {"chr": "ACGT" * 4},
            "outcomes": outcomes,
        },
    }


def test_run_history_filters_and_compare() -> None:
    _get_app()
    session = SessionModel()
    compared: list[tuple[dict, dict]] = []
    widget = RunHistoryWidget(session, on_compare=lambda a, b: compared.append((a, b)))
    snapshots = [
        _snapshot_payload("CRISPR", 1, "chr1.fa", 1, 2),
        _snapshot_payload("PRIME", 2, "chr2.fa", 0, 3),
        _snapshot_payload("CRISPR", 3, "demo", 2, 4),
    ]
    for snap in snapshots:
        widget._add_entry(snap)

    # Apply intended filter
    widget._intended_only.setChecked(True)
    assert widget._tree.topLevelItemCount() == 2, "Only runs with intended outcomes remain (plus pinned)"

    # Pin a run without intended outcomes and ensure it remains visible when filtered
    widget._intended_only.setChecked(False)
    widget._tree.setCurrentItem(widget._tree.topLevelItem(1))
    widget._toggle_pin()
    widget._intended_only.setChecked(True)
    labels = [widget._tree.topLevelItem(i).text(0) for i in range(widget._tree.topLevelItemCount())]
    assert "2" in labels, "Pinned run should remain visible even if it does not match filters"

    # Focus search via shortcut helper and filter genome
    widget._clear_filters()
    widget._search_box.setText("demo")
    widget._refilter()
    assert widget._tree.topLevelItemCount() >= 1

    # Compare two runs
    widget._clear_filters()
    widget._tree.topLevelItem(0).setSelected(True)
    widget._tree.topLevelItem(1).setSelected(True)
    widget._compare_selected()
    assert compared, "Compare callback should receive selected snapshots"
