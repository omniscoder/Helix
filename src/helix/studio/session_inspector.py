from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Optional

from PySide6.QtWidgets import QLabel, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget

from .session import ExperimentState, SessionModel


class SessionInspector(QWidget):
    def __init__(self, session: SessionModel, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._session = session

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(4)

        self._status = QLabel("Awaiting simulation run.", self)
        layout.addWidget(self._status)

        self._tree = QTreeWidget(self)
        self._tree.setHeaderLabels(["Field", "Value"])
        self._tree.setRootIsDecorated(True)
        layout.addWidget(self._tree, stretch=1)

        self._session.stateChanged.connect(self._on_state_changed)
        self._on_state_changed(self._session.state)

    def _on_state_changed(self, state: ExperimentState) -> None:
        self._tree.clear()
        self._status.setText("Session updated.")

        genome_item = QTreeWidgetItem(["Genome", self._describe_genome(state.genome)])
        self._tree.addTopLevelItem(genome_item)
        QTreeWidgetItem(self._tree, ["Genome source", self._describe_genome_source(state)])
        QTreeWidgetItem(self._tree, ["Genome hash", state.genome_hash or "—"])

        guide_info = self._guide_summary(state.config)
        QTreeWidgetItem(self._tree, ["Guide", guide_info])

        peg_info = self._peg_summary(state.peg)
        QTreeWidgetItem(self._tree, ["pegRNA", peg_info])

        editor_info = self._editor_summary(state.editor)
        QTreeWidgetItem(self._tree, ["Editor", editor_info])

        run_config = state.config.get("run_config") if isinstance(state.config, Mapping) else {}
        QTreeWidgetItem(self._tree, ["Last run", self._format_run_config(run_config)])
        QTreeWidgetItem(self._tree, ["Run kind", state.run_kind.name])
        QTreeWidgetItem(self._tree, ["Viz dirty", "yes" if state.viz_dirty else "no"])

        outcomes_count = len(state.outcomes)
        QTreeWidgetItem(self._tree, ["Outcomes", f"{outcomes_count} recorded"])

        self._tree.expandAll()

    def _describe_genome(self, genome: Optional[Any]) -> str:
        if genome is None or not getattr(genome, "sequences", None):
            return "—"
        seqs = getattr(genome, "sequences", {})
        lengths = [len(seq) for seq in seqs.values()]
        total = sum(lengths)
        chroms = len(lengths)
        return f"{total} bp across {chroms} contig(s)"

    def _guide_summary(self, config: Mapping[str, Any] | None) -> str:
        if not isinstance(config, Mapping):
            return "—"
        guide = config.get("guide")
        if isinstance(guide, Mapping):
            label = guide.get("id") or guide.get("name") or "guide"
            start = guide.get("start")
            end = guide.get("end")
            strand = guide.get("strand", "+")
            if start is not None and end is not None:
                return f"{label}: {start}-{end} ({strand})"
            return str(label)
        return "—"

    def _peg_summary(self, peg: Optional[Any]) -> str:
        if peg is None:
            return "—"
        name = getattr(peg, "name", None) or "pegRNA"
        spacer = getattr(peg, "spacer", "")
        return f"{name} (spacer {len(spacer)} bp)"

    def _editor_summary(self, editor: Optional[Any]) -> str:
        if editor is None:
            return "—"
        name = getattr(editor, "name", "Editor")
        offset = getattr(editor, "nick_to_edit_offset", 0)
        return f"{name} (nick offset {offset})"

    def _format_run_config(self, config: Mapping[str, Any] | None) -> str:
        if not isinstance(config, Mapping) or not config:
            return "—"
        parts: list[str] = []
        for key in ("mode", "draws", "seed", "pam_profile", "priors_profile"):
            if key in config and config[key] is not None:
                parts.append(f"{key}={config[key]}")
        return ", ".join(parts) if parts else "—"

    def _describe_genome_source(self, state: ExperimentState) -> str:
        source = state.genome_source or "inline"
        uri = state.genome_uri
        if source == "file" and uri:
            return f"file: {uri}"
        if uri:
            return f"{source}: {uri}"
        return source
