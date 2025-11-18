from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable, Mapping, List

from PySide6.QtCore import Qt
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLineEdit,
    QFileDialog,
    QMessageBox,
    QPushButton,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
    QComboBox,
    QCheckBox,
    QLabel,
    QInputDialog,
)

from .session import RunKind, SessionModel
from .run_metrics import (
    build_report,
    render_markdown_report,
    summarize_snapshot,
    ensure_state_payload,
    compute_run_metrics,
    format_snapshot_text,
)


class RunHistoryWidget(QWidget):
    def __init__(
        self,
        session: SessionModel,
        parent: QWidget | None = None,
        *,
        on_compare: Callable[[Mapping[str, Any], Mapping[str, Any]], None] | None = None,
    ) -> None:
        super().__init__(parent)
        self._session = session
        self._on_compare = on_compare
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)

        filter_row = QHBoxLayout()
        filter_row.addWidget(QLabel("Run kind:", self))
        self._kind_filter = QComboBox(self)
        self._kind_filter.addItem("All", userData="ALL")
        for kind in RunKind:
            if kind == RunKind.NONE:
                continue
            self._kind_filter.addItem(kind.name.title(), userData=kind.name)
        self._kind_filter.currentIndexChanged.connect(self._refilter)
        filter_row.addWidget(self._kind_filter)

        filter_row.addWidget(QLabel("Genome contains:", self))
        self._search_box = QLineEdit(self)
        self._search_box.setPlaceholderText("e.g. chr19, demo...")
        self._search_box.textChanged.connect(self._refilter)
        filter_row.addWidget(self._search_box)

        self._intended_only = QCheckBox("Has intended", self)
        self._intended_only.stateChanged.connect(self._refilter)
        filter_row.addWidget(self._intended_only)

        filter_row.addWidget(QLabel("Profile:", self))
        self._profile_filter = QComboBox(self)
        self._profile_filter.addItem("All runs", userData=None)
        self._profile_filter.currentIndexChanged.connect(self._refilter)
        filter_row.addWidget(self._profile_filter)
        filter_row.addStretch(1)
        layout.addLayout(filter_row)

        self._tree = QTreeWidget(self)
        self._tree.setSelectionMode(QTreeWidget.SelectionMode.ExtendedSelection)
        self._tree.setHeaderLabels(
            ["Run #", "Pinned", "Kind", "Guide", "Genome", "Label", "Timestamp", "Draws", "Outcomes"]
        )
        self._tree.itemDoubleClicked.connect(self._restore_item)
        self._tree.currentItemChanged.connect(lambda *_: self._update_buttons())
        layout.addWidget(self._tree, stretch=1)

        btn_row = QHBoxLayout()
        self._restore_btn = QPushButton("Restore Selected", self)
        self._restore_btn.clicked.connect(self._restore_selected)
        btn_row.addWidget(self._restore_btn)

        self._pin_btn = QPushButton("Pin", self)
        self._pin_btn.clicked.connect(self._toggle_pin)
        btn_row.addWidget(self._pin_btn)

        self._compare_btn = QPushButton("Compare…", self)
        self._compare_btn.clicked.connect(self._compare_selected)
        btn_row.addWidget(self._compare_btn)

        self._report_btn = QPushButton("Export Report…", self)
        self._report_btn.clicked.connect(self._export_report)
        btn_row.addWidget(self._report_btn)

        self._label_btn = QPushButton("Set Label…", self)
        self._label_btn.clicked.connect(self._set_label)
        btn_row.addWidget(self._label_btn)

        self._export_btn = QPushButton("Export…", self)
        self._export_btn.clicked.connect(self._export_selected)
        btn_row.addWidget(self._export_btn)

        btn_row.addStretch(1)
        layout.addLayout(btn_row)

        hint = QLabel("Search Ctrl+F · Intended Ctrl+Shift+I · Cycle Ctrl+Shift+K · Clear Ctrl+Shift+F", self)
        hint.setStyleSheet("color: #7a8099; font-size: 11px;")
        layout.addWidget(hint)

        self._snapshots: list[dict[str, Any]] = []
        self._pinned: set[int] = set()
        self._session.runRecorded.connect(self._add_entry)
        self._populate_tree()
        self._refresh_profiles()
        self._update_buttons()
        self._setup_shortcuts()
        self._announce_filter_state()

    def _add_entry(self, payload: Mapping[str, Any]) -> None:
        snapshot = dict(payload)
        self._snapshots.append(snapshot)
        self._populate_tree(keep_selection=True)
        self._refresh_profiles()

    def _populate_tree(self, keep_selection: bool = False) -> None:
        previous_index = self._current_index() if keep_selection else None
        self._tree.clear()
        for idx, snapshot in enumerate(self._snapshots):
            metrics = compute_run_metrics(snapshot)
            summary = metrics.summary
            if not self._matches_filters(summary, idx):
                continue
            item = QTreeWidgetItem(
                [
                    summary["run_id"],
                    "📌" if idx in self._pinned else "",
                    summary["run_kind"],
                    summary["guide"],
                    summary["genome"],
                    summary["label"],
                    summary["timestamp"],
                    summary["draws"],
                    summary["outcomes"],
                ]
            )
            item.setData(0, Qt.UserRole, idx)
            tooltip = format_snapshot_text(metrics)
            if metrics.is_pcr:
                pcr_lines = [
                    "PCR run",
                    f"Amplicon: {metrics.pcr_amplicon_length} bp",
                ]
                if (
                    metrics.pcr_amplicon_start is not None
                    and metrics.pcr_amplicon_end is not None
                    and metrics.pcr_amplicon_end > metrics.pcr_amplicon_start
                ):
                    pcr_lines.append(f"Region: {metrics.pcr_amplicon_start}–{metrics.pcr_amplicon_end}")
                if metrics.pcr_cycles:
                    pcr_lines.append(f"Cycles: {metrics.pcr_cycles}")
                if metrics.pcr_forward_primer:
                    pcr_lines.append(f"Forward primer: {metrics.pcr_forward_primer}")
                if metrics.pcr_reverse_primer:
                    pcr_lines.append(f"Reverse primer: {metrics.pcr_reverse_primer}")
                if metrics.pcr_final_mass_ng is not None:
                    pcr_lines.append(f"Final mass: {metrics.pcr_final_mass_ng:.2f} ng")
                if metrics.pcr_final_mutation_rate is not None:
                    pcr_lines.append(f"Mutation rate: {metrics.pcr_final_mutation_rate:.4f}")
                tooltip = tooltip + "\n\n" + "\n".join(pcr_lines)
            item.setToolTip(0, tooltip)
            self._tree.addTopLevelItem(item)
            if previous_index is not None and previous_index == idx:
                self._tree.setCurrentItem(item)
        self._update_buttons()

    def _matches_filters(self, summary: Mapping[str, Any], idx: int) -> bool:
        if idx in self._pinned:
            return True
        kind_filter = self._kind_filter.currentData()
        if kind_filter and kind_filter != "ALL" and summary["run_kind_key"] != kind_filter:
            return False
        if self._intended_only.isChecked() and not summary.get("has_intended"):
            return False
        profile_key = self._profile_filter.currentData()
        if profile_key and summary.get("profile_key") != profile_key:
            return False
        query = self._search_box.text().strip().lower()
        if query:
            haystack = f"{summary['genome']} {summary['guide']} {summary['label']}".lower()
            if query not in haystack:
                return False
        return True

    def _refilter(self) -> None:
        self._populate_tree(keep_selection=True)
        self._update_buttons()
        self._announce_filter_state()

    def _refresh_profiles(self) -> None:
        current = self._profile_filter.currentData()
        profiles: dict[tuple[str, str, str], dict[str, Any]] = {}
        for snapshot in self._snapshots:
            metrics = compute_run_metrics(snapshot)
            summary = metrics.summary
            key = summary.get("profile_key")
            label = summary.get("profile_label")
            if not key or not label:
                continue
            entry = profiles.setdefault(key, {"label": label, "count": 0, "pinned": 0})
            entry["count"] += 1
        for idx in self._pinned:
            if 0 <= idx < len(self._snapshots):
                metrics = compute_run_metrics(self._snapshots[idx])
                key = metrics.summary.get("profile_key")
                if key in profiles:
                    profiles[key]["pinned"] += 1

        self._profile_filter.blockSignals(True)
        self._profile_filter.clear()
        self._profile_filter.addItem("All runs", userData=None)
        self._profile_filter.setItemData(0, "Show every recorded run", Qt.ItemDataRole.ToolTipRole)
        for key, data in sorted(profiles.items(), key=lambda item: item[1]["label"]):
            label = data["label"]
            count = data["count"]
            pinned = data["pinned"]
            run_word = "run" if count == 1 else "runs"
            display = f"{label} — {count} {run_word}"
            if pinned:
                display += f" ({pinned} pinned)"
            index = self._profile_filter.count()
            self._profile_filter.addItem(display, userData=key)
            tooltip_parts = [f"{count} {run_word}"]
            if pinned:
                tooltip_parts.append(f"{pinned} pinned")
            self._profile_filter.setItemData(
                index,
                " • ".join(tooltip_parts),
                Qt.ItemDataRole.ToolTipRole,
            )
        if current:
            idx = self._profile_filter.findData(current)
            if idx >= 0:
                self._profile_filter.setCurrentIndex(idx)
        self._profile_filter.blockSignals(False)

    def _restore_selected(self) -> None:
        item = self._tree.currentItem()
        if not item:
            QMessageBox.information(self, "Run History", "Select a run to restore.")
            return
        self._restore_item(item)

    def _restore_item(self, item: QTreeWidgetItem) -> None:
        idx = self._current_index()
        if idx is None:
            return
        try:
            snapshot = self._snapshots[idx]
        except IndexError:
            return
        state_payload = ensure_state_payload(snapshot)
        if not state_payload:
            QMessageBox.warning(self, "Run History", "Stored run payload is invalid.")
            return
        restored = dict(state_payload)
        restored["viz_dirty"] = True
        self._session.load_payload(restored)
        self._session.set_status(f"Restored run #{restored.get('run_id', '?')}")

    def _toggle_pin(self) -> None:
        item = self._tree.currentItem()
        if not item:
            return
        idx = self._current_index()
        if idx is None:
            return
        if idx in self._pinned:
            self._pinned.remove(idx)
            item.setText(1, "")
            self._pin_btn.setText("Pin")
        else:
            self._pinned.add(idx)
            item.setText(1, "📌")
            self._pin_btn.setText("Unpin")

    def _export_selected(self) -> None:
        item = self._tree.currentItem()
        if not item:
            QMessageBox.information(self, "Run History", "Select a run to export.")
            return
        idx = self._current_index()
        if idx is None:
            return
        try:
            snapshot = self._snapshots[idx]
        except IndexError:
            return
        path_str, _ = QFileDialog.getSaveFileName(
            self,
            "Export Run Snapshot",
            filter="Helix Session (*.helix.json);;JSON (*.json)",
        )
        if not path_str:
            return
        Path(path_str).write_text(json.dumps(snapshot, indent=2), encoding="utf-8")
        self._session.set_status(f"Exported run snapshot to {path_str}")

    def _export_report(self) -> None:
        idx = self._current_index()
        if idx is None:
            QMessageBox.information(self, "Run History", "Select a run before exporting a report.")
            return
        snapshot = self._snapshots[idx]
        report = build_report(snapshot)
        path_str, selected_filter = QFileDialog.getSaveFileName(
            self,
            "Export Run Report",
            filter="Markdown (*.md);;JSON (*.json)",
        )
        if not path_str:
            return
        path = Path(path_str)
        use_json = selected_filter.startswith("JSON") or path.suffix.lower() == ".json"
        payload = json.dumps(report, indent=2) if use_json else render_markdown_report(report)
        path.write_text(payload, encoding="utf-8")
        self._session.set_status(f"Exported run report to {path}")

    def _set_label(self) -> None:
        idx = self._current_index()
        if idx is None:
            QMessageBox.information(self, "Run History", "Select a run to label.")
            return
        snapshot = self._snapshots[idx]
        state_payload = ensure_state_payload(snapshot)
        if state_payload is None:
            QMessageBox.warning(self, "Run History", "Cannot label this snapshot (invalid state).")
            return
        current = state_payload.get("label", "")
        text, ok = QInputDialog.getText(self, "Set Run Label", "Label:", text=str(current))
        if not ok:
            return
        state_payload["label"] = text.strip()
        self._populate_tree(keep_selection=True)
        self._session.set_status(f"Updated label for run #{state_payload.get('run_id', '?')}")
        self._refresh_profiles()

    def _update_buttons(self) -> None:
        has_selection = self._tree.currentItem() is not None
        selected_indices = self._selected_indices()
        self._restore_btn.setEnabled(has_selection)
        self._pin_btn.setEnabled(has_selection)
        self._export_btn.setEnabled(has_selection)
        self._compare_btn.setEnabled(len(selected_indices) == 2)
        self._label_btn.setEnabled(has_selection)
        self._report_btn.setEnabled(has_selection)
        if has_selection:
            idx = selected_indices[0] if selected_indices else None
            if idx is not None and idx in self._pinned:
                self._pin_btn.setText("Unpin")
            else:
                self._pin_btn.setText("Pin")
        else:
            self._pin_btn.setText("Pin")

    def _setup_shortcuts(self) -> None:
        self._shortcuts: list[QShortcut] = []
        search_sc = QShortcut(QKeySequence("Ctrl+F"), self)
        search_sc.activated.connect(self._focus_search)
        self._shortcuts.append(search_sc)

        intended_sc = QShortcut(QKeySequence("Ctrl+Shift+I"), self)
        intended_sc.activated.connect(lambda: self._intended_only.setChecked(not self._intended_only.isChecked()))
        self._shortcuts.append(intended_sc)

        kind_sc = QShortcut(QKeySequence("Ctrl+Shift+K"), self)
        kind_sc.activated.connect(self._cycle_kind_filter)
        self._shortcuts.append(kind_sc)

        clear_sc = QShortcut(QKeySequence("Ctrl+Shift+F"), self)
        clear_sc.activated.connect(self._clear_filters)
        self._shortcuts.append(clear_sc)

    def _current_index(self) -> int | None:
        item = self._tree.currentItem()
        if not item:
            return None
        index = item.data(0, Qt.UserRole)
        if index is None:
            return None
        try:
            return int(index)
        except (TypeError, ValueError):
            return None

    def _selected_indices(self) -> list[int]:
        indices: list[int] = []
        for item in self._tree.selectedItems():
            data = item.data(0, Qt.UserRole)
            try:
                indices.append(int(data))
            except (TypeError, ValueError):
                continue
        if not indices and self._tree.currentItem() is not None:
            idx = self._current_index()
            if idx is not None:
                indices.append(idx)
        return sorted(set(indices))

    def _focus_search(self) -> None:
        self._search_box.setFocus()
        self._search_box.selectAll()

    def _cycle_kind_filter(self) -> None:
        idx = self._kind_filter.currentIndex()
        idx = (idx + 1) % self._kind_filter.count()
        self._kind_filter.setCurrentIndex(idx)

    def _clear_filters(self) -> None:
        self._search_box.clear()
        self._intended_only.setChecked(False)
        self._kind_filter.setCurrentIndex(0)
        self._profile_filter.setCurrentIndex(0)
        self._announce_filter_state()

    def _compare_selected(self) -> None:
        if self._on_compare is None:
            QMessageBox.information(self, "Run History", "Comparison view not available in this build.")
            return
        indices = self._selected_indices()
        if len(indices) != 2:
            QMessageBox.information(self, "Run History", "Select exactly two runs to compare.")
            return
        try:
            lhs = self._snapshots[indices[0]]
            rhs = self._snapshots[indices[1]]
        except IndexError:
            return
        self._on_compare(lhs, rhs)

    def select_run_by_position(self, position: int) -> None:
        if position < 0 or position >= self._tree.topLevelItemCount():
            return
        item = self._tree.topLevelItem(position)
        if item:
            self._tree.setCurrentItem(item)
            self._tree.scrollToItem(item)

    def focusInEvent(self, event) -> None:  # type: ignore[override]
        super().focusInEvent(event)
        self._announce_filter_state()

    def _announce_filter_state(self) -> None:
        visible = self._tree.topLevelItemCount()
        total = len(self._snapshots)
        profile = self._profile_filter.currentText() or "All runs"
        if not self._profile_filter.currentData():
            profile = "All runs"
        kind = self._kind_filter.currentText() if self._kind_filter.currentIndex() > 0 else "All kinds"
        intended = "on" if self._intended_only.isChecked() else "any"
        query = self._search_box.text().strip()
        search_fragment = f'"{query}"' if query else "—"
        pinned_in_view = sum(1 for idx in self._pinned if idx < len(self._snapshots))
        pinned_fragment = f" · ★{pinned_in_view}" if pinned_in_view else ""
        message = (
            f"Runs {visible}/{total} · Profile: {profile} · Kind: {kind} · "
            f"Intended: {intended} · Search: {search_fragment}{pinned_fragment} · "
            "Shortcuts: Ctrl+F search · Ctrl+Shift+I intended · Ctrl+Shift+K kind · "
            "Ctrl+Shift+F clear"
        )
        self._session.set_status(message)
