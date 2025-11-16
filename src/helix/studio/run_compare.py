from __future__ import annotations

from typing import Mapping

from PySide6.QtWidgets import QHBoxLayout, QLabel, QTextEdit, QVBoxLayout, QWidget

from .run_metrics import (
    classify_run_delta,
    compute_run_metrics,
    diff_summary_text,
    format_snapshot_text,
)


class RunCompareWidget(QWidget):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(6)
        layout.addWidget(QLabel("Select two runs in history and click Compare to see details.", self))
        self._verdict_label = QLabel("", self)
        self._verdict_label.setWordWrap(True)
        layout.addWidget(self._verdict_label)

        body = QHBoxLayout()
        self._left = QTextEdit(self)
        self._left.setReadOnly(True)
        self._right = QTextEdit(self)
        self._right.setReadOnly(True)
        body.addWidget(self._left)
        body.addWidget(self._right)
        layout.addLayout(body)

        layout.addWidget(QLabel("Diff summary", self))
        self._diff = QTextEdit(self)
        self._diff.setReadOnly(True)
        layout.addWidget(self._diff)

    def compare_snapshots(self, lhs: Mapping[str, object], rhs: Mapping[str, object]) -> None:
        lhs_metrics = compute_run_metrics(lhs)
        rhs_metrics = compute_run_metrics(rhs)
        self._left.setPlainText(format_snapshot_text(lhs_metrics))
        self._right.setPlainText(format_snapshot_text(rhs_metrics))
        verdict = classify_run_delta(lhs_metrics, rhs_metrics)
        self._verdict_label.setText(f"Verdict: {verdict.label} – {verdict.details}")
        self._diff.setPlainText(diff_summary_text(lhs_metrics, rhs_metrics))
