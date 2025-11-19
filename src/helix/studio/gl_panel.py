from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Optional

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QStackedLayout,
    QVBoxLayout,
    QWidget,
)

from helix.gui.modern.builders import (
    build_crispr_orbitals_3d_spec,
    build_crispr_rail_3d_spec,
    build_crispr_viz_spec,
    build_prime_scaffold_3d_spec,
    build_prime_viz_spec,
)
from helix.gui.modern.dag_adapters import crispr_dag_to_viz_spec, crispr_dag_to_workflow
from helix.gui.modern.qt import HelixModernWidget
from helix.gui.railbands.widget import RailBandWidget
from helix.gui.modern.spec import EditVisualizationSpec, load_viz_spec
from helix.studio.session import ExperimentState, RunKind, SessionModel
from helix.studio.run_metrics import RunMetrics, classify_run_quality, compute_run_metrics


def _workflow_from_serializable(data: Mapping[str, Any]) -> dict[str, Any]:
    buffers: dict[str, Any] = {}
    scalar_keys = {
        "x_extent",
        "y_extent",
        "captured_mass",
        "total_mass",
        "cut_index",
        "heat_y",
        "heat_height",
        "slice_width",
    }
    for key, value in data.items():
        if key in scalar_keys:
            try:
                buffers[key] = float(value)
            except Exception:
                buffers[key] = 0.0
        else:
            buffers[key] = np.array(value, dtype="f4")
    return buffers

class GLViewport(QWidget):
    def __init__(self, session: SessionModel, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._session = session
        self._last_error_run_id: Optional[int] = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        # Local viewport controls (rail view + animation).
        controls = QWidget(self)
        controls_layout = QHBoxLayout(controls)
        controls_layout.setContentsMargins(8, 4, 8, 4)
        controls_layout.setSpacing(12)

        self._animate_toggle = QCheckBox("Animate", controls)
        self._animate_toggle.setChecked(True)
        controls_layout.addWidget(self._animate_toggle)

        controls_layout.addWidget(QLabel("Rail view", controls))
        self._rail_mode_combo = QComboBox(controls)
        self._rail_mode_combo.addItem("Flat rails", "flat")
        self._rail_mode_combo.addItem("3D helix", "helix")
        self._rail_mode_combo.setCurrentIndex(0)
        controls_layout.addWidget(self._rail_mode_combo)

        controls_layout.addStretch(1)

        self._stack_host = QWidget(self)
        self._stack = QStackedLayout(self._stack_host)

        self._placeholder = QLabel(
            "Viewport\n\nNo simulation yet.\nRun a CRISPR or Prime simulation to populate DAG, rails, and workflow views.",
            self,
        )
        self._placeholder.setAlignment(Qt.AlignCenter)
        self._placeholder.setStyleSheet("color: #9aa4c2; font-size: 15px;")
        self._stack.addWidget(self._placeholder)

        self._gl = HelixModernWidget()
        self._gl_container = QWidget.createWindowContainer(self._gl, self)
        self._stack.addWidget(self._gl_container)
        self._rail_widget = RailBandWidget(self)
        self._stack.addWidget(self._rail_widget)
        self._stack.setCurrentWidget(self._placeholder)

        self._overlay = QLabel(self._stack_host)
        self._overlay.setObjectName("ViewportOverlay")
        self._overlay.setStyleSheet(
            "#ViewportOverlay { background-color: rgba(8, 12, 24, 190); color: #e0e4ff; padding: 6px 10px; border-radius: 4px; font-size: 12px; }"
        )
        self._overlay.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        self._overlay.setAttribute(Qt.WA_TransparentForMouseEvents)
        self._overlay.hide()

        layout.addWidget(controls)
        layout.addWidget(self._stack_host)

        # Wire controls once the GL widget exists.
        self._animate_toggle.toggled.connect(self._on_animate_toggled)
        self._rail_mode_combo.currentIndexChanged.connect(self._on_rail_mode_changed)

        self._session.stateChanged.connect(self._on_state_changed)
        self._on_state_changed(self._session.state)

    def _on_state_changed(self, state: ExperimentState) -> None:
        if not state.viz_dirty:
            return
        spec = self._extract_spec(state)
        if spec is None and state.outcomes:
            spec = self._build_spec_from_state(state)
        config = state.config if isinstance(state.config, Mapping) else {}
        if spec is not None:
            settings = config.get("sim_settings") or {}
            viz_mode = str((settings.get("viz_mode") or "rail_3d")).lower()
            metadata = spec.metadata or {}
            workflow_meta = config.get("workflow_view") or metadata.get("workflow_view")

            # If we have a CRISPR edit DAG but no workflow metadata yet,
            # derive a 2.5D workflow view directly from the DAG.
            if workflow_meta is None and isinstance(state.dag_from_runtime, Mapping):
                try:
                    wf = crispr_dag_to_workflow(state.dag_from_runtime)
                except Exception:
                    wf = None
                if wf is not None:
                    workflow_meta = wf
                    config["workflow_view"] = wf

            if viz_mode == "rail_2d":
                rail_meta = metadata.get("rail_bands")
                if isinstance(rail_meta, Mapping) and rail_meta.get("bands"):
                    self._rail_widget.set_spec(spec)
                    self._gl.set_workflow_view(None)
                    self._stack.setCurrentWidget(self._rail_widget)
                    self._overlay.hide()
                    self._last_error_run_id = None
                    return

            self._stack.setCurrentWidget(self._gl_container)
            self._gl.set_spec(spec)
            if viz_mode == "workflow_2p5d" and isinstance(workflow_meta, Mapping):
                self._gl.set_workflow_view(_workflow_from_serializable(workflow_meta))
            else:
                self._gl.set_workflow_view(None)
            self._gl.update()
            metrics = None
            try:
                metrics = compute_run_metrics(self._session.to_payload())
            except Exception:
                metrics = None
            self._refresh_overlay(metrics)
            self._last_error_run_id = None
        else:
            self._gl.set_workflow_view(None)
            if state.dag_from_runtime:
                self._show_placeholder("Viewport\n\nDAG available – awaiting visualization.")
            else:
                self._show_placeholder(
                    "Viewport\n\nNo simulation yet.\nRun a CRISPR or Prime simulation to populate DAG, rails, and workflow views."
                )
            self._overlay.hide()
            if (state.viz_spec_payload or state.dag_from_runtime or state.outcomes) and (
                state.run_id != self._last_error_run_id
            ):
                self._session.error("Failed to build visualization spec from session; see logs.")
                self._last_error_run_id = state.run_id
        self._session.mark_viz_clean()

    def set_dag(self, dag) -> None:
        # self._gl.set_dag(dag)
        pass

    def _extract_spec(self, state: ExperimentState) -> Optional[EditVisualizationSpec]:
        config = state.config if isinstance(state.config, Mapping) else {}
        candidate = state.viz_spec_payload or config.get("viz_spec_payload") or config.get("viz_spec")
        if candidate is not None:
            try:
                return load_viz_spec(candidate)
            except Exception:
                self._session.log("Failed to load viz spec from session config.")
        dag_payload = state.dag_from_runtime
        if isinstance(dag_payload, Mapping):
            try:
                return crispr_dag_to_viz_spec(dag_payload)
            except Exception:
                self._session.log("Failed to adapt DAG payload into visualization spec.")
        return None

    def _build_spec_from_state(self, state: ExperimentState) -> Optional[EditVisualizationSpec]:
        config = state.config if isinstance(state.config, Mapping) else {}
        payload = config.get("sim_payload")
        if not isinstance(payload, Mapping):
            return None
        viz_mode = str(((config.get("sim_settings") or {}).get("viz_mode") or "rail_3d")).lower()
        sim_type = state.run_kind.name if state.run_kind != RunKind.NONE else str(config.get("sim_type") or "").upper()
        if sim_type == "CRISPR":
            sequence = self._sequence_from_state(state) or self._sequence_from_payload(payload)
            guide = config.get("guide") or payload.get("guide")
            if not sequence or not isinstance(guide, Mapping):
                return None
            if viz_mode == "orbitals_3d":
                spec = build_crispr_orbitals_3d_spec(sequence, guide, payload)
            else:
                spec = build_crispr_rail_3d_spec(sequence, guide, payload)
            if spec is None:
                spec = build_crispr_viz_spec(sequence, guide, payload)
            return spec
        if sim_type == "PRIME":
            genome = state.genome
            peg = state.peg
            editor = state.editor
            if genome is None or peg is None or editor is None:
                return None
            if viz_mode == "prime_scaffold_3d":
                return build_prime_scaffold_3d_spec(genome, peg, editor, payload)
            return build_prime_viz_spec(genome, peg, editor, payload)
        return None

    def _sequence_from_state(self, state: ExperimentState) -> Optional[str]:
        genome = state.genome
        sequences = getattr(genome, "sequences", None) if genome is not None else None
        if isinstance(sequences, Mapping):
            for seq in sequences.values():
                if isinstance(seq, str) and seq:
                    return seq
        return None

    def _sequence_from_payload(self, payload: Mapping[str, Any]) -> Optional[str]:
        site = payload.get("site")
        if isinstance(site, Mapping):
            seq = site.get("sequence")
            if isinstance(seq, str):
                return seq
        seq = payload.get("sequence")
        if isinstance(seq, str):
            return seq
        return None

    def _show_placeholder(self, message: str) -> None:
        self._placeholder.setText(message)
        self._stack.setCurrentWidget(self._placeholder)
        self._overlay.hide()

    def _refresh_overlay(self, snapshot_metrics: RunMetrics | None) -> None:
        lines = build_overlay_lines(snapshot_metrics)
        if not lines:
            self._overlay.hide()
            return
        self._overlay.setText("\n".join(lines))
        self._overlay.adjustSize()
        self._overlay.move(12, 12)
        self._overlay.show()

    def _on_animate_toggled(self, enabled: bool) -> None:
        if self._gl is not None:
            self._gl.set_animate(enabled)

    def _on_rail_mode_changed(self, index: int) -> None:
        if self._gl is None:
            return
        mode = self._rail_mode_combo.itemData(index) or "flat"
        self._gl.set_rail_mode(str(mode))

def build_overlay_lines(metrics: RunMetrics | None) -> list[str]:
    if metrics is None:
        return []
    summary = metrics.summary
    draw_label = summary.get("draws") or "—"
    run_label = f"Run #{summary.get('run_id', '—')} [{summary.get('run_kind', '—')}]"
    genome_label = summary.get("genome", "—")
    lines: list[str] = [f"{run_label} | Draws: {draw_label} | Genome: {genome_label}"]
    dag = metrics.dag_stats
    if any(dag.values()):
        lines.append(
            "DAG nodes: {nodes} · branches: {branches} · leaves: {leaves} · depth: {depth}".format(
                nodes=dag.get("nodes", 0),
                branches=dag.get("branch_nodes", 0),
                leaves=dag.get("leaves", 0),
                depth=dag.get("depth", 0),
            )
        )
    total_outcomes = metrics.intended_count + metrics.off_target_count + metrics.neutral_count
    lines.append(
        f"Outcomes: {total_outcomes} (Intended: {metrics.intended_count}, Off-target: {metrics.off_target_count}, Neutral: {metrics.neutral_count})"
    )
    lines.append(
        f"Prob mass → Intended: {metrics.intended_mass:.2f}  Off-target: {metrics.off_target_mass:.2f}  Neutral: {metrics.neutral_mass:.2f}"
    )
    perf_line = _build_perf_line(metrics.perf)
    if perf_line:
        lines.append(perf_line)
    if metrics.is_pcr:
        region = ""
        if (
            metrics.pcr_amplicon_start is not None
            and metrics.pcr_amplicon_end is not None
            and metrics.pcr_amplicon_end > metrics.pcr_amplicon_start
        ):
            region = f" ({metrics.pcr_amplicon_start}–{metrics.pcr_amplicon_end})"
        lines.append(
            "PCR: {length} bp{region} · {cycles} cycles · {mass:.2f} ng · mut {mut:.4f}".format(
                length=metrics.pcr_amplicon_length or 0,
                region=region,
                cycles=metrics.pcr_cycles or "?",
                mass=metrics.pcr_final_mass_ng or 0.0,
                mut=metrics.pcr_final_mutation_rate or 0.0,
            )
        )
    quality_line = _build_quality_line(metrics)
    if quality_line:
        lines.append(quality_line)
    return lines


def _build_perf_line(perf: Mapping[str, float]) -> str:
    if not perf:
        return ""
    parts: list[str] = []
    sim_ms = perf.get("sim_ms")
    if sim_ms is not None:
        parts.append(f"Sim {sim_ms:.0f} ms")
    viz_ms = perf.get("viz_ms")
    if viz_ms is not None:
        parts.append(f"Viz {viz_ms:.0f} ms")
    if not parts:
        return ""
    return "Perf: " + " | ".join(parts)


def _build_quality_line(metrics: RunMetrics) -> str:
    verdict = classify_run_quality(metrics)
    if not verdict.label:
        return ""
    return f"Quality: {verdict.label} – {verdict.details}"
