"""Interactive PyQt Sim Builder for CRISPR/prime editing workflows."""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass, fields
import os
from functools import partial
from pathlib import Path
from typing import Any, Callable, Optional

import numpy as np
from PySide6.QtCore import QByteArray, QObject, QSettings, Qt, Signal, QThread, QTimer
from PySide6.QtGui import QCloseEvent, QColor, QPalette, QSurfaceFormat
from PySide6.QtWidgets import (
    QApplication,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QPlainTextEdit,
    QSizePolicy,
    QSpinBox,
    QSplitter,
    QStackedLayout,
    QTabWidget,
    QVBoxLayout,
    QWidget,
    QComboBox,
    QDoubleSpinBox,
    QCompleter,
)

from helix import bioinformatics
from helix.crispr import guide as guide_module
from helix.crispr.model import CasSystem, CasSystemType, DigitalGenome, PAMRule, GuideRNA
from helix.crispr.pam import build_crispr_pam_mask, build_prime_pam_mask
from helix.crispr.simulate import resolve_crispr_priors, simulate_cut_repair
from helix.crispr.dag_api import build_crispr_edit_dag
from helix.gui.modern.qt import HelixModernWidget
from helix.gui.modern.builders import (
    build_crispr_viz_spec,
    build_prime_viz_spec,
    build_crispr_rail_3d_spec,
    build_crispr_orbitals_3d_spec,
    build_prime_scaffold_3d_spec,
)
from helix.gui.modern.dag_adapters import crispr_dag_to_viz_spec
from helix.gui.modern.builders_2p5d import build_workflow2p5d_crispr, build_workflow2p5d_prime
from helix.gui.railbands.widget import RailBandWidget
from helix.gui.modern.spec import EditVisualizationSpec
from helix.gui.theme import apply_helix_theme

LOGGER = logging.getLogger(__name__)
from helix.prime.model import PegRNA, PrimeEditor
from helix.prime.priors import resolve_prime_priors
from helix.prime.simulator import simulate_prime_edit


def _read_fasta(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    sequence_lines = [line for line in lines if not line.startswith(">")]
    return "".join(sequence_lines)


def _workflow_to_serializable(data: dict[str, Any]) -> dict[str, Any]:
    serializable: dict[str, Any] = {}
    for key, value in data.items():
        if isinstance(value, np.ndarray):
            serializable[key] = value.tolist()
        else:
            serializable[key] = value
    return serializable


def _workflow_from_serializable(data: dict[str, Any]) -> dict[str, np.ndarray | float]:
    buffers: dict[str, np.ndarray | float] = {}
    for key, value in data.items():
        if key in ("x_extent", "y_extent", "captured_mass", "total_mass", "cut_index", "heat_y", "heat_height", "slice_width"):
            buffers[key] = float(value)
        else:
            buffers[key] = np.array(value, dtype="f4")
    return buffers


class SequenceInput(QWidget):
    """Shared sequence editor with load-from-file support."""

    sequenceChanged = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        controls = QHBoxLayout()
        self._load_btn = QPushButton("Load FASTA…", self)
        self._clear_btn = QPushButton("Clear", self)
        controls.addWidget(self._load_btn)
        controls.addWidget(self._clear_btn)
        layout.addLayout(controls)
        self._text = QPlainTextEdit(self)
        self._text.setPlaceholderText("Paste sequence or FASTA content here…")
        layout.addWidget(self._text)
        self._length_label = QLabel("Length: 0 bp", self)
        layout.addWidget(self._length_label)

        self._load_btn.clicked.connect(self._load_from_file)
        self._clear_btn.clicked.connect(self.clear)
        self._text.textChanged.connect(self._on_text_changed)

    def sequence(self) -> str:
        raw = self._text.toPlainText().strip()
        if not raw:
            return ""
        lines = [line.strip() for line in raw.splitlines() if line.strip()]
        seq = "".join(line for line in lines if not line.startswith(">"))
        return bioinformatics.normalize_sequence(seq)

    def current_sequence(self) -> str:
        """Expose a naming-convention-friendly accessor for consumers."""
        return self.sequence()

    def set_sequence(self, text: str) -> None:
        self._text.setPlainText(text)

    def clear(self) -> None:
        self._text.clear()

    def _on_text_changed(self) -> None:
        seq = self.sequence()
        self._length_label.setText(f"Length: {len(seq)} bp")
        self.sequenceChanged.emit(seq)

    def _load_from_file(self) -> None:
        path_str, _ = QFileDialog.getOpenFileName(self, "Select FASTA", filter="FASTA (*.fa *.fasta *.fna);;All files (*.*)")
        if not path_str:
            return
        path = Path(path_str)
        try:
            seq = _read_fasta(path)
        except Exception as exc:  # pragma: no cover - QFileDialog path errors are UI-bound
            QMessageBox.critical(self, "Failed to load FASTA", str(exc))
            return
        self._text.setPlainText(seq)


class ViewerPlaceholder(QWidget):
    """Simple fallback canvas shown when OpenGL hasn't rendered yet."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setAutoFillBackground(True)
        palette = QPalette(self.palette())
        palette.setColor(QPalette.Window, QColor("#050910"))
        self.setPalette(palette)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addStretch(1)
        self._label = QLabel("Viewer initializing…", self)
        self._label.setAlignment(Qt.AlignCenter)
        self._label.setStyleSheet("color: #6f7899; font-size: 16px;")
        layout.addWidget(self._label)
        layout.addStretch(1)

    def set_message(self, message: str) -> None:
        self._label.setText(message)


@dataclass
class SimSettings:
    """Shared simulation configuration passed across tabs."""

    pam_profile: str = "SpCas9_NGG"
    pam_softness: float = 1.0
    draws: int = 1000
    seed: Optional[int] = None
    priors_profile: str = "default_indel"
    viz_mode: str = "rail_3d"


class SimSettingsModel(QObject):
    """Wraps SimSettings with change notifications."""

    changed = Signal()

    def __init__(self, initial: SimSettings | None = None, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._settings = initial or SimSettings()

    @property
    def value(self) -> SimSettings:
        return self._settings

    def update(self, **kwargs: Any) -> None:
        updated = False
        field_names = {field.name for field in fields(SimSettings)}
        for key, val in kwargs.items():
            if key in field_names:
                current = getattr(self._settings, key)
                if current != val:
                    setattr(self._settings, key, val)
                    updated = True
        if updated:
            self.changed.emit()

    def apply_dict(self, data: dict[str, Any]) -> None:
        filtered = {k: v for k, v in data.items() if hasattr(self._settings, k)}
        if filtered:
            self.update(**filtered)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self._settings)


PAM_PROFILE_PRESETS: dict[str, dict[str, str]] = {
    "SpCas9_NGG": {"pam": "SpCas9-NGG"},
    "SpG_NGN": {"pattern": "NGN", "orientation": "3prime", "notes": "SpG relaxed SpCas9"},
    "SpRY_NNN": {"pam": "SpRY-NNN"},
    "SaCas9_NNGRRT": {"pam": "SaCas9-NNGRRT"},
    "Cas12a_TTTN": {"pattern": "TTTN", "orientation": "5prime", "notes": "Cas12a/Cpf1"},
}
DEFAULT_PAM_PROFILES: tuple[str, ...] = tuple(PAM_PROFILE_PRESETS)
PRIOR_PROFILE_CHOICES: tuple[str, ...] = ("default_indel", "scarless", "indel_heavy", "hdr_only")
VIZ_MODE_CHOICES: tuple[str, ...] = ("rail_3d", "rail_2d", "orbitals_3d", "prime_scaffold_3d", "workflow_2p5d")
VIZ_MODE_LABELS: dict[str, str] = {
    "rail_3d": "Rail 3D (linear genome)",
    "rail_2d": "Rail 2D (bands)",
    "orbitals_3d": "Orbitals 3D (looped)",
    "prime_scaffold_3d": "Prime scaffold 3D",
    "workflow_2p5d": "Workflow 2.5D (timeline)",
}


class PamConfigPanel(QWidget):
    """Shared configuration form that feeds the SimSettingsModel."""

    settingsChanged = Signal(dict)

    def __init__(self, settings_model: SimSettingsModel, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._settings_model = settings_model
        self._syncing = False

        layout = QFormLayout(self)

        self.pam_profile_combo = HelixComboBox(self, searchable=True)
        self.pam_profile_combo.setEditable(True)
        for key in DEFAULT_PAM_PROFILES:
            preset = PAM_PROFILE_PRESETS.get(key, {})
            self.pam_profile_combo.addItem(key)
            notes = preset.get("notes")
            if notes:
                idx = self.pam_profile_combo.count() - 1
                self.pam_profile_combo.setItemData(idx, notes, Qt.ToolTipRole)
        self._configure_combo(self.pam_profile_combo)
        layout.addRow("PAM profile:", self.pam_profile_combo)

        self.pam_softness_spin = QDoubleSpinBox(self)
        self.pam_softness_spin.setRange(0.0, 1.0)
        self.pam_softness_spin.setSingleStep(0.1)
        layout.addRow("PAM softness:", self.pam_softness_spin)

        self.draws_spin = QSpinBox(self)
        self.draws_spin.setRange(1, 1_000_000)
        layout.addRow("Draws:", self.draws_spin)

        self.seed_spin = QSpinBox(self)
        self.seed_spin.setRange(-1, 2**31 - 1)
        self.seed_spin.setSpecialValueText("Random each run")
        self.seed_spin.setToolTip("Use -1 for random seed per run.")
        layout.addRow("Seed:", self.seed_spin)

        self.priors_combo = HelixComboBox(self, searchable=False)
        self.priors_combo.addItems(PRIOR_PROFILE_CHOICES)
        self._configure_combo(self.priors_combo)
        layout.addRow("Priors profile:", self.priors_combo)

        self.viz_mode_combo = HelixComboBox(self, searchable=False)
        for key in VIZ_MODE_CHOICES:
            label = VIZ_MODE_LABELS.get(key, key)
            self.viz_mode_combo.addItem(label, userData=key)
        self._configure_combo(self.viz_mode_combo)
        layout.addRow("Viz mode:", self.viz_mode_combo)

        self._settings_model.changed.connect(self._load_from_model)
        self.pam_profile_combo.currentTextChanged.connect(self._on_changed)
        self.pam_softness_spin.valueChanged.connect(self._on_changed)
        self.draws_spin.valueChanged.connect(self._on_changed)
        self.seed_spin.valueChanged.connect(self._on_changed)
        self.priors_combo.currentTextChanged.connect(self._on_changed)
        self.viz_mode_combo.currentIndexChanged.connect(self._on_changed)

        self._load_from_model()

    def _load_from_model(self) -> None:
        self._syncing = True
        try:
            settings = self._settings_model.value
            self._ensure_combo_text(self.pam_profile_combo, settings.pam_profile or "SpCas9_NGG")
            self.pam_softness_spin.setValue(settings.pam_softness)
            self.draws_spin.setValue(max(1, settings.draws))
            self.seed_spin.setValue(-1 if settings.seed is None else settings.seed)
            self._ensure_combo_text(self.priors_combo, settings.priors_profile or PRIOR_PROFILE_CHOICES[0])
            self._set_combo_by_data(self.viz_mode_combo, settings.viz_mode or "rail_3d")
        finally:
            self._syncing = False

    def _ensure_combo_text(self, combo: QComboBox, text: str) -> None:
        if combo.findText(text) < 0:
            combo.addItem(text)
        combo.setCurrentText(text)

    def _set_combo_by_data(self, combo: QComboBox, value: str) -> None:
        idx = combo.findData(value)
        if idx >= 0:
            combo.setCurrentIndex(idx)

    def _configure_combo(self, combo: QComboBox) -> None:
        combo.setMinimumWidth(220)
        combo.setMinimumContentsLength(16)
        combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        view = combo.view()
        if view is not None:
            view.setMinimumWidth(260)
        combo.setMaxVisibleItems(16)
        combo.setInsertPolicy(QComboBox.NoInsert)

    def _on_changed(self, *_args: object) -> None:
        if self._syncing:
            return
        payload = {
            "pam_profile": self.pam_profile_combo.currentText().strip(),
            "pam_softness": float(self.pam_softness_spin.value()),
            "draws": int(self.draws_spin.value()),
            "seed": None if self.seed_spin.value() < 0 else int(self.seed_spin.value()),
            "priors_profile": self.priors_combo.currentText().strip(),
            "viz_mode": str(self.viz_mode_combo.currentData() or self.viz_mode_combo.currentText().strip()),
        }
        self._settings_model.update(**payload)
        self.settingsChanged.emit(payload)


class LogPanel(QPlainTextEdit):
    """Simple log viewer."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setReadOnly(True)

    def log(self, message: str) -> None:
        self.appendPlainText(message)


class SimJobWorker(QObject):
    """Background worker that executes simulations off the UI thread."""

    finished = Signal(object, object)  # payload, spec
    failed = Signal(str)

    def __init__(self, task: Callable[[], tuple[dict, Optional[EditVisualizationSpec]]], parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._task = task

    def run(self) -> None:
        try:
            payload, spec = self._task()
            if spec is None:
                raise ValueError("Simulation succeeded but no visualization spec was generated.")
            self.finished.emit(payload, spec)
        except Exception as exc:  # pragma: no cover - errors reported via UI
            self.failed.emit(str(exc))


class CrisprTab(QWidget):
    """CRISPR simulation form."""

    specReady = Signal(object)
    simPayloadReady = Signal(object)
    logMessage = Signal(str)
    runCompleted = Signal(str, object, object)  # kind, payload, spec

    def __init__(self, sequence_input: SequenceInput, sim_settings: SimSettingsModel, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._sequence_input = sequence_input
        self._sim_settings = sim_settings
        self._guides: list[dict] = []
        self._last_payload: Optional[dict] = None
        self._last_spec: Optional[EditVisualizationSpec] = None
        self._active_thread: QThread | None = None
        self._active_worker: SimJobWorker | None = None

        layout = QVBoxLayout(self)
        form = QFormLayout()
        self._guide_len = QSpinBox(self)
        self._guide_len.setRange(10, 40)
        self._guide_len.setValue(20)
        form.addRow("Guide length", self._guide_len)
        layout.addLayout(form)

        self._find_btn = QPushButton("Find Guides", self)
        layout.addWidget(self._find_btn)
        self._guide_list = QListWidget(self)
        layout.addWidget(self._guide_list)

        btn_row = QHBoxLayout()
        self._simulate_btn = QPushButton("Simulate CRISPR", self)
        self._save_sim_btn = QPushButton("Save .sim…", self)
        self._save_spec_btn = QPushButton("Save Viz Spec…", self)
        self._save_sim_btn.setEnabled(False)
        self._save_spec_btn.setEnabled(False)
        btn_row.addWidget(self._simulate_btn)
        btn_row.addWidget(self._save_sim_btn)
        btn_row.addWidget(self._save_spec_btn)
        layout.addLayout(btn_row)
        layout.addStretch(1)

        self._find_btn.clicked.connect(self._find_guides)
        self._simulate_btn.clicked.connect(self._simulate)
        self._save_sim_btn.clicked.connect(partial(self._export_json, "crispr.sim json", lambda: self._last_payload))
        self._save_spec_btn.clicked.connect(partial(self._export_spec, lambda: self._last_spec))

    def _sequence_or_warn(self) -> Optional[str]:
        sequence = self._sequence_input.current_sequence()
        if not sequence:
            self.logMessage.emit("Provide a sequence (paste or load FASTA) before running CRISPR sims.")
            return None
        return sequence

    def _find_guides(self) -> None:
        sequence = self._sequence_or_warn()
        if not sequence:
            return
        pam_profile = self._sim_settings.value.pam_profile or "SpCas9_NGG"
        pam = self._pam_from_profile(pam_profile)
        guides = guide_module.find_guides(sequence, pam, self._guide_len.value(), strand="both")
        self._guides = guides
        self._guide_list.clear()
        for guide in guides:
            label = f"{guide['id']} – {guide['start']}-{guide['end']} ({guide['strand']}) GC={guide['gc_content']}"
            item = QListWidgetItem(label)
            item.setData(Qt.UserRole, guide)
            self._guide_list.addItem(item)
        if not guides:
            self.logMessage.emit("No guides found with the current parameters.")
        else:
            self.logMessage.emit(f"Located {len(guides)} guides for PAM {pam_profile}. Select one and click Simulate.")

    def _selected_guide(self) -> Optional[dict]:
        item = self._guide_list.currentItem()
        if not item:
            return None
        return item.data(Qt.UserRole)

    def _simulate(self) -> None:
        sequence = self._sequence_or_warn()
        if not sequence:
            return

        guide = self._selected_guide()
        if not guide:
            self.logMessage.emit("Select a guide before running the simulation.")
            return

        shared = self._sim_settings.value
        draws = max(1, shared.draws)
        seed = shared.seed
        pam_profile = shared.pam_profile or "SpCas9_NGG"
        pam_softness = shared.pam_softness
        priors_profile = shared.priors_profile or "default_indel"
        seed_label = "auto" if seed is None else str(seed)
        raw_dag_flag = (os.environ.get("HELIX_DAG_MODE") or "").strip().lower()
        use_dag_mode = raw_dag_flag in {"1", "true", "yes", "on"}
        mode_label = "DAG" if use_dag_mode else "legacy"
        self.logMessage.emit(
            f"Running CRISPR simulation [{mode_label} mode] "
            f"(draws={draws}, seed={seed_label}, pam={pam_profile}, softness={pam_softness}, priors={priors_profile})…"
        )

        # Capture only primitives and immutable data into the closure.
        def task():
            if use_dag_mode:
                # DAG-first path: build a CRISPR edit DAG and adapt it to a viz spec.
                norm_seq = bioinformatics.normalize_sequence(sequence)
                genome = DigitalGenome(sequences={"chrGui": norm_seq})
                pam_cfg = guide_module.get_pam(pam_profile) if pam_profile else {"pattern": "NGG"}
                pattern = pam_cfg.get("pattern") or "NGG"
                cas = CasSystem(
                    name=f"GUI-{pam_profile}",
                    system_type=CasSystemType.CAS9,
                    pam_rules=[PAMRule(pattern=pattern)],
                    cut_offset=3,
                    max_mismatches=3,
                    weight_mismatch_penalty=1.0,
                    weight_pam_penalty=2.0,
                )
                guide_seq = guide.get("sequence") or norm_seq
                guide_obj = GuideRNA(
                    sequence=bioinformatics.normalize_sequence(str(guide_seq)),
                    pam=pattern,
                    name=str(guide.get("id") or guide.get("name") or "guide"),
                    metadata={"strand": str(guide.get("strand", "+"))},
                )
                dag = build_crispr_edit_dag(
                    genome,
                    cas,
                    guide_obj,
                    rng_seed=seed or 0,
                    max_depth=2,
                    min_prob=1e-4,
                    max_sites=50,
                    use_gpu=False,
                    frame_consumer=None,
                )
                # Reuse CLI payload shape for adapters and downstream tooling.
                from helix.cli import _edit_dag_to_payload  # local import to avoid cycles

                dag_payload = _edit_dag_to_payload(
                    dag,
                    artifact="helix.crispr.edit_dag.v1.1",
                    metadata={
                        "source": "gui.crispr_dag_mode",
                        "guide": guide,
                        "site_sequence": norm_seq,
                    },
                )
                spec = crispr_dag_to_viz_spec(dag_payload)
                if spec is None:
                    raise ValueError("Failed to build CRISPR visualization spec from DAG.")
                return dag_payload, spec

            # Legacy sim-payload path (default).
            priors = resolve_crispr_priors(priors_profile)
            pam_mask = build_crispr_pam_mask(sequence, guide, pam_profile, pam_softness)
            payload = simulate_cut_repair(
                site_seq=sequence,
                guide=guide,
                priors=priors,
                draws=draws,
                seed=seed,
                emit_sequence=True,
                pam_mask=pam_mask,
            )

            viz_mode = (shared.viz_mode or "rail_3d").lower()
            if viz_mode == "orbitals_3d":
                spec = build_crispr_orbitals_3d_spec(sequence, guide, payload)
            else:
                spec = build_crispr_rail_3d_spec(sequence, guide, payload)
            if spec is None:
                spec = build_crispr_viz_spec(sequence, guide, payload)
            if spec is None:
                raise ValueError("Failed to build CRISPR visualization spec.")

            meta = spec.metadata or {}
            workflow = build_workflow2p5d_crispr(sequence, guide, payload)
            meta["workflow_view"] = _workflow_to_serializable(workflow)
            spec.metadata = meta

            return payload, spec

        self._start_worker(task)

    def _start_worker(self, task: Callable[[], tuple[dict, Optional[EditVisualizationSpec]]]) -> None:
        if self._active_thread is not None:
            self.logMessage.emit("Simulation already running. Please wait for it to finish.")
            return

        self._simulate_btn.setEnabled(False)

        worker = SimJobWorker(task)
        thread = QThread(self)

        # Keep references so Python can't collect them early.
        self._active_thread = thread
        self._active_worker = worker

        worker.moveToThread(thread)

        # Success/failure handling
        worker.finished.connect(self._on_worker_success, Qt.QueuedConnection)
        worker.failed.connect(self._on_worker_failure, Qt.QueuedConnection)

        # Clean up worker + thread
        worker.finished.connect(worker.deleteLater)
        worker.failed.connect(worker.deleteLater)

        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)

        thread.finished.connect(self._on_thread_finished)
        thread.finished.connect(thread.deleteLater)

        thread.started.connect(worker.run)

        thread.start()


    def _on_thread_finished(self) -> None:
        self._active_thread = None
        self._active_worker = None
        self._simulate_btn.setEnabled(True)
        self.logMessage.emit("Simulation finished.")

    def _on_worker_success(self, payload: dict, spec: EditVisualizationSpec) -> None:
        self._last_payload = payload
        self._last_spec = spec
        self._save_sim_btn.setEnabled(True)
        self._save_spec_btn.setEnabled(True)
        self.simPayloadReady.emit(payload)
        self.specReady.emit(spec)
        self.runCompleted.emit("CRISPR", payload, spec)
        self.logMessage.emit("CRISPR simulation complete.")
        outcomes = payload.get("outcomes") or []
        top = outcomes[0] if outcomes else None
        LOGGER.debug("CRISPR top outcome: %r", top)
        LOGGER.debug("CRISPR priors (keys): %r", list(payload.get("priors", {}).keys()))
        LOGGER.debug("CRISPR spec.edit_type: %s", getattr(spec, "edit_type", "<none>"))

    def _on_worker_failure(self, message: str) -> None:
        QMessageBox.critical(self, "CRISPR simulation failed", message)
        self.logMessage.emit(f"CRISPR simulation failed: {message}")

    def shutdown(self) -> None:
        if self._active_thread is not None:
            self._active_thread.quit()
            self._active_thread.wait()
            self._active_thread = None
        self._simulate_btn.setEnabled(True)

    def _export_json(self, title: str, payload_getter: Callable[[], Optional[dict]]) -> None:
        payload = payload_getter()
        if not payload:
            self.logMessage.emit(f"No {title} available yet.")
            return
        path_str, _ = QFileDialog.getSaveFileName(self, f"Save {title}", filter="JSON (*.json)")
        if not path_str:
            return
        Path(path_str).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        self.logMessage.emit(f"Saved {title} to {path_str}.")

    def _export_spec(self, spec_getter: Callable[[], Optional[EditVisualizationSpec]]) -> None:
        spec = spec_getter()
        if not spec:
            self.logMessage.emit("No viz spec available yet.")
            return
        payload = spec.to_payload()
        path_str, _ = QFileDialog.getSaveFileName(self, "Save Viz Spec", filter="JSON (*.json)")
        if not path_str:
            return
        Path(path_str).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        self.logMessage.emit(f"Saved viz spec to {path_str}.")

    def save_state(self, settings: QSettings) -> None:
        settings.setValue("crispr/guide_len", self._guide_len.value())

    def load_state(self, settings: QSettings) -> None:
        guide_len = settings.value("crispr/guide_len")
        if guide_len is not None:
            self._guide_len.setValue(int(float(guide_len)))

    def _pam_from_profile(self, profile: str) -> dict[str, str]:
        profile = profile or "SpCas9_NGG"
        preset = PAM_PROFILE_PRESETS.get(profile)
        if preset:
            pam_name = preset.get("pam")
            if pam_name:
                try:
                    return guide_module.get_pam(pam_name)
                except Exception:  # pragma: no cover - fallback to literal definitions
                    LOGGER.warning("Unknown PAM preset '%s'. Falling back to literal pattern.", pam_name)
            pattern = preset.get("pattern")
            if pattern:
                return {"pattern": pattern, "orientation": preset.get("orientation", "3prime")}
        normalized = profile.replace("_", "-")
        for token in (profile, normalized):
            if token:
                try:
                    return guide_module.get_pam(token)
                except Exception:
                    continue
        pattern = profile.replace("-", "").replace("_", "")
        if not pattern:
            pattern = "NGG"
        LOGGER.info("Treating PAM profile '%s' as literal pattern '%s'.", profile, pattern)
        return {"pattern": pattern.upper(), "orientation": "3prime"}


class RunHistoryPanel(QWidget):
    """Displays completed simulations with preview/export helpers."""

    previewRequested = Signal(object)
    logMessage = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Run History", self))
        self._list = QListWidget(self)
        layout.addWidget(self._list)
        btn_row = QHBoxLayout()
        self._preview_btn = QPushButton("Preview", self)
        self._save_sim_btn = QPushButton("Save .sim…", self)
        self._save_spec_btn = QPushButton("Save Viz Spec…", self)
        self._save_dag_btn = QPushButton("Save DAG…", self)
        self._save_dag_btn.setEnabled(False)
        btn_row.addWidget(self._preview_btn)
        btn_row.addWidget(self._save_sim_btn)
        btn_row.addWidget(self._save_spec_btn)
        btn_row.addWidget(self._save_dag_btn)
        layout.addLayout(btn_row)

        self._preview_btn.clicked.connect(self._preview_selected)
        self._save_sim_btn.clicked.connect(self._save_sim)
        self._save_spec_btn.clicked.connect(self._save_spec)
        self._save_dag_btn.clicked.connect(self._save_dag)
        self._list.itemDoubleClicked.connect(lambda _: self._preview_selected())
        self._list.currentItemChanged.connect(lambda _current, _previous: self._update_save_buttons())

    def add_entry(self, kind: str, spec: EditVisualizationSpec, payload: dict) -> None:
        label = f"[{kind}] {spec.edit_type}"
        item = QListWidgetItem(label)
        is_dag = isinstance(payload, dict) and isinstance(payload.get("artifact"), str)
        item.setData(Qt.UserRole, {"spec": spec, "payload": payload, "kind": kind, "is_dag": is_dag})
        self._list.addItem(item)
        self._list.setCurrentItem(item)
        self._update_save_buttons()

    def _current_entry(self) -> Optional[dict]:
        item = self._list.currentItem()
        if not item:
            return None
        return item.data(Qt.UserRole)

    def _preview_selected(self) -> None:
        entry = self._current_entry()
        if not entry:
            self.logMessage.emit("Select a run from the history to preview.")
            return
        spec: EditVisualizationSpec = entry["spec"]
        self.previewRequested.emit(spec)

    def _update_save_buttons(self) -> None:
        entry = self._current_entry()
        is_dag = bool(entry and entry.get("is_dag"))
        self._save_dag_btn.setEnabled(is_dag)

    def _save_sim(self) -> None:
        entry = self._current_entry()
        if not entry:
            self.logMessage.emit("Select a run before saving a simulation payload.")
            return
        payload = entry["payload"]
        path_str, _ = QFileDialog.getSaveFileName(self, "Save Simulation JSON", filter="JSON (*.json)")
        if not path_str:
            return
        Path(path_str).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        self.logMessage.emit(f"Saved simulation payload to {path_str}.")

    def _save_spec(self) -> None:
        entry = self._current_entry()
        if not entry:
            self.logMessage.emit("Select a run before saving a viz spec.")
            return
        spec: EditVisualizationSpec = entry["spec"]
        path_str, _ = QFileDialog.getSaveFileName(self, "Save Viz Spec", filter="JSON (*.json)")
        if not path_str:
            return
        Path(path_str).write_text(json.dumps(spec.to_payload(), indent=2) + "\n", encoding="utf-8")
        self.logMessage.emit(f"Saved viz spec to {path_str}.")

    def _save_dag(self) -> None:
        entry = self._current_entry()
        if not entry:
            self.logMessage.emit("Select a run before saving an edit DAG payload.")
            return
        payload = entry["payload"]
        if not isinstance(payload, dict) or "artifact" not in payload:
            self.logMessage.emit("Selected run does not contain an edit DAG payload.")
            return
        path_str, _ = QFileDialog.getSaveFileName(
            self,
            "Save Edit DAG",
            filter="Edit DAG (*.edit_dag.json);;JSON (*.json)",
        )
        if not path_str:
            return
        Path(path_str).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        self.logMessage.emit(f"Saved edit DAG payload to {path_str}.")


class PrimeTab(QWidget):
    """Prime editing form."""

    specReady = Signal(object)
    simPayloadReady = Signal(object)
    logMessage = Signal(str)
    runCompleted = Signal(str, object, object)

    def __init__(self, sequence_input: SequenceInput, sim_settings: SimSettingsModel, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._sequence_input = sequence_input
        self._sim_settings = sim_settings
        self._last_payload: Optional[dict] = None
        self._last_spec: Optional[EditVisualizationSpec] = None
        self._active_thread: QThread | None = None
        self._active_worker: SimJobWorker | None = None

        layout = QVBoxLayout(self)
        form = QFormLayout()
        self._peg_name = QLineEdit("peg-demo", self)
        self._spacer = QLineEdit(self)
        self._pbs = QLineEdit(self)
        self._rtt = QLineEdit(self)
        form.addRow("pegRNA name", self._peg_name)
        form.addRow("Spacer", self._spacer)
        form.addRow("PBS", self._pbs)
        form.addRow("RTT", self._rtt)
        layout.addLayout(form)

        editor_form = QFormLayout()
        self._editor_name = QLineEdit("PrimeEditor", self)
        self._nick_offset = QSpinBox(self)
        self._nick_offset.setRange(-50, 50)
        self._nick_offset.setValue(0)
        self._efficiency = QDoubleSpinBox(self)
        self._efficiency.setRange(0.0, 2.0)
        self._efficiency.setSingleStep(0.1)
        self._efficiency.setValue(0.8)
        self._indel_bias = QDoubleSpinBox(self)
        self._indel_bias.setRange(0.0, 1.0)
        self._indel_bias.setSingleStep(0.05)
        self._indel_bias.setValue(0.2)
        self._mismatch_tol = QSpinBox(self)
        self._mismatch_tol.setRange(0, 20)
        self._mismatch_tol.setValue(3)
        self._max_outcomes = QSpinBox(self)
        self._max_outcomes.setRange(1, 64)
        self._max_outcomes.setValue(8)

        editor_form.addRow("Editor name", self._editor_name)
        editor_form.addRow("Nick offset", self._nick_offset)
        editor_form.addRow("Efficiency scale", self._efficiency)
        editor_form.addRow("Indel bias", self._indel_bias)
        editor_form.addRow("Mismatch tolerance", self._mismatch_tol)
        editor_form.addRow("Max outcomes", self._max_outcomes)
        layout.addLayout(editor_form)

        btn_row = QHBoxLayout()
        self._simulate_btn = QPushButton("Simulate Prime Edit", self)
        self._save_sim_btn = QPushButton("Save prime.edit_sim…", self)
        self._save_spec_btn = QPushButton("Save Viz Spec…", self)
        self._save_sim_btn.setEnabled(False)
        self._save_spec_btn.setEnabled(False)
        btn_row.addWidget(self._simulate_btn)
        btn_row.addWidget(self._save_sim_btn)
        btn_row.addWidget(self._save_spec_btn)
        layout.addLayout(btn_row)
        layout.addStretch(1)

        self._simulate_btn.clicked.connect(self._simulate)
        self._save_sim_btn.clicked.connect(partial(self._export_json, "prime.edit_sim json", lambda: self._last_payload))
        self._save_spec_btn.clicked.connect(partial(self._export_spec, lambda: self._last_spec))

    def _sequence_or_warn(self) -> Optional[str]:
        sequence = self._sequence_input.current_sequence()
        if not sequence:
            self.logMessage.emit("Provide a sequence before running prime simulations.")
            return None
        return sequence

    def _simulate(self) -> None:
        sequence = self._sequence_or_warn()
        if not sequence:
            return
        spacer = self._spacer.text().strip()
        pbs = self._pbs.text().strip()
        rtt = self._rtt.text().strip()
        if not spacer or not pbs or not rtt:
            self.logMessage.emit("Spacer, PBS, and RTT fields are required for prime simulations.")
            return

        peg = PegRNA(spacer=spacer, pbs=pbs, rtt=rtt, name=self._peg_name.text().strip() or None)
        cas = CasSystem(
            name="SimBuilder-Cas9",
            system_type=CasSystemType.CAS9,
            pam_rules=[PAMRule(pattern="NGG", description="SpCas9")],
            cut_offset=3,
        )
        editor = PrimeEditor(
            name=self._editor_name.text().strip() or "PrimeEditor",
            cas=cas,
            nick_to_edit_offset=self._nick_offset.value(),
            efficiency_scale=self._efficiency.value(),
            indel_bias=self._indel_bias.value(),
            mismatch_tolerance=self._mismatch_tol.value(),
        )
        genome = DigitalGenome({"chr": sequence})
        shared = self._sim_settings.value
        draws = max(1, shared.draws)
        seed = shared.seed
        pam_profile = shared.pam_profile or "SpCas9_NGG"
        pam_softness = shared.pam_softness
        priors_profile = shared.priors_profile or "default_indel"
        seed_label = "auto" if seed is None else str(seed)
        self.logMessage.emit(
            "Running prime editing simulation "
            f"(draws={draws}, seed={seed_label}, priors={priors_profile}, PAM={pam_profile})…"
        )

        def task() -> tuple[dict, Optional[EditVisualizationSpec]]:
            priors = resolve_prime_priors(priors_profile)
            pam_mask = build_prime_pam_mask(sequence, peg, pam_profile, pam_softness)
            payload = simulate_prime_edit(
                site_seq=sequence,
                peg=peg,
                priors=priors,
                draws=draws,
                seed=seed,
                emit_sequence=True,
                pam_mask=pam_mask,
            )
            viz_mode = (shared.viz_mode or "").lower()
            if viz_mode == "prime_scaffold_3d":
                spec = build_prime_scaffold_3d_spec(genome, peg, editor, payload)
            else:
                spec = build_prime_viz_spec(genome, peg, editor, payload)
            if spec is None:
                raise ValueError("Failed to build prime viz spec.")

            workflow = build_workflow2p5d_prime(sequence, payload)
            meta = spec.metadata or {}
            meta["workflow_view"] = _workflow_to_serializable(workflow)
            spec.metadata = meta

            return payload, spec

        self._start_worker(task)

    def _start_worker(self, task: Callable[[], tuple[dict, Optional[EditVisualizationSpec]]]) -> None:
        if self._active_thread is not None:
            self.logMessage.emit("Simulation already running. Please wait for it to finish.")
            return
        self._simulate_btn.setEnabled(False)
        def safe_task():
            try:
                return task()
            except Exception as exc:
                LOGGER.exception("Prime simulation failed:")
                raise

        worker = SimJobWorker(safe_task)
        thread = QThread(self)
        self._active_thread = thread
        self._active_worker = worker

        worker.moveToThread(thread)
        worker.finished.connect(self._on_worker_success, Qt.QueuedConnection)
        worker.failed.connect(self._on_worker_failure, Qt.QueuedConnection)
        worker.finished.connect(worker.deleteLater)
        worker.failed.connect(worker.deleteLater)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(self._on_thread_finished)
        thread.finished.connect(thread.deleteLater)
        thread.started.connect(worker.run)
        thread.start()

    def _on_thread_finished(self) -> None:
        self._active_thread = None
        self._active_worker = None
        self._simulate_btn.setEnabled(True)

    def _on_worker_success(self, payload: dict, spec: EditVisualizationSpec) -> None:
        self._last_payload = payload
        self._last_spec = spec
        self._save_sim_btn.setEnabled(True)
        self._save_spec_btn.setEnabled(True)
        self.simPayloadReady.emit(payload)
        self.specReady.emit(spec)
        self.runCompleted.emit("PRIME", payload, spec)
        self.logMessage.emit("Prime simulation complete.")
        top = payload.get("outcomes", [{}])
        top_label = top[0].get("label") if top and isinstance(top[0], dict) else None
        LOGGER.debug(
            "Prime sim outcome: top=%s priors=%s",
            top_label,
            list(payload.get("priors", {}).keys()),
        )

    def _on_worker_failure(self, message: str) -> None:
        QMessageBox.critical(self, "Prime simulation failed", message)
        self.logMessage.emit(f"Prime simulation failed: {message}")

    def shutdown(self) -> None:
        if self._active_thread is not None:
            self._active_thread.quit()
            self._active_thread.wait()
            self._active_thread = None
        self._simulate_btn.setEnabled(True)

    def _export_json(self, title: str, payload_getter: Callable[[], Optional[dict]]) -> None:
        payload = payload_getter()
        if not payload:
            self.logMessage.emit(f"No {title} available yet.")
            return
        path_str, _ = QFileDialog.getSaveFileName(self, f"Save {title}", filter="JSON (*.json)")
        if not path_str:
            return
        Path(path_str).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        self.logMessage.emit(f"Saved {title} to {path_str}.")

    def _export_spec(self, spec_getter: Callable[[], Optional[EditVisualizationSpec]]) -> None:
        spec = spec_getter()
        if not spec:
            self.logMessage.emit("No viz spec available yet.")
            return
        payload = spec.to_payload()
        path_str, _ = QFileDialog.getSaveFileName(self, "Save Viz Spec", filter="JSON (*.json)")
        if not path_str:
            return
        Path(path_str).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        self.logMessage.emit(f"Saved viz spec to {path_str}.")

    def save_state(self, settings: QSettings) -> None:
        settings.setValue("prime/peg_name", self._peg_name.text())
        settings.setValue("prime/spacer", self._spacer.text())
        settings.setValue("prime/pbs", self._pbs.text())
        settings.setValue("prime/rtt", self._rtt.text())
        settings.setValue("prime/editor_name", self._editor_name.text())
        settings.setValue("prime/nick_offset", self._nick_offset.value())
        settings.setValue("prime/efficiency", self._efficiency.value())
        settings.setValue("prime/indel_bias", self._indel_bias.value())
        settings.setValue("prime/mismatch_tol", self._mismatch_tol.value())
        settings.setValue("prime/max_outcomes", self._max_outcomes.value())

    def load_state(self, settings: QSettings) -> None:
        peg_name = settings.value("prime/peg_name")
        if peg_name is not None:
            self._peg_name.setText(str(peg_name))
        spacer = settings.value("prime/spacer")
        if spacer is not None:
            self._spacer.setText(str(spacer))
        pbs = settings.value("prime/pbs")
        if pbs is not None:
            self._pbs.setText(str(pbs))
        rtt = settings.value("prime/rtt")
        if rtt is not None:
            self._rtt.setText(str(rtt))
        editor_name = settings.value("prime/editor_name")
        if editor_name is not None:
            self._editor_name.setText(str(editor_name))
        nick = settings.value("prime/nick_offset")
        if nick is not None:
            self._nick_offset.setValue(int(float(nick)))
        efficiency = settings.value("prime/efficiency")
        if efficiency is not None:
            self._efficiency.setValue(float(efficiency))
        indel = settings.value("prime/indel_bias")
        if indel is not None:
            self._indel_bias.setValue(float(indel))
        mismatch = settings.value("prime/mismatch_tol")
        if mismatch is not None:
            self._mismatch_tol.setValue(int(float(mismatch)))
        max_outcomes = settings.value("prime/max_outcomes")
        if max_outcomes is not None:
            self._max_outcomes.setValue(int(float(max_outcomes)))


class SimBuilderWindow(QMainWindow):
    """Main window combining forms, logs, and the ModernGL viewer."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Helix Sim Builder")
        central = QWidget(self)
        self.setCentralWidget(central)
        root_layout = QHBoxLayout(central)
        self.sim_settings = SimSettingsModel(parent=self)

        left = QWidget(self)
        left_layout = QVBoxLayout(left)
        self.sequence_input = SequenceInput(self)
        left_layout.addWidget(self.sequence_input)
        self.pam_config = PamConfigPanel(self.sim_settings, self)
        left_layout.addWidget(self.pam_config)

        self.tabs = QTabWidget(self)
        self.crispr_tab = CrisprTab(self.sequence_input, self.sim_settings, self)
        self.prime_tab = PrimeTab(self.sequence_input, self.sim_settings, self)
        self.tabs.addTab(self.crispr_tab, "CRISPR Sim")
        self.tabs.addTab(self.prime_tab, "Prime Sim")
        left_layout.addWidget(self.tabs)

        self.log_panel = LogPanel(self)
        self.log_panel.setMaximumHeight(150)
        left_layout.addWidget(self.log_panel)

        self.gl_window = HelixModernWidget()
        gl_container = QWidget.createWindowContainer(self.gl_window, self)
        gl_container.setMinimumWidth(600)
        gl_container.setMinimumHeight(320)
        gl_container.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.viewer = gl_container
        self.gl_window.frameDrawn.connect(self._on_viewer_frame)

        self.rail_widget = RailBandWidget(self)
        self.viewer_placeholder = ViewerPlaceholder(self)
        self.viewer_container = QWidget(self)
        self.viewer_container.setAutoFillBackground(False)
        self.viewer_stack = QStackedLayout(self.viewer_container)
        self.viewer_stack.addWidget(self.viewer_placeholder)
        self.viewer_stack.addWidget(self.rail_widget)
        self.viewer_stack.addWidget(self.viewer)
        self.viewer_stack.setCurrentWidget(self.viewer_placeholder)
        self._active_view_widget: QWidget = self.viewer_placeholder
        right = QWidget(self)
        right_layout = QVBoxLayout(right)
        self.right_splitter = QSplitter(Qt.Vertical, right)
        self.right_splitter.addWidget(self.viewer_container)
        self.history_panel = RunHistoryPanel(self)
        self.right_splitter.addWidget(self.history_panel)
        self.right_splitter.setStretchFactor(0, 3)
        self.right_splitter.setStretchFactor(1, 2)
        self.right_splitter.setCollapsible(0, False)
        right_layout.addWidget(self.right_splitter)
        controls = QHBoxLayout()
        controls.addStretch(1)
        self.viewer_reset_btn = QPushButton("Reset Camera", self)
        self.viewer_reset_btn.setToolTip("Return the 3D view to its default orbit.")
        self.viewer_reset_btn.clicked.connect(self._reset_viewer_camera)
        controls.addWidget(self.viewer_reset_btn)
        controls.addStretch(1)
        right_layout.addLayout(controls)

        self.main_splitter = QSplitter(Qt.Horizontal, central)
        self.main_splitter.addWidget(left)
        self.main_splitter.addWidget(right)
        self.main_splitter.setStretchFactor(0, 1)
        self.main_splitter.setStretchFactor(1, 2)
        self.main_splitter.setCollapsible(0, False)
        self.main_splitter.setCollapsible(1, False)
        root_layout.addWidget(self.main_splitter)

        self._settings = QSettings("Helix", "SimBuilder")
        self._load_settings()

        self.crispr_tab.specReady.connect(self._display_spec)
        self.prime_tab.specReady.connect(self._display_spec)
        self.history_panel.previewRequested.connect(self._display_spec)
        self.crispr_tab.logMessage.connect(self._log)
        self.prime_tab.logMessage.connect(self._log)
        self.history_panel.logMessage.connect(self._log)
        self.crispr_tab.runCompleted.connect(self._record_run)
        self.prime_tab.runCompleted.connect(self._record_run)
        self._show_viewer_canvas(False, "Viewer initializing…")

        # DEBUG bootstrap: load a trivial CRISPR spec so the viewer renders something on startup.
        try:
            fake_payload = {
                "schema": {"kind": "crispr.sim", "spec_version": "1.0"},
                "meta": {},
                "site": {"length": 40, "sequence": "ACGT" * 10},
                "guide": {"id": "debug", "start": 10, "end": 30, "strand": "+", "gc_content": 0.5},
                "priors": {},
                "draws": 1,
                "outcomes": [
                    {
                        "label": "fake_del",
                        "count": 1,
                        "probability": 1.0,
                        "diff": {
                            "kind": "deletion",
                            "start": 15,
                            "end": 18,
                        },
                    }
                ],
            }
            mode = (self.sim_settings.value.viz_mode or "rail_3d").lower()
            seq = fake_payload["site"]["sequence"]
            guide = fake_payload["guide"]
            if mode == "orbitals_3d":
                fake_spec = build_crispr_orbitals_3d_spec(seq, guide, fake_payload)
            else:
                fake_spec = build_crispr_rail_3d_spec(seq, guide, fake_payload)
            if fake_spec is None:
                fake_spec = build_crispr_viz_spec(seq, guide, fake_payload)
            if fake_spec and mode == "workflow_2p5d":
                workflow = build_workflow2p5d_crispr(seq, guide, fake_payload)
                meta = fake_spec.metadata or {}
                meta["workflow_view"] = _workflow_to_serializable(workflow)
                fake_spec.metadata = meta
            def _load_debug_spec(spec=fake_spec) -> None:
                self._display_spec(spec)
                self._log("Loaded DEBUG visualization: " + spec.edit_type)
            QTimer.singleShot(0, _load_debug_spec)
        except Exception as exc:
            self._show_viewer_canvas(False, f"Failed to load debug viz: {exc}")
            self._log(f"Failed to load debug spec: {exc}")

    def _display_spec(self, spec: EditVisualizationSpec) -> None:
        if spec is None:
            self._show_viewer_canvas(False, "No visualization available for this run.")
            return
        metadata = spec.metadata or {}
        viz_mode = (self.sim_settings.value.viz_mode or "rail_3d").lower()
        workflow_meta = metadata.get("workflow_view")
        rail_meta = metadata.get("rail_bands")

        # Always load the spec into the GL widget so mode switches stay in sync.
        self.gl_window.set_spec(spec)

        if viz_mode == "rail_2d":
            if rail_meta and rail_meta.get("bands"):
                self.rail_widget.set_spec(spec)
                self.gl_window.set_workflow_view(None)
                self._show_viewer_canvas(True, widget=self.rail_widget)
                self._log(f"Loaded visualization: {spec.edit_type} (rail 2D)")
                return
            else:
                self._log("Rail 2D view unavailable for this spec. Falling back to 3D.")

        workflow_buffers = None
        if viz_mode == "workflow_2p5d":
            if workflow_meta:
                workflow_buffers = _workflow_from_serializable(workflow_meta)
            else:
                self._log("Workflow view unavailable for this spec. Showing 3D view instead.")

        self.gl_window.set_workflow_view(workflow_buffers)  # type: ignore[arg-type]
        self._show_viewer_canvas(True, widget=self.viewer)
        self.gl_window.update()
        self._log(f"Loaded visualization: {spec.edit_type}")

    def _record_run(self, kind: str, payload: dict, spec: EditVisualizationSpec) -> None:
        self.history_panel.add_entry(kind, spec, payload)
        self._log(f"Recorded {kind} run: {spec.edit_type}")

    def _log(self, message: str) -> None:
        self.log_panel.log(message)
        self._console_log(message)

    def _show_viewer_canvas(self, ready: bool, message: str | None = None, widget: QWidget | None = None) -> None:
        if ready:
            target = widget or self.viewer
            self.viewer_stack.setCurrentWidget(target)
            self._active_view_widget = target
        else:
            if message:
                self.viewer_placeholder.set_message(message)
            self.viewer_stack.setCurrentWidget(self.viewer_placeholder)
            self._active_view_widget = self.viewer_placeholder

    def _on_viewer_frame(self) -> None:
        if not hasattr(self, "_gl_logged"):
            self._gl_logged = True
            try:
                ctx = self.gl_window.context()
                if ctx is not None:
                    fmt = ctx.format()
                    profile_map = {
                        QSurfaceFormat.NoProfile: "none",
                        QSurfaceFormat.CoreProfile: "core",
                        QSurfaceFormat.CompatibilityProfile: "compat",
                    }
                    profile_name = profile_map.get(fmt.profile(), str(fmt.profile()))
                    self._log(
                        "GL context: "
                        f"{fmt.majorVersion()}.{fmt.minorVersion()} "
                        f"profile={profile_name} "
                        f"samples={fmt.samples()} "
                        f"depth={fmt.depthBufferSize()} "
                        f"swap={fmt.swapBehavior()}"
                    )
                else:
                    self._log("GL context unavailable (no QOpenGLContext).")
            except Exception as exc:
                self._log(f"GL context query failed: {exc}")
        if getattr(self, "_active_view_widget", None) is self.viewer:
            self._show_viewer_canvas(True, widget=self.viewer)

    def _reset_viewer_camera(self) -> None:
        self.gl_window.reset_camera()
        self._log("Viewer camera reset.")

    def _load_settings(self) -> None:
        geometry = self._settings.value("geometry", type=QByteArray)
        if geometry and isinstance(geometry, QByteArray):
            self.restoreGeometry(geometry)
        else:
            self.resize(1400, 900)
        window_state = self._settings.value("windowState", type=QByteArray)
        if window_state and isinstance(window_state, QByteArray):
            self.restoreState(window_state)
        splitter_state = self._settings.value("main_splitter_state", type=QByteArray)
        if not splitter_state or not isinstance(splitter_state, QByteArray) or not self.main_splitter.restoreState(splitter_state):
            self.main_splitter.setSizes([450, 950])
        right_splitter_state = self._settings.value("right_splitter_state", type=QByteArray)
        if not right_splitter_state or not isinstance(right_splitter_state, QByteArray) or not self.right_splitter.restoreState(right_splitter_state):
            self.right_splitter.setSizes([600, 300])
        sizes = self.right_splitter.sizes()
        if not sizes or len(sizes) < 2 or sizes[0] < 100:
            self.right_splitter.setSizes([600, 300])
        main_sizes = self.main_splitter.sizes()
        if not main_sizes or len(main_sizes) < 2 or main_sizes[0] < 200:
            self.main_splitter.setSizes([450, 950])
        sim_settings_blob = self._settings.value("sim_settings")
        if sim_settings_blob:
            try:
                if isinstance(sim_settings_blob, (bytes, bytearray)):
                    raw = bytes(sim_settings_blob).decode("utf-8")
                else:
                    raw = str(sim_settings_blob)
                payload = json.loads(raw)
            except Exception as exc:
                LOGGER.warning("Failed to restore sim settings: %s", exc)
            else:
                if isinstance(payload, dict):
                    self.sim_settings.apply_dict(payload)
        sequence_text = self._settings.value("sequence_text")
        if sequence_text:
            self.sequence_input.set_sequence(str(sequence_text))
        self.crispr_tab.load_state(self._settings)
        self.prime_tab.load_state(self._settings)

    def _save_settings(self) -> None:
        self._settings.setValue("geometry", self.saveGeometry())
        self._settings.setValue("windowState", self.saveState())
        self._settings.setValue("main_splitter_state", self.main_splitter.saveState())
        self._settings.setValue("right_splitter_state", self.right_splitter.saveState())
        self._settings.setValue("sequence_text", self.sequence_input.sequence())
        self._settings.setValue("sim_settings", json.dumps(self.sim_settings.to_dict()))
        self.crispr_tab.save_state(self._settings)
        self.prime_tab.save_state(self._settings)

    def closeEvent(self, event: QCloseEvent) -> None:
        self.crispr_tab.shutdown()
        self.prime_tab.shutdown()
        self._save_settings()
        self._settings.sync()
        super().closeEvent(event)

    @staticmethod
    def _console_log(message: str) -> None:
        print(f"[SimBuilder] {message}", flush=True)


def run_sim_builder() -> None:
    """Launch the Sim Builder UI."""

    app = QApplication.instance()
    owns_app = False
    if app is None:
        app = QApplication([])
        owns_app = True
    fmt = QSurfaceFormat()
    fmt.setRenderableType(QSurfaceFormat.OpenGL)
    fmt.setProfile(QSurfaceFormat.CoreProfile)
    fmt.setVersion(3, 3)
    fmt.setSwapBehavior(QSurfaceFormat.DoubleBuffer)
    fmt.setDepthBufferSize(24)
    fmt.setStencilBufferSize(8)
    QSurfaceFormat.setDefaultFormat(fmt)
    apply_helix_theme(app)
    window = SimBuilderWindow()
    window.show()
    if owns_app:
        app.exec()
class HelixComboBox(QComboBox):
    """ComboBox tuned for Sim Builder interactions."""

    def __init__(self, parent: QWidget | None = None, *, searchable: bool = True) -> None:
        super().__init__(parent)
        self._searchable = searchable
        self._init_completer()

    def _init_completer(self) -> None:
        if not self._searchable:
            return
        completer = QCompleter(self.model(), self)
        completer.setCompletionMode(QCompleter.PopupCompletion)
        completer.setFilterMode(Qt.MatchContains)
        completer.setCaseSensitivity(Qt.CaseInsensitive)
        self.setCompleter(completer)

    def wheelEvent(self, event) -> None:  # type: ignore[override]
        view = self.view()
        if view is not None and view.isVisible():
            super().wheelEvent(event)
        else:
            event.ignore()

    def showPopup(self) -> None:  # type: ignore[override]
        super().showPopup()
        view = self.view()
        if view is not None:
            width = max(self.width(), view.sizeHintForColumn(0) + 40)
            view.setMinimumWidth(width)
