from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Any, Mapping, Optional

from PySide6.QtCore import QObject, Signal

from helix.crispr.model import CasSystem, CasSystemType, DigitalGenome, PAMRule
from helix.prime.model import PegRNA, PrimeEditor


class RunKind(Enum):
    NONE = auto()
    CRISPR = auto()
    PRIME = auto()
    PCR = auto()

SESSION_SNAPSHOT_VERSION = 1


@dataclass
class ExperimentState:
    genome: Optional[DigitalGenome] = None
    peg: Optional[PegRNA] = None
    editor: Optional[PrimeEditor] = None
    genome_source: str = "inline"
    genome_uri: Optional[str] = None
    genome_hash: Optional[str] = None

    dag_from_runtime: Optional[Any] = None
    viz_spec_payload: Optional[dict[str, Any]] = None
    outcomes: list[dict[str, Any]] = field(default_factory=list)
    config: dict[str, Any] = field(default_factory=dict)

    run_kind: RunKind = RunKind.NONE
    viz_dirty: bool = True
    run_id: int = 0

    @property
    def dag(self) -> Optional[Any]:
        return self.dag_from_runtime


class SessionModel(QObject):
    stateChanged = Signal(ExperimentState)
    logMessage = Signal(str)
    statusChanged = Signal(str)
    errorOccurred = Signal(str)
    runRecorded = Signal(dict)

    def __init__(self, parent: Optional[QObject] = None) -> None:
        super().__init__(parent)
        self._state = ExperimentState()
        self._run_counter = 0

    @property
    def state(self) -> ExperimentState:
        return self._state

    def update(self, **kwargs: Any) -> None:
        changed = False
        record_snapshot = bool(kwargs.get("viz_dirty"))
        target_kind = kwargs.get("run_kind")
        for key, value in kwargs.items():
            if not hasattr(self._state, key):
                continue
            current = getattr(self._state, key)
            if current is value:
                continue
            setattr(self._state, key, value)
            changed = True
        if changed:
            if record_snapshot:
                kind = target_kind or self._state.run_kind
                if kind and kind != RunKind.NONE:
                    self._run_counter += 1
                    self._state.run_id = self._run_counter
                    snapshot = self.to_payload(timestamp=True)
                    self.stateChanged.emit(self._state)
                    self.runRecorded.emit(snapshot)
                    return
            self.stateChanged.emit(self._state)

    def append_outcome(self, outcome: dict[str, Any]) -> None:
        self._state.outcomes.append(outcome)
        self._state.viz_dirty = True
        self.stateChanged.emit(self._state)

    def log(self, msg: str) -> None:
        self.logMessage.emit(str(msg))

    def mark_viz_clean(self) -> None:
        if self._state.viz_dirty:
            self._state.viz_dirty = False

    def set_status(self, message: str) -> None:
        self.statusChanged.emit(str(message))
        self.log(str(message))

    def error(self, message: str) -> None:
        self.errorOccurred.emit(str(message))
        self.log(str(message))

    def to_payload(self, *, timestamp: bool = False) -> dict[str, Any]:
        state_payload = {
            "run_kind": self._state.run_kind.name,
            "run_id": self._state.run_id,
            "viz_dirty": self._state.viz_dirty,
            "config": deepcopy(self._state.config),
            "outcomes": deepcopy(self._state.outcomes),
            "dag_from_runtime": deepcopy(self._state.dag_from_runtime),
            "viz_spec_payload": deepcopy(self._state.viz_spec_payload),
            "genome_source": self._state.genome_source,
            "genome_uri": self._state.genome_uri,
            "genome_hash": self._state.genome_hash,
        }
        snapshot_timestamp: Optional[str] = None
        if timestamp:
            from datetime import datetime, timezone

            snapshot_timestamp = datetime.now(timezone.utc).isoformat()
        if self._state.genome is not None:
            state_payload["genome"] = deepcopy(self._state.genome.sequences)
        if self._state.peg is not None:
            state_payload["peg"] = _serialize_peg(self._state.peg)
        if self._state.editor is not None:
            state_payload["editor"] = _serialize_editor(self._state.editor)
        payload = {
            "version": SESSION_SNAPSHOT_VERSION,
            "state": state_payload,
        }
        if snapshot_timestamp:
            payload["timestamp"] = snapshot_timestamp
        return payload

    def load_payload(self, payload: Mapping[str, Any]) -> None:
        if "state" in payload:
            version = payload.get("version", 0)
            state_payload = payload.get("state", {})
        else:
            version = 0
            state_payload = payload
        state = ExperimentState()
        state.run_kind = RunKind[state_payload.get("run_kind", "NONE")]
        state.run_id = int(state_payload.get("run_id", 0))
        state.viz_dirty = bool(state_payload.get("viz_dirty", True))
        state.config = deepcopy(state_payload.get("config", {}))
        state.outcomes = deepcopy(state_payload.get("outcomes", []))
        state.dag_from_runtime = deepcopy(state_payload.get("dag_from_runtime"))
        state.viz_spec_payload = deepcopy(state_payload.get("viz_spec_payload"))
        state.genome_source = str(state_payload.get("genome_source", "inline"))
        state.genome_uri = state_payload.get("genome_uri")
        state.genome_hash = state_payload.get("genome_hash")
        genome_data = state_payload.get("genome")
        if isinstance(genome_data, Mapping):
            state.genome = DigitalGenome(sequences=dict(genome_data))
        peg_data = state_payload.get("peg")
        if isinstance(peg_data, Mapping):
            state.peg = _deserialize_peg(peg_data)
        editor_data = state_payload.get("editor")
        if isinstance(editor_data, Mapping):
            state.editor = _deserialize_editor(editor_data)
        snapshot_timestamp = payload.get("timestamp")
        if version == 0 and "timestamp" in state_payload:
            snapshot_timestamp = state_payload["timestamp"]
        if snapshot_timestamp:
            state.config.setdefault("metadata", {})["timestamp"] = snapshot_timestamp
        self._state = state
        self._run_counter = max(self._run_counter, state.run_id)
        self.stateChanged.emit(self._state)


def _serialize_peg(peg: PegRNA) -> dict[str, Any]:
    return {
        "spacer": peg.spacer,
        "pbs": peg.pbs,
        "rtt": peg.rtt,
        "name": peg.name,
        "metadata": dict(peg.metadata),
    }


def _deserialize_peg(data: Mapping[str, Any]) -> PegRNA:
    return PegRNA(
        spacer=str(data.get("spacer", "")),
        pbs=str(data.get("pbs", "")),
        rtt=str(data.get("rtt", "")),
        name=data.get("name"),
        metadata=dict(data.get("metadata", {})),
    )


def _serialize_editor(editor: PrimeEditor) -> dict[str, Any]:
    return {
        "name": editor.name,
        "nick_to_edit_offset": editor.nick_to_edit_offset,
        "efficiency_scale": editor.efficiency_scale,
        "indel_bias": editor.indel_bias,
        "mismatch_tolerance": editor.mismatch_tolerance,
        "flap_balance": editor.flap_balance,
        "reanneal_bias": editor.reanneal_bias,
        "metadata": dict(editor.metadata),
        "cas": _serialize_cas(editor.cas),
    }


def _deserialize_editor(data: Mapping[str, Any]) -> PrimeEditor:
    cas_data = data.get("cas") or {}
    cas = _deserialize_cas(cas_data) if isinstance(cas_data, Mapping) else CasSystem(
        name="Imported",
        system_type=CasSystemType.CAS9,
        pam_rules=[PAMRule(pattern="NGG")],
        cut_offset=3,
    )
    return PrimeEditor(
        name=str(data.get("name", "PrimeEditor")),
        cas=cas,
        nick_to_edit_offset=int(data.get("nick_to_edit_offset", 0)),
        efficiency_scale=float(data.get("efficiency_scale", 1.0)),
        indel_bias=float(data.get("indel_bias", 0.0)),
        mismatch_tolerance=int(data.get("mismatch_tolerance", 3)),
        flap_balance=float(data.get("flap_balance", 0.5)),
        reanneal_bias=float(data.get("reanneal_bias", 0.1)),
        metadata=dict(data.get("metadata", {})),
    )


def _serialize_cas(cas: CasSystem) -> dict[str, Any]:
    return {
        "name": cas.name,
        "system_type": cas.system_type.value,
        "pam_rules": [{"pattern": rule.pattern, "description": rule.description} for rule in cas.pam_rules],
        "cut_offset": cas.cut_offset,
        "max_mismatches": cas.max_mismatches,
        "weight_mismatch_penalty": cas.weight_mismatch_penalty,
        "weight_pam_penalty": cas.weight_pam_penalty,
    }


def _deserialize_cas(data: Mapping[str, Any]) -> CasSystem:
    rules = [
        PAMRule(pattern=str(entry.get("pattern", "")), description=str(entry.get("description", "")))
        for entry in data.get("pam_rules", [])
        if isinstance(entry, Mapping)
    ]
    system_type = data.get("system_type", CasSystemType.CAS9.value)
    return CasSystem(
        name=str(data.get("name", "CasSystem")),
        system_type=CasSystemType(system_type),
        pam_rules=rules or [PAMRule(pattern="NGG")],
        cut_offset=int(data.get("cut_offset", 3)),
        max_mismatches=int(data.get("max_mismatches", 3)),
        weight_mismatch_penalty=float(data.get("weight_mismatch_penalty", 1.0)),
        weight_pam_penalty=float(data.get("weight_pam_penalty", 2.0)),
    )
