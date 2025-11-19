"""Pure Python core of the PyQt + ModernGL helix visualizer."""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, Iterable, List, Mapping

import numpy as np

from .spec import EditEvent, EditEventType, EditVisualizationSpec


@dataclass(slots=True)
class HelixGeometry:
    """Vertex bundles describing the primary helix rails."""

    rail_vertices: np.ndarray  # (N, 3) float32 positions
    rail_ids: np.ndarray  # (N,) float32 0 for strand A, 1 for strand B
    param_coords: np.ndarray  # (N,) float32 [0, 1] along helix
    base_centers: np.ndarray  # (B, 3) float32
    base_indices: np.ndarray  # (B,) int32 indexes into sequence


@dataclass(slots=True)
class RepairArcGeometry:
    """A packed vertex buffer representing independent repair arcs."""

    vertices: np.ndarray  # (M, 3) float32
    offsets: np.ndarray  # (A + 1,) int32 prefix sums for arc segments
    event_types: list[EditEventType]
    event_times: list[float]
    weights: np.ndarray | None = None
    kinds: np.ndarray | None = None


class AnimationPhase(str, Enum):
    """Named beats aligned with the user's storyboard."""

    APPROACH = "approach"
    RECOGNITION = "recognition"
    CUT = "cut"
    PRIME = "prime"
    REPAIR = "repair"
    RESOLVE = "resolve"


@dataclass(slots=True)
class AnimationState:
    """Mutable animation context shared by the Qt shell and shader uniforms."""

    time: float = 0.0
    phase: AnimationPhase = AnimationPhase.APPROACH
    loop_duration: float = 12.0


class HelixVizEngine:
    """Constructs helix/base/arc geometry and tracks the animation timeline."""

    def __init__(
        self,
        spec: EditVisualizationSpec,
        *,
        bases_per_turn: float = 10.0,
        segments_per_base: int = 5,
        radius: float = 0.38,
        pitch: float = 0.08,
    ) -> None:
        self.spec = spec
        self.bases_per_turn = bases_per_turn
        self.segments_per_base = max(1, segments_per_base)
        self.radius = radius
        self.pitch = pitch
        # Loop duration controls how long ribbons and phase beats stay
        # visible in realtime. Keep a reasonably long loop so edits
        # linger, with a small tail after the last event.
        self.loop_duration = max(12.0, spec.duration + 2.0)
        self.state = AnimationState(time=0.0, phase=AnimationPhase.APPROACH, loop_duration=self.loop_duration)
        self._base_positions: Dict[int, np.ndarray] = {}
        custom_geom = (spec.metadata or {}).get("custom_geometry") if spec.metadata else None
        if isinstance(custom_geom, Mapping):
            self.helix_geometry = self._geometry_from_custom(custom_geom)
            self.repair_geometry = self._repair_geometry_from_custom(custom_geom)
        else:
            self.helix_geometry = self._build_helix_geometry()
            self.repair_geometry = self._build_repair_geometry()
        self._phase_thresholds = self._build_phase_thresholds()

    # ------------------------------------------------------------------
    # Geometry builders
    def _build_helix_geometry(self) -> HelixGeometry:
        sequence = self.spec.sequence
        base_count = len(sequence)
        theta_step = 2.0 * math.pi / self.bases_per_turn
        total_samples = base_count * self.segments_per_base
        # Parameter along helix from 0..1 to keep consistent camera scale.
        samples = np.linspace(0.0, base_count - 1, total_samples, dtype=np.float32)
        thetas = samples * theta_step
        z = (samples / max(1, base_count - 1)) * 2.0 - 1.0

        def _rail_vertices(phase_offset: float) -> np.ndarray:
            theta = thetas + phase_offset
            x = self.radius * np.cos(theta)
            y = self.radius * np.sin(theta)
            return np.stack([x, y, z], axis=1).astype("f4")

        rail_a = _rail_vertices(0.0)
        rail_b = _rail_vertices(math.pi)
        rail_vertices = np.concatenate([rail_a, rail_b], axis=0)
        rail_ids = np.concatenate([
            np.zeros(len(rail_a), dtype="f4"),
            np.ones(len(rail_b), dtype="f4"),
        ])
        param_coords = np.concatenate([samples, samples], axis=0)
        if param_coords.size > 0:
            span = float(np.ptp(param_coords))
            param_coords = (param_coords - param_coords.min()) / max(1e-6, span)

        # Cache per-base centers for arc anchoring.
        base_centers = []
        base_indices = []
        for idx in range(base_count):
            theta = theta_step * idx
            x = self.radius * math.cos(theta)
            y = self.radius * math.sin(theta)
            zc = float(idx) / max(1, base_count - 1) * 2.0 - 1.0
            pos = np.array([x, y, zc], dtype="f4")
            self._base_positions[idx] = pos
            base_centers.append(pos)
            base_indices.append(idx)
        base_centers_arr = np.stack(base_centers, axis=0) if base_centers else np.zeros((0, 3), dtype="f4")
        base_idx_arr = np.array(base_indices, dtype="i4") if base_indices else np.zeros((0,), dtype="i4")
        return HelixGeometry(
            rail_vertices=rail_vertices,
            rail_ids=rail_ids,
            param_coords=param_coords.astype("f4"),
            base_centers=base_centers_arr,
            base_indices=base_idx_arr,
        )

    def _build_repair_geometry(self) -> RepairArcGeometry:
        arc_vertices: List[np.ndarray] = []
        offsets: List[int] = [0]
        event_types: List[EditEventType] = []
        event_times: List[float] = []
        for event in self.spec.events:
            anchor_idx = event.index
            if anchor_idx is None:
                anchor_idx = event.start
            if anchor_idx is None:
                continue
            if anchor_idx not in self._base_positions:
                continue
            if event.type not in {
                EditEventType.RECOGNITION,
                EditEventType.RT_SYNTHESIS,
                EditEventType.FLAP_RESOLUTION,
                EditEventType.REPAIR_COMPLETE,
            }:
                continue
            arc = self._arc_for_base(anchor_idx, event)
            arc_vertices.append(arc)
            offsets.append(offsets[-1] + len(arc))
            event_types.append(event.type)
            event_times.append(event.t)
        if arc_vertices:
            vertex_blob = np.vstack(arc_vertices).astype("f4")
        else:
            vertex_blob = np.zeros((0, 3), dtype="f4")
        return RepairArcGeometry(
            vertices=vertex_blob,
            offsets=np.array(offsets, dtype="i4"),
            event_types=event_types,
            event_times=event_times,
        )

    def _arc_for_base(self, base_idx: int, event: EditEvent) -> np.ndarray:
        base_pos = self._base_positions[base_idx]
        radial = base_pos.copy()
        radial[2] = 0.0
        norm = np.linalg.norm(radial)
        if norm < 1e-6:
            radial = np.array([1.0, 0.0, 0.0], dtype="f4")
        else:
            radial /= norm
        arc_height = 0.25 + 0.05 * (event.length or 1)
        control = base_pos + radial * arc_height + np.array([0.0, 0.0, 0.15], dtype="f4")
        samples = np.linspace(0.0, 1.0, 32, dtype="f4")
        curve = ((1 - samples)[:, None] ** 2) * base_pos + 2 * ((1 - samples)[:, None]) * (samples[:, None]) * control + (samples[:, None] ** 2) * base_pos
        return curve.astype("f4")

    def _geometry_from_custom(self, payload: Mapping[str, Any]) -> HelixGeometry:
        rv_data = payload.get("rail_vertices")
        rv = np.array(rv_data if rv_data is not None else [], dtype="f4")
        if rv.ndim == 1:
            rv = rv.reshape((-1, 3))
        rail_ids_data = payload.get("rail_ids")
        rail_ids = np.array(rail_ids_data if rail_ids_data is not None else [], dtype="f4")
        if rail_ids.size != rv.shape[0]:
            rail_ids = np.resize(rail_ids, (rv.shape[0],))
        param_data = payload.get("param_coords")
        param = np.array(param_data if param_data is not None else [], dtype="f4")
        if param.size > 0:
            span = float(np.ptp(param))
            param = (param - param.min()) / max(1e-6, span)
        else:
            param = np.zeros((rv.shape[0],), dtype="f4")
        if param.size != rv.shape[0]:
            param = np.resize(param, (rv.shape[0],))
        base_centers = np.zeros((0, 3), dtype="f4")
        base_indices = np.zeros((0,), dtype="i4")
        return HelixGeometry(
            rail_vertices=rv,
            rail_ids=rail_ids,
            param_coords=param,
            base_centers=base_centers,
            base_indices=base_indices,
        )

    def _repair_geometry_from_custom(self, payload: Mapping[str, Any]) -> RepairArcGeometry:
        vert_data = payload.get("arc_vertices")
        vertices = np.array(vert_data if vert_data is not None else [], dtype="f4")
        if vertices.ndim == 1:
            vertices = vertices.reshape((-1, 3))
        offsets_data = payload.get("arc_offsets")
        if offsets_data is None:
            offsets_arr = np.array([0], dtype="i4")
        else:
            offsets_arr = np.array(offsets_data, dtype="i4")
        weights_data = payload.get("arc_weights")
        kinds_data = payload.get("arc_kinds")
        weights = np.array(weights_data if weights_data is not None else [], dtype="f4")
        kinds = np.array(kinds_data if kinds_data is not None else [], dtype="f4")
        return RepairArcGeometry(
            vertices=vertices,
            offsets=offsets_arr,
            event_types=[],
            event_times=[],
            weights=weights if weights.size else None,
            kinds=kinds if kinds.size else None,
        )

    # ------------------------------------------------------------------
    # Animation timeline helpers
    def _build_phase_thresholds(self) -> List[tuple[float, AnimationPhase]]:
        thresholds: List[tuple[float, AnimationPhase]] = [(0.0, AnimationPhase.APPROACH)]
        def _record(first_event: EditEvent | None, phase: AnimationPhase) -> None:
            if first_event:
                thresholds.append((first_event.t, phase))

        events = self.spec.events
        def _first_of(types: Iterable[EditEventType]) -> EditEvent | None:
            type_set = set(types)
            for ev in events:
                if ev.type in type_set:
                    return ev
            return None

        _record(_first_of([EditEventType.RECOGNITION]), AnimationPhase.RECOGNITION)
        _record(_first_of([EditEventType.NICK_PRIMARY, EditEventType.CUT]), AnimationPhase.CUT)
        _record(_first_of([EditEventType.PRIME_INIT, EditEventType.RT_SYNTHESIS]), AnimationPhase.PRIME)
        _record(_first_of([EditEventType.FLAP_RESOLUTION]), AnimationPhase.REPAIR)
        last_event = events[-1] if events else None
        if last_event:
            thresholds.append((last_event.t, AnimationPhase.RESOLVE))
        thresholds = sorted({pair for pair in thresholds}, key=lambda x: x[0])
        if thresholds[-1][1] != AnimationPhase.RESOLVE:
            thresholds.append((self.loop_duration, AnimationPhase.RESOLVE))
        return thresholds

    def phase_for_time(self, t: float) -> AnimationPhase:
        loop_t = t % self.loop_duration
        current = AnimationPhase.APPROACH
        for time_marker, phase in self._phase_thresholds:
            if loop_t >= time_marker:
                current = phase
            else:
                break
        return current

    def set_time(self, t: float) -> AnimationState:
        self.state.time = t % self.loop_duration
        self.state.phase = self.phase_for_time(self.state.time)
        return self.state

    def advance(self, dt: float) -> AnimationState:
        return self.set_time(self.state.time + dt)

    # ------------------------------------------------------------------
    # Export helpers
    def active_events_within(self, t: float, window: float = 0.3) -> list[EditEvent]:
        loop_t = t % self.loop_duration
        start = max(0.0, loop_t - window)
        end = min(self.loop_duration, loop_t + window)
        hits = [ev for ev in self.spec.events if start <= ev.t <= end]
        if hits:
            return hits
        if loop_t + window > self.loop_duration:
            wrap_hits = [ev for ev in self.spec.events if 0.0 <= ev.t <= (loop_t + window - self.loop_duration)]
            return wrap_hits
        return []

    def buffer_payload(self) -> Dict[str, np.ndarray]:
        """Return ModernGL-friendly ndarray bundles for the Qt shell."""

        return {
            "rail_vertices": self.helix_geometry.rail_vertices,
            "rail_ids": self.helix_geometry.rail_ids,
            "param_coords": self.helix_geometry.param_coords,
            "base_centers": self.helix_geometry.base_centers,
            "base_indices": self.helix_geometry.base_indices,
            "arc_vertices": self.repair_geometry.vertices,
            "arc_offsets": self.repair_geometry.offsets,
            "arc_weights": self.repair_geometry.weights if self.repair_geometry.weights is not None else np.zeros((0,), dtype="f4"),
            "arc_kinds": self.repair_geometry.kinds if self.repair_geometry.kinds is not None else np.zeros((0,), dtype="f4"),
        }

    def export_timeline(self) -> dict[str, Any]:
        """Structured dict that Unity (or other DCC tools) can ingest."""

        return {
            "duration": self.loop_duration,
            "phase_markers": [
                {"time": t, "phase": phase.value}
                for t, phase in self._phase_thresholds
            ],
            "events": [
                {
                    "time": event.t,
                    "type": event.type.value,
                    "index": event.index,
                    "start": event.start,
                    "length": event.length,
                    "metadata": dict(event.metadata),
                }
                for event in self.spec.events
            ],
        }
