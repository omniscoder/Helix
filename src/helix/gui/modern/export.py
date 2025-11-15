"""Unity export helpers for the ModernGL helix visualizer."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Iterable, Sequence

from .engine import HelixVizEngine
from .spec import EditVisualizationSpec, load_viz_specs


def _camera_track(duration: float, samples: int, radius: float = 2.2) -> list[dict[str, object]]:
    frames: list[dict[str, object]] = []
    for i in range(max(2, samples)):
        t = duration * (i / (samples - 1)) if samples > 1 else 0.0
        angle = t * 0.25
        frames.append(
            {
                "t": round(t, 4),
                "position": [round(radius * math.cos(angle), 4), round(radius * 0.35, 4), round(radius * math.sin(angle), 4)],
                "target": [0.0, 0.0, 0.0],
            }
        )
    return frames


def build_unity_payload(engine: HelixVizEngine, *, camera_samples: int = 120) -> dict[str, object]:
    """Convert the current engine to a Unity-friendly JSON payload."""

    timeline = engine.export_timeline()
    payload = {
        "kind": "helix.viz.timeline.v1",
        "duration": timeline["duration"],
        "camera_orbit": _camera_track(engine.loop_duration, camera_samples),
        "phase_markers": timeline["phase_markers"],
        "events": timeline["events"],
        "metadata": engine.spec.metadata,
    }
    return payload


def export_unity_json(engine: HelixVizEngine, path: Path, *, camera_samples: int = 120) -> Path:
    """Write the Unity payload to disk, creating parent folders as needed."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = build_unity_payload(engine, camera_samples=camera_samples)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def export_unity_bundle(
    specs: Sequence[EditVisualizationSpec] | Iterable[EditVisualizationSpec],
    output_dir: Path,
    *,
    camera_samples: int = 120,
) -> list[Path]:
    """Export multiple specs to JSON files for Unity ingestion."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for idx, spec in enumerate(specs):
        engine = HelixVizEngine(spec)
        label = spec.metadata.get("label") or spec.metadata.get("name") or spec.edit_type or f"edit_{idx+1}"
        slug = re.sub(r"[^a-z0-9_]+", "-", label.lower()).strip("-") or f"edit_{idx+1}"
        filename = output_dir / f"{slug}.json"
        export_unity_json(engine, filename, camera_samples=camera_samples)
        paths.append(filename)
    return paths


def export_unity_from_json(source: str | Path, output_dir: Path, *, camera_samples: int = 120) -> list[Path]:
    """Load JSON viz specs (from the CRISPR sim) and export Unity timelines."""

    specs = load_viz_specs(source)
    return export_unity_bundle(specs, output_dir=output_dir, camera_samples=camera_samples)
