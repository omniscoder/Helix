from __future__ import annotations

from collections.abc import Mapping
from typing import Any, List


def summarize_workflow_meta(meta: Mapping[str, Any]) -> List[str]:
    lines: list[str] = []
    inst_events = meta.get("inst_events")
    if isinstance(inst_events, list):
        lines.append(f"Bars (inst_events): {len(inst_events)}")
    heat_values = meta.get("heat_values")
    if isinstance(heat_values, list):
        lines.append(f"Heat samples: {len(heat_values)}")
    captured = meta.get("captured_mass")
    total = meta.get("total_mass")
    if captured is not None and total is not None:
        try:
            pct = (float(captured) / float(total)) * 100.0 if float(total) else 0.0
            lines.append(f"Captured mass: {pct:.2f}% of distribution")
        except Exception:
            pass
    extent = meta.get("x_extent")
    if extent is not None:
        try:
            lines.append(f"Timeline extent: {float(extent):.1f} units")
        except Exception:
            pass
    return lines
