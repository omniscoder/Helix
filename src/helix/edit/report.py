"""HTML report generator for edit DAG artifacts."""
from __future__ import annotations

import base64
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from .dag import EditDAG


@dataclass
class DagMetrics:
    total_nodes: int
    total_edges: int
    top_outcomes: Iterable[Dict[str, object]]
    entropy: float
    max_prob: float
    total_prob: float


def compute_metrics(dag: EditDAG, *, max_outcomes: int = 5) -> DagMetrics:
    terminals = dag.terminal_nodes()
    raw_probs = [node.log_prob for node in terminals]
    max_prob = max(raw_probs) if raw_probs else 0.0
    total_linear = sum(pow(2.718281828, p) for p in raw_probs) or 1.0
    normalized = [pow(2.718281828, p) / total_linear for p in raw_probs]
    import math

    entropy = -sum(p * math.log(p) for p in normalized if p > 0)
    ranked = sorted(zip(terminals, normalized), key=lambda item: item[1], reverse=True)
    top = [
        {
            "node_id": node.id,
            "probability": round(prob, 4),
            "stage": node.metadata.get("stage"),
        }
        for node, prob in ranked[:max_outcomes]
    ]
    return DagMetrics(
        total_nodes=len(dag.nodes),
        total_edges=len(dag.edges),
        top_outcomes=top,
        entropy=entropy,
        max_prob=max_prob,
        total_prob=total_linear,
    )


def _encode_png(path: Path) -> str:
    data = path.read_bytes()
    return "data:image/png;base64," + base64.b64encode(data).decode()


def render_html_report(dag: EditDAG, *, png_path: Optional[Path], title: str = "Helix Edit DAG Report") -> str:
    metrics = compute_metrics(dag)
    png_data = _encode_png(png_path) if png_path and png_path.exists() else None
    top_rows = "\n".join(
        f"<tr><td>{row['node_id']}</td><td>{row['stage']}</td><td>{row['probability']}</td></tr>"
        for row in metrics.top_outcomes
    )
    image_html = f"<img src='{png_data}' alt='Edit DAG' />" if png_data else "<em>No image provided.</em>"
    return f"""
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{title}</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif; background:#0b0e14; color:#f0f4ff; }}
    table {{ width:100%; border-collapse: collapse; margin-top:1rem; }}
    th, td {{ border-bottom:1px solid #1f2430; padding:0.5rem; text-align:left; }}
    section {{ margin-bottom:2rem; }}
  </style>
</head>
<body>
  <h1>{title}</h1>
  <section>
    <h2>Summary</h2>
    <ul>
      <li>Total nodes: {metrics.total_nodes}</li>
      <li>Total edges: {metrics.total_edges}</li>
      <li>Entropy: {metrics.entropy:.3f}</li>
      <li>Max branch probability (log space): {metrics.max_prob:.3f}</li>
    </ul>
  </section>
  <section>
    <h2>Visualization</h2>
    {image_html}
  </section>
  <section>
    <h2>Top Outcomes</h2>
    <table>
      <thead><tr><th>Node</th><th>Stage</th><th>Probability</th></tr></thead>
      <tbody>
        {top_rows}
      </tbody>
    </table>
  </section>
</body>
</html>
"""
