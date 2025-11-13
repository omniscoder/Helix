"""Edit DAG data structures used by Helix simulators."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

from helix.genome.digital import DigitalGenome, DigitalGenomeView

from .events import EditEvent


@dataclass
class EditNode:
    """Single node within an edit DAG (a genome view + log probability)."""

    id: str
    genome_view: DigitalGenomeView
    log_prob: float
    metadata: Dict[str, object] = field(default_factory=dict)


@dataclass
class EditEdge:
    """Directed edge describing an edit event between two nodes."""

    source: str
    target: str
    rule_name: str
    event: EditEvent
    metadata: Dict[str, object] = field(default_factory=dict)


@dataclass
class EditDAG:
    """Full edit DAG artifact (nodes + edges + root identifier)."""

    nodes: Dict[str, EditNode]
    edges: List[EditEdge]
    root_id: str

    def terminal_nodes(self) -> List[EditNode]:
        """Return terminal nodes (no outgoing edges)."""
        outgoing = {edge.source for edge in self.edges}
        return [node for node_id, node in self.nodes.items() if node_id not in outgoing]


def dag_from_payload(payload: Dict[str, object]) -> EditDAG:
    """Reconstruct an EditDAG from a serialized payload (artifact JSON)."""
    nodes: Dict[str, EditNode] = {}
    for node_id, node_entry in payload.get("nodes", {}).items():
        sequences = node_entry.get("sequences", {}) if isinstance(node_entry, dict) else {}
        genome = DigitalGenome(sequences=dict(sequences))
        view = genome.view()
        metadata = node_entry.get("metadata", {}) if isinstance(node_entry, dict) else {}
        log_prob = float(node_entry.get("log_prob", 0.0))
        nodes[node_id] = EditNode(
            id=node_id,
            genome_view=view,
            log_prob=log_prob,
            metadata=metadata,
        )

    edges: List[EditEdge] = []
    for edge_entry in payload.get("edges", []):
        event_data = edge_entry.get("event", {})
        event = EditEvent(
            chrom=event_data.get("chrom", ""),
            start=int(event_data.get("start", 0)),
            end=int(event_data.get("end", 0)),
            replacement=event_data.get("replacement", ""),
            metadata=event_data.get("metadata", {}),
        )
        edges.append(
            EditEdge(
                source=edge_entry.get("source"),
                target=edge_entry.get("target"),
                rule_name=edge_entry.get("rule", ""),
                event=event,
                metadata=edge_entry.get("metadata", {}),
            )
        )

    root_id = payload.get("root_id") or next(iter(nodes))
    return EditDAG(nodes=nodes, edges=edges, root_id=root_id)
