"""Edit DAG data structures used by Helix simulators."""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

from helix.genome.digital import DigitalGenome, DigitalGenomeView

from .events import EditEvent


@dataclass
class EditNode:
    """Single node within an edit DAG (a genome view + log probability)."""

    id: str
    genome_view: DigitalGenomeView
    log_prob: float
    metadata: Dict[str, object] = field(default_factory=dict)
    parents: Sequence[str] = field(default_factory=tuple)
    seq_hashes: Dict[str, str] = field(default_factory=dict)
    diffs: Sequence[EditEvent] = field(default_factory=tuple)


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


@dataclass
class EditDAGFrame:
    """
    Incremental snapshot of new nodes/edges introduced during DAG expansion.
    """

    step: int
    new_nodes: Dict[str, EditNode]
    new_edges: List[EditEdge]


def _compute_seq_hashes(view: DigitalGenomeView) -> Dict[str, str]:
    hashes: Dict[str, str] = {}
    for chrom, sequence in view.materialize_all().items():
        hashes[chrom] = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
    return hashes


def dag_from_payload(payload: Dict[str, object]) -> EditDAG:
    """Reconstruct an EditDAG from a serialized payload (artifact JSON)."""
    nodes: Dict[str, EditNode] = {}
    for node_id, node_entry in payload.get("nodes", {}).items():
        metadata = node_entry.get("metadata", {}) if isinstance(node_entry, dict) else {}
        log_prob = float(node_entry.get("log_prob", 0.0))
        seq_hashes = node_entry.get("seq_hashes") or {}
        parent_ids = node_entry.get("parent_ids")
        if parent_ids is None:
            parent_ids = node_entry.get("parents")
        parents = tuple(parent_ids or [])
        sequences = node_entry.get("sequences", {}) if isinstance(node_entry, dict) else {}
        genome = DigitalGenome(sequences=dict(sequences))
        view = genome.view()
        diffs_raw = node_entry.get("diffs", [])
        diff_events = tuple(
            EditEvent(
                chrom=entry.get("chrom", ""),
                start=int(entry.get("start", 0)),
                end=int(entry.get("end", 0)),
                replacement=entry.get("replacement", ""),
                metadata=entry.get("metadata", {}),
            )
            for entry in diffs_raw
        )
        if not sequences and diff_events:
            for event in diff_events:
                view = view.apply(event)
        nodes[node_id] = EditNode(
            id=node_id,
            genome_view=view,
            log_prob=log_prob,
            metadata=metadata,
            parents=parents,
            seq_hashes=dict(seq_hashes),
            diffs=diff_events,
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
