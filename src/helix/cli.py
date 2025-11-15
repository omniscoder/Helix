"""Unified Helix CLI that wraps DNA, peptide, RNA, protein, and triage helpers."""
from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import json
import math
import random
import re
import sys
import queue
import threading
import time
import uuid
from collections import OrderedDict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import TYPE_CHECKING, Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, TextIO
from types import SimpleNamespace

import numpy as np

from . import bioinformatics, cyclospectrum, triage
from .crispr import guide as crispr_guide
from .crispr import pam as crispr_pam
from .crispr import score
from .crispr import simulate as crispr_simulate
from .crispr.model import (
    CasSystem,
    CasSystemType,
    DigitalGenome as CrisprDigitalGenome,
    GuideRNA,
    PAMRule,
    TargetSite,
)
from .crispr.simulator import CutEvent, simulate_cuts
from .io import read_fasta
from .string import fm as string_fm
from .string import edit as string_edit
from .seed import minimizers as seed_minimizers
from .seed import syncmers as seed_syncmers
from .seed.extend import SeedMatch, extend_alignment
from .graphs import (
    build_dbg as graph_build_dbg,
    clean_dbg as graph_clean_dbg,
    serialize_graph as graph_serialize,
    deserialize_graph as graph_deserialize,
    export_graphml as graph_export_graphml,
    build_colored_dbg,
)
from .sketch import compute_minhash, mash_distance, compute_hll, union_hll
from .motif import discover_motifs
from .prime.model import PegRNA, PrimeEditOutcome, PrimeEditor
from .prime.simulator import locate_prime_target_site, simulate_prime_edit
from .prime.priors import resolve_prime_priors
from .crispr.pam import build_prime_pam_mask
from .genome.digital import DigitalGenome as CoreDigitalGenome
from .crispr.dag_api import build_crispr_edit_dag
from .prime.dag_api import build_prime_edit_dag
from .pcr.model import Primer, PrimerPair, PCRConfig
from .pcr.dag_api import pcr_edit_dag
from .edit.dag import EditDAG, EditDAGFrame, EditEdge, EditNode, dag_from_payload
from .edit.events import EditEvent
from .edit.report import render_html_report
from .experiment.spec import (
    CrisprExperimentSpec,
    ExperimentSpecError,
    PrimeExperimentSpec,
    load_experiment_spec,
)
from .dsl import load_hgx, build_graph_from_hgx
from .core.scheduler import Island, LiveScheduler
from .live import StateReducer, EventBus, HotSwapManager
from .live.metadata import describe_graph_nodes
from .live.realtime import RealtimeHook, RealtimeServer
from .live.cli_utils import (
    compose_run_id,
    current_command_str,
    dump_json,
    finalize_live_run,
    maybe_create_realtime_queue,
    prepare_live_bundle_dir,
    print_plan_summary,
    sanitize_run_name,
)
from .live.livelab import LiveSession, LiveDevShell
from .slices import resolve_slice, SLICE_REGISTRY

if TYPE_CHECKING:
    from .gui.modern.spec import EditVisualizationSpec

_CRISPR_EXPERIMENT_TEMPLATE = """\
kind: helix.crispr.experiment.v1
name: example_crispr_experiment
description: "Describe the goal of this CRISPR experiment."

genome:
  fasta: path/to/reference.fa
  region: chrDemo:1-500

cas:
  config: path/to/cas9.json

guide:
  sequence: ACGTACGTACGTACGTACGT
  name: GUIDE_NAME

simulation:
  max_depth: 2
  min_prob: 1.0e-4
  max_sites: 20
  seed: 0
"""

_PRIME_EXPERIMENT_TEMPLATE = """\
kind: helix.prime.experiment.v1
name: example_prime_experiment
description: "Describe the goal of this prime editing experiment."

genome:
  fasta: path/to/reference.fa
  region: chrDemo:1-500

editor:
  config: path/to/prime_editor.json

peg:
  config: path/to/peg_config.json

simulation:
  max_depth: 3
  min_prob: 1.0e-4
  seed: 0
"""
from .schema import (
    SchemaError,
    SPEC_VERSION,
    describe_schema,
    diff_manifests,
    format_manifest_diff,
    load_manifest,
    manifest,
    validate_viz_payload,
)
from .veribiota import dag_payload_to_lean, dag_payloads_to_lean, module_name_from_path
from .veribiota.checks import build_lean_check, validate_lean_check, LeanCheckError
from . import __version__ as HELIX_VERSION

try:  # optional viz imports (matplotlib)
    from .viz import rna as viz_rna
    from .viz import (
        plot_alignment_ribbon,
        plot_distance_heatmap,
        plot_minimizer_density,
        plot_motif_logo,
        plot_seed_chain,
        plot_rna_dotplot,
        render_crispr_track,
    )
    VIZ_AVAILABLE = True
except ImportError:  # pragma: no cover - viz extra not installed
    viz_rna = None  # type: ignore
    plot_alignment_ribbon = None  # type: ignore
    plot_distance_heatmap = None  # type: ignore
    plot_minimizer_density = None  # type: ignore
    plot_motif_logo = None  # type: ignore
    plot_seed_chain = None  # type: ignore
    plot_rna_dotplot = None  # type: ignore
    VIZ_AVAILABLE = False
from .rna import mfe_dotbracket, partition_posteriors, mea_structure, centroid_structure

try:
    from . import protein as protein_module

    PROTEIN_AVAILABLE = getattr(protein_module, "BIOPYTHON_AVAILABLE", True)
except ImportError:  # pragma: no cover - protein extras optional
    protein_module = None
    PROTEIN_AVAILABLE = False


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _load_sequence_arg(sequence: str | None, path: Path | None, *, default: str | None = None) -> str:
    if sequence and path:
        raise ValueError("Provide either an inline sequence or --input path, not both.")
    if path:
        return _read_text(path)
    if sequence:
        return sequence
    if default is not None:
        return default
    raise ValueError("Missing sequence data; provide a positional sequence or use --input.")


@dataclass
class TranscriptDef:
    id: str
    name: str
    strand: str
    cds_start: int
    cds_end: int
    exons: List[tuple[int, int]]
    chrom: Optional[str] = None


_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_RC_MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _normalize_exon_ranges(exons: Sequence[Sequence[int | float | str]]) -> List[tuple[int, int]]:
    normalized: List[tuple[int, int]] = []
    for exon in exons:
        if len(exon) < 2:
            continue
        try:
            start = int(exon[0])
            end = int(exon[1])
        except (TypeError, ValueError):
            continue
        if start <= 0 or end <= 0:
            continue
        if end < start:
            start, end = end, start
        normalized.append((start, end))
    return normalized


def _load_transcripts_from_json(path: Path) -> List[TranscriptDef]:
    try:
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise SystemExit(f"Transcript JSON '{path}' not found.") from exc
    except json.JSONDecodeError as exc:
        raise SystemExit(f"Transcript JSON '{path}' is not valid JSON: {exc}") from exc
    entries = raw
    if isinstance(raw, dict) and "transcripts" in raw:
        entries = raw["transcripts"]
    if not isinstance(entries, list):
        raise SystemExit(f"Transcript JSON '{path}' must contain a list of transcripts.")
    transcripts: List[TranscriptDef] = []
    for idx, entry in enumerate(entries):
        if not isinstance(entry, Mapping):
            continue
        name = str(entry.get("name") or entry.get("id") or f"tx_{idx + 1}")
        strand = (entry.get("strand") or "+").strip() or "+"
        try:
            cds_start = int(entry.get("cds_start") or entry.get("cdsStart"))
            cds_end = int(entry.get("cds_end") or entry.get("cdsEnd"))
        except (TypeError, ValueError):
            continue
        exons_raw = entry.get("exons") or []
        exons = _normalize_exon_ranges(exons_raw) or [(cds_start, cds_end)]
        if cds_start <= 0 or cds_end <= 0 or cds_end < cds_start:
            continue
        chrom = entry.get("chrom") or entry.get("chromosome")
        transcripts.append(
            TranscriptDef(
                id=str(entry.get("id") or name),
                name=name,
                strand=strand,
                cds_start=cds_start,
                cds_end=cds_end,
                exons=exons,
                chrom=str(chrom) if chrom else None,
            )
        )
    if not transcripts:
        raise SystemExit(f"No valid transcripts found in {path}.")
    return transcripts


def _select_transcript(transcripts: List[TranscriptDef], selector: Optional[str]) -> TranscriptDef:
    if not selector:
        return transcripts[0]
    selector_lower = selector.lower()
    for transcript in transcripts:
        if transcript.id.lower() == selector_lower or transcript.name.lower() == selector_lower:
            return transcript
    available = ", ".join(t.name for t in transcripts)
    raise SystemExit(f"Transcript '{selector}' not found. Available: {available}.")


def _build_coding_map_for_sequence(transcript: TranscriptDef, seq: str) -> List[int]:
    blocks = sorted(transcript.exons, key=lambda item: item[0])
    coding: List[int] = []
    if transcript.strand == "+":
        for start1, end1 in blocks:
            start = max(1, start1)
            end = min(len(seq), end1)
            for pos in range(start - 1, end):
                coding.append(pos)
    else:
        for start1, end1 in reversed(blocks):
            start = max(1, start1)
            end = min(len(seq), end1)
            for pos in range(end - 1, start - 2, -1):
                coding.append(pos)
    return coding


def _build_coding_maps(transcript: TranscriptDef, sequences: Mapping[str, str]) -> Dict[str, List[int]]:
    coding_maps: Dict[str, List[int]] = {}
    if transcript.chrom:
        seq = sequences.get(transcript.chrom)
        if seq:
            cmap = _build_coding_map_for_sequence(transcript, seq)
            coding_maps[transcript.chrom] = cmap
        return coding_maps
    for chrom, seq in sequences.items():
        coding_maps[chrom] = _build_coding_map_for_sequence(transcript, seq)
    return coding_maps


def _protein_impact_for_event(event: EditEvent, ref_seq: str, strand: str, coding_map: Optional[List[int]]) -> Optional[str]:
    if coding_map is None:
        return None
    length = max(0, int(event.end) - int(event.start))
    ins = len(event.replacement or "")
    if not (length == 1 and ins == 1):
        return None
    pos0 = int(event.start)
    try:
        t_idx = coding_map.index(pos0)
    except ValueError:
        return None
    frame = t_idx % 3
    if t_idx - frame < 0 or t_idx - frame + 2 >= len(coding_map):
        return None
    codon_positions = [coding_map[t_idx - frame + i] for i in range(3)]
    bases = [ref_seq[p] if 0 <= p < len(ref_seq) else "N" for p in codon_positions]
    if strand == "-":
        ref_codon = "".join(_RC_MAP.get(b, b) for b in reversed(bases))
        idx = frame
        alt = list(ref_codon)
        alt[idx] = _RC_MAP.get((event.replacement or "N")[0], (event.replacement or "N")[0])
        alt_codon = "".join(alt)
    else:
        ref_codon = "".join(bases)
        idx = frame
        alt = list(ref_codon)
        alt[idx] = (event.replacement or "N")[0]
        alt_codon = "".join(alt)
    aa0 = _CODON_TABLE.get(ref_codon)
    aa1 = _CODON_TABLE.get(alt_codon)
    if not aa0 or not aa1:
        return None
    if aa1 == "*":
        return "nonsense"
    if aa0 == aa1:
        return "silent"
    return "missense"


def _annotate_dag_with_transcript(dag: EditDAG, sequences: Mapping[str, str], transcript: TranscriptDef) -> bool:
    coding_maps = _build_coding_maps(transcript, sequences)
    if not coding_maps:
        return False
    updated = False
    for node in dag.nodes.values():
        if not node.diffs:
            continue
        event = node.diffs[-1]
        seq = sequences.get(event.chrom)
        cmap = coding_maps.get(event.chrom) or next(iter(coding_maps.values()), None)
        if seq is None or not cmap:
            continue
        impact = _protein_impact_for_event(event, seq, transcript.strand, cmap)
        if impact:
            node.metadata = dict(node.metadata)
            node.metadata["protein_impact"] = impact
            updated = True
    return updated


def _load_transcript_from_args(json_path: Optional[Path], transcript_name: Optional[str]) -> Optional[TranscriptDef]:
    if not json_path:
        return None
    transcripts = _load_transcripts_from_json(json_path)
    return _select_transcript(transcripts, transcript_name)


def _parse_spectrum(text: str | None) -> List[int]:
    if not text:
        return []
    tokens = text.replace(",", " ").split()
    if not tokens:
        return []
    return [int(token) for token in tokens]


def _payload_hash(payload: dict) -> str:
    normalized = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(normalized).hexdigest()


def _validate_payload_or_exit(kind: str, payload: dict) -> dict:
    try:
        return validate_viz_payload(kind, payload)
    except SchemaError as exc:
        raise SystemExit(str(exc))


def _stamp_spec_version(payload: dict, *, to_meta: bool = True) -> dict:
    if to_meta:
        meta = payload.setdefault("meta", {})
        meta.setdefault("spec_version", SPEC_VERSION)
    else:
        payload.setdefault("spec_version", SPEC_VERSION)
    return payload


def _input_meta(payload: dict) -> Dict[str, str]:
    return {"input_sha256": _payload_hash(payload)}


def _schema_sample(name: str) -> Optional[dict]:
    entry = VIZ_DEMO_PAYLOADS.get(name)
    if not entry:
        return None
    return json.loads(json.dumps(entry["data"]))


def _load_dag_viz_exporter():
    """
    Import the edit DAG visualization helper only when needed.

    The dependency chain pulls in optional packages (networkx, matplotlib), so
    we defer the import until a visualization command is invoked.
    """
    try:
        from .visualization.dag_viz import save_edit_dag_png
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise SystemExit(
            "Edit DAG visualization requires optional dependencies 'networkx' and 'matplotlib'. "
            "Install them via 'pip install helix[viz]' or 'pip install networkx matplotlib'."
        ) from exc
    return save_edit_dag_png


def _load_dag_compare_exporter():
    try:
        from .visualization.dag_viz import save_edit_dag_compare_png
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise SystemExit(
            "Edit DAG visualization requires optional dependencies 'networkx' and 'matplotlib'. "
            "Install them via 'pip install helix[viz]' or 'pip install networkx matplotlib'."
        ) from exc
    return save_edit_dag_compare_png


def _compute_dag_diff_summary(dag_a: EditDAG, dag_b: EditDAG, *, label_a: str, label_b: str) -> Dict[str, Any]:
    nodes_a = set(dag_a.nodes.keys())
    nodes_b = set(dag_b.nodes.keys())
    unique_a = sorted(nodes_a - nodes_b)
    unique_b = sorted(nodes_b - nodes_a)
    shared = sorted(nodes_a & nodes_b)

    terminals_a = {node.id for node in dag_a.terminal_nodes()}
    terminals_b = {node.id for node in dag_b.terminal_nodes()}

    return {
        "label_a": label_a,
        "label_b": label_b,
        "nodes": {
            "unique_to_a": unique_a,
            "unique_to_b": unique_b,
            "shared": shared,
        },
        "terminal_nodes": {
            "unique_to_a": sorted(terminals_a - terminals_b),
            "unique_to_b": sorted(terminals_b - terminals_a),
            "shared": sorted(terminals_a & terminals_b),
        },
        "top_outcomes": {
            label_a: _top_outcomes_from_dag(dag_a, topk=5),
            label_b: _top_outcomes_from_dag(dag_b, topk=5),
        },
    }


def _load_dag_animator():
    """
    Import the GIF animation helper lazily because it depends on viz extras.
    """
    try:
        from .visualization.animate import animate_edit_dag
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise SystemExit(
            "Edit DAG animation requires optional dependencies 'networkx', 'matplotlib', and 'imageio'. "
            "Install them via 'pip install helix[viz]' or 'pip install networkx matplotlib imageio'."
        ) from exc
    return animate_edit_dag


def _print_schema_help(kind: str, sample_name: str | None = None) -> None:
    print(describe_schema(kind))
    if sample_name:
        sample = _schema_sample(sample_name)
        if sample:
            print("\nSample payload:")
            print(json.dumps(sample, indent=2))


VIZ_DEMO_PAYLOADS: Dict[str, Dict[str, object]] = {
    "minimizers": {
        "kind": "viz_minimizers",
        "data": {"sequence_length": 500, "minimizers": [5, [25, "AAA", 1], {"pos": 120}]},
    },
    "seed-chain": {
        "kind": "viz_seed_chain",
        "data": {
            "ref_length": 400,
            "qry_length": 380,
            "chains": [
                [{"ref_start": 10, "ref_end": 40, "qry_start": 12, "qry_end": 42}],
                [{"ref_start": 80, "ref_end": 120, "qry_start": 78, "qry_end": 118}],
            ],
        },
    },
    "rna-dotplot": {
        "kind": "viz_rna_dotplot",
        "data": {"posterior": [[0.0, 0.6, 0.0], [0.6, 0.0, 0.4], [0.0, 0.4, 0.0]]},
    },
    "alignment-ribbon": {
        "kind": "viz_alignment_ribbon",
        "data": {
            "ref_length": 300,
            "qry_length": 290,
            "ref_start": 20,
            "qry_start": 18,
            "cigar": "30M2I10M3D15M",
            "metadata": {"name": "demo_read"},
        },
    },
    "distance-heatmap": {
        "kind": "viz_distance_heatmap",
        "data": {"labels": ["A", "B", "C"], "matrix": [[0.0, 0.05, 0.1], [0.05, 0.0, 0.08], [0.1, 0.08, 0.0]]},
    },
    "motif-logo": {
        "kind": "viz_motif_logo",
        "data": {
            "alphabet": ["A", "C", "G", "T"],
            "pwm": [
                [0.25, 0.25, 0.25, 0.25],
                [0.05, 0.05, 0.85, 0.05],
                [0.6, 0.1, 0.1, 0.2],
            ],
            "background": [0.25, 0.25, 0.25, 0.25],
        },
    },
}


def _default_viz_spec_path(save: Path | None, provided: Path | None) -> Path | None:
    if provided:
        return provided
    if save:
        return save.with_name(f"{save.stem}.viz.json")
    return None


def _file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()

def _write_provenance(
    image_path: Optional[Path],
    *,
    schema_kind: str,
    spec: Dict[str, Any],
    input_sha: Optional[str],
    command: str,
    viz_spec_path: Optional[Path],
) -> None:
    if not image_path:
        return
    img_path = Path(image_path)
    if not img_path.exists():
        return
    image_sha = _file_sha256(img_path)
    if viz_spec_path and Path(viz_spec_path).exists():
        viz_spec_sha = _file_sha256(Path(viz_spec_path))
    else:
        viz_spec_sha = _payload_hash(spec)
    spec_version = spec.get("spec_version") or spec.get("meta", {}).get("spec_version")
    provenance = {
        "schema_kind": schema_kind,
        "spec_version": spec_version,
        "input_sha256": input_sha,
        "viz_spec_sha256": viz_spec_sha,
        "image_sha256": image_sha,
        "helix_version": HELIX_VERSION,
        "command": command,
    }
    prov_path = img_path.with_name(img_path.stem + ".provenance.json")
    prov_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")


def _make_realtime_command_handler(scheduler: LiveScheduler, *, default_field_node: str = "egf_field"):
    graph = scheduler.graph

    def _handler(command: Dict[str, Any]) -> None:
        if command.get("kind") == "live_control":
            _dispatch_control_doc(command)
            return
        action = command.get("command")
        if action == "pause":
            scheduler.pause()
        elif action == "resume":
            scheduler.resume()
        elif action == "toggle":
            scheduler.toggle_pause()
        elif action == "set_input":
            node = command.get("node", default_field_node)
            values = command.get("values") or {}
            if node and isinstance(values, Mapping):
                scheduler.update_external(str(node), {k: float(v) for k, v in values.items() if isinstance(v, (int, float))})
        elif action == "switch_variant":
            variant = command.get("variant")
            if isinstance(variant, str):
                scheduler.apply_hot_swap("variant", {"variant": variant})
        elif action == "set_egf":
            value = command.get("value")
            if isinstance(value, (int, float)):
                scheduler.update_external(default_field_node, {"control": float(value)})

    def _dispatch_control_doc(control: Mapping[str, Any]) -> None:
        control_type = control.get("type")
        if control_type == "pause":
            scheduler.pause()
        elif control_type == "resume":
            scheduler.resume()
        elif control_type == "toggle":
            scheduler.toggle_pause()
        elif control_type == "set_variant":
            variant = control.get("variant")
            if isinstance(variant, str):
                scheduler.apply_hot_swap("variant", {"variant": variant})
        elif control_type == "set_param":
            target = control.get("target")
            param = control.get("param")
            value = control.get("value")
            if isinstance(target, str) and isinstance(param, str) and isinstance(value, (int, float)):
                try:
                    node = graph.node(target)
                except KeyError:
                    return
                params = getattr(node, "params", None)
                if isinstance(params, dict):
                    params[param] = float(value)
        elif control_type == "set_field":
            target_node = control.get("target") or default_field_node
            value = control.get("value")
            if isinstance(target_node, str) and isinstance(value, (int, float)):
                scheduler.update_external(target_node, {"control": float(value)})

    return _handler


def _parse_dt_overrides(entries: Optional[Sequence[str]]) -> Dict[str, float]:
    overrides: Dict[str, float] = {}
    for raw in entries or []:
        if "=" not in raw:
            raise SystemExit(f"Invalid dt override '{raw}'. Expected format node=0.05")
        name, value = raw.split("=", 1)
        name = name.strip()
        if not name:
            raise SystemExit(f"Invalid dt override '{raw}'. Node name missing.")
        try:
            overrides[name] = float(value)
        except ValueError as exc:
            raise SystemExit(f"Invalid dt override '{raw}'. Expected a numeric dt.") from exc
    return overrides


def _parse_node_inputs(entries: Optional[Sequence[str]]) -> Dict[str, Dict[str, float]]:
    inputs: Dict[str, Dict[str, float]] = {}
    for raw in entries or []:
        if "=" not in raw or "." not in raw.split("=", 1)[0]:
            raise SystemExit(f"Invalid input '{raw}'. Expected node.port=value")
        key, value = raw.split("=", 1)
        node, port = key.split(".", 1)
        node = node.strip()
        port = port.strip()
        if not node or not port:
            raise SystemExit(f"Invalid input '{raw}'. Node or port missing.")
        try:
            numeric = float(value)
        except ValueError as exc:
            raise SystemExit(f"Invalid input '{raw}'. Expected numeric value.") from exc
        inputs.setdefault(node, {})[port] = numeric
    return inputs


class _TimeSeriesSignal:
    def __init__(self, points: Sequence[Mapping[str, Any]]):
        if not points:
            raise SystemExit("Time-series inputs must contain at least one point.")
        normalized = []
        for idx, point in enumerate(points):
            if "t" not in point or "value" not in point:
                raise SystemExit(f"Time-series point #{idx} missing 't' or 'value'.")
            try:
                t_val = float(point["t"])
                v_val = float(point["value"])
            except (TypeError, ValueError) as exc:
                raise SystemExit(f"Invalid time-series point #{idx}: expected numeric t/value.") from exc
            normalized.append((t_val, v_val))
        normalized.sort(key=lambda item: item[0])
        self.points = normalized

    def __call__(self, t: float) -> float:
        value = self.points[0][1]
        for time_point, candidate in self.points:
            if t < time_point:
                break
            value = candidate
        return value


def _load_input_series(path: Optional[Path]) -> Dict[str, Dict[str, Callable[[float], float]]]:
    if not path:
        return {}
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    signals: Dict[str, Dict[str, Callable[[float], float]]] = {}
    for key, value in raw.items():
        if "." not in key:
            raise SystemExit(f"Invalid input key '{key}' (expected node.port).")
        node, port = key.split(".", 1)
        node = node.strip()
        port = port.strip()
        if not node or not port:
            raise SystemExit(f"Invalid input key '{key}' (missing node or port).")
        if isinstance(value, (int, float)):
            fn = (lambda constant: (lambda _t: float(constant)))(value)
        elif isinstance(value, list):
            fn = _TimeSeriesSignal(value)
        else:
            raise SystemExit(
                f"Unsupported input spec for '{key}'. Use a number or a list of {{'t','value'}} points."
            )
        signals.setdefault(node, {})[port] = fn
    return signals


def _build_input_provider(
    constants: Mapping[str, Mapping[str, float]],
    series: Mapping[str, Mapping[str, Callable[[float], float]]],
) -> Callable[[str, float], Mapping[str, float]]:
    registry: Dict[str, Dict[str, Callable[[float], float]]] = {}

    def register(node: str, port: str, fn: Callable[[float], float]) -> None:
        registry.setdefault(node, {})[port] = fn

    # file-provided series first
    for node, ports in series.items():
        for port, fn in ports.items():
            register(node, port, fn)

    # CLI constants override file-defined entries
    for node, ports in constants.items():
        for port, value in ports.items():
            register(node, port, lambda _t, v=value: v)

    def provider(node_name: str, t: float) -> Mapping[str, float]:
        node_ports = registry.get(node_name, {})
        return {port: fn(t) for port, fn in node_ports.items()}

    return provider


def _read_json_optional(path: Path) -> Any:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _ensure_viz_available(feature: str = "visualization") -> None:
    if not VIZ_AVAILABLE:
        raise SystemExit(
            f"{feature} requires the 'viz' extra. Install with `pip install \"veri-helix[viz]\"`."
        )


def _ensure_realtime_available() -> None:
    try:
        import glfw  # type: ignore  # noqa: F401
        import moderngl  # type: ignore  # noqa: F401
    except ImportError as exc:  # pragma: no cover - optional deps
        raise SystemExit(
            "Realtime visualization requires the 'realtime' extra (glfw + moderngl). "
            "Install with `pip install \"veri-helix[realtime]\"`."
        ) from exc


def _ensure_modern_gui_available() -> None:
    try:
        from PySide6 import QtWidgets  # type: ignore  # noqa: F401
        import moderngl  # type: ignore  # noqa: F401
    except ImportError as exc:  # pragma: no cover - optional deps
        raise SystemExit(
            "Modern CRISPR viz requires PySide6 + moderngl. Install with "
            "`pip install \"veri-helix[gui,realtime]\"`."
        ) from exc


CAS_PRESET_CONFIGS: Dict[str, Dict[str, Any]] = {
    "cas9": {
        "name": "SpCas9",
        "system_type": "cas9",
        "pam_rules": [{"pattern": "NGG", "description": "SpCas9 canonical PAM"}],
        "cut_offset": 3,
        "max_mismatches": 3,
        "weight_mismatch_penalty": 1.0,
        "weight_pam_penalty": 2.0,
    },
    "cas12a": {
        "name": "LbCas12a",
        "system_type": "cas12a",
        "pam_rules": [{"pattern": "TTTV", "description": "LbCas12a canonical PAM"}],
        "cut_offset": -18,
        "max_mismatches": 4,
        "weight_mismatch_penalty": 1.0,
        "weight_pam_penalty": 2.5,
    },
}


def _load_json_config(path: Path) -> Dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _cas_system_from_config(config: Mapping[str, Any]) -> CasSystem:
    try:
        name = str(config["name"])
        cut_offset = int(config["cut_offset"])
    except KeyError as exc:  # pragma: no cover - defensive, CLI validated via tests
        raise ValueError(f"CasSystem config missing required field: {exc}") from exc
    system_type_value = config.get("system_type") or config.get("type")
    if not system_type_value:
        raise ValueError("CasSystem config requires a 'system_type' (or legacy 'type') field.")
    system_type = CasSystemType(str(system_type_value).lower())
    pam_entries = config.get("pam_rules")
    if not pam_entries:
        pam_pattern = config.get("pam_pattern")
        if pam_pattern:
            pam_entries = [{"pattern": pam_pattern, "description": config.get("pam_description", "")}]
        else:
            pam_entries = []
    if not pam_entries:
        raise ValueError("CasSystem config requires at least one PAM rule.")
    pam_rules: List[PAMRule] = []
    for entry in pam_entries:
        if isinstance(entry, str):
            pam_rules.append(PAMRule(pattern=entry))
            continue
        pattern = entry.get("pattern")
        if not pattern:
            raise ValueError("PAM rule entries require a 'pattern' value.")
        pam_rules.append(PAMRule(pattern=str(pattern), description=str(entry.get("description", ""))))
    return CasSystem(
        name=name,
        system_type=system_type,
        pam_rules=pam_rules,
        cut_offset=cut_offset,
        max_mismatches=int(config.get("max_mismatches", 3)),
        weight_mismatch_penalty=float(config.get("weight_mismatch_penalty", 1.0)),
        weight_pam_penalty=float(config.get("weight_pam_penalty", 2.0)),
    )


def _resolve_cas_system(
    preset_name: Optional[str],
    config_path: Optional[Path] = None,
    *,
    inline_config: Optional[Mapping[str, Any]] = None,
) -> CasSystem:
    if inline_config:
        return _cas_system_from_config(inline_config)
    if config_path:
        return _cas_system_from_config(_load_json_config(config_path))
    if preset_name:
        preset_cfg = CAS_PRESET_CONFIGS.get(preset_name.lower())
        if not preset_cfg:
            raise SystemExit(f"Unknown CasSystem preset '{preset_name}'.")
        return _cas_system_from_config(preset_cfg)
    raise SystemExit("Provide --cas or --cas-config to select a CasSystem.")


def _load_digital_genome(path: Path) -> CrisprDigitalGenome:
    try:
        core = CoreDigitalGenome.from_fasta(path)
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc
    return CrisprDigitalGenome(sequences=dict(core.sequences))


def _digital_genome_summary(genome: CrisprDigitalGenome, *, source: Optional[Path] = None) -> Dict[str, Any]:
    chromosomes = [
        {"name": chrom, "length": len(sequence)}
        for chrom, sequence in genome.sequences.items()
    ]
    total_length = sum(entry["length"] for entry in chromosomes)
    summary: Dict[str, Any] = {
        "chromosomes": chromosomes,
        "total_length": total_length,
        "count": len(chromosomes),
    }
    if source is not None:
        summary["source"] = str(source)
    return summary


def _core_genome_from_legacy(genome: CrisprDigitalGenome) -> CoreDigitalGenome:
    return CoreDigitalGenome(sequences=dict(genome.sequences))


def _genome_from_experiment(fasta: Path, region: Optional[str]) -> CrisprDigitalGenome:
    core = CoreDigitalGenome.from_fasta(fasta)
    if region:
        core = core.slice_region(region)
    return CrisprDigitalGenome(sequences=dict(core.sequences))


def _run_experiment_to_dag(
    config_path: Path,
    *,
    use_gpu_override: Optional[bool] = None,
    frame_consumer: Optional[Callable[[EditDAGFrame], None]] = None,
) -> tuple[EditDAG, Dict[str, Any], str]:
    try:
        spec = load_experiment_spec(config_path)
    except ExperimentSpecError as exc:
        raise SystemExit(str(exc)) from exc

    if isinstance(spec, CrisprExperimentSpec):
        genome = _genome_from_experiment(spec.genome.fasta, spec.genome.region)
        cas = _cas_system_from_config(_load_json_config(spec.cas_config))
        guide = GuideRNA(
            sequence=bioinformatics.normalize_sequence(spec.guide_sequence),
            pam=spec.guide_pam,
            name=spec.guide_name,
        )
        use_gpu = use_gpu_override if use_gpu_override is not None else spec.simulation.use_gpu
        dag = build_crispr_edit_dag(
            genome,
            cas,
            guide,
            rng_seed=spec.simulation.seed,
            max_depth=spec.simulation.max_depth,
            min_prob=spec.simulation.min_prob,
            max_sites=spec.simulation.max_sites,
            use_gpu=use_gpu,
            frame_consumer=frame_consumer,
        )
        metadata = {
            "experiment_name": spec.name,
            "experiment_description": spec.description,
            "experiment_kind": spec.kind,
            "experiment_config": str(config_path),
            "genome_fasta": str(spec.genome.fasta),
            "genome_region": spec.genome.region,
            "cas_config": str(spec.cas_config),
            "guide_name": spec.guide_name,
            "physics_backend": "gpu" if use_gpu else "cpu",
        }
        artifact = "helix.crispr.edit_dag.v1.1"
        return dag, metadata, artifact

    if isinstance(spec, PrimeExperimentSpec):
        genome = _genome_from_experiment(spec.genome.fasta, spec.genome.region)
        editor = _load_prime_editor_from_config(spec.editor_config)
        if spec.peg_config:
            peg = _load_peg_from_json(spec.peg_config)
        elif spec.peg_inline:
            peg_data = spec.peg_inline
            peg = PegRNA(
                spacer=bioinformatics.normalize_sequence(peg_data["spacer"]),
                pbs=bioinformatics.normalize_sequence(peg_data["pbs"]),
                rtt=bioinformatics.normalize_sequence(peg_data["rtt"]),
                name=peg_data.get("name"),
            )
        else:  # pragma: no cover - defensive
            raise SystemExit("Prime experiment missing peg definition.")
        dag = build_prime_edit_dag(
            genome,
            editor,
            peg,
            rng_seed=spec.simulation.seed,
            max_depth=spec.simulation.max_depth,
            min_prob=spec.simulation.min_prob,
            frame_consumer=frame_consumer,
        )
        metadata = {
            "experiment_name": spec.name,
            "experiment_description": spec.description,
            "experiment_kind": spec.kind,
            "experiment_config": str(config_path),
            "genome_fasta": str(spec.genome.fasta),
            "genome_region": spec.genome.region,
            "editor_config": str(spec.editor_config),
            "peg_name": peg.name,
        }
        artifact = "helix.prime.edit_dag.v1.1"
        if use_gpu_override:
            print("Warning: --use-gpu is not currently available for prime experiments; running on CPU.")
        return dag, metadata, artifact

    raise SystemExit("Unsupported experiment specification.")


def _serialize_pam_rule(rule: PAMRule) -> Dict[str, Any]:
    return {"pattern": rule.pattern, "description": rule.description}


def _serialize_cas_system(cas: CasSystem) -> Dict[str, Any]:
    return {
        "name": cas.name,
        "system_type": cas.system_type.value,
        "cut_offset": cas.cut_offset,
        "max_mismatches": cas.max_mismatches,
        "weight_mismatch_penalty": cas.weight_mismatch_penalty,
        "weight_pam_penalty": cas.weight_pam_penalty,
        "pam_rules": [_serialize_pam_rule(rule) for rule in cas.pam_rules],
    }


def _serialize_guide(guide: GuideRNA) -> Dict[str, Any]:
    payload = {
        "sequence": guide.sequence,
        "pam": guide.pam,
        "name": guide.name,
    }
    if guide.metadata:
        payload["metadata"] = guide.metadata
    return payload


def _serialize_target_site(site: TargetSite) -> Dict[str, Any]:
    data: Dict[str, Any] = {
        "chrom": site.chrom,
        "start": site.start,
        "end": site.end,
        "strand": site.strand,
        "sequence": site.sequence,
    }
    if site.on_target_score is not None:
        data["on_target_score"] = site.on_target_score
    if site.off_target_score is not None:
        data["off_target_score"] = site.off_target_score
    if site.pam_match_score is not None:
        data["pam_match_score"] = site.pam_match_score
    return data


def _serialize_cut_event(event: CutEvent) -> Dict[str, Any]:
    return {
        "site": _serialize_target_site(event.site),
        "cut_position": event.cut_position,
        "guide": _serialize_guide(event.guide),
        "cas": _serialize_cas_system(event.cas),
        "score": event.score,
    }


def _serialize_prime_editor(editor: PrimeEditor) -> Dict[str, Any]:
    return {
        "name": editor.name,
        "cas": _serialize_cas_system(editor.cas),
        "nick_to_edit_offset": editor.nick_to_edit_offset,
        "efficiency_scale": editor.efficiency_scale,
        "indel_bias": editor.indel_bias,
        "mismatch_tolerance": editor.mismatch_tolerance,
        "flap_balance": editor.flap_balance,
        "reanneal_bias": editor.reanneal_bias,
    }


def _serialize_peg(peg: PegRNA) -> Dict[str, Any]:
    data: Dict[str, Any] = {
        "spacer": peg.spacer,
        "pbs": peg.pbs,
        "rtt": peg.rtt,
    }
    if peg.name:
        data["name"] = peg.name
    if peg.metadata:
        data["metadata"] = peg.metadata
    return data


def _serialize_prime_outcome(outcome: PrimeEditOutcome) -> Dict[str, Any]:
    return {
        "site": _serialize_target_site(outcome.site),
        "edited_sequence": outcome.edited_sequence,
        "logit_score": outcome.logit_score,
        "description": outcome.description,
        "stage": outcome.stage,
        "metadata": outcome.metadata,
    }


def _write_json_output(payload: Dict[str, Any], output_path: Optional[Path]) -> None:
    text = json.dumps(payload, indent=2)
    if output_path:
        out_path = Path(output_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)


def _write_viz_spec_output(spec: "EditVisualizationSpec", output_path: Path) -> None:
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(spec.to_payload(), indent=2)
    out_path.write_text(text + "\n", encoding="utf-8")


def _build_guide(
    sequence: str,
    *,
    name: Optional[str],
    pam: Optional[str],
    metadata: Optional[Mapping[str, str]] = None,
) -> GuideRNA:
    normalized = bioinformatics.normalize_sequence(sequence)
    guide = GuideRNA(sequence=normalized, pam=pam, name=name)
    if metadata:
        guide.metadata.update({key: str(value) for key, value in metadata.items()})
    return guide


def _normalize_region_hint(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    trimmed = value.strip()
    if not trimmed or trimmed == ".":
        return None
    return trimmed


def _load_guides_from_tsv(path: Path) -> List[GuideRNA]:
    guides: List[GuideRNA] = []
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "sequence" not in reader.fieldnames:
            raise SystemExit("Guides TSV must contain at least a 'sequence' column.")
        for idx, row in enumerate(reader, start=1):
            sequence = (row.get("sequence") or "").strip()
            if not sequence:
                raise SystemExit(f"Guide row {idx} missing 'sequence'.")
            name = (row.get("name") or "").strip() or None
            region_hint = _normalize_region_hint(row.get("region"))
            metadata = {"region": region_hint} if region_hint else None
            guides.append(_build_guide(sequence, name=name, pam=None, metadata=metadata))
    if not guides:
        raise SystemExit(f"No guides found in {path}.")
    return guides


def _load_guides_from_json(path: Path) -> List[GuideRNA]:
    data = _load_json_config(path)
    entries: List[Mapping[str, Any]]
    if isinstance(data, list):
        entries = data  # type: ignore[assignment]
    else:
        maybe_guides = data.get("guides")
        if not isinstance(maybe_guides, list):
            raise SystemExit("Guide JSON must be a list or contain a top-level 'guides' list.")
        entries = maybe_guides  # type: ignore[assignment]
    guides: List[GuideRNA] = []
    for idx, entry in enumerate(entries, start=1):
        sequence = str(entry.get("sequence", "")).strip()
        if not sequence:
            raise SystemExit(f"Guide entry {idx} missing 'sequence'.")
        name_value = entry.get("name")
        name = str(name_value).strip() if name_value else None
        region_hint = _normalize_region_hint(entry.get("region"))
        metadata = {"region": region_hint} if region_hint else None
        guides.append(_build_guide(sequence, name=name, pam=None, metadata=metadata))
    if not guides:
        raise SystemExit(f"No guides found in {path}.")
    return guides


def _load_guides_file(path: Path) -> List[GuideRNA]:
    suffix = path.suffix.lower()
    if suffix in (".tsv", ".txt"):
        return _load_guides_from_tsv(path)
    if suffix == ".json":
        return _load_guides_from_json(path)
    raise SystemExit("Guide files must be .tsv or .json.")


def _safe_guide_label(name: str) -> str:
    label = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    return label or "guide"


def _load_peg_from_args(args: argparse.Namespace) -> PegRNA:
    config: Dict[str, Any] = {}
    if getattr(args, "peg_config", None):
        config = _load_json_config(args.peg_config)

    def _value(key: str, attr: str) -> Optional[str]:
        arg_value = getattr(args, attr, None)
        if arg_value:
            return arg_value
        return config.get(key)

    spacer = _value("spacer", "peg_spacer")
    pbs = _value("pbs", "peg_pbs")
    rtt = _value("rtt", "peg_rtt")
    if not spacer or not pbs or not rtt:
        raise SystemExit("Provide peg spacer/PBS/RTT via --peg-config or CLI flags.")
    metadata = config.get("metadata") or {}
    if metadata and not isinstance(metadata, dict):
        raise SystemExit("peg metadata must be a JSON object.")
    name = getattr(args, "peg_name", None) or config.get("name")
    return PegRNA(
        spacer=bioinformatics.normalize_sequence(spacer),
        pbs=bioinformatics.normalize_sequence(pbs),
        rtt=bioinformatics.normalize_sequence(rtt),
        name=name,
        metadata=dict(metadata),
    )


def _load_peg_from_json(path: Path) -> PegRNA:
    data = _load_json_config(path)
    spacer = str(data.get("spacer", "")).strip()
    pbs = str(data.get("pbs", "")).strip()
    rtt = str(data.get("rtt", "")).strip()
    if not spacer or not pbs or not rtt:
        raise SystemExit(f"PegRNA config at {path} must include spacer, pbs, and rtt.")
    metadata = data.get("metadata") or {}
    if metadata and not isinstance(metadata, dict):
        raise SystemExit("peg metadata must be a JSON object.")
    name = data.get("name")
    return PegRNA(
        spacer=bioinformatics.normalize_sequence(spacer),
        pbs=bioinformatics.normalize_sequence(pbs),
        rtt=bioinformatics.normalize_sequence(rtt),
        name=name,
        metadata=dict(metadata),
    )


def _load_pegs_from_tsv(path: Path) -> List[PegRNA]:
    pegs: List[PegRNA] = []
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"spacer", "pbs", "rtt"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            missing = required - set(reader.fieldnames or [])
            raise SystemExit(f"Peg TSV must include columns: {', '.join(sorted(missing))}.")
        for idx, row in enumerate(reader, start=1):
            spacer = (row.get("spacer") or "").strip()
            pbs = (row.get("pbs") or "").strip()
            rtt = (row.get("rtt") or "").strip()
            if not spacer or not pbs or not rtt:
                raise SystemExit(f"Peg row {idx} missing spacer/pbs/rtt.")
            name = (row.get("name") or "").strip() or None
            meta: Dict[str, str] = {}
            region_hint = _normalize_region_hint(row.get("region"))
            if region_hint:
                meta["region"] = region_hint
            pegs.append(
                PegRNA(
                    spacer=bioinformatics.normalize_sequence(spacer),
                    pbs=bioinformatics.normalize_sequence(pbs),
                    rtt=bioinformatics.normalize_sequence(rtt),
                    name=name,
                    metadata=meta,
                )
            )
    if not pegs:
        raise SystemExit(f"No pegRNAs found in {path}.")
    return pegs


def _load_pegs_from_json(path: Path) -> List[PegRNA]:
    data = _load_json_config(path)
    entries: List[Mapping[str, Any]]
    if isinstance(data, list):
        entries = data  # type: ignore[assignment]
    else:
        maybe_pegs = data.get("pegs")
        if not isinstance(maybe_pegs, list):
            raise SystemExit("Peg JSON must be a list or contain a top-level 'pegs' list.")
        entries = maybe_pegs  # type: ignore[assignment]
    pegs: List[PegRNA] = []
    for idx, entry in enumerate(entries, start=1):
        spacer = str(entry.get("spacer", "")).strip()
        pbs = str(entry.get("pbs", "")).strip()
        rtt = str(entry.get("rtt", "")).strip()
        if not spacer or not pbs or not rtt:
            raise SystemExit(f"Peg entry {idx} missing spacer/pbs/rtt.")
        metadata = entry.get("metadata") or {}
        if metadata and not isinstance(metadata, dict):
            raise SystemExit("Peg metadata must be a JSON object.")
        region_hint = _normalize_region_hint(entry.get("region"))
        if region_hint:
            metadata = dict(metadata)
            metadata["region"] = region_hint
        pegs.append(
            PegRNA(
                spacer=bioinformatics.normalize_sequence(spacer),
                pbs=bioinformatics.normalize_sequence(pbs),
                rtt=bioinformatics.normalize_sequence(rtt),
                name=entry.get("name"),
                metadata=dict(metadata),
            )
        )
    if not pegs:
        raise SystemExit(f"No pegRNAs found in {path}.")
    return pegs


def _load_pegs_file(path: Path) -> List[PegRNA]:
    suffix = path.suffix.lower()
    if suffix in (".tsv", ".txt"):
        return _load_pegs_from_tsv(path)
    if suffix == ".json":
        return _load_pegs_from_json(path)
    raise SystemExit("Peg files must be .tsv or .json.")


def _load_prime_editor_from_args(args: argparse.Namespace) -> PrimeEditor:
    config: Dict[str, Any] = {}
    if getattr(args, "editor_config", None):
        config = _load_json_config(args.editor_config)
    inline_cas_config = config.get("cas")
    cas = _resolve_cas_system(
        preset_name=config.get("cas_preset") or getattr(args, "cas", None),
        config_path=getattr(args, "cas_config", None),
        inline_config=inline_cas_config,
    )
    name = getattr(args, "editor_name", None) or config.get("name") or cas.name

    def _override(attr: str, key: str, default: Any) -> Any:
        arg_value = getattr(args, attr, None)
        if arg_value is not None:
            return arg_value
        return config.get(key, default)

    nick_offset = int(_override("nick_offset", "nick_to_edit_offset", 0))
    efficiency_scale = float(_override("efficiency_scale", "efficiency_scale", 1.0))
    indel_bias = float(_override("indel_bias", "indel_bias", 0.0))
    mismatch_tolerance = int(_override("mismatch_tolerance", "mismatch_tolerance", 3))
    flap_balance = float(_override("flap_balance", "flap_balance", 0.5))
    reanneal_bias = float(_override("reanneal_bias", "reanneal_bias", 0.1))
    metadata = config.get("metadata") or {}
    if metadata and not isinstance(metadata, dict):
        raise SystemExit("Prime editor metadata must be a JSON object.")
    return PrimeEditor(
        name=name,
        cas=cas,
        nick_to_edit_offset=nick_offset,
        efficiency_scale=efficiency_scale,
        indel_bias=indel_bias,
        mismatch_tolerance=mismatch_tolerance,
        flap_balance=flap_balance,
        reanneal_bias=reanneal_bias,
        metadata=dict(metadata),
    )


def _load_prime_editor_from_config(path: Path) -> PrimeEditor:
    args = SimpleNamespace(
        editor_config=path,
        editor_name=None,
        cas=None,
        cas_config=None,
        nick_offset=None,
        efficiency_scale=None,
        indel_bias=None,
        mismatch_tolerance=None,
        flap_balance=None,
        reanneal_bias=None,
    )
    return _load_prime_editor_from_args(args)


def _load_primer_pair_from_args(args: argparse.Namespace) -> PrimerPair:
    if not getattr(args, "primer_config", None):
        raise SystemExit("Provide --primer-config pointing to a primer pair JSON file.")
    config = _load_json_config(args.primer_config)
    forward_cfg = config.get("forward")
    reverse_cfg = config.get("reverse")
    if not isinstance(forward_cfg, Mapping) or not isinstance(reverse_cfg, Mapping):
        raise SystemExit("Primer config must contain 'forward' and 'reverse' objects.")

    def _primer(entry: Mapping[str, Any], fallback_name: str) -> Primer:
        sequence = entry.get("sequence")
        if not sequence:
            raise SystemExit(f"Primer '{fallback_name}' requires a 'sequence' field.")
        metadata = entry.get("metadata") or {}
        if metadata and not isinstance(metadata, dict):
            raise SystemExit(f"Metadata for primer '{fallback_name}' must be a JSON object.")
        return Primer(
            name=str(entry.get("name") or fallback_name),
            sequence=bioinformatics.normalize_sequence(sequence),
            max_mismatches=int(entry.get("max_mismatches", 2)),
            metadata=dict(metadata),
        )

    forward = _primer(forward_cfg, "forward")
    reverse = _primer(reverse_cfg, "reverse")
    metadata = config.get("metadata") or {}
    if metadata and not isinstance(metadata, dict):
        raise SystemExit("Primer pair metadata must be a JSON object.")
    pair_name = str(config.get("name") or f"{forward.name}/{reverse.name}")
    return PrimerPair(name=pair_name, forward=forward, reverse=reverse, metadata=dict(metadata))


def _load_pcr_config_from_args(args: argparse.Namespace) -> PCRConfig:
    data: Dict[str, Any] = {}
    if getattr(args, "pcr_config", None):
        data = _load_json_config(args.pcr_config)

    def _value(attr: str, key: str, default: Any) -> Any:
        arg_value = getattr(args, attr, None)
        if arg_value is not None:
            return arg_value
        return data.get(key, default)

    return PCRConfig(
        cycles=int(_value("cycles", "cycles", 10)),
        per_cycle_efficiency=float(_value("per_cycle_efficiency", "per_cycle_efficiency", 0.9)),
        error_rate=float(_value("error_rate", "error_rate", 1e-3)),
        max_amplicon_length=int(_value("max_amplicon_length", "max_amplicon_length", 2000)),
        min_amplicon_length=int(_value("min_amplicon_length", "min_amplicon_length", 50)),
        max_amplicons=int(_value("max_amplicons", "max_amplicons", 256)),
    )


def _edit_dag_to_payload(dag: EditDAG, *, artifact: str, metadata: Dict[str, Any]) -> Dict[str, Any]:
    metadata = dict(metadata)
    metadata.setdefault("schema_version", "1.1")
    metadata.setdefault("rule_version", "1.0")
    metadata.setdefault("created_at", datetime.now(timezone.utc).isoformat())
    nodes_payload: Dict[str, Any] = {}
    for node_id, node in dag.nodes.items():
        nodes_payload[node_id] = _edit_node_to_payload(node)
    edges_payload = []
    for edge in dag.edges:
        edges_payload.append(_edit_edge_to_payload(edge))
    return {
        "artifact": artifact,
        "version": "1.1",
        "schema_version": "1.1",
        "meta": metadata,
        "nodes": nodes_payload,
        "edges": edges_payload,
        "root_id": dag.root_id,
    }


def _edit_node_to_payload(node: EditNode) -> Dict[str, Any]:
    entry: Dict[str, Any] = {
        "log_prob": node.log_prob,
        "metadata": node.metadata,
        "parent_ids": list(node.parents),
        "seq_hashes": node.seq_hashes,
    }
    if node.diffs:
        entry["diffs"] = [
            {
                "chrom": event.chrom,
                "start": event.start,
                "end": event.end,
                "replacement": event.replacement,
                "metadata": event.metadata,
            }
            for event in node.diffs
        ]
    entry["sequences"] = node.genome_view.materialize_all()
    return entry


def _edit_edge_to_payload(edge: EditEdge) -> Dict[str, Any]:
    return {
        "source": edge.source,
        "target": edge.target,
        "rule": edge.rule_name,
        "event": {
            "chrom": edge.event.chrom,
            "start": edge.event.start,
            "end": edge.event.end,
            "replacement": edge.event.replacement,
            "metadata": edge.event.metadata,
        },
        "metadata": edge.metadata,
    }


def _frame_to_payload(frame: EditDAGFrame) -> Dict[str, Any]:
    return {
        "kind": "helix.edit_dag.frame.v1",
        "step": frame.step,
        "new_nodes": {node_id: _edit_node_to_payload(node) for node_id, node in frame.new_nodes.items()},
        "new_edges": [_edit_edge_to_payload(edge) for edge in frame.new_edges],
    }


def _frames_to_artifact_payload(frames: Iterable[Dict[str, Any]]) -> Dict[str, Any]:
    nodes: "OrderedDict[str, Dict[str, Any]]" = OrderedDict()
    edges: List[Dict[str, Any]] = []
    artifact = "helix.crispr.edit_dag.v1.1"
    for frame in frames:
        mechanism = (frame.get("meta") or {}).get("mechanism")
        if mechanism == "prime":
            artifact = "helix.prime.edit_dag.v1.1"
        for node_id, node_payload in (frame.get("new_nodes") or {}).items():
            nodes[node_id] = node_payload
        for edge_payload in frame.get("new_edges") or []:
            edges.append(edge_payload)
    if not nodes:
        raise ValueError("Frame stream contained no nodes.")
    root_id = next(iter(nodes))
    return {
        "artifact": artifact,
        "version": "1.1",
        "nodes": nodes,
        "edges": edges,
        "root_id": root_id,
    }


def _read_frame_jsonl(path: Path) -> List[Dict[str, Any]]:
    frames: List[Dict[str, Any]] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            frames.append(json.loads(line))
    return frames


def _open_frame_stream(path_value: str) -> tuple[TextIO, bool]:
    if path_value == "-":
        return sys.stdout, False
    frame_path = Path(path_value)
    frame_path.parent.mkdir(parents=True, exist_ok=True)
    handle = frame_path.open("w", encoding="utf-8")
    return handle, True


def _top_outcomes_from_dag(dag: EditDAG, *, topk: int) -> List[Dict[str, object]]:
    terminals = dag.terminal_nodes()
    if not terminals:
        return []
    weights = [math.exp(node.log_prob) for node in terminals]
    total = sum(weights) or 1.0
    probs = [w / total for w in weights]
    ranked = sorted(zip(terminals, probs), key=lambda pair: pair[1], reverse=True)
    results = []
    for node, prob in ranked[:topk]:
        results.append(
            {
                "node_id": node.id,
                "probability": round(prob, 6),
                "stage": node.metadata.get("stage"),
            }
        )
    return results


def command_dna(args: argparse.Namespace) -> None:
    raw = _load_sequence_arg(args.sequence, args.input, default=bioinformatics.seq)
    genome = bioinformatics.normalize_sequence(raw)
    print(f"Sequence length: {len(genome)} nt")
    print(f"GC content: {bioinformatics.gc_content(genome) * 100:.2f}%")

    if args.window > 0 and len(genome) >= args.window:
        windows = bioinformatics.windowed_gc_content(genome, args.window, args.step)
        if windows:
            richest = max(windows, key=lambda win: win.gc_fraction)
            poorest = min(windows, key=lambda win: win.gc_fraction)
            print(
                f"GC window extremes ({args.window} nt): "
                f"max={richest.gc_fraction*100:.2f}% [{richest.start}-{richest.end}), "
                f"min={poorest.gc_fraction*100:.2f}% [{poorest.start}-{poorest.end})"
            )
    else:
        print("GC window summary skipped (window disabled or longer than the sequence).")

    clusters = bioinformatics.find_kmers_with_differences(genome, args.k, args.max_diff)
    sorted_clusters = sorted(clusters.items(), key=lambda item: item[1]["count"], reverse=True)
    if not sorted_clusters:
        print("No k-mer clusters detected with the current parameters.")
    else:
        print(f"\nTop {min(args.top, len(sorted_clusters))} clusters (k={args.k}, max_diff={args.max_diff}):")
        for canonical, info in sorted_clusters[: args.top]:
            patterns = ",".join(info["patterns"])
            positions = ",".join(map(str, info["positions"]))
            print(f"{canonical}\tcount={info['count']}\tpatterns=[{patterns}]\tpositions=[{positions}]")


def command_spectrum(args: argparse.Namespace) -> None:
    spectrum = _parse_spectrum(args.spectrum)
    if args.spectrum_file:
        spectrum.extend(_parse_spectrum(_read_text(args.spectrum_file)))

    if args.peptide:
        theoretical = cyclospectrum.theoretical_spectrum(args.peptide, cyclic=not args.linear)
        mode = "cyclic" if not args.linear else "linear"
        print(f"{mode.title()} spectrum for {args.peptide}:")
        print(" ".join(str(mass) for mass in theoretical))
        if spectrum:
            score = cyclospectrum.score_peptide(args.peptide, spectrum, cyclic=not args.linear)
            print(f"Score vs provided spectrum: {score}")

    if spectrum:
        hits = cyclospectrum.leaderboard_cyclopeptide_sequencing(
            spectrum,
            leaderboard_size=args.leaderboard,
        )
        if hits:
            print(f"\nLeaderboard candidates (top {len(hits)}):")
            for peptide, score in hits:
                print(f"{peptide}\tscore={score}")
        else:
            print("No leaderboard candidates matched the spectrum.")
    elif not args.peptide:
        raise SystemExit("Provide at least --peptide or --spectrum/--spectrum-file.")


def command_rna_mfe(args: argparse.Namespace) -> None:
    records = read_fasta(args.fasta)
    if not records:
        raise SystemExit(f"No sequences found in {args.fasta}")
    header, seq = records[0]
    result = mfe_dotbracket(seq)
    payload = {
        "sequence_id": header,
        "dotbracket": result["dotbracket"],
        "energy": result["energy"],
        "pairs": result["pairs"],
    }
    text_json = json.dumps(payload, indent=2)
    print(text_json)
    if args.json:
        args.json.write_text(text_json + "\n", encoding="utf-8")
    if args.dotbracket:
        Path(args.dotbracket).write_text(result["dotbracket"] + "\n", encoding="utf-8")
        print(f"Dot-bracket saved to {args.dotbracket}")



def command_rna_ensemble(args: argparse.Namespace) -> None:
    records = read_fasta(args.fasta)
    if not records:
        raise SystemExit(f"No sequences found in {args.fasta}")
    header, seq = records[0]
    ensemble = partition_posteriors(seq)
    posterior = ensemble["P"]
    mea = mea_structure(seq, posterior, gamma=args.gamma)
    centroid = centroid_structure(seq, posterior)
    payload = {
        "sequence_id": header,
        "partition_function": ensemble["Q"],
        "entropy": ensemble["entropy"],
        "p_unpaired": ensemble["p_unpaired"],
        "mea_structure": mea,
        "centroid_structure": centroid,
        "gamma": args.gamma,
    }
    text_json = json.dumps(payload, indent=2)
    print(text_json)
    if args.json:
        args.json.write_text(text_json + "\n", encoding="utf-8")
    if args.dotplot:
        _ensure_viz_available("RNA ensemble dot-plot")
        dot_path = str(args.dotplot)
        spec_path = _default_viz_spec_path(args.dotplot, getattr(args, "save_viz_spec", None))
        dot_payload = _validate_payload_or_exit("viz_rna_dotplot", {"posterior": posterior})
        extra_meta = _input_meta(dot_payload)
        _, spec = plot_rna_dotplot(
            posterior=posterior,
            save=dot_path,
            save_viz_spec=str(spec_path) if spec_path else None,
            extra_meta=extra_meta,
        )
        _write_provenance(
            Path(dot_path),
            schema_kind="viz_rna_dotplot",
            spec=spec,
            input_sha=extra_meta.get("input_sha256"),
            command=current_command_str(),
            viz_spec_path=spec_path,
        )
    if args.arc:
        _ensure_viz_available("RNA ensemble arc plotting")
        viz_rna.plot_arc(mea["dotbracket"], Path(args.arc), title=f"MEA structure ({header})")
    if args.entropy_plot:
        _ensure_viz_available("RNA ensemble entropy plotting")
        viz_rna.plot_entropy(ensemble["entropy"], Path(args.entropy_plot), title=f"Entropy ({header})")


def command_crispr_find_guides(args: argparse.Namespace) -> None:
    raw_sequence = _load_sequence_arg(args.sequence, args.fasta)
    normalized_sequence = bioinformatics.normalize_sequence(raw_sequence)
    pam_def = crispr_pam.get_pam(args.pam)
    window = tuple(args.window) if args.window else None
    try:
        guides = crispr_guide.find_guides(
            normalized_sequence,
            pam_def,
            args.guide_len,
            strand=args.strand,
            window=window,
        )
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    if not args.emit_sequences:
        for guide in guides:
            guide["sequence"] = None

    payload = {
        "schema": {"kind": "crispr.guides", "spec_version": SPEC_VERSION},
        "meta": {"helix_version": HELIX_VERSION, "timestamp": datetime.now(timezone.utc).isoformat()},
        "input_sha256": _sequence_sha256(normalized_sequence),
        "sequence_length": len(normalized_sequence),
        "pam": pam_def,
        "params": {"guide_length": args.guide_len, "strand": args.strand},
        "guides": guides,
    }
    if window:
        payload["params"]["window"] = list(window)

    validated = validate_viz_payload("crispr.guides", payload)
    if not args.emit_sequences:
        for guide in validated.get("guides", []):
            guide.setdefault("sequence", None)
    text = json.dumps(validated, indent=2)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text + "\n", encoding="utf-8")
        print(f"Guide JSON saved to {args.json} ({len(guides)} guides).")
    else:
        print(text)
        print(f"Guides discovered: {len(guides)}")
    if not guides:
        print("No guides found with the current parameters.")


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def command_crispr_offtargets(args: argparse.Namespace) -> None:
    raw_sequence = _load_sequence_arg(args.genome, args.fasta)
    genome = bioinformatics.normalize_sequence(raw_sequence)
    guides_payload = _read_json(args.guides)
    pam_def = crispr_pam.get_pam(args.pam)
    params = {"max_mismatches": args.max_mm, "max_gaps": args.max_gap}

    guide_hits = []
    for guide in guides_payload.get("guides", []):
        hits = score.enumerate_off_targets(
            genome,
            guide,
            pam_def,
            max_mm=args.max_mm,
            max_gap=args.max_gap,
        )
        guide_hits.append({"guide_id": guide.get("id"), "hits": hits})

    payload = {
        "schema": {"kind": "crispr.offtargets", "spec_version": SPEC_VERSION},
        "meta": {"helix_version": HELIX_VERSION, "timestamp": datetime.now(timezone.utc).isoformat()},
        "input_sha256": _sequence_sha256(genome),
        "pam": pam_def,
        "params": params,
        "guides": guide_hits,
    }

    validated = validate_viz_payload("crispr.offtargets", payload)
    text = json.dumps(validated, indent=2)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text + "\n", encoding="utf-8")
        print(f"Off-target JSON saved to {args.json}.")
    else:
        print(text)
    if not guide_hits:
        print("No guides provided; off-target search skipped.")


def command_crispr_score(args: argparse.Namespace) -> None:
    guides_payload = _read_json(args.guides)
    hits_payload = _read_json(args.hits)
    weights = score.load_weights(args.weights)
    guide_lookup = {guide.get("id"): guide for guide in guides_payload.get("guides", []) if guide.get("id")}

    on_params = weights.get("on_target", {})
    off_params = weights.get("off_target", {})

    for entry in hits_payload.get("guides", []):
        guide_id = entry.get("guide_id")
        guide_info = guide_lookup.get(guide_id, {})
        entry["on_target_score"] = score.score_on_target(guide_info, on_params)
        entry["hits"] = score.score_off_targets(entry.get("hits", []), off_params)

    hits_payload.setdefault("meta", {})["weights"] = {"path": str(args.weights) if args.weights else None}
    validated = validate_viz_payload("crispr.offtargets", hits_payload)
    text = json.dumps(validated, indent=2)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text + "\n", encoding="utf-8")
        print(f"Scored off-target JSON saved to {args.json}.")
    else:
        print(text)


def command_crispr_simulate(args: argparse.Namespace) -> None:
    raw_sequence = _load_sequence_arg(args.sequence, args.fasta)
    guides_payload = _read_json(args.guides)
    guide = next((g for g in guides_payload.get("guides", []) if g.get("id") == args.guide_id), None)
    if not guide:
        raise SystemExit(f"Guide '{args.guide_id}' not found in {args.guides}.")
    priors = _read_json(args.priors) if args.priors else None
    result = crispr_simulate.simulate_cut_repair(
        raw_sequence,
        guide,
        priors,
        draws=args.draws,
        seed=args.seed,
        emit_sequence=args.emit_sequences,
    )
    result.setdefault("meta", {})["helix_version"] = HELIX_VERSION
    validated = validate_viz_payload("crispr.sim", result)
    text = json.dumps(validated, indent=2)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text + "\n", encoding="utf-8")
        print(f"Simulation JSON saved to {args.json}.")
    else:
        print(text)

    if args.viz_spec:
        from .gui.modern.builders import build_crispr_viz_spec

        spec = build_crispr_viz_spec(raw_sequence, guide, validated)
        if spec is None:
            print("Unable to derive a visualization spec from crispr.sim output.", file=sys.stderr)
        else:
            _write_viz_spec_output(spec, args.viz_spec)


def command_crispr_genome_sim(args: argparse.Namespace) -> None:
    genome = _load_digital_genome(args.genome)
    cas_system = _resolve_cas_system(args.cas, args.cas_config)
    guide = _build_guide(args.guide_sequence, name=args.guide_name, pam=args.guide_pam)
    try:
        events = simulate_cuts(
            genome,
            cas_system,
            guide,
            max_events=args.max_events,
        ) or []
    except NotImplementedError as exc:
        raise SystemExit(str(exc))

    params: Dict[str, Any] = {}
    if args.max_events is not None:
        params["max_events"] = args.max_events

    payload = {
        "schema": {"kind": "crispr.cut_events", "spec_version": SPEC_VERSION},
        "meta": {
            "helix_version": HELIX_VERSION,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "command": current_command_str(),
        },
        "cas": _serialize_cas_system(cas_system),
        "guide": _serialize_guide(guide),
        "genome": _digital_genome_summary(genome, source=args.genome),
        "events": [_serialize_cut_event(event) for event in events],
    }
    if params:
        payload["params"] = params
    validated = validate_viz_payload("crispr.cut_events", payload)
    _write_json_output(validated, args.json)
    if not events:
        print("No candidate cut events were produced with the current parameters.")


def command_crispr_dag(args: argparse.Namespace) -> None:
    try:
        base_core = CoreDigitalGenome.from_fasta(args.genome)
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc
    cas_system = _resolve_cas_system(args.cas, args.cas_config)
    guides: List[GuideRNA]
    if args.guides_file:
        guides = _load_guides_file(args.guides_file)
    else:
        if not args.guide_sequence:
            raise SystemExit("Provide --guide-sequence or --guides-file.")
        region_meta = {}
        region_hint = _normalize_region_hint(args.region)
        if region_hint:
            region_meta["region"] = region_hint
        guides = [
            _build_guide(
                args.guide_sequence,
                name=args.guide_name,
                pam=args.guide_pam,
                metadata=region_meta or None,
            )
        ]

    global_region = _normalize_region_hint(args.region)
    out_single = args.out or args.json
    multiple_guides = len(guides) > 1

    if multiple_guides and not args.out_dir:
        raise SystemExit("Multi-guide runs require --out-dir.")
    if args.out_dir:
        args.out_dir.mkdir(parents=True, exist_ok=True)

    frame_stream = None
    close_stream = False
    if args.frames:
        frame_stream, close_stream = _open_frame_stream(args.frames)

    for idx, guide in enumerate(guides, start=1):
        guide_region = _normalize_region_hint(guide.metadata.get("region")) if guide.metadata else None
        region = guide_region or global_region
        working_core = base_core.slice_region(region) if region else base_core
        genome = CrisprDigitalGenome(sequences=dict(working_core.sequences))
        frame_consumer = None
        if frame_stream:
            label = guide.name or f"guide_{idx}"

            def frame_consumer(frame: EditDAGFrame, handle=frame_stream, guide_label=label):
                payload = _frame_to_payload(frame)
                payload.setdefault("meta", {})["guide_name"] = guide_label
                payload["meta"]["mechanism"] = "crispr"
                handle.write(json.dumps(payload) + "\n")
                handle.flush()

        dag = build_crispr_edit_dag(
            genome,
            cas_system,
            guide,
            rng_seed=args.seed or 0,
            max_depth=args.max_depth,
            min_prob=args.min_prob,
            max_sites=args.max_sites,
            use_gpu=bool(getattr(args, "use_gpu", False)),
            frame_consumer=frame_consumer,
        )
        transcript = _load_transcript_from_args(getattr(args, "coding_json", None), getattr(args, "coding_transcript", None))
        if transcript:
            if _annotate_dag_with_transcript(dag, genome.sequences, transcript):
                metadata.setdefault("protein_annotation", {})
                metadata["protein_annotation"].update(
                    {
                        "transcript": transcript.name,
                        "transcript_id": transcript.id,
                        "coding_source": str(args.coding_json),
                    }
                )
        metadata: Dict[str, Any] = {
            "helix_version": HELIX_VERSION,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "genome_source": str(args.genome),
            "guide_name": guide.name or f"guide_{idx}",
            "guide_sequence": guide.sequence,
        }
        if region:
            metadata["region"] = region
        if args.guides_file:
            metadata["guides_file"] = str(args.guides_file)
        metadata["physics_backend"] = "gpu" if getattr(args, "use_gpu", False) else "cpu"

        payload = _edit_dag_to_payload(
            dag,
            artifact="helix.crispr.edit_dag.v1.1",
            metadata=metadata,
        )

        if args.out_dir:
            safe_name = _safe_guide_label(guide.name or f"guide_{idx}")
            filename = f"crispr_{idx:03d}_{safe_name}.edit_dag.json"
            out_path = args.out_dir / filename
        else:
            out_path = out_single

        _write_json_output(payload, out_path)
        if out_path:
            print(f"CRISPR edit DAG saved to {out_path}.")
        if not dag.edges:
            label = guide.name or f"guide_{idx}"
            print(f"No edit events generated for guide '{label}'; DAG contains only the root state.")
    if frame_stream and close_stream:
        frame_stream.close()


def command_prime_simulate(args: argparse.Namespace) -> None:
    genome = _load_digital_genome(args.genome)
    peg = _load_peg_from_args(args)
    editor = _load_prime_editor_from_args(args)
    site = locate_prime_target_site(genome, peg)
    if site is None:
        raise SystemExit("Unable to locate a target site for the provided pegRNA.")
    chrom_seq = bioinformatics.normalize_sequence(genome.sequences.get(site.chrom, ""))
    if not chrom_seq:
        raise SystemExit(f"Chromosome '{site.chrom}' contained no sequence data.")
    flank = max(len(peg.pbs or ""), len(peg.rtt or ""), 20)
    window_start = max(0, site.start - flank)
    window_end = min(len(chrom_seq), site.end + flank)
    site_seq = chrom_seq[window_start:window_end]

    priors = resolve_prime_priors(args.priors_profile)
    pam_mask = build_prime_pam_mask(site_seq, peg, args.pam_profile, args.pam_softness)
    payload = simulate_prime_edit(
        site_seq=site_seq,
        peg=peg,
        priors=priors,
        draws=args.draws,
        seed=args.seed,
        emit_sequence=True,
        pam_mask=pam_mask,
    )
    payload.setdefault("meta", {}).update(
        {
            "helix_version": HELIX_VERSION,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "command": current_command_str(),
        }
    )
    payload["params"] = {
        "draws": args.draws,
        "seed": args.seed,
        "priors_profile": args.priors_profile,
        "pam_profile": args.pam_profile,
        "pam_softness": args.pam_softness,
    }
    payload["editor"] = _serialize_prime_editor(editor)
    payload["peg"] = _serialize_peg(peg)
    payload["genome"] = _digital_genome_summary(genome, source=args.genome)
    validated = validate_viz_payload("prime.edit_sim", payload)
    _write_json_output(validated, args.json)

    if args.viz_spec:
        from .gui.modern.builders import build_prime_viz_spec

        spec = build_prime_viz_spec(genome, peg, editor, validated)
        if spec is None:
            print("Unable to derive a visualization spec from prime.edit_sim output.", file=sys.stderr)
        else:
            _write_viz_spec_output(spec, args.viz_spec)


def command_prime_dag(args: argparse.Namespace) -> None:
    try:
        base_core = CoreDigitalGenome.from_fasta(args.genome)
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc
    editor = _load_prime_editor_from_args(args)
    if args.pegs_file:
        pegs = _load_pegs_file(args.pegs_file)
    elif args.peg_config:
        pegs = [_load_peg_from_json(args.peg_config)]
    else:
        peg = _load_peg_from_args(args)
        pegs = [peg]

    global_region = _normalize_region_hint(args.region)
    out_single = args.out or args.json
    multiple_pegs = len(pegs) > 1
    if multiple_pegs and not args.out_dir:
        raise SystemExit("Multi-peg runs require --out-dir.")
    if args.out_dir:
        args.out_dir.mkdir(parents=True, exist_ok=True)

    frame_stream = None
    close_stream = False
    if args.frames:
        frame_stream, close_stream = _open_frame_stream(args.frames)

    for idx, peg in enumerate(pegs, start=1):
        peg_region = _normalize_region_hint(peg.metadata.get("region")) if peg.metadata else None
        region = peg_region or global_region
        working_core = base_core.slice_region(region) if region else base_core
        genome = CrisprDigitalGenome(sequences=dict(working_core.sequences))
        frame_consumer = None
        if frame_stream:
            label = peg.name or f"peg_{idx}"

            def frame_consumer(frame: EditDAGFrame, handle=frame_stream, peg_label=label):
                payload = _frame_to_payload(frame)
                payload.setdefault("meta", {})["peg_name"] = peg_label
                payload["meta"]["mechanism"] = "prime"
                handle.write(json.dumps(payload) + "\n")
                handle.flush()

        dag = build_prime_edit_dag(
            genome,
            editor,
            peg,
            rng_seed=args.seed or 0,
            max_depth=args.max_depth,
            min_prob=args.min_prob,
            frame_consumer=frame_consumer,
        )
        transcript = _load_transcript_from_args(getattr(args, "coding_json", None), getattr(args, "coding_transcript", None))
        if transcript:
            if _annotate_dag_with_transcript(dag, genome.sequences, transcript):
                metadata.setdefault("protein_annotation", {})
                metadata["protein_annotation"].update(
                    {
                        "transcript": transcript.name,
                        "transcript_id": transcript.id,
                        "coding_source": str(args.coding_json),
                    }
                )
        metadata: Dict[str, Any] = {
            "helix_version": HELIX_VERSION,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "genome_source": str(args.genome),
            "peg_name": peg.name or f"peg_{idx}",
        }
        if region:
            metadata["region"] = region
        if args.pegs_file:
            metadata["pegs_file"] = str(args.pegs_file)
        payload = _edit_dag_to_payload(
            dag,
            artifact="helix.prime.edit_dag.v1.1",
            metadata=metadata,
        )
        if args.out_dir:
            safe_name = _safe_guide_label(peg.name or f"peg_{idx}")
            out_path = args.out_dir / f"prime_{idx:03d}_{safe_name}.edit_dag.json"
        else:
            out_path = out_single
        _write_json_output(payload, out_path)
        if out_path:
            print(f"Prime edit DAG saved to {out_path}.")
        if not dag.edges:
            label = peg.name or f"peg_{idx}"
            print(f"No prime edit events generated for peg '{label}'; DAG contains only the root state.")
    if frame_stream and close_stream:
        frame_stream.close()


def command_pcr_dag(args: argparse.Namespace) -> None:
    genome = _load_digital_genome(args.genome)
    primer_pair = _load_primer_pair_from_args(args)
    config = _load_pcr_config_from_args(args)
    dag = pcr_edit_dag(
        genome,
        primer_pair,
        config,
        rng_seed=args.seed or 0,
        min_prob=args.min_prob,
    )
    metadata = {
        "helix_version": HELIX_VERSION,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "genome_source": str(args.genome),
        "primer_pair": primer_pair.name,
        "cycles": config.cycles,
    }
    payload = _edit_dag_to_payload(dag, artifact="helix.pcr.amplicon_dag.v1", metadata=metadata)
    _write_json_output(payload, args.out)
    if not dag.edges:
        print("No PCR amplification products were generated (no binding sites detected).")


def command_edit_dag_viz(args: argparse.Namespace) -> None:
    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    dag = dag_from_payload(payload)
    save_edit_dag_png = _load_dag_viz_exporter()
    save_edit_dag_png(
        dag,
        str(args.out),
        layout=args.layout,
        min_prob_filter=args.min_prob,
        max_time_filter=args.max_time,
    )
    print(f"Edit DAG visualization saved to {args.out}.")


def command_edit_dag_animate(args: argparse.Namespace) -> None:
    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    dag = dag_from_payload(payload)
    animate_edit_dag = _load_dag_animator()
    animate_edit_dag(
        dag,
        str(args.out),
        fps=args.fps,
        layout=args.layout,
        figsize=(args.width, args.height),
        background_color=args.bg_color,
        prob_cmap=args.prob_cmap,
        stage_cmap=args.stage_cmap,
        highlight_color=args.highlight_color,
        dpi=args.dpi,
        min_prob=args.min_prob,
        max_time_filter=args.max_time,
    )
    print(f"Edit DAG animation saved to {args.out} at {args.fps} fps.")


def command_edit_dag_compare(args: argparse.Namespace) -> None:
    payload_a = json.loads(Path(args.a).read_text(encoding="utf-8"))
    payload_b = json.loads(Path(args.b).read_text(encoding="utf-8"))
    dag_a = dag_from_payload(payload_a)
    dag_b = dag_from_payload(payload_b)
    save_compare = _load_dag_compare_exporter()
    save_compare(
        dag_a,
        dag_b,
        str(args.out),
        label_a=args.label_a,
        label_b=args.label_b,
        layout=args.layout,
        min_prob_filter=args.min_prob,
        max_time_filter=args.max_time,
    )
    print(f"Edit DAG comparison saved to {args.out}.")
    if args.summary:
        summary = _compute_dag_diff_summary(dag_a, dag_b, label_a=args.label_a, label_b=args.label_b)
        Path(args.summary).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        print(f"Diff summary saved to {args.summary}.")


def command_edit_dag_report(args: argparse.Namespace) -> None:
    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    dag = dag_from_payload(payload)
    html = render_html_report(dag, png_path=args.png, title=args.title)
    Path(args.out).write_text(html, encoding="utf-8")
    print(f"Edit DAG report saved to {args.out}.")


def command_experiment_new(args: argparse.Namespace) -> None:
    template = _CRISPR_EXPERIMENT_TEMPLATE if args.type == "crispr" else _PRIME_EXPERIMENT_TEMPLATE
    out_path = Path(args.out)
    if out_path.exists() and not args.force:
        raise SystemExit(f"Refusing to overwrite existing file: {out_path}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(template, encoding="utf-8")
    print(f"Wrote experiment template to {out_path}.")


def command_experiment_run(args: argparse.Namespace) -> None:
    if args.out is None:
        raise SystemExit("--out is required for 'experiment run'.")
    frame_stream = None
    close_stream = False
    frame_consumer = None
    if args.frames:
        frame_stream, close_stream = _open_frame_stream(args.frames)

        def frame_consumer(frame: EditDAGFrame, handle=frame_stream):
            payload = _frame_to_payload(frame)
            payload.setdefault("meta", {})["experiment_config"] = str(args.config)
            handle.write(json.dumps(payload) + "\n")
            handle.flush()

    try:
        dag, metadata, artifact = _run_experiment_to_dag(
            args.config,
            use_gpu_override=getattr(args, "use_gpu", None),
            frame_consumer=frame_consumer,
        )
    finally:
        if frame_stream and close_stream:
            frame_stream.close()
    clean_meta = {k: v for k, v in metadata.items() if v is not None}
    payload = _edit_dag_to_payload(dag, artifact=artifact, metadata=clean_meta)
    _write_json_output(payload, args.out)
    print(f"Experiment DAG saved to {args.out}.")


def command_experiment_viz(args: argparse.Namespace) -> None:
    dag, _, _ = _run_experiment_to_dag(args.config, use_gpu_override=getattr(args, "use_gpu", None))
    save_edit_dag_png = _load_dag_viz_exporter()
    save_edit_dag_png(
        dag,
        str(args.out),
        layout=args.layout,
        min_prob_filter=args.min_prob,
        max_time_filter=args.max_time,
    )
    print(f"Experiment visualization saved to {args.out}.")


def command_experiment_report(args: argparse.Namespace) -> None:
    dag, _, _ = _run_experiment_to_dag(args.config, use_gpu_override=getattr(args, "use_gpu", None))
    html = render_html_report(dag, png_path=args.png, title=args.title)
    Path(args.out).write_text(html, encoding="utf-8")
    print(f"Experiment report saved to {args.out}.")


def command_experiment_validate(args: argparse.Namespace) -> None:
    dag, _, _ = _run_experiment_to_dag(args.config, use_gpu_override=getattr(args, "use_gpu", None))
    print(
        f"Experiment '{args.config}' is valid (nodes={len(dag.nodes)}, edges={len(dag.edges)})."
    )


def command_gui(_args: argparse.Namespace) -> None:
    try:
        from helix.gui.main import run_gui
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise SystemExit(
            "PySide6 is required for the Helix GUI. Install it via 'pip install veri-helix[gui]'."
        ) from exc
    run_gui()


def _random_dna(rng: random.Random, length: int) -> str:
    bases = "ACGT"
    return "".join(rng.choice(bases) for _ in range(length))


def _random_crispr_dag(rng: random.Random) -> EditDAG:
    genome_seq = _random_dna(rng, 120)
    genome = CrisprDigitalGenome({"chr": genome_seq})
    start = rng.randrange(0, len(genome_seq) - 21)
    guide_seq = genome_seq[start : start + 20]
    guide = GuideRNA(sequence=guide_seq)
    cas = CasSystem(
        name="dataset-cas9",
        system_type=CasSystemType.CAS9,
        pam_rules=[PAMRule(pattern="NGG")],
        cut_offset=3,
        max_mismatches=3,
    )
    return build_crispr_edit_dag(
        genome,
        cas,
        guide,
        rng_seed=rng.randint(0, 10_000),
        max_depth=1,
        max_sites=5,
    )


def _random_prime_dag(rng: random.Random) -> EditDAG:
    genome_seq = _random_dna(rng, 140)
    genome = CrisprDigitalGenome({"chr": genome_seq})
    start = rng.randrange(0, len(genome_seq) - 25)
    spacer = genome_seq[start : start + 20]
    pbs = genome_seq[start : start + 5]
    rtt = genome_seq[start + 5 : start + 12]
    peg = PegRNA(spacer=spacer, pbs=pbs, rtt=rtt)
    cas = CasSystem(
        name="dataset-cas9",
        system_type=CasSystemType.CAS9,
        pam_rules=[PAMRule(pattern="NGG")],
        cut_offset=3,
        max_mismatches=3,
    )
    editor = PrimeEditor(name="dataset-prime", cas=cas, nick_to_edit_offset=1)
    return build_prime_edit_dag(
        genome,
        editor,
        peg,
        rng_seed=rng.randint(0, 10_000),
        max_depth=2,
    )


def command_edit_dag_generate_dataset(args: argparse.Namespace) -> None:
    rng = random.Random(args.seed)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    next_id = 0
    with out_path.open("w", encoding="utf-8") as handle:
        for _ in range(args.n):
            if args.mechanism == "crispr":
                dag = _random_crispr_dag(rng)
                artifact_name = "helix.crispr.edit_dag.v1.1"
            else:
                dag = _random_prime_dag(rng)
                artifact_name = "helix.prime.edit_dag.v1.1"
            metadata = {
                "helix_version": HELIX_VERSION,
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "mechanism": args.mechanism,
            }
            payload = _edit_dag_to_payload(dag, artifact=artifact_name, metadata=metadata)
            record = {
                "id": next_id,
                "mechanism": args.mechanism,
                "node_count": len(dag.nodes),
                "edge_count": len(dag.edges),
                "top_outcomes": _top_outcomes_from_dag(dag, topk=args.topk),
                "artifact": payload,
            }
            handle.write(json.dumps(record) + "\n")
            next_id += 1
        for frame_path in args.frames_input or []:
            frames = _read_frame_jsonl(frame_path)
            payload = _frames_to_artifact_payload(frames)
            payload.setdefault("meta", {})["frame_source"] = str(frame_path)
            dag = dag_from_payload(payload)
            mechanism = "prime" if payload.get("artifact", "").startswith("helix.prime") else "crispr"
            record = {
                "id": next_id,
                "mechanism": mechanism,
                "node_count": len(dag.nodes),
                "edge_count": len(dag.edges),
                "top_outcomes": _top_outcomes_from_dag(dag, topk=args.topk),
                "artifact": payload,
                "frame_source": str(frame_path),
            }
            handle.write(json.dumps(record) + "\n")
            next_id += 1
    print(f"Edit DAG dataset written to {args.out}.")


def command_protein(args: argparse.Namespace) -> None:
    if not PROTEIN_AVAILABLE:
        raise SystemExit("Biopython is required for protein helpers. Install it via 'pip install biopython'.")
    raw = _load_sequence_arg(args.sequence, args.input)
    summary = protein_module.summarize_sequence(raw)
    print(f"Length: {summary.length}")
    print(f"Molecular weight: {summary.molecular_weight:.2f} Da")
    print(f"GRAVY: {summary.gravy:.3f}")
    print(f"Aromaticity: {summary.aromaticity:.3f}")
    print(f"Instability index: {summary.instability_index:.2f}")
    print(f"Charge @ pH 7.0: {summary.charge_at_pH7:.2f}")

    if summary.length >= args.window:
        windows = protein_module.hydropathy_profile(
            summary.sequence,
            window=args.window,
            step=args.step,
            scale=args.scale,
        )
        sorted_windows = sorted(windows, key=lambda win: win.score, reverse=True)
        print(f"\nTop {min(args.top, len(sorted_windows))} hydrophobic windows:")
        for window in sorted_windows[: args.top]:
            print(f"{window.start:>4}-{window.end:<4}\tscore={window.score:.3f}")
    else:
        print("Hydropathy profile skipped (sequence shorter than requested window).")


def _orf_to_dict(orf) -> dict:
    return {
        "start": orf.start,
        "end": orf.end,
        "strand": orf.strand,
        "frame": orf.frame,
        "length_nt": orf.length_nt(),
        "length_aa": orf.length_aa(),
        "peptide": orf.peptide,
    }


def _triage_to_dict(report: triage.TriageReport) -> dict:
    return {
        "sequence": report.sequence,
        "skew": report.skew,
        "clusters": [
            {
                "canonical": cluster.canonical,
                "count": cluster.count,
                "patterns": list(cluster.patterns),
                "positions": list(cluster.positions),
            }
            for cluster in report.clusters
        ],
        "orfs": [_orf_to_dict(orf) for orf in report.orfs],
    }


def command_triage(args: argparse.Namespace) -> None:
    raw = _load_sequence_arg(args.sequence, args.input)
    report = triage.compute_triage_report(
        raw,
        k=args.k,
        max_diff=args.max_diff,
        min_orf_length=args.min_orf_length,
    )
    print(f"Sequence length: {len(report.sequence)} nt")
    print(f"Skew span: min={min(report.skew)} max={max(report.skew)}")
    print(f"Detected {len(report.clusters)} k-mer clusters and {len(report.orfs)} ORFs >= {args.min_orf_length} nt.")
    if report.clusters:
        print("\nTop clusters:")
        for cluster in report.clusters[: args.top]:
            patterns = ",".join(cluster.patterns)
            print(f"{cluster.canonical}\tcount={cluster.count}\tpatterns=[{patterns}]")
    if report.orfs:
        print("\nTop ORFs:")
        for orf in report.orfs[: args.top]:
            print(
                f"start={orf.start} end={orf.end} strand={orf.strand} frame={orf.frame} "
                f"length_nt={orf.length_nt()} length_aa={orf.length_aa()}"
            )

    if args.json:
        payload = _triage_to_dict(report)
        args.json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"\nJSON report saved to {args.json}")


def command_string_search(args: argparse.Namespace) -> None:
    records = read_fasta(args.input)
    if not records:
        raise SystemExit(f"No sequences found in {args.input}")

    pattern = args.pattern.upper()
    results = []
    for idx, (header, sequence) in enumerate(records):
        label = header or f"seq_{idx}"
        seq = sequence.upper()
        if args.k == 0:
            fm_index = string_fm.build_fm(seq)
            hits = string_fm.search(fm_index, pattern)
            payload = {
                "sequence_id": label,
                "mode": "exact",
                "hits": hits,
            }
        else:
            matches = string_edit.myers_search(pattern, seq, args.k)
            payload = {
                "sequence_id": label,
                "mode": f"myers_k_{args.k}",
                "pattern": pattern,
                "matches": matches,
            }
        results.append(payload)

    output = {
        "meta": {
            "pattern": pattern,
            "k": args.k,
            "sequence_count": len(records),
        },
        "results": results,
    }

    text = json.dumps(output, indent=2)
    print(text)
    if args.json:
        args.json.write_text(text + "\n", encoding="utf-8")


def command_seed_index(args: argparse.Namespace) -> None:
    records = read_fasta(args.input)
    if not records:
        raise SystemExit(f"No sequences found in {args.input}")

    all_results = []
    for header, seq in records:
        if args.method == "minimizer":
            seeds = seed_minimizers(seq, args.k, args.window)
        else:
            seeds = seed_syncmers(seq, args.k, args.sync)
        payload = {
            "sequence_id": header or "seq",
            "length": len(seq),
            "method": args.method,
            "seed_count": len(seeds),
            "seeds": [{"pos": pos, "kmer": kmer, "hash": h} for pos, kmer, h in seeds],
        }
        all_results.append(payload)
        if args.plot:
            _ensure_viz_available("helix seed index plotting")
            from .viz import seed as viz_seed  # local import to defer matplotlib

            if len(records) == 1:
                output = args.plot
            else:
                output = args.plot.with_name(f"{args.plot.stem}_{len(all_results)}{args.plot.suffix}")
            viz_seed.plot_density(seeds, len(seq), output, title=f"{payload['sequence_id']} ({args.method})")

    data = {"meta": {"method": args.method, "k": args.k}, "results": all_results}
    text = json.dumps(data, indent=2)
    print(text)
    if args.json:
        args.json.write_text(text + "\n", encoding="utf-8")


def command_seed_map(args: argparse.Namespace) -> None:
    ref_records = read_fasta(args.ref)
    if not ref_records:
        raise SystemExit(f"No reference sequences in {args.ref}")
    ref_header, ref_seq = ref_records[0]
    ref_seeds = seed_minimizers(ref_seq, args.k, args.window)
    index = {}
    for pos, _, h in ref_seeds:
        index.setdefault(h, []).append(pos)

    read_records = read_fasta(args.reads)
    results = []
    for header, seq in read_records:
        read_seeds = seed_minimizers(seq, args.k, args.window)
        matches = []
        for pos, _, h in read_seeds:
            if h not in index:
                continue
            for ref_pos in index[h][: args.max_matches]:
                seed = SeedMatch(ref_pos=ref_pos, read_pos=pos, length=args.k)
                aln = extend_alignment(seed, ref_seq, seq, band=args.band, xdrop=args.xdrop)
                matches.append({"seed_ref": ref_pos, "seed_read": pos, "alignment": aln})
        results.append(
            {
                "read_id": header or "read",
                "read_length": len(seq),
                "seed_hits": len(matches),
                "alignments": matches,
            }
        )

    payload = {
        "meta": {
            "reference": ref_header,
            "ref_length": len(ref_seq),
            "k": args.k,
            "window": args.window,
            "band": args.band,
            "xdrop": args.xdrop,
        },
        "results": results,
    }
    payload = _stamp_spec_version(payload)
    payload = _validate_payload_or_exit("viz_alignment_ribbon", payload)
    text = json.dumps(payload, indent=2)
    print(text)
    if args.json:
        args.json.write_text(text + "\n", encoding="utf-8")


def _load_reads_from_paths(paths: List[Path]) -> List[str]:
    reads: List[str] = []
    for path in paths:
        records = read_fasta(path)
        for _, seq in records:
            reads.append(seq)
    return reads


def command_dbg_build(args: argparse.Namespace) -> None:
    reads = _load_reads_from_paths(args.reads)
    graph = graph_build_dbg(reads, args.k)
    payload = graph_serialize(graph)
    args.graph.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Graph saved to {args.graph} (nodes={len(graph.nodes)})")
    if args.graphml:
        args.graphml.write_text(graph_export_graphml(graph), encoding="utf-8")
        print(f"GraphML saved to {args.graphml}")


def command_dbg_clean(args: argparse.Namespace) -> None:
    graph_json = json.loads(args.graph.read_text(encoding="utf-8"))
    graph = graph_deserialize(graph_json)
    graph_clean_dbg(graph, tips=not args.no_tips, bubbles=not args.no_bubbles, tip_length=args.tip_length)
    payload = graph_serialize(graph)
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Cleaned graph written to {args.out}")


def command_dbg_color(args: argparse.Namespace) -> None:
    if args.labels and len(args.labels) != len(args.reads):
        raise SystemExit("Number of labels must match number of read files.")
    labels = args.labels or [path.stem for path in args.reads]
    reads_by_sample = {}
    for label, path in zip(labels, args.reads):
        reads_by_sample[label] = [seq for _, seq in read_fasta(path)]
    colored = build_colored_dbg(reads_by_sample, args.k)
    presence = {node: sorted(samples) for node, samples in colored.presence.items() if samples}
    payload = {
        "k": colored.graph.k,
        "graph": graph_serialize(colored.graph),
        "samples": colored.samples,
        "presence": presence,
    }
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Colored graph written to {args.out}")


def command_sketch_build(args: argparse.Namespace) -> None:
    records = read_fasta(args.fasta)
    if not records:
        raise SystemExit(f"No sequences found in {args.fasta}")
    header, seq = records[0]
    if args.method == "minhash":
        sketch = compute_minhash(seq, k=args.k, sketch_size=args.size)
        payload = {
            "sequence_id": header,
            "method": "minhash",
            "sketch": sketch.to_dict(),
        }
    else:
        sketch = compute_hll(seq, k=args.k, p=args.precision)
        payload = {
            "sequence_id": header,
            "method": "hll",
            "sketch": sketch.to_dict(),
        }
    text = json.dumps(payload, indent=2)
    print(text)
    if args.json:
        args.json.write_text(text + "\n", encoding="utf-8")


def command_sketch_compare(args: argparse.Namespace) -> None:
    records_a = read_fasta(args.fasta_a)
    records_b = read_fasta(args.fasta_b)
    if not records_a or not records_b:
        raise SystemExit("Both FASTA inputs must contain sequences.")
    if args.method == "minhash":
        sketch_a = compute_minhash(records_a[0][1], k=args.k, sketch_size=args.size)
        sketch_b = compute_minhash(records_b[0][1], k=args.k, sketch_size=args.size)
        distance = mash_distance(sketch_a, sketch_b)
        labels = [records_a[0][0] or "seq_a", records_b[0][0] or "seq_b"]
        matrix = [
            [0.0, float(distance)],
            [float(distance), 0.0],
        ]
        payload = {
            "method": "minhash",
            "distance": distance,
            "sketch_a": sketch_a.to_dict(),
            "sketch_b": sketch_b.to_dict(),
            "labels": labels,
            "matrix": matrix,
        }
    else:
        sketch_a = compute_hll(records_a[0][1], k=args.k, p=args.precision)
        sketch_b = compute_hll(records_b[0][1], k=args.k, p=args.precision)
        union = union_hll(sketch_a, sketch_b)
        est_a = sketch_a.estimate()
        est_b = sketch_b.estimate()
        est_union = union.estimate()
        inter = max(0.0, est_a + est_b - est_union)
        jaccard = inter / est_union if est_union else 0.0
        labels = [records_a[0][0] or "seq_a", records_b[0][0] or "seq_b"]
        distance = 1.0 - jaccard
        matrix = [
            [0.0, float(distance)],
            [float(distance), 0.0],
        ]
        payload = {
            "method": "hll",
            "jaccard": jaccard,
            "cardinality_a": est_a,
            "cardinality_b": est_b,
            "cardinality_union": est_union,
            "labels": labels,
            "matrix": matrix,
        }
    payload = _stamp_spec_version(payload, to_meta=False)
    _validate_payload_or_exit("viz_distance_heatmap", payload)
    text = json.dumps(payload, indent=2)
    print(text)
    if args.json:
        args.json.write_text(text + "\n", encoding="utf-8")


def command_motif_find(args: argparse.Namespace) -> None:
    records = read_fasta(args.fasta)
    if not records:
        raise SystemExit(f"No sequences found in {args.fasta}")
    sequences = [seq for _, seq in records]
    kwargs = {"iterations": args.iterations}
    if args.solver == "steme":
        kwargs["restarts"] = args.restarts
    if args.solver == "online":
        kwargs["learning_rate"] = args.learning_rate
        kwargs["passes"] = args.passes
    result = discover_motifs(sequences, width=args.width, solver=args.solver, **kwargs)
    payload = result.as_json()
    payload = _stamp_spec_version(payload, to_meta=False)
    payload = _validate_payload_or_exit("viz_motif_logo", payload)
    text = json.dumps(payload, indent=2)
    print(text)
    if args.json:
        args.json.write_text(text + "\n", encoding="utf-8")
    if args.plot:
        spec_path = _default_viz_spec_path(args.plot, getattr(args, "plot_viz_spec", None))
        plot_motif_logo(
            pwm=result.pwm,
            title=f"Motif consensus {payload['consensus']}",
            save=str(args.plot),
            save_viz_spec=str(spec_path) if spec_path else None,
        )


def command_veribiota_export(args: argparse.Namespace) -> None:
    payload = _read_json(args.input)
    lean_text = dag_payload_to_lean(
        payload,
        dag_name=args.dag_name,
        module_name=args.module_name,
        import_module=args.lean_import,
        include_eval=not args.skip_eval,
        include_theorem=not args.skip_theorem,
    )
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(lean_text + "\n", encoding="utf-8")
        print(f"Lean export saved to {args.out}.")
    else:
        print(lean_text)


def command_veribiota_lean_check(args: argparse.Namespace) -> None:
    payload = _read_json(args.input)
    dag_name = args.dag_name or args.input.stem
    summary = build_lean_check(payload, dag_name)
    if not args.skip_validate:
        validate_lean_check(summary, prob_tolerance=args.prob_tol)
    text = json.dumps(summary, indent=2)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
        print(f"Lean-check summary saved to {args.out}.")
    else:
        print(text)


def command_veribiota_preflight(args: argparse.Namespace) -> None:
    checks = args.checks
    payloads = args.payloads or []
    if payloads and len(payloads) != len(checks):
        raise SystemExit("--payloads length must match --checks when provided.")
    for idx, check_path in enumerate(checks):
        check = _read_json(check_path)
        try:
            validate_lean_check(check, prob_tolerance=args.prob_tol)
        except LeanCheckError as exc:
            raise SystemExit(f"{check_path}: {exc}") from exc
        if payloads:
            payload = _read_json(payloads[idx])
            dag_name = check.get("dag_name") or payload.get("dag_name") or check_path.stem
            reference = build_lean_check(payload, dag_name)
            for node_id, summary in reference["nodes"].items():
                expected = summary.get("seq_hashes") or {}
                actual = (check.get("nodes") or {}).get(node_id, {}).get("seq_hashes") or {}
                if expected != actual:
                    raise SystemExit(
                        f"{check_path}: sequence hashes for node '{node_id}' did not match payload."
                    )
    print(f"Validated {len(checks)} lean-check file(s).")


def command_veribiota_export_dags(args: argparse.Namespace) -> None:
    inputs: List[Path] = args.inputs
    if args.dag_names and len(args.dag_names) != len(inputs):
        raise SystemExit("--dag-names length must match --inputs.")
    dag_names = args.dag_names or [path.stem for path in inputs]
    payloads = [_read_json(path) for path in inputs]
    lean_text = dag_payloads_to_lean(
        payloads,
        dag_names,
        module_name=args.module_name,
        import_module=args.lean_import,
        include_eval=args.eval,
        include_theorem=not args.skip_theorem,
        list_name=args.list_name,
        theorem_name=args.theorem_name,
    )
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(lean_text + "\n", encoding="utf-8")
        print(f"Lean export saved to {args.out}.")
    else:
        print(lean_text)


def command_veribiota_export_suite(args: argparse.Namespace) -> None:
    inputs: List[Path] = args.inputs
    if not inputs:
        raise SystemExit("Provide at least one DAG JSON via --inputs.")
    if args.dag_names and len(args.dag_names) != len(inputs):
        raise SystemExit("--dag-names length must match --inputs when provided.")
    dag_names = args.dag_names or [path.stem for path in inputs]
    payloads = [_read_json(path) for path in inputs]
    module_path = Path(args.module_path)
    if module_path.suffix != ".lean":
        module_path = module_path.with_suffix(".lean")
    module_name = args.module_name or module_name_from_path(module_path)
    lean_text = dag_payloads_to_lean(
        payloads,
        dag_names,
        module_name=module_name,
        import_module=args.lean_import,
        include_eval=args.eval,
        include_theorem=not args.skip_theorem,
        list_name=args.list_name,
        theorem_name=args.theorem_name,
    )
    root_path = args.veribiota_root
    target_path = root_path / module_path
    target_path.parent.mkdir(parents=True, exist_ok=True)
    target_path.write_text(lean_text + "\n", encoding="utf-8")
    print(f"Lean suite exported to {target_path}.")


def command_workflows(args: argparse.Namespace) -> None:
    from helix_workflows import run_workflow_config

    if getattr(args, "as_json", False) and not getattr(args, "with_schema", False):
        raise SystemExit("--as-json requires --with-schema.")
    results = run_workflow_config(
        args.config,
        output_dir=args.output_dir,
        selected=args.name,
    )
    schema_json: List[Dict[str, Any]] = []
    for result in results:
        print(f"Workflow '{result.name}' completed. Logs at {result.output_dir}")
        for step in result.steps:
            print(f"  - {step.command} -> {step.output_path or 'stdout captured'}")
        if getattr(args, "with_schema", False):
            rows: List[Dict[str, Any]] = []
            for idx, step in enumerate(result.steps, start=1):
                rows.append(
                    {
                        "step": idx,
                        "command": step.command,
                        "schema_kind": step.schema_kind,
                        "spec_version": step.schema_version,
                        "input_sha256": step.schema_hash,
                        "status": "ok" if step.schema_kind else "n/a",
                    }
                )
            if getattr(args, "as_json", False):
                schema_json.append({"workflow": result.name, "steps": rows})
            else:
                print("  Schema provenance:")
                header = f"{'Step':<4} {'Command':<15} {'Schema':<25} {'Spec':<6} {'SHA256':<64} {'Status':<6}"
                print("    " + header)
                for row in rows:
                    row_str = (
                        f"{row['step']:<4} {row['command']:<15} "
                        f"{(row['schema_kind'] or '-'):<25} {(row['spec_version'] or '-'):<6} "
                        f"{(row['input_sha256'] or '-')[:64]:<64} {row['status']:<6}"
                    )
                    print("    " + row_str)
    if getattr(args, "with_schema", False) and getattr(args, "as_json", False):
        print(json.dumps(schema_json, indent=2))


def _require_matplotlib():
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise SystemExit(f"matplotlib is required for visualization commands ({exc}).")
    return plt


def command_viz_triage(args: argparse.Namespace) -> None:
    plt = _require_matplotlib()
    data = json.loads(args.json.read_text(encoding="utf-8"))
    skew = data.get("skew", [])
    clusters = data.get("clusters", [])
    orfs = data.get("orfs", [])

    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=False)
    axes[0].plot(skew)
    axes[0].set_title("GC Skew")
    axes[0].set_xlabel("Position")
    axes[0].set_ylabel("Cumulative skew")

    subset_orfs = orfs[: args.top]
    if subset_orfs:
        y_pos = range(len(subset_orfs))
        lengths = [entry["length_nt"] for entry in subset_orfs]
        labels = [f"{entry['strand']}:{entry['frame']}" for entry in subset_orfs]
        axes[1].barh(list(y_pos), lengths)
        axes[1].set_yticks(list(y_pos))
        axes[1].set_yticklabels(labels)
        axes[1].set_xlabel("Length (nt)")
        axes[1].set_title("Top ORFs")
    else:
        axes[1].text(0.5, 0.5, "No ORFs", ha="center", va="center")

    subset_clusters = clusters[: args.top]
    if subset_clusters:
        axes[2].bar([c["canonical"] for c in subset_clusters], [c["count"] for c in subset_clusters])
        axes[2].tick_params(axis="x", rotation=45)
        axes[2].set_title("Top k-mer clusters")
        axes[2].set_ylabel("Count")
    else:
        axes[2].text(0.5, 0.5, "No clusters", ha="center", va="center")

    fig.tight_layout()
    fig.savefig(args.output)
    if args.show:  # pragma: no cover - interactive path
        plt.show()
    plt.close(fig)
    print(f"Triage visualization saved to {args.output}")


def command_viz_hydropathy(args: argparse.Namespace) -> None:
    if not PROTEIN_AVAILABLE:
        raise SystemExit("Biopython is required for hydropathy visualization (pip install biopython).")
    plt = _require_matplotlib()
    raw = _load_sequence_arg(args.sequence, args.input)
    windows = protein_module.hydropathy_profile(raw, window=args.window, step=args.step, scale=args.scale)
    if not windows:
        raise SystemExit("Sequence is shorter than the requested window size.")

    xs = [window.start for window in windows]
    ys = [window.score for window in windows]
    plt.figure(figsize=(10, 4))
    plt.plot(xs, ys, marker="o")
    plt.axhline(0, color="black", linewidth=0.5)
    plt.title(f"Hydropathy profile (window={args.window}, scale={args.scale})")
    plt.xlabel("Position")
    plt.ylabel("Score")
    plt.tight_layout()
    plt.savefig(args.output)
    if args.show:  # pragma: no cover - interactive path
        plt.show()
    plt.close()
    print(f"Hydropathy chart saved to {args.output}")


def command_viz_minimizers(args: argparse.Namespace) -> None:
    if getattr(args, "schema", False):
        _print_schema_help("viz_minimizers", "minimizers")
        return
    _ensure_viz_available("helix viz minimizers")
    if not args.input:
        raise SystemExit("--input is required unless --schema is provided.")
    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    payload = _validate_payload_or_exit("viz_minimizers", payload)
    extra_meta = _input_meta(payload)
    seq_len = int(payload["sequence_length"])
    minimizers = payload.get("minimizers", [])
    save_path = args.save
    save = str(save_path) if save_path else None
    spec_path = _default_viz_spec_path(save_path, args.save_viz_spec)
    _, spec = plot_minimizer_density(
        sequence_length=seq_len,
        minimizers=minimizers,
        bin_count=args.bins,
        save=save,
        save_viz_spec=str(spec_path) if spec_path else None,
        extra_meta=extra_meta,
    )
    if save_path:
        _write_provenance(
            save_path,
            schema_kind="viz_minimizers",
            spec=spec,
            input_sha=extra_meta.get("input_sha256"),
            command=current_command_str(),
            viz_spec_path=spec_path,
        )


def command_viz_seed_chain(args: argparse.Namespace) -> None:
    if getattr(args, "schema", False):
        _print_schema_help("viz_seed_chain", "seed-chain")
        return
    _ensure_viz_available("helix viz seed-chain")
    if not args.input:
        raise SystemExit("--input is required unless --schema is provided.")
    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    payload = _validate_payload_or_exit("viz_seed_chain", payload)
    extra_meta = _input_meta(payload)
    save_path = args.save
    save = str(save_path) if save_path else None
    spec_path = _default_viz_spec_path(save_path, args.save_viz_spec)
    _, spec = plot_seed_chain(
        ref_length=int(payload["ref_length"]),
        qry_length=int(payload["qry_length"]),
        chains=payload.get("chains", []),
        save=save,
        save_viz_spec=str(spec_path) if spec_path else None,
        extra_meta=extra_meta,
    )
    if save_path:
        _write_provenance(
            save_path,
            schema_kind="viz_seed_chain",
            spec=spec,
            input_sha=extra_meta.get("input_sha256"),
            command=current_command_str(),
            viz_spec_path=spec_path,
        )


def command_viz_rna_dotplot(args: argparse.Namespace) -> None:
    if getattr(args, "schema", False):
        _print_schema_help("viz_rna_dotplot", "rna-dotplot")
        return
    _ensure_viz_available("helix viz rna-dotplot")
    if not args.input:
        raise SystemExit("--input is required unless --schema is provided.")
    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    payload = _validate_payload_or_exit("viz_rna_dotplot", payload)
    extra_meta = _input_meta(payload)
    save_path = args.save
    save = str(save_path) if save_path else None
    spec_path = _default_viz_spec_path(save_path, args.save_viz_spec)
    _, spec = plot_rna_dotplot(
        posterior=payload["posterior"],
        vmin=args.vmin,
        vmax=args.vmax,
        save=save,
        save_viz_spec=str(spec_path) if spec_path else None,
        extra_meta=extra_meta,
    )
    if save_path:
        _write_provenance(
            save_path,
            schema_kind="viz_rna_dotplot",
            spec=spec,
            input_sha=extra_meta.get("input_sha256"),
            command=current_command_str(),
            viz_spec_path=spec_path,
        )


def command_viz_alignment_ribbon(args: argparse.Namespace) -> None:
    if getattr(args, "schema", False):
        _print_schema_help("viz_alignment_ribbon", "alignment-ribbon")
        return
    _ensure_viz_available("helix viz alignment-ribbon")
    if not args.input:
        raise SystemExit("--input is required unless --schema is provided.")
    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    payload = _validate_payload_or_exit("viz_alignment_ribbon", payload)
    extra_meta = _input_meta(payload)
    alignment: Dict[str, Any]
    ref_length: int | None = None
    qry_length: int | None = None
    metadata: Dict[str, Any] | None = None
    title: str | None = args.title

    if "results" in payload:
        results = payload.get("results", [])
        if not results:
            raise SystemExit("No alignment results in payload.")
        target = None
        if args.read_id:
            for entry in results:
                if entry.get("read_id") == args.read_id:
                    target = entry
                    break
            if target is None:
                raise SystemExit(f"Read '{args.read_id}' not found in payload.")
        else:
            target = results[0]
        alignments = target.get("alignments", [])
        if not alignments:
            raise SystemExit(f"No alignments for read '{target.get('read_id', 'read')}'.")
        if args.alignment_index < 0 or args.alignment_index >= len(alignments):
            raise SystemExit("alignment-index out of range for selected read.")
        entry = alignments[args.alignment_index]
        alignment = entry.get("alignment", entry)
        ref_length = payload.get("meta", {}).get("ref_length", alignment.get("ref_end"))
        qry_length = target.get("read_length", alignment.get("read_end"))
        metadata = {
            "read_id": target.get("read_id"),
            "seed_hits": target.get("seed_hits"),
        }
        if entry.get("metadata"):
            metadata.update(entry["metadata"])
        if payload.get("meta", {}).get("reference"):
            metadata["reference"] = payload["meta"]["reference"]
        if not title:
            title = f"{target.get('read_id', 'read')} vs {payload.get('meta', {}).get('reference', 'reference')}"
    else:
        alignment = payload.get("alignment", payload)
        ref_length = payload.get("ref_length", alignment.get("ref_length") or alignment.get("ref_end"))
        qry_length = payload.get("qry_length", alignment.get("qry_length") or alignment.get("read_end"))
        metadata = payload.get("metadata", alignment.get("metadata"))
        if not title:
            title = payload.get("metadata", {}).get("name") or "Alignment ribbon"

    if ref_length is None or qry_length is None:
        raise SystemExit("ref_length and qry_length must be provided.")
    save_path = args.save
    save = str(save_path) if save_path else None
    spec_path = _default_viz_spec_path(save_path, args.save_viz_spec)
    _, spec = plot_alignment_ribbon(
        ref_length=int(ref_length or 0),
        qry_length=int(qry_length or 0),
        alignment=alignment,
        metadata=metadata,
        title=title,
        save=save,
        save_viz_spec=str(spec_path) if spec_path else None,
        extra_meta=extra_meta,
    )
    if save_path:
        _write_provenance(
            save_path,
            schema_kind="viz_alignment_ribbon",
            spec=spec,
            input_sha=extra_meta.get("input_sha256"),
            command=current_command_str(),
            viz_spec_path=spec_path,
        )


def command_viz_distance_heatmap(args: argparse.Namespace) -> None:
    if getattr(args, "schema", False):
        _print_schema_help("viz_distance_heatmap", "distance-heatmap")
        return
    _ensure_viz_available("helix viz distance-heatmap")
    if not args.input:
        raise SystemExit("--input is required unless --schema is provided.")
    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    payload = _validate_payload_or_exit("viz_distance_heatmap", payload)
    extra_meta = _input_meta(payload)
    if "matrix" not in payload or "labels" not in payload:
        raise SystemExit("Distance payload must include 'matrix' and 'labels'.")
    save_path = args.save
    save = str(save_path) if save_path else None
    spec_path = _default_viz_spec_path(save_path, args.save_viz_spec)
    _, spec = plot_distance_heatmap(
        matrix=payload["matrix"],
        labels=payload["labels"],
        method=payload.get("method", "minhash"),
        save=save,
        save_viz_spec=str(spec_path) if spec_path else None,
        extra_meta=extra_meta,
    )
    if save_path:
        _write_provenance(
            save_path,
            schema_kind="viz_distance_heatmap",
            spec=spec,
            input_sha=extra_meta.get("input_sha256"),
            command=current_command_str(),
            viz_spec_path=spec_path,
        )


def command_viz_motif_logo(args: argparse.Namespace) -> None:
    if getattr(args, "schema", False):
        _print_schema_help("viz_motif_logo", "motif-logo")
        return
    _ensure_viz_available("helix viz motif-logo")
    if not args.input:
        raise SystemExit("--input is required unless --schema is provided.")
    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    payload = _validate_payload_or_exit("viz_motif_logo", payload)
    extra_meta = _input_meta(payload)
    if "pwm" not in payload:
        raise SystemExit("Motif payload must include 'pwm'.")
    save_path = args.save
    save = str(save_path) if save_path else None
    spec_path = _default_viz_spec_path(save_path, args.save_viz_spec)
    _, spec = plot_motif_logo(
        pwm=payload["pwm"],
        title=args.title or payload.get("consensus", "Motif logo"),
        alphabet=payload.get("alphabet"),
        background=payload.get("background"),
        save=save,
        save_viz_spec=str(spec_path) if spec_path else None,
        extra_meta=extra_meta,
    )
    if save_path:
        _write_provenance(
            save_path,
            schema_kind="viz_motif_logo",
            spec=spec,
            input_sha=extra_meta.get("input_sha256"),
            command=current_command_str(),
            viz_spec_path=spec_path,
        )


def _load_modern_specs_from_paths(paths: Sequence[Path] | None) -> list[EditVisualizationSpec]:
    from .gui.modern.spec import load_viz_specs

    specs: list[EditVisualizationSpec] = []
    if not paths:
        return specs
    for path in paths:
        specs.extend(load_viz_specs(path))
    return specs


def command_viz_modern(args: argparse.Namespace) -> None:
    _ensure_modern_gui_available()
    from .gui.modern.qt import run_modern_viz

    specs = _load_modern_specs_from_paths(args.spec)
    if specs:
        run_modern_viz(specs)
    else:
        run_modern_viz(None)


def command_viz_modern_export(args: argparse.Namespace) -> None:
    from .gui.modern.export import export_unity_bundle

    specs = _load_modern_specs_from_paths(args.spec)
    if not specs:
        raise SystemExit("Provide at least one --spec path containing edit events.")
    output_paths = export_unity_bundle(specs, output_dir=args.output_dir, camera_samples=args.camera_samples)
    for path in output_paths:
        print(f"Wrote {path}")


def command_viz_sim_builder(args: argparse.Namespace) -> None:
    _ensure_modern_gui_available()
    from .gui.simbuilder import run_sim_builder

    run_sim_builder()


def command_viz_crispr_track(args: argparse.Namespace) -> None:
    _ensure_viz_available("helix viz crispr-track")
    payload = _validate_payload_or_exit("crispr.sim", _read_json(args.input))
    extra_meta = _input_meta(payload)
    save_path = args.save
    save = str(save_path) if save_path else None
    spec_path = _default_viz_spec_path(save_path, args.save_viz_spec)
    _, spec = render_crispr_track(
        payload,
        save=save,
        show=args.show,
        extra_meta=extra_meta,
        save_viz_spec=str(spec_path) if spec_path else None,
    )
    if save_path:
        _write_provenance(
            save_path,
            schema_kind="viz_crispr_track",
            spec=spec,
            input_sha=extra_meta.get("input_sha256"),
            command=current_command_str(),
            viz_spec_path=spec_path,
        )


def command_viz_schema(args: argparse.Namespace) -> None:
    kind = args.kind
    try:
        description = describe_schema(kind)
    except SchemaError as exc:
        raise SystemExit(str(exc))
    print(description)


def command_schema_manifest(args: argparse.Namespace) -> None:
    data = manifest()
    text = json.dumps(data, indent=2)
    if args.out:
        args.out.write_text(text + "\n", encoding="utf-8")
        print(f"Schema manifest written to {args.out}")
    else:
        print(text)


def command_schema_diff(args: argparse.Namespace) -> None:
    base_manifest = load_manifest(args.base)
    if args.target:
        target_manifest = load_manifest(args.target)
    else:
        target_manifest = manifest()
    diff = diff_manifests(base_manifest, target_manifest)
    fmt = args.format
    if fmt == "json":
        print(json.dumps(diff, indent=2))
    else:
        print(format_manifest_diff(diff, fmt="table"))


def command_demo_viz(args: argparse.Namespace) -> None:
    _ensure_viz_available("helix demo viz")
    output_dir: Path = args.output
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, spec in VIZ_DEMO_PAYLOADS.items():
        data = json.loads(json.dumps(spec["data"]))  # shallow copy
        kind = spec["kind"]
        data = _validate_payload_or_exit(kind, data)
        payload_path = output_dir / f"{name}.json"
        payload_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
        png_path = output_dir / f"{name}.png"
        viz_path = output_dir / f"{name}.viz.json"
        extra_meta = _input_meta(data)
        if name == "minimizers":
            _, spec = plot_minimizer_density(
                sequence_length=data["sequence_length"],
                minimizers=data["minimizers"],
                save=str(png_path),
                save_viz_spec=str(viz_path),
                extra_meta=extra_meta,
            )
            _write_provenance(
                png_path,
                schema_kind="viz_minimizers",
                spec=spec,
                input_sha=extra_meta.get("input_sha256"),
                command=current_command_str(),
                viz_spec_path=viz_path,
            )
        elif name == "seed-chain":
            _, spec = plot_seed_chain(
                ref_length=data["ref_length"],
                qry_length=data["qry_length"],
                chains=data["chains"],
                save=str(png_path),
                save_viz_spec=str(viz_path),
                extra_meta=extra_meta,
            )
            _write_provenance(
                png_path,
                schema_kind="viz_seed_chain",
                spec=spec,
                input_sha=extra_meta.get("input_sha256"),
                command=current_command_str(),
                viz_spec_path=viz_path,
            )
        elif name == "rna-dotplot":
            _, spec = plot_rna_dotplot(
                posterior=data["posterior"],
                save=str(png_path),
                save_viz_spec=str(viz_path),
                extra_meta=extra_meta,
            )
            _write_provenance(
                png_path,
                schema_kind="viz_rna_dotplot",
                spec=spec,
                input_sha=extra_meta.get("input_sha256"),
                command=current_command_str(),
                viz_spec_path=viz_path,
            )
        elif name == "alignment-ribbon":
            _, spec = plot_alignment_ribbon(
                ref_length=data["ref_length"],
                qry_length=data["qry_length"],
                alignment=data,
                metadata=data.get("metadata"),
                save=str(png_path),
                save_viz_spec=str(viz_path),
                extra_meta=extra_meta,
            )
            _write_provenance(
                png_path,
                schema_kind="viz_alignment_ribbon",
                spec=spec,
                input_sha=extra_meta.get("input_sha256"),
                command=current_command_str(),
                viz_spec_path=viz_path,
            )
        elif name == "distance-heatmap":
            _, spec = plot_distance_heatmap(
                matrix=data["matrix"],
                labels=data["labels"],
                method=data.get("method", "demo"),
                save=str(png_path),
                save_viz_spec=str(viz_path),
                extra_meta=extra_meta,
            )
            _write_provenance(
                png_path,
                schema_kind="viz_distance_heatmap",
                spec=spec,
                input_sha=extra_meta.get("input_sha256"),
                command=current_command_str(),
                viz_spec_path=viz_path,
            )
        elif name == "motif-logo":
            _, spec = plot_motif_logo(
                pwm=data["pwm"],
                alphabet=data.get("alphabet"),
                background=data.get("background"),
                save=str(png_path),
                save_viz_spec=str(viz_path),
                extra_meta=extra_meta,
            )
            _write_provenance(
                png_path,
                schema_kind="viz_motif_logo",
                spec=spec,
                input_sha=extra_meta.get("input_sha256"),
                command=current_command_str(),
                viz_spec_path=viz_path,
            )
    print(f"Demo visualizations written to {output_dir}")


def command_live_run(args: argparse.Namespace) -> None:
    random.seed(args.seed)
    slice_path = None
    if getattr(args, "slice", None):
        try:
            slice_path = resolve_slice(args.slice)
        except KeyError as exc:
            raise SystemExit(str(exc)) from exc
    target = args.config or slice_path or args.hgx or args.target
    if target is None:
        raise SystemExit("Provide a .hgx model (--hgx) or an orchestrator config path.")
    target_path = Path(target)
    suffix = target_path.suffix.lower()
    if suffix in {".yaml", ".yml"}:
        _run_config_live(args, target_path)
        return
    _run_hgx_live(args, target_path)


def _run_hgx_live(args: argparse.Namespace, hgx_path: Path) -> None:
    model = load_hgx(hgx_path)
    graph = build_graph_from_hgx(model)
    dt_overrides = _parse_dt_overrides(args.dt_override)
    const_inputs = _parse_node_inputs(args.input)
    series_inputs = _load_input_series(args.inputs_json)
    components = graph.strongly_connected_components()
    if not components:
        components = [(node.name,) for node in graph.iter_nodes()]
    islands: List[Island] = []
    for idx, comp in enumerate(components):
        nodes = [graph.node(name) for name in comp]
        node_dts = [dt_overrides.get(name, args.default_dt) for name in comp]
        dt = min(node_dts) if node_dts else args.default_dt
        island = Island(name=f"island_{idx}", nodes=nodes, dt=dt)
        islands.append(island)

    realtime_hook = maybe_create_realtime_queue(getattr(args, "realtime", False), getattr(args, "realtime_endpoint", None))
    realtime_queue = realtime_hook.queue if realtime_hook else None
    node_meta = describe_graph_nodes(graph)
    reducer = StateReducer(hz=args.hz, realtime_queue=realtime_queue, node_meta=node_meta)
    event_bus = EventBus()
    hot_swap_manager = HotSwapManager()
    external_inputs = _build_input_provider(const_inputs, series_inputs)
    scheduler = LiveScheduler(
        graph=graph,
        islands=islands,
        sync_dt=args.sync_dt,
        state_reducer=reducer,
        external_inputs=external_inputs,
        event_bus=event_bus,
        hot_swap_manager=hot_swap_manager,
    )
    run_id = compose_run_id(model.name or "model", getattr(args, "variant", None), args.seed)
    scheduler.runtime_meta = {
        "model": model.name,
        "slice": getattr(args, "slice", None),
        "variant": args.variant,
        "run_id": run_id,
        "variants": [args.variant] if args.variant else [],
    }
    if realtime_hook and realtime_hook.server:
        realtime_hook.server.set_command_handler(_make_realtime_command_handler(scheduler))
    duration = args.duration if args.duration is not None else 1.0
    hz = getattr(reducer, "hz", args.hz)
    input_series_file = str(args.inputs_json) if args.inputs_json else None
    finalize_live_run(
        args,
        scheduler=scheduler,
        spec_path=hgx_path,
        model_name=model.name,
        duration=duration,
        sync_dt=args.sync_dt,
        dt_overrides=dt_overrides,
        const_inputs=const_inputs,
        input_series_file=input_series_file,
        hz=hz,
        reducer=reducer,
        event_bus=event_bus,
        realtime_hook=realtime_hook,
        meta_extra={"target_type": "hgx"},
    )


def _run_config_live(args: argparse.Namespace, config_path: Path) -> None:
    try:
        from helixtasks.egfr_grb2_orchestrator import build_from_config
    except ImportError as exc:  # pragma: no cover - optional repo scripts
        raise SystemExit(
            "EGFR/GRB2 orchestrator is unavailable. Ensure `helixtasks` is importable."
        ) from exc

    variant = args.variant or "wt"
    graph, scheduler, cfg = build_from_config(config_path, variant)
    sim_cfg = cfg.get("sim", {})
    duration = args.duration if args.duration is not None else float(sim_cfg.get("duration", 60.0))
    sync_dt = float(sim_cfg.get("sync_dt", args.sync_dt))
    realtime_hook = maybe_create_realtime_queue(getattr(args, "realtime", False), getattr(args, "realtime_endpoint", None))
    realtime_queue = realtime_hook.queue if realtime_hook else None
    node_meta = describe_graph_nodes(graph)
    reducer = StateReducer(hz=args.hz, realtime_queue=realtime_queue, node_meta=node_meta)
    scheduler.state_reducer = reducer
    hz = getattr(reducer, "hz", args.hz)
    model_base = cfg.get("name", config_path.stem)
    model_name = f"{model_base}-{variant}"
    run_id = compose_run_id(model_name, variant, args.seed)
    meta_extra = {
        "target_type": "config",
        "config_path": str(config_path),
        "variant": variant,
        "model_label": cfg.get("name"),
        "description": cfg.get("description"),
        "slice": getattr(args, "slice", None),
        "run_id": run_id,
        "variants": sorted(cfg.get("variants", {}).keys()),
    }
    scheduler.runtime_meta = dict(meta_extra)
    if realtime_hook and realtime_hook.server:
        realtime_hook.server.set_command_handler(_make_realtime_command_handler(scheduler))
    finalize_live_run(
        args,
        scheduler=scheduler,
        spec_path=config_path,
        model_name=model_name,
        duration=duration,
        sync_dt=sync_dt,
        dt_overrides={},
        const_inputs={},
        input_series_file=None,
        hz=hz,
        reducer=reducer,
        event_bus=None,
        realtime_hook=realtime_hook,
        meta_extra=meta_extra,
    )


def _format_kv_pairs(values: Mapping[str, float], limit: int) -> str:
    items = list(values.items())[:limit]
    return ", ".join(f"{k}={v:.4g}" for k, v in items)




def command_live_dev(args: argparse.Namespace) -> None:
    if args.hgx:
        session = LiveSession.from_hgx(Path(args.hgx), sync_dt=args.sync_dt, default_dt=args.default_dt)
    else:
        session = LiveSession(name=args.name, sync_dt=args.sync_dt, default_dt=args.default_dt)
    shell = LiveDevShell(session)
    shell.loop()


def command_live_inspect(args: argparse.Namespace) -> None:
    bundle = Path(args.bundle)
    meta = _read_json(bundle / "meta.json")
    snapshots = _read_json_optional(bundle / "snapshots.json") or []
    deltas = _read_json_optional(bundle / "deltas.json") or []
    events = _read_json_optional(bundle / "events.json") or []

    model = meta.get("model", "<unknown>")
    duration = meta.get("duration")
    sync_dt = meta.get("sync_dt")
    seed = meta.get("seed")
    islands = meta.get("islands", [])

    model_path = bundle / "model.hgx"
    model_hash = _file_sha256(model_path) if model_path.exists() else None
    graph = None
    graph_nodes = {}
    node_kind_map: Dict[str, str] = {}
    node_edges = 0
    model_name = meta.get("model_label") or model
    if model_path.exists() and meta.get("target_type") != "config":
        try:
            hgx_model = load_hgx(model_path)
            model_name = hgx_model.name or model
            graph = build_graph_from_hgx(hgx_model)
            graph_nodes = graph.nodes
            node_kind_map = {name: node.kind for name, node in graph.nodes.items()}
            node_edges = len(graph.edges)
        except Exception:
            graph = None

    island_lookup: Dict[str, str] = {}
    for island in islands:
        for node_name in island.get("nodes", []):
            island_lookup[node_name] = island.get("name")

    class _MetricAggregator:
        def __init__(self) -> None:
            self.count = 0
            self.mean = 0.0
            self.m2 = 0.0
            self.min = None
            self.max = None

        def update(self, value: float) -> None:
            if self.count == 0:
                self.mean = float(value)
                self.min = float(value)
                self.max = float(value)
                self.count = 1
                self.m2 = 0.0
                return
            self.count += 1
            delta = value - self.mean
            self.mean += delta / self.count
            delta2 = value - self.mean
            self.m2 += delta * delta2
            self.min = min(self.min, value) if self.min is not None else float(value)
            self.max = max(self.max, value) if self.max is not None else float(value)

        def as_dict(self) -> Dict[str, float]:
            var = self.m2 / (self.count - 1) if self.count > 1 else 0.0
            std = var ** 0.5
            return {
                "count": self.count,
                "min": float(self.min) if self.min is not None else None,
                "max": float(self.max) if self.max is not None else None,
                "mean": float(self.mean),
                "var": float(var),
                "std": float(std),
            }

    collect_metrics = not args.no_metrics
    node_metric_map: Dict[str, Dict[str, _MetricAggregator]] = {}
    if collect_metrics:
        for snapshot in snapshots:
            for node_name, ports in snapshot.items():
                for port_name, value in ports.items():
                    if not isinstance(value, (int, float)):
                        continue
                    node_metric_map.setdefault(node_name, {}).setdefault(port_name, _MetricAggregator()).update(
                        float(value)
                    )

    metric_filters = set(args.metric or [])
    node_summaries: Dict[str, Any] = {}
    node_iterable = node_kind_map.keys() if node_kind_map else node_metric_map.keys()
    for node_name in node_iterable:
        metrics = node_metric_map.get(node_name, {})
        filtered_metrics = {}
        if collect_metrics:
            if metric_filters:
                for port, agg in metrics.items():
                    if port in metric_filters:
                        filtered_metrics[port] = agg.as_dict()
            else:
                filtered_metrics = {port: agg.as_dict() for port, agg in metrics.items()}
        node_entry: Dict[str, Any] = {
            "kind": node_kind_map.get(node_name),
            "island": island_lookup.get(node_name),
            "metrics": filtered_metrics if collect_metrics else {},
        }
        node_summaries[node_name] = node_entry

    summary = {
        "schema_version": 1,
        "bundle": str(bundle),
        "model": {
            "name": model_name,
            "path": str(model_path) if model_path.exists() else None,
            "hash": model_hash,
            "nodes": len(graph_nodes) if graph_nodes else len(node_summaries),
            "edges": node_edges,
        },
        "runtime": {
            "duration": duration,
            "sync_dt": sync_dt,
            "seed": seed,
            "islands": len(islands),
            "snapshots": len(snapshots),
            "events": len(events),
            "wall_time_sec": meta.get("wall_time_sec"),
        },
        "nodes": node_summaries,
    }

    if args.json:
        dump_json(args.json, summary)
        return

    print(f"Bundle: {summary['bundle']}")
    print(
        f"  Model: {summary['model']['name']}  Duration: {summary['runtime']['duration']}  "
        f"sync_dt: {summary['runtime']['sync_dt']}  snapshots: {summary['runtime']['snapshots']}"
    )
    wall = summary["runtime"].get("wall_time_sec")
    wall_str = f"  Wall: {wall:.3f}s" if isinstance(wall, (int, float)) else ""
    print(
        f"  Seed: {summary['runtime']['seed']}  Islands: {summary['runtime']['islands']}  "
        f"Events: {summary['runtime']['events']}{wall_str}"
    )
    if islands:
        for island in islands[: args.max_islands]:
            nodes = ", ".join(island.get("nodes", [])[: args.max_nodes])
            if len(island.get("nodes", [])) > args.max_nodes:
                nodes += ", …"
            print(f"    - {island['name']} dt={island.get('dt')} nodes: {nodes}")
    if snapshots:
        first = snapshots[0]
        last = snapshots[-1]
        print("  Snapshot[0]:")
        for node, values in list(first.items())[: args.head]:
            print(f"    {node}: {_format_kv_pairs(values, args.metrics)}")
        if len(first) > args.head:
            print("    …")
        print("  Snapshot[-1]:")
        for node, values in list(last.items())[: args.head]:
            print(f"    {node}: {_format_kv_pairs(values, args.metrics)}")
        if len(last) > args.head:
            print("    …")
    if deltas:
        added = sum(len(delta.get("added", {})) for delta in deltas)
        removed = sum(len(delta.get("removed", {})) for delta in deltas)
        updated = sum(len(delta.get("updated", {})) for delta in deltas)
        print(f"  Deltas: added={added} updated={updated} removed={removed} (across {len(deltas)} samples)")


def command_live_viz(args: argparse.Namespace) -> None:
    _ensure_realtime_available()
    bundle = Path(args.bundle) if args.bundle else None
    endpoint = args.endpoint if not bundle else None
    from .live.viz.app import LiveVizApp  # noqa: WPS433

    app = LiveVizApp(endpoint=endpoint, bundle=bundle, target_hz=args.hz)
    app.run()


def command_live_demo(args: argparse.Namespace) -> None:
    _ensure_realtime_available()
    endpoint = args.realtime_endpoint or "tcp://127.0.0.1:8765"
    hook = maybe_create_realtime_queue(True, endpoint)
    if not hook:
        raise SystemExit("Unable to start realtime server.")
    rig = _RealtimeDemo(hz=args.hz)
    if hook.server:
        hook.server.set_command_handler(rig.handle_control)
    duration = args.duration or 0.0
    try:
        rig.run(hook.queue, duration=duration)
    except KeyboardInterrupt:
        pass
    finally:
        hook.queue.put(None)
        if hook.server:
            hook.server.close()


class _RealtimeDemo:
    """Synthetic EGF + agent streamer for testing the realtime viz."""

    def __init__(self, *, field_size: Tuple[int, int] = (128, 128), hz: float = 30.0, agents: int = 200):
        self.nx, self.ny = field_size
        self.hz = hz
        self.dt = 1.0 / max(hz, 1.0)
        self.t = 0.0
        self.paused = False
        self.grb2_scale = 1.0
        self.variant = "wt"
        self.variants = ["wt", "ko", "hig"]
        self.rng = np.random.default_rng(17)
        self._agents: Dict[int, Dict[str, float]] = {
            idx: {
                "id": idx,
                "x": float(self.rng.uniform(0.0, self.nx)),
                "y": float(self.rng.uniform(0.0, self.ny)),
                "perk": float(self.rng.uniform(0.0, 1.0)),
            }
            for idx in range(agents)
        }
        self._full_sent = False

    def run(self, queue: "queue.Queue", *, duration: float) -> None:
        deadline = time.time() + duration if duration > 0 else None
        try:
            while deadline is None or time.time() < deadline:
                start = time.time()
                delta = self._build_delta()
                queue.put(delta)
                sleep_for = max(0.0, self.dt - (time.time() - start))
                time.sleep(sleep_for)
        except KeyboardInterrupt:
            raise

    def _build_delta(self) -> Dict[str, Any]:
        if not self.paused:
            self._step_agents()
            self.t += self.dt
        field = self._make_field(self.t)
        tile = self._encode_tile(field)
        agents = self._agents_payload()
        metrics = {
            "grb2_scale": self.grb2_scale,
            "cell_count": len(self._agents),
            "time": self.t,
        }
        return {
            "schema_version": 1,
            "kind": "live_delta",
            "slice": "egfr_demo",
            "variant": self.variant,
            "run_id": "egfr_demo",
            "t": self.t,
            "dt": 0.0 if self.paused else self.dt,
            "metrics": metrics,
            "fields": {
                "egf": {
                    "id": 0,
                    "nx": self.nx,
                    "ny": self.ny,
                    "min": float(field.min()),
                    "max": float(field.max()),
                    "dtype": "f32",
                    "tiles": [tile],
                }
            },
            "agents": agents,
        }

    def _make_field(self, t: float) -> np.ndarray:
        x = np.linspace(0.0, 1.0, self.nx, dtype="float32")[None, :]
        y = np.linspace(0.0, 1.0, self.ny, dtype="float32")[:, None]
        wave = np.sin(0.5 * t)
        gradient = np.clip(x + y + 0.2 * wave * self.grb2_scale, 0.0, 2.0)
        gradient /= gradient.max() if gradient.max() else 1.0
        return gradient.astype("float32")

    def _encode_tile(self, arr: np.ndarray) -> Dict[str, Any]:
        contiguous = np.ascontiguousarray(arr)
        return {
            "x0": 0,
            "y0": 0,
            "nx": contiguous.shape[1],
            "ny": contiguous.shape[0],
            "encoding": "raw_f32_le",
            "data": base64.b64encode(contiguous.tobytes()).decode("ascii"),
        }

    def _agents_payload(self) -> Dict[str, Any]:
        if not self._agents:
            return {"full": True, "add": [], "update": [], "remove": []}
        if not self._full_sent:
            add = [dict(agent) for agent in self._agents.values()]
            updates: List[Dict[str, float]] = []
            self._full_sent = True
            full = True
        else:
            add = []
            updates = [
                {"id": agent_id, "x": agent["x"], "y": agent["y"], "perk": agent["perk"]}
                for agent_id, agent in self._agents.items()
            ]
            full = False
        return {
            "full": full,
            "add": add,
            "update": updates,
            "remove": [],
        }

    def _step_agents(self) -> None:
        for agent in self._agents.values():
            agent["x"] = float((agent["x"] + self.rng.normal(scale=1.5)) % self.nx)
            agent["y"] = float((agent["y"] + self.rng.normal(scale=1.5)) % self.ny)
            agent["perk"] = float(np.clip(agent["perk"] + self.rng.normal(scale=0.04), 0.0, 1.0))

    def handle_control(self, payload: Dict[str, Any]) -> None:
        if payload.get("kind") == "live_control":
            control_type = payload.get("type")
            if control_type == "set_param" and payload.get("param") == "grb2_scale":
                value = payload.get("value")
                if isinstance(value, (int, float)):
                    self.grb2_scale = float(value)
            elif control_type == "set_variant":
                variant = payload.get("variant")
                if isinstance(variant, str) and variant in self.variants:
                    self.variant = variant
            elif control_type == "pause":
                self.paused = True
            elif control_type == "resume":
                self.paused = False
            elif control_type == "toggle":
                self.paused = not self.paused
        else:
            command = payload.get("command")
            if command == "set_egf":
                value = payload.get("value")
                if isinstance(value, (int, float)):
                    self.grb2_scale = float(value)
            elif command == "pause":
                self.paused = True
            elif command == "resume":
                self.paused = False
            elif command == "toggle":
                self.paused = not self.paused


def command_live_plan(args: argparse.Namespace) -> None:
    slice_path = None
    if getattr(args, "slice", None):
        try:
            slice_path = resolve_slice(args.slice)
        except KeyError as exc:
            raise SystemExit(str(exc)) from exc
    target = args.config or slice_path or args.hgx or args.target
    if target is None:
        raise SystemExit("Provide a .hgx model (--hgx) or an orchestrator config path.")
    path = Path(target)
    if path.suffix.lower() in {".yaml", ".yml"} and not args.hgx:
        _plan_config_live(args, path)
    else:
        _plan_hgx_live(args, path)


def _plan_hgx_live(args: argparse.Namespace, hgx_path: Path) -> None:
    model = load_hgx(hgx_path)
    graph = build_graph_from_hgx(model)
    dt_overrides = _parse_dt_overrides(getattr(args, "dt_override", []))
    default_dt = args.default_dt
    sync_dt = args.sync_dt
    components = graph.strongly_connected_components()
    if not components:
        components = [(node.name,) for node in graph.iter_nodes()]
    islands = []
    for idx, comp in enumerate(components):
        comp_dt = min(dt_overrides.get(name, default_dt) for name in comp)
        islands.append(
            {
                "name": f"island_{idx}",
                "dt": comp_dt,
                "nodes": list(comp),
            }
        )
    print_plan_summary(
        title=f"Live plan for {hgx_path}",
        model_name=model.name or hgx_path.stem,
        sync_dt=sync_dt,
        default_dt=default_dt,
        overrides=dt_overrides,
        islands=islands,
        node_count=len(list(graph.iter_nodes())),
        edge_count=len(graph.edges),
    )


def _plan_config_live(args: argparse.Namespace, config_path: Path) -> None:
    try:
        from helixtasks.egfr_grb2_orchestrator import build_from_config
    except ImportError as exc:  # pragma: no cover - optional repo scripts
        raise SystemExit(
            "EGFR/GRB2 orchestrator is unavailable. Ensure `helixtasks` is importable."
        ) from exc
    variant = args.variant or "wt"
    graph, scheduler, cfg = build_from_config(config_path, variant)
    sim_cfg = cfg.get("sim", {})
    sync_dt = float(sim_cfg.get("sync_dt", args.sync_dt))
    default_dt = float(sim_cfg.get("default_dt", args.default_dt))
    islands = []
    for idx, island in enumerate(getattr(scheduler, "islands", [])):
        nodes = [getattr(node, "name", "") for node in getattr(island, "nodes", []) if getattr(node, "name", None)]
        islands.append(
            {
                "name": getattr(island, "name", f"island_{idx}"),
                "dt": getattr(island, "dt", None),
                "nodes": nodes,
            }
        )
    print_plan_summary(
        title=f"Live plan for {config_path}",
        model_name=cfg.get("name") or config_path.stem,
        sync_dt=sync_dt,
        default_dt=default_dt,
        overrides={},
        islands=islands,
        node_count=len(graph.nodes),
        edge_count=len(graph.edges),
    )


def build_parser() -> argparse.ArgumentParser:
    description = (
        "Helix unified CLI for computational bioinformatics workflows.\n\n"
        "Simulation only: Helix operates purely in silico on digital sequences/datasets and does not "
        "control lab equipment or prescribe wet-lab procedures."
    )
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    dna = subparsers.add_parser("dna", help="Summarize GC, windows, and k-mer hotspots.")
    dna.add_argument("--sequence", help="Inline DNA string (defaults to the bundled sample).")
    dna.add_argument("--input", type=Path, help="Path to a FASTA/text file.")
    dna.add_argument("--window", type=int, default=200, help="GC content window size (default: 200).")
    dna.add_argument("--step", type=int, default=50, help="GC window stride (default: 50).")
    dna.add_argument("--k", type=int, default=5, help="k-mer size (default: 5).")
    dna.add_argument("--max-diff", type=int, default=1, help="Maximum mismatches when clustering (default: 1).")
    dna.add_argument("--top", type=int, default=10, help="Print this many top clusters (default: 10).")
    dna.set_defaults(func=command_dna)

    spectrum = subparsers.add_parser("spectrum", help="Compute theoretical spectra or run leaderboard sequencing.")
    spectrum.add_argument("--peptide", help="Peptide sequence to analyse.")
    spectrum.add_argument("--linear", action="store_true", help="Use the linear spectrum instead of cyclic.")
    spectrum.add_argument("--spectrum", help="Comma/space-separated experimental masses.")
    spectrum.add_argument("--spectrum-file", type=Path, help="File containing experimental masses.")
    spectrum.add_argument("--leaderboard", type=int, default=5, help="Leaderboard size (default: 5).")
    spectrum.set_defaults(func=command_spectrum)

    rna = subparsers.add_parser("rna", help="RNA folding + ensemble helpers.")
    rna_sub = rna.add_subparsers(dest="rna_command", required=True)

    rna_mfe = rna_sub.add_parser("mfe", help="Zuker-style MFE folding.")
    rna_mfe.add_argument("--fasta", type=Path, required=True, help="FASTA file containing a single sequence.")
    rna_mfe.add_argument("--json", type=Path, help="Optional JSON output path.")
    rna_mfe.add_argument("--dotbracket", type=Path, help="Optional file to save dot-bracket.")
    rna_mfe.set_defaults(func=command_rna_mfe)

    rna_ensemble = rna_sub.add_parser("ensemble", help="Partition function + MEA/centroid structures.")
    rna_ensemble.add_argument("--fasta", type=Path, required=True, help="FASTA file containing a single sequence.")
    rna_ensemble.add_argument("--gamma", type=float, default=1.0, help="MEA gamma parameter (default: 1.0).")
    rna_ensemble.add_argument("--json", type=Path, help="Optional JSON output path.")
    rna_ensemble.add_argument("--dotplot", type=Path, help="Optional dot-plot path (requires matplotlib).")
    rna_ensemble.add_argument("--arc", type=Path, help="Optional arc diagram path (requires matplotlib).")
    rna_ensemble.add_argument("--entropy-plot", type=Path, help="Optional entropy plot path (requires matplotlib).")
    rna_ensemble.add_argument("--save-viz-spec", type=Path, help="Optional viz-spec JSON for dot-plot.")
    rna_ensemble.set_defaults(func=command_rna_ensemble)

    live = subparsers.add_parser("live", help="Run LiveGraph HGX models.")
    live_sub = live.add_subparsers(dest="live_command", required=True)
    live_run = live_sub.add_parser("run", help="Execute a model using the LiveScheduler.")
    live_run.add_argument(
        "target",
        nargs="?",
        type=Path,
        help="Path to a .hgx LiveGraph model or an orchestrator config YAML.",
    )
    live_run.add_argument("--hgx", type=Path, help="Path to the .hgx model (legacy flag).")
    live_run.add_argument("--config", type=Path, help="Path to an orchestrator config YAML.")
    live_run.add_argument(
        "--slice",
        type=str,
        help=f"Named slice declared in the registry (available: {', '.join(sorted(SLICE_REGISTRY))}).",
    )
    live_run.add_argument("--variant", type=str, help="Variant name defined by the orchestrator config.")
    live_run.add_argument("--duration", type=float, default=None, help="Simulation time horizon in minutes.")
    live_run.add_argument("--sync-dt", type=float, default=0.25, help="Synchronization cadence Δt (default: 0.25).")
    live_run.add_argument("--default-dt", type=float, default=0.1, help="Default island Δt (default: 0.1).")
    live_run.add_argument(
        "--dt",
        dest="dt_override",
        action="append",
        default=[],
        help="Override per-node dt (format: node=0.05). Repeat for multiple nodes.",
    )
    live_run.add_argument(
        "--input",
        action="append",
        default=[],
        help="Inject constant signals (format: node.port=value). Repeat as needed.",
    )
    live_run.add_argument(
        "--inputs-json",
        type=Path,
        help="Path to a JSON file describing constant/time-series inputs.",
    )
    live_run.add_argument("--hz", type=float, default=60.0, help="StateReducer target Hz (default: 60).")
    live_run.add_argument("--wall-ms", type=float, default=5.0, help="Wall-clock pacing budget in ms (default: 5).")
    live_run.add_argument("--seed", type=int, default=0, help="Deterministic seed (default: 0).")
    live_run.add_argument(
        "--bundle",
        type=Path,
        help="Directory for the COMBINE-style bundle (default: live_runs/<model>-<timestamp>).",
    )
    live_run.add_argument(
        "--realtime",
        action="store_true",
        help="Stream state deltas to a realtime queue for visualization clients.",
    )
    live_run.add_argument(
        "--realtime-endpoint",
        type=str,
        help="Endpoint for realtime clients (tcp://host:port or ipc://path). Default: tcp://127.0.0.1:8765.",
    )
    live_run.set_defaults(func=command_live_run)
    live_inspect = live_sub.add_parser("inspect", help="Inspect a LiveGraph bundle.")
    live_inspect.add_argument("--bundle", type=Path, required=True, help="Bundle directory to inspect.")
    live_inspect.add_argument("--head", type=int, default=3, help="How many nodes to show per snapshot (default: 3).")
    live_inspect.add_argument(
        "--metrics", type=int, default=3, help="How many port values to show per node (default: 3)."
    )
    live_inspect.add_argument(
        "--max-islands", type=int, default=4, help="Maximum number of islands to list (default: 4)."
    )
    live_inspect.add_argument(
        "--max-nodes", type=int, default=5, help="Maximum number of node names to show per island (default: 5)."
    )
    live_inspect.add_argument(
        "--no-metrics",
        action="store_true",
        help="Skip aggregating per-node metrics (useful for huge bundles). Metrics will be empty objects.",
    )
    live_inspect.add_argument(
        "--metric",
        action="append",
        help="Restrict metrics to these keys (repeatable). Useful for CI filters.",
    )
    live_inspect.add_argument("--json", type=Path, help="Optional JSON path for structured summaries.")
    live_inspect.set_defaults(func=command_live_inspect)
    live_viz = live_sub.add_parser("viz", help="Launch the realtime renderer.")
    live_viz.add_argument(
        "--endpoint",
        type=str,
        default="tcp://127.0.0.1:8765",
        help="Realtime endpoint to attach (default: tcp://127.0.0.1:8765).",
    )
    live_viz.add_argument(
        "--bundle",
        type=Path,
        help="Replay an existing bundle instead of attaching to a live endpoint.",
    )
    live_viz.add_argument("--hz", type=float, default=60.0, help="Target redraw rate (default: 60 Hz).")
    live_viz.set_defaults(func=command_live_viz)
    live_demo = live_sub.add_parser("demo", help="Stream a synthetic EGF field + agents for testing the viz stack.")
    live_demo.add_argument(
        "--duration",
        type=float,
        default=60.0,
        help="How long to stream (seconds). Use 0 for an endless loop.",
    )
    live_demo.add_argument("--hz", type=float, default=30.0, help="Demo frame rate (default: 30 Hz).")
    live_demo.add_argument(
        "--realtime-endpoint",
        type=str,
        help="Endpoint for the realtime server (default: tcp://127.0.0.1:8765).",
    )
    live_demo.set_defaults(func=command_live_demo)
    live_plan = live_sub.add_parser("plan", help="Inspect the island layout/dt table without running the model.")
    live_plan.add_argument(
        "target",
        nargs="?",
        type=Path,
        help="Path to a .hgx LiveGraph model or an orchestrator config YAML.",
    )
    live_plan.add_argument("--hgx", type=Path, help="Path to the .hgx model (legacy flag).")
    live_plan.add_argument("--config", type=Path, help="Path to an orchestrator config YAML.")
    live_plan.add_argument(
        "--slice",
        type=str,
        help=f"Named slice declared in the registry (available: {', '.join(sorted(SLICE_REGISTRY))}).",
    )
    live_plan.add_argument("--variant", type=str, help="Variant name defined by the orchestrator config.")
    live_plan.add_argument("--sync-dt", type=float, default=0.25, help="Synchronization cadence Δt (default: 0.25).")
    live_plan.add_argument("--default-dt", type=float, default=0.1, help="Default island Δt (default: 0.1).")
    live_plan.add_argument(
        "--dt",
        dest="dt_override",
        action="append",
        default=[],
        help="Override per-node dt (format: node=0.05). Repeat for multiple nodes.",
    )
    live_plan.set_defaults(func=command_live_plan)
    live_dev = live_sub.add_parser("dev", help="Interactive LiveLab session (REPL).")
    live_dev.add_argument(
        "--hgx",
        type=Path,
        help="Optional HGX graph to load as the starting point.",
    )
    live_dev.add_argument("--name", type=str, default="LiveSession", help="Session name (default: LiveSession).")
    live_dev.add_argument("--sync-dt", type=float, default=0.25, help="Synchronization Δt (default: 0.25).")
    live_dev.add_argument("--default-dt", type=float, default=0.1, help="Default island Δt (default: 0.1).")
    live_dev.set_defaults(func=command_live_dev)

    crispr = subparsers.add_parser("crispr", help="CRISPR design helpers.")
    crispr_sub = crispr.add_subparsers(dest="crispr_command", required=True)

    crispr_find = crispr_sub.add_parser("find-guides", help="Discover CRISPR guide candidates near PAM matches.")
    crispr_find.add_argument("--sequence", help="Inline target sequence.")
    crispr_find.add_argument("--fasta", type=Path, help="Path to a FASTA file containing the target sequence.")
    crispr_find.add_argument("--pam", default="SpCas9-NGG", help="Named PAM definition (default: SpCas9-NGG).")
    crispr_find.add_argument("--guide-len", type=int, default=20, help="Guide length (default: 20).")
    crispr_find.add_argument(
        "--strand",
        choices=["+", "-", "both"],
        default="both",
        help="Strand(s) to search (default: both).",
    )
    crispr_find.add_argument(
        "--window",
        type=int,
        nargs=2,
        metavar=("START", "END"),
        help="Optional 0-based window [START, END) to limit guide discovery.",
    )
    crispr_find.add_argument("--json", type=Path, help="Optional path to write the crispr.guides JSON payload.")
    crispr_find.add_argument(
        "--emit-sequences",
        action="store_true",
        help="Include guide sequences in the JSON output (default masks sequences).",
    )
    crispr_find.set_defaults(func=command_crispr_find_guides)

    crispr_off = crispr_sub.add_parser("offtargets", help="Enumerate CRISPR off-target sites for each guide.")
    crispr_off.add_argument("--genome", help="Inline genomic sequence to search.")
    crispr_off.add_argument("--fasta", type=Path, help="FASTA file containing the genome/contig to search.")
    crispr_off.add_argument("--guides", type=Path, required=True, help="Path to a crispr.guides JSON file.")
    crispr_off.add_argument("--pam", default="SpCas9-NGG", help="Named PAM definition (default: SpCas9-NGG).")
    crispr_off.add_argument("--max-mm", type=int, default=3, help="Maximum mismatches allowed (default: 3).")
    crispr_off.add_argument(
        "--max-gap",
        type=int,
        default=0,
        help="Maximum gaps allowed (default: 0; currently gaps>0 unsupported).",
    )
    crispr_off.add_argument("--json", type=Path, help="Optional path to write crispr.offtargets JSON.")
    crispr_off.set_defaults(func=command_crispr_offtargets)

    crispr_score = crispr_sub.add_parser("score", help="Apply on/off-target scoring using a weights plugin.")
    crispr_score.add_argument("--guides", type=Path, required=True, help="Path to crispr.guides JSON.")
    crispr_score.add_argument("--hits", type=Path, required=True, help="Path to crispr.offtargets JSON.")
    crispr_score.add_argument(
        "--weights",
        type=Path,
        help="Optional JSON file describing scoring weights (position/mismatch penalties).",
    )
    crispr_score.add_argument("--json", type=Path, help="Optional output JSON path for scored payload.")
    crispr_score.set_defaults(func=command_crispr_score)

    crispr_sim = crispr_sub.add_parser("simulate", help="Simulate CRISPR cut/repair outcomes.")
    crispr_sim.add_argument("--sequence", help="Inline site sequence.")
    crispr_sim.add_argument("--fasta", type=Path, help="FASTA file containing the site sequence.")
    crispr_sim.add_argument("--guides", type=Path, required=True, help="Path to crispr.guides JSON.")
    crispr_sim.add_argument("--guide-id", required=True, help="Guide identifier to simulate.")
    crispr_sim.add_argument("--priors", type=Path, help="JSON file describing outcome priors/weights.")
    crispr_sim.add_argument("--draws", type=int, default=1000, help="Number of Monte Carlo draws (default: 1000).")
    crispr_sim.add_argument("--seed", type=int, help="Optional RNG seed for reproducibility.")
    crispr_sim.add_argument("--json", type=Path, help="Optional path to write crispr.sim JSON.")
    crispr_sim.add_argument(
        "--emit-sequences",
        action="store_true",
        help="Include raw site/guide sequences in the JSON (masked by default).",
    )
    crispr_sim.add_argument(
        "--viz-spec",
        type=Path,
        help="Optional path to write a ModernGL viz spec for helix viz modern.",
    )
    crispr_sim.set_defaults(func=command_crispr_simulate)

    crispr_genome = crispr_sub.add_parser("genome-sim", help="Simulate digital cut events across a genome.")
    crispr_genome.add_argument("--genome", type=Path, required=True, help="Genome FASTA file to scan.")
    crispr_genome.add_argument("--guide-sequence", required=True, help="Guide RNA sequence (5'->3').")
    crispr_genome.add_argument("--guide-name", help="Optional guide identifier.")
    crispr_genome.add_argument("--guide-pam", help="Optional PAM label stored with the guide metadata.")
    crispr_genome.add_argument(
        "--cas",
        default="cas9",
        help=f"Cas preset to use (options: {', '.join(sorted(CAS_PRESET_CONFIGS))}; default: cas9).",
    )
    crispr_genome.add_argument("--cas-config", type=Path, help="JSON file describing a custom CasSystem.")
    crispr_genome.add_argument("--max-events", type=int, help="Limit the number of events returned.")
    crispr_genome.add_argument("--json", type=Path, help="Optional output JSON path (defaults to stdout).")
    crispr_genome.set_defaults(func=command_crispr_genome_sim)

    crispr_dag = crispr_sub.add_parser("dag", help="Construct a CRISPR edit DAG.")
    crispr_dag.add_argument("--genome", type=Path, required=True, help="Genome FASTA file.")
    crispr_dag.add_argument("--guide-sequence", help="Guide RNA sequence (5'->3').")
    crispr_dag.add_argument("--guide-name", help="Optional guide identifier.")
    crispr_dag.add_argument("--guide-pam", help="Optional PAM label stored with the guide metadata.")
    crispr_dag.add_argument("--guides-file", type=Path, help="TSV/JSON file describing multiple guides.")
    crispr_dag.add_argument(
        "--region",
        help="Restrict simulations to a genomic slice (chrom:start-end). Overridden by guide-specific regions.",
    )
    crispr_dag.add_argument(
        "--cas",
        default="cas9",
        help=f"Cas preset to use (options: {', '.join(sorted(CAS_PRESET_CONFIGS))}; default: cas9).",
    )
    crispr_dag.add_argument("--cas-config", type=Path, help="JSON file describing a custom CasSystem.")
    crispr_dag.add_argument("--max-sites", type=int, default=50, help="Maximum candidate sites per guide (default: 50).")
    crispr_dag.add_argument("--max-depth", type=int, default=2, help="Maximum DAG depth (default: 2).")
    crispr_dag.add_argument("--min-prob", type=float, default=1e-4, help="Minimum branch probability (default: 1e-4).")
    crispr_dag.add_argument("--seed", type=int, default=0, help="Random seed for stochastic rules (default: 0).")
    crispr_dag.add_argument("--out", type=Path, help="Output path for single-guide simulations.")
    crispr_dag.add_argument("--out-dir", type=Path, help="Directory for per-guide DAG outputs.")
    crispr_dag.add_argument("--json", type=Path, help="Deprecated alias for --out (defaults to stdout).")
    crispr_dag.add_argument("--use-gpu", action="store_true", help="Use GPU-accelerated scoring (Numba/CUDA).")
    crispr_dag.add_argument(
        "--frames",
        type=str,
        help="Optional JSONL path to stream DAG frames (use '-' for stdout).",
    )
    crispr_dag.add_argument(
        "--coding-json",
        type=Path,
        help="Transcript JSON for protein-impact annotations (same schema as the GUI transcript loader).",
    )
    crispr_dag.add_argument(
        "--coding-transcript",
        help="Transcript ID or name to use from --coding-json (defaults to the first entry).",
    )
    crispr_dag.set_defaults(func=command_crispr_dag)

    prime_cmd = subparsers.add_parser("prime", help="Prime editing helpers.")
    prime_sub = prime_cmd.add_subparsers(dest="prime_command", required=True)

    prime_sim = prime_sub.add_parser("simulate", help="Simulate prime editing outcomes.")
    prime_sim.add_argument("--genome", type=Path, required=True, help="Genome FASTA file containing the target site(s).")
    prime_sim.add_argument("--peg-config", type=Path, help="JSON file describing pegRNA components.")
    prime_sim.add_argument("--peg-spacer", help="pegRNA spacer sequence (5'->3').")
    prime_sim.add_argument("--peg-pbs", help="pegRNA primer binding site sequence.")
    prime_sim.add_argument("--peg-rtt", help="pegRNA reverse transcription template sequence.")
    prime_sim.add_argument("--peg-name", help="Optional pegRNA identifier.")
    prime_sim.add_argument("--editor-config", type=Path, help="JSON file describing a PrimeEditor.")
    prime_sim.add_argument("--editor-name", help="Override the editor name.")
    prime_sim.add_argument(
        "--cas",
        default="cas9",
        help=f"Cas preset for inline editor definitions (options: {', '.join(sorted(CAS_PRESET_CONFIGS))}; default: cas9).",
    )
    prime_sim.add_argument("--cas-config", type=Path, help="JSON file describing a custom CasSystem for the editor.")
    prime_sim.add_argument("--nick-offset", type=int, help="Override nick-to-edit offset.")
    prime_sim.add_argument("--efficiency-scale", type=float, help="Override efficiency scale.")
    prime_sim.add_argument("--indel-bias", type=float, help="Override indel bias.")
    prime_sim.add_argument("--mismatch-tolerance", type=int, help="Override mismatch tolerance.")
    prime_sim.add_argument("--flap-balance", type=float, help="Relative weight for left vs right flap resolution (0-1).")
    prime_sim.add_argument("--reanneal-bias", type=float, help="Probability mass for reanneal/no-edit branches.")
    prime_sim.add_argument(
        "--max-outcomes",
        type=int,
        default=16,
        help="Deprecated; retained for compatibility. Use --draws to control sampling size.",
    )
    prime_sim.add_argument("--draws", type=int, default=1000, help="Number of Monte Carlo draws (default: 1000).")
    prime_sim.add_argument("--seed", type=int, help="Optional RNG seed.")
    prime_sim.add_argument(
        "--priors-profile",
        default="default_indel",
        help="Prime priors profile to use (default: default_indel).",
    )
    prime_sim.add_argument(
        "--pam-profile",
        default="SpCas9_NGG",
        help="Named PAM profile or preset used to build the PAM mask (default: SpCas9_NGG).",
    )
    prime_sim.add_argument(
        "--pam-softness",
        type=float,
        default=1.0,
        help="Penalty applied to off-PAM positions (0.0=disallow, 1.0=ignore PAM; default: 1.0).",
    )
    prime_sim.add_argument("--json", type=Path, help="Optional output JSON path (defaults to stdout).")
    prime_sim.add_argument(
        "--viz-spec",
        type=Path,
        help="Optional path to write a ModernGL viz spec for helix viz modern.",
    )
    prime_sim.set_defaults(func=command_prime_simulate)

    prime_dag = prime_sub.add_parser("dag", help="Construct a prime editing edit DAG.")
    prime_dag.add_argument("--genome", type=Path, required=True, help="Genome FASTA file.")
    prime_dag.add_argument("--peg-config", type=Path, help="JSON file describing pegRNA components.")
    prime_dag.add_argument("--peg-spacer", help="pegRNA spacer sequence (5'->3').")
    prime_dag.add_argument("--peg-pbs", help="pegRNA primer binding site sequence.")
    prime_dag.add_argument("--peg-rtt", help="pegRNA reverse transcription template sequence.")
    prime_dag.add_argument("--peg-name", help="Optional pegRNA identifier.")
    prime_dag.add_argument("--editor-config", type=Path, help="JSON file describing a PrimeEditor.")
    prime_dag.add_argument("--editor-name", help="Override the editor name.")
    prime_dag.add_argument(
        "--cas",
        default="cas9",
        help=f"Cas preset for inline editor definitions (options: {', '.join(sorted(CAS_PRESET_CONFIGS))}; default: cas9).",
    )
    prime_dag.add_argument("--cas-config", type=Path, help="JSON file describing a custom CasSystem for the editor.")
    prime_dag.add_argument("--nick-offset", type=int, help="Override nick-to-edit offset.")
    prime_dag.add_argument("--efficiency-scale", type=float, help="Override efficiency scale.")
    prime_dag.add_argument("--indel-bias", type=float, help="Override indel bias.")
    prime_dag.add_argument("--mismatch-tolerance", type=int, help="Override mismatch tolerance.")
    prime_dag.add_argument("--flap-balance", type=float, help="Relative weight for left vs right flap resolution (0-1).")
    prime_dag.add_argument("--reanneal-bias", type=float, help="Probability mass for reanneal/no-edit branches.")
    prime_dag.add_argument("--max-depth", type=int, default=2, help="Maximum DAG depth (default: 2).")
    prime_dag.add_argument("--min-prob", type=float, default=1e-4, help="Minimum branch probability (default: 1e-4).")
    prime_dag.add_argument("--seed", type=int, default=0, help="Random seed (default: 0).")
    prime_dag.add_argument("--pegs-file", type=Path, help="TSV/JSON file describing multiple pegRNAs.")
    prime_dag.add_argument("--region", help="Restrict simulations to chrom:start-end.")
    prime_dag.add_argument("--out", type=Path, help="Output path for single-peg runs.")
    prime_dag.add_argument("--out-dir", type=Path, help="Directory for per-peg outputs.")
    prime_dag.add_argument("--json", type=Path, help="Deprecated alias for --out (defaults to stdout).")
    prime_dag.add_argument(
        "--frames",
        type=str,
        help="Optional JSONL path to stream DAG frames (use '-' for stdout).",
    )
    prime_dag.add_argument(
        "--coding-json",
        type=Path,
        help="Transcript JSON for protein-impact annotations (same schema as the GUI transcript loader).",
    )
    prime_dag.add_argument(
        "--coding-transcript",
        help="Transcript ID or name to use from --coding-json (defaults to the first entry).",
    )
    prime_dag.set_defaults(func=command_prime_dag)

    pcr_cmd = subparsers.add_parser("pcr", help="PCR amplification helpers.")
    pcr_sub = pcr_cmd.add_subparsers(dest="pcr_command", required=True)

    pcr_dag = pcr_sub.add_parser("dag", help="Construct a PCR amplicon Edit DAG.")
    pcr_dag.add_argument("--genome", type=Path, required=True, help="Genome FASTA file.")
    pcr_dag.add_argument("--primer-config", type=Path, required=True, help="Primer pair JSON configuration.")
    pcr_dag.add_argument("--pcr-config", type=Path, help="PCR configuration JSON file.")
    pcr_dag.add_argument("--cycles", type=int, help="Override PCR cycle count.")
    pcr_dag.add_argument(
        "--per-cycle-efficiency",
        type=float,
        help="Override per-cycle amplification efficiency (0-1).",
    )
    pcr_dag.add_argument("--error-rate", type=float, help="Override per-base error rate.")
    pcr_dag.add_argument("--min-amplicon-length", type=int, help="Minimum amplicon length filter.")
    pcr_dag.add_argument("--max-amplicon-length", type=int, help="Maximum amplicon length filter.")
    pcr_dag.add_argument("--max-amplicons", type=int, help="Maximum amplicon branches to keep.")
    pcr_dag.add_argument("--seed", type=int, default=0, help="Random seed (default: 0).")
    pcr_dag.add_argument(
        "--min-prob",
        type=float,
        default=1e-6,
        help="Minimum branch probability (default: 1e-6).",
    )
    pcr_dag.add_argument(
        "--out",
        type=Path,
        help="JSON output path for the PCR DAG (default: stdout).",
    )
    pcr_dag.set_defaults(func=command_pcr_dag)

    edit_dag_cmd = subparsers.add_parser("edit-dag", help="Edit DAG utilities.")
    edit_dag_sub = edit_dag_cmd.add_subparsers(dest="edit_dag_command", required=True)

    edit_dag_viz = edit_dag_sub.add_parser("viz", help="Render an edit DAG artifact as a PNG.")
    edit_dag_viz.add_argument("--input", type=Path, required=True, help="Path to an edit DAG JSON artifact.")
    edit_dag_viz.add_argument("--out", type=Path, required=True, help="PNG output path.")
    edit_dag_viz.add_argument(
        "--layout",
        choices=["timeline", "spring", "kamada-kawai", "spectral", "circular", "shell", "planar", "random"],
        default="timeline",
        help="Graph layout algorithm (default: timeline).",
    )
    edit_dag_viz.add_argument(
        "--min-prob",
        type=float,
        default=0.0,
        help="Hide nodes with normalized probability below this threshold (default: 0).",
    )
    edit_dag_viz.add_argument(
        "--max-time",
        type=int,
        help="Hide nodes with time_step greater than this value.",
    )
    edit_dag_viz.set_defaults(func=command_edit_dag_viz)
    edit_dag_report = edit_dag_sub.add_parser("report", help="Generate an HTML report for an edit DAG.")
    edit_dag_report.add_argument("--input", type=Path, required=True, help="Path to an edit DAG JSON artifact.")
    edit_dag_report.add_argument("--out", type=Path, required=True, help="HTML output path.")
    edit_dag_report.add_argument("--png", type=Path, help="Optional pre-rendered PNG to embed.")
    edit_dag_report.add_argument("--title", default="Helix Edit DAG Report", help="Title for the report.")
    edit_dag_report.set_defaults(func=command_edit_dag_report)
    edit_dag_anim = edit_dag_sub.add_parser("animate", help="Render an edit DAG animation (GIF).")
    edit_dag_anim.add_argument("--input", type=Path, required=True, help="Path to an edit DAG JSON artifact.")
    edit_dag_anim.add_argument("--out", type=Path, required=True, help="GIF output path.")
    edit_dag_anim.add_argument("--fps", type=int, default=3, help="Frames per second for the animation (default: 3).")
    edit_dag_anim.add_argument(
        "--layout",
        choices=["timeline", "spring", "kamada-kawai", "spectral", "circular", "shell", "planar", "random"],
        default="timeline",
        help="Graph layout algorithm (default: timeline).",
    )
    edit_dag_anim.add_argument("--width", type=float, default=8.0, help="Figure width in inches (default: 8).")
    edit_dag_anim.add_argument("--height", type=float, default=6.0, help="Figure height in inches (default: 6).")
    edit_dag_anim.add_argument(
        "--bg-color",
        default="#0b0e14",
        help="Background color for frames (default: #0b0e14).",
    )
    edit_dag_anim.add_argument("--prob-cmap", default="viridis", help="Matplotlib colormap for probabilities.")
    edit_dag_anim.add_argument("--stage-cmap", default="tab10", help="Matplotlib colormap for stage outlines.")
    edit_dag_anim.add_argument(
        "--highlight-color",
        default="#f5a97f",
        help="Accent color for new nodes/edges (default: #f5a97f).",
    )
    edit_dag_anim.add_argument("--dpi", type=int, default=150, help="Figure DPI (default: 150).")
    edit_dag_anim.add_argument(
        "--min-prob",
        type=float,
        default=0.0,
        help="Hide nodes with normalized probability below this threshold before animating.",
    )
    edit_dag_anim.add_argument(
        "--max-time",
        type=int,
        help="Hide nodes with time_step greater than this value before animating.",
    )
    edit_dag_anim.set_defaults(func=command_edit_dag_animate)
    edit_dag_compare = edit_dag_sub.add_parser("compare", help="Compare two edit DAG artifacts.")
    edit_dag_compare.add_argument("--a", type=Path, required=True, help="Path to the first edit DAG JSON artifact.")
    edit_dag_compare.add_argument("--b", type=Path, required=True, help="Path to the second edit DAG JSON artifact.")
    edit_dag_compare.add_argument("--out", type=Path, required=True, help="PNG output path.")
    edit_dag_compare.add_argument(
        "--layout",
        choices=["timeline", "spring", "kamada-kawai", "spectral", "circular", "shell", "planar", "random"],
        default="timeline",
        help="Graph layout algorithm (default: timeline).",
    )
    edit_dag_compare.add_argument(
        "--min-prob",
        type=float,
        default=0.0,
        help="Hide nodes with normalized probability below this threshold.",
    )
    edit_dag_compare.add_argument("--max-time", type=int, help="Hide nodes with time_step greater than this value.")
    edit_dag_compare.add_argument("--label-a", default="A", help="Legend label for the first DAG.")
    edit_dag_compare.add_argument("--label-b", default="B", help="Legend label for the second DAG.")
    edit_dag_compare.add_argument("--summary", type=Path, help="Optional JSON file to write diff summary.")
    edit_dag_compare.set_defaults(func=command_edit_dag_compare)
    edit_dag_dataset = edit_dag_sub.add_parser(
        "generate-dataset", help="Synthesize edit DAG records for downstream ML workflows."
    )
    edit_dag_dataset.add_argument("--mechanism", choices=["crispr", "prime"], default="crispr")
    edit_dag_dataset.add_argument("--n", type=int, required=True, help="Number of records to generate.")
    edit_dag_dataset.add_argument("--topk", type=int, default=5, help="How many outcomes to summarize per record.")
    edit_dag_dataset.add_argument("--out", type=Path, required=True, help="Output JSONL path.")
    edit_dag_dataset.add_argument("--seed", type=int, default=7, help="Random seed (default: 7).")
    edit_dag_dataset.add_argument(
        "--frames-input",
        type=Path,
        action="append",
        help="Optional JSONL frame files to convert into dataset records.",
    )
    edit_dag_dataset.set_defaults(func=command_edit_dag_generate_dataset)

    experiment_cmd = subparsers.add_parser("experiment", help="Experiment utilities.")
    experiment_sub = experiment_cmd.add_subparsers(dest="experiment_command", required=True)

    experiment_new = experiment_sub.add_parser("new", help="Create a starter experiment config.")
    experiment_new.add_argument("--type", choices=("crispr", "prime"), required=True, help="Experiment type.")
    experiment_new.add_argument("--out", type=Path, required=True, help="Output path for the template.")
    experiment_new.add_argument("--force", action="store_true", help="Overwrite existing files.")
    experiment_new.set_defaults(func=command_experiment_new)

    experiment_run = experiment_sub.add_parser("run", help="Execute an experiment config and write a DAG.")
    experiment_run.add_argument("--config", type=Path, required=True, help="Experiment config (.helix.yml).")
    experiment_run.add_argument("--out", type=Path, required=True, help="Output DAG JSON path.")
    experiment_run.add_argument("--use-gpu", action="store_true", help="Use GPU backend when available.")
    experiment_run.add_argument(
        "--frames",
        type=str,
        help="Optional JSONL path to stream DAG frames (use '-' for stdout).",
    )
    experiment_run.set_defaults(func=command_experiment_run)

    experiment_viz = experiment_sub.add_parser("viz", help="Run an experiment and render a PNG.")
    experiment_viz.add_argument("--config", type=Path, required=True, help="Experiment config (.helix.yml).")
    experiment_viz.add_argument("--out", type=Path, required=True, help="PNG output path.")
    experiment_viz.add_argument("--layout", default="cose", help="Graph layout (default: cose).")
    experiment_viz.add_argument("--min-prob", type=float, default=0.0, help="Minimum branch probability filter.")
    experiment_viz.add_argument("--max-time", type=int, help="Maximum time_step to display.")
    experiment_viz.add_argument("--use-gpu", action="store_true", help="Use GPU backend when available.")
    experiment_viz.set_defaults(func=command_experiment_viz)

    experiment_report = experiment_sub.add_parser("report", help="Run an experiment and generate an HTML report.")
    experiment_report.add_argument("--config", type=Path, required=True, help="Experiment config (.helix.yml).")
    experiment_report.add_argument("--out", type=Path, required=True, help="HTML output path.")
    experiment_report.add_argument("--png", type=Path, help="Optional PNG to embed in the report.")
    experiment_report.add_argument("--title", default="Helix Experiment Report", help="Report title.")
    experiment_report.add_argument("--use-gpu", action="store_true", help="Use GPU backend when available.")
    experiment_report.set_defaults(func=command_experiment_report)

    experiment_validate = experiment_sub.add_parser("validate", help="Validate an experiment config.")
    experiment_validate.add_argument("--config", type=Path, required=True, help="Experiment config (.helix.yml).")
    experiment_validate.add_argument("--use-gpu", action="store_true", help="Use GPU backend when available.")
    experiment_validate.set_defaults(func=command_experiment_validate)

    gui_cmd = subparsers.add_parser("gui", help="Launch the PySide6 realtime simulator.")
    gui_cmd.set_defaults(func=command_gui)

    protein = subparsers.add_parser("protein", help="Summarize protein sequences (requires Biopython).")
    protein.add_argument("--sequence", help="Inline amino-acid string.")
    protein.add_argument("--input", type=Path, help="FASTA/text file containing a protein sequence.")
    protein.add_argument("--window", type=int, default=9, help="Hydropathy window size (default: 9).")
    protein.add_argument("--step", type=int, default=1, help="Hydropathy step size (default: 1).")
    protein.add_argument("--scale", default="kd", help="Hydropathy scale key (default: kd).")
    protein.add_argument("--top", type=int, default=5, help="How many windows to print (default: 5).")
    protein.set_defaults(func=command_protein)

    triage_cmd = subparsers.add_parser("triage", help="Run the combined GC/k-mer/ORF triage report.")
    triage_cmd.add_argument("--sequence", help="Inline DNA/RNA sequence.")
    triage_cmd.add_argument("--input", type=Path, help="Path to a FASTA/text file.")
    triage_cmd.add_argument("--k", type=int, default=5, help="k-mer length (default: 5).")
    triage_cmd.add_argument("--max-diff", type=int, default=1, help="Allowed mismatches for k-mer clustering (default: 1).")
    triage_cmd.add_argument("--min-orf-length", type=int, default=90, help="Minimum ORF length in nucleotides (default: 90).")
    triage_cmd.add_argument("--top", type=int, default=5, help="Show this many clusters/ORFs (default: 5).")
    triage_cmd.add_argument("--json", type=Path, help="Optional path to write the entire report as JSON.")
    triage_cmd.set_defaults(func=command_triage)

    workflows_cmd = subparsers.add_parser("workflows", help="Run YAML-defined workflows.")
    workflows_cmd.add_argument("--config", type=Path, required=True, help="Path to a workflow YAML file.")
    workflows_cmd.add_argument(
        "--output-dir",
        type=Path,
        default=Path("workflow_runs"),
        help="Directory for workflow logs/output (default: workflow_runs).",
    )
    workflows_cmd.add_argument("--name", help="Optional workflow name to run.")
    workflows_cmd.add_argument(
        "--with-schema",
        action="store_true",
        help="Print a schema provenance summary for each step after execution.",
    )
    workflows_cmd.add_argument(
        "--as-json",
        action="store_true",
        help="Emit schema provenance as JSON (requires --with-schema).",
    )
    workflows_cmd.set_defaults(func=command_workflows)

    viz_cmd = subparsers.add_parser("viz", help="Visualization helpers.")
    viz_subparsers = viz_cmd.add_subparsers(dest="viz_command", required=True)

    viz_triage = viz_subparsers.add_parser("triage", help="Plot a triage JSON payload.")
    viz_triage.add_argument("--json", type=Path, required=True, help="Path to a triage JSON file.")
    viz_triage.add_argument(
        "--output",
        type=Path,
        default=Path("triage_viz.png"),
        help="Output image path (default: triage_viz.png).",
    )
    viz_triage.add_argument("--top", type=int, default=5, help="Top N clusters/ORFs to visualize (default: 5).")
    viz_triage.add_argument("--show", action="store_true", help="Display interactively.")
    viz_triage.set_defaults(func=command_viz_triage)

    viz_hydro = viz_subparsers.add_parser("hydropathy", help="Plot a hydropathy profile for a protein.")
    viz_hydro.add_argument("--sequence", help="Inline amino-acid string.")
    viz_hydro.add_argument("--input", type=Path, help="FASTA/text file containing a protein sequence.")
    viz_hydro.add_argument("--window", type=int, default=9, help="Window size (default: 9).")
    viz_hydro.add_argument("--step", type=int, default=1, help="Step size (default: 1).")
    viz_hydro.add_argument("--scale", default="kd", help="Hydropathy scale (default: kd).")
    viz_hydro.add_argument(
        "--output",
        type=Path,
        default=Path("hydropathy.png"),
        help="Output image path (default: hydropathy.png).",
    )
    viz_hydro.add_argument("--show", action="store_true", help="Display interactively.")
    viz_hydro.set_defaults(func=command_viz_hydropathy)

    viz_crispr = viz_subparsers.add_parser("crispr-track", help="Plot a CRISPR guide/outcome track.")
    viz_crispr.add_argument("--input", type=Path, required=True, help="Path to a crispr.sim JSON payload.")
    viz_crispr.add_argument("--save", type=Path, help="Optional output image path.")
    viz_crispr.add_argument("--save-viz-spec", type=Path, help="Optional viz-spec JSON output path.")
    viz_crispr.add_argument("--show", action="store_true", help="Display interactively.")
    viz_crispr.set_defaults(func=command_viz_crispr_track)

    viz_min = viz_subparsers.add_parser("minimizers", help="Plot minimizer density from a JSON payload.")
    viz_min.add_argument("--input", type=Path, help="JSON with sequence_length and minimizers.")
    viz_min.add_argument("--bins", type=int, default=200, help="Number of bins (default: 200).")
    viz_min.add_argument("--save", type=Path, help="Optional output image path.")
    viz_min.add_argument("--save-viz-spec", type=Path, help="Optional viz-spec JSON output.")
    viz_min.add_argument("--schema", action="store_true", help="Print schema/sample and exit.")
    viz_min.set_defaults(func=command_viz_minimizers)

    viz_seed = viz_subparsers.add_parser("seed-chain", help="Plot chained seed anchors.")
    viz_seed.add_argument("--input", type=Path, help="JSON with ref_length, qry_length, chains.")
    viz_seed.add_argument("--save", type=Path, help="Optional output image path.")
    viz_seed.add_argument("--save-viz-spec", type=Path, help="Optional viz-spec JSON output.")
    viz_seed.add_argument("--schema", action="store_true", help="Print schema/sample and exit.")
    viz_seed.set_defaults(func=command_viz_seed_chain)

    viz_rna_plot = viz_subparsers.add_parser("rna-dotplot", help="Plot RNA pairing posterior from JSON.")
    viz_rna_plot.add_argument("--input", type=Path, help="JSON with posterior matrix.")
    viz_rna_plot.add_argument("--vmin", type=float, default=0.0)
    viz_rna_plot.add_argument("--vmax", type=float, default=1.0)
    viz_rna_plot.add_argument("--save", type=Path, help="Optional output image path.")
    viz_rna_plot.add_argument("--save-viz-spec", type=Path, help="Optional viz-spec JSON output.")
    viz_rna_plot.add_argument("--schema", action="store_true", help="Print schema/sample and exit.")
    viz_rna_plot.set_defaults(func=command_viz_rna_dotplot)

    viz_schema_cmd = viz_subparsers.add_parser("schema", help="Describe viz schemas.")
    viz_schema_cmd.add_argument("--kind", help="Optional schema key (e.g., viz_minimizers).")
    viz_schema_cmd.set_defaults(func=command_viz_schema)

    viz_align = viz_subparsers.add_parser("alignment-ribbon", help="Plot an alignment ribbon from mapping JSON.")
    viz_align.add_argument("--input", type=Path, help="JSON output from 'helix seed map'.")
    viz_align.add_argument("--read-id", help="Read ID to plot (defaults to the first entry).")
    viz_align.add_argument("--alignment-index", type=int, default=0, help="Alignment index for the selected read.")
    viz_align.add_argument("--title", help="Override plot title.")
    viz_align.add_argument("--save", type=Path, help="Optional output image path.")
    viz_align.add_argument("--save-viz-spec", type=Path, help="Optional viz-spec JSON output.")
    viz_align.add_argument("--schema", action="store_true", help="Print schema/sample and exit.")
    viz_align.set_defaults(func=command_viz_alignment_ribbon)

    viz_dist = viz_subparsers.add_parser("distance-heatmap", help="Plot a distance matrix heatmap.")
    viz_dist.add_argument("--input", type=Path, help="JSON with 'matrix' and 'labels'.")
    viz_dist.add_argument("--save", type=Path, help="Optional output image path.")
    viz_dist.add_argument("--save-viz-spec", type=Path, help="Optional viz-spec JSON output.")
    viz_dist.add_argument("--schema", action="store_true", help="Print schema/sample and exit.")
    viz_dist.set_defaults(func=command_viz_distance_heatmap)

    viz_motif = viz_subparsers.add_parser("motif-logo", help="Plot a motif logo from a PWM JSON payload.")
    viz_motif.add_argument("--input", type=Path, help="JSON containing a 'pwm' entry.")
    viz_motif.add_argument("--title", help="Optional plot title override.")
    viz_motif.add_argument("--save", type=Path, help="Optional output image path.")
    viz_motif.add_argument("--save-viz-spec", type=Path, help="Optional viz-spec JSON output.")
    viz_motif.add_argument("--schema", action="store_true", help="Print schema/sample and exit.")
    viz_motif.set_defaults(func=command_viz_motif_logo)

    viz_modern = viz_subparsers.add_parser(
        "modern",
        help="Launch the PyQt + ModernGL CRISPR/prime-editing viewer.",
    )
    viz_modern.add_argument(
        "--spec",
        type=Path,
        action="append",
        help="Path to a visualization spec JSON (repeat for multiple edits).",
    )
    viz_modern.set_defaults(func=command_viz_modern)

    viz_modern_export = viz_subparsers.add_parser(
        "modern-export",
        help="Export CRISPR/prime edit timelines for Unity/Three.js.",
    )
    viz_modern_export.add_argument(
        "--spec",
        type=Path,
        action="append",
        required=True,
        help="Visualization spec JSON path(s) or bundles.",
    )
    viz_modern_export.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory to write timeline JSON files.",
    )
    viz_modern_export.add_argument(
        "--camera-samples",
        type=int,
        default=120,
        help="Number of orbit keyframes to bake into the export (default: 120).",
    )
    viz_modern_export.set_defaults(func=command_viz_modern_export)

    viz_sim_builder = viz_subparsers.add_parser(
        "sim-builder",
        help="Launch the interactive CRISPR/prime simulation UI.",
    )
    viz_sim_builder.set_defaults(func=command_viz_sim_builder)

    demo_cmd = subparsers.add_parser("demo", help="Demo helpers for Helix.")
    demo_sub = demo_cmd.add_subparsers(dest="demo_command", required=True)
    demo_viz = demo_sub.add_parser("viz", help="Render sample visualization payloads.")
    demo_viz.add_argument(
        "--output",
        type=Path,
        default=Path("demo_viz"),
        help="Directory to write demo JSON/PNG artifacts (default: demo_viz).",
    )
    demo_viz.set_defaults(func=command_demo_viz)

    schema_cmd = subparsers.add_parser("schema", help="Schema utilities.")
    schema_sub = schema_cmd.add_subparsers(dest="schema_command", required=True)
    schema_manifest_cmd = schema_sub.add_parser("manifest", help="Export schema manifest JSON.")
    schema_manifest_cmd.add_argument("--out", type=Path, help="Optional output path for the manifest JSON.")
    schema_manifest_cmd.set_defaults(func=command_schema_manifest)

    schema_diff_cmd = schema_sub.add_parser("diff", help="Diff schema manifests.")
    schema_diff_cmd.add_argument("--base", required=True, type=Path, help="Path to the base manifest JSON.")
    schema_diff_cmd.add_argument("--target", type=Path, help="Optional target manifest (defaults to current).")
    schema_diff_cmd.add_argument("--format", choices=["table", "json"], default="table")
    schema_diff_cmd.set_defaults(func=command_schema_diff)

    string_cmd = subparsers.add_parser("string", help="String / sequence search helpers.")
    string_sub = string_cmd.add_subparsers(dest="string_command", required=True)
    string_search = string_sub.add_parser("search", help="Exact or <=k edit-distance search.")
    string_search.add_argument("input", type=Path, help="FASTA or raw text file containing sequence(s).")
    string_search.add_argument("--pattern", required=True, help="Pattern to search for.")
    string_search.add_argument("--k", type=int, default=0, help="Maximum edit distance (default: 0).")
    string_search.add_argument("--json", type=Path, help="Optional path to write the JSON output.")
    string_search.set_defaults(func=command_string_search)

    seed_cmd = subparsers.add_parser("seed", help="Seed extraction and mapping helpers.")
    seed_sub = seed_cmd.add_subparsers(dest="seed_command", required=True)

    seed_index = seed_sub.add_parser("index", help="Compute minimizers or syncmers for a sequence.")
    seed_index.add_argument("input", type=Path, help="FASTA file.")
    seed_index.add_argument("--method", choices=["minimizer", "syncmer"], default="minimizer")
    seed_index.add_argument("--k", type=int, default=15, help="k-mer length.")
    seed_index.add_argument("--window", type=int, default=10, help="Window size for minimizers.")
    seed_index.add_argument("--sync", type=int, default=5, help="s-mer size for syncmers.")
    seed_index.add_argument("--json", type=Path, help="Optional JSON output file.")
    seed_index.add_argument("--plot", type=Path, help="Optional density plot path (requires matplotlib).")
    seed_index.set_defaults(func=command_seed_index)

    seed_map = seed_sub.add_parser("map", help="Seed-and-extend mapping (toy).")
    seed_map.add_argument("--ref", type=Path, required=True, help="Reference FASTA.")
    seed_map.add_argument("--reads", type=Path, required=True, help="Reads FASTA.")
    seed_map.add_argument("--k", type=int, default=15, help="k-mer length.")
    seed_map.add_argument("--window", type=int, default=10, help="Window size for minimizers.")
    seed_map.add_argument("--band", type=int, default=64, help="Band size for extension.")
    seed_map.add_argument("--xdrop", type=int, default=10, help="X-drop threshold.")
    seed_map.add_argument("--max-matches", type=int, default=3, help="Cap per-seed matches to avoid blowups.")
    seed_map.add_argument("--json", type=Path, help="Optional output path.")
    seed_map.set_defaults(func=command_seed_map)

    dbg_cmd = subparsers.add_parser("dbg", help="De Bruijn graph helpers.")
    dbg_sub = dbg_cmd.add_subparsers(dest="dbg_command", required=True)

    dbg_build = dbg_sub.add_parser("build", help="Build a DBG from reads.")
    dbg_build.add_argument("--reads", type=Path, nargs="+", required=True, help="FASTA/FASTQ files.")
    dbg_build.add_argument("--k", type=int, required=True, help="k-mer size.")
    dbg_build.add_argument("--graph", type=Path, required=True, help="Output JSON graph path.")
    dbg_build.add_argument("--graphml", type=Path, help="Optional GraphML output path.")
    dbg_build.set_defaults(func=command_dbg_build)

    dbg_clean = dbg_sub.add_parser("clean", help="Clean a DBG JSON (tips/bubbles).")
    dbg_clean.add_argument("--graph", type=Path, required=True, help="Input JSON graph path.")
    dbg_clean.add_argument("--out", type=Path, required=True, help="Output JSON path.")
    dbg_clean.add_argument("--tip-length", type=int, default=2, help="Tip length threshold (default: 2).")
    dbg_clean.add_argument("--no-tips", action="store_true", help="Disable tip removal.")
    dbg_clean.add_argument("--no-bubbles", action="store_true", help="Disable bubble removal.")
    dbg_clean.set_defaults(func=command_dbg_clean)

    dbg_color = dbg_sub.add_parser("color", help="Build a colored DBG from labeled read sets.")
    dbg_color.add_argument("--reads", type=Path, nargs="+", required=True, help="FASTA files per sample.")
    dbg_color.add_argument("--labels", nargs="+", help="Optional sample labels (defaults to filename stems).")
    dbg_color.add_argument("--k", type=int, required=True, help="k-mer size.")
    dbg_color.add_argument("--out", type=Path, required=True, help="Output colored graph JSON.")
    dbg_color.set_defaults(func=command_dbg_color)

    sketch_cmd = subparsers.add_parser("sketch", help="Sketch-based genome similarity helpers.")
    sketch_sub = sketch_cmd.add_subparsers(dest="sketch_command", required=True)
    sketch_build = sketch_sub.add_parser("build", help="Compute a MinHash or HLL sketch for a FASTA sequence.")
    sketch_build.add_argument("--method", choices=["minhash", "hll"], default="minhash")
    sketch_build.add_argument("--fasta", type=Path, required=True, help="Input FASTA.")
    sketch_build.add_argument("--k", type=int, default=21, help="k-mer size (default: 21).")
    sketch_build.add_argument("--size", type=int, default=1000, help="Sketch size (minhash only).")
    sketch_build.add_argument("--precision", type=int, default=10, help="HLL precision p (default: 10).")
    sketch_build.add_argument("--json", type=Path, help="Optional output path.")
    sketch_build.set_defaults(func=command_sketch_build)

    sketch_compare = sketch_sub.add_parser("compare", help="Compare two sequences using MinHash or HLL.")
    sketch_compare.add_argument("--method", choices=["minhash", "hll"], default="minhash")
    sketch_compare.add_argument("--fasta-a", type=Path, required=True)
    sketch_compare.add_argument("--fasta-b", type=Path, required=True)
    sketch_compare.add_argument("--k", type=int, default=21)
    sketch_compare.add_argument("--size", type=int, default=1000, help="Sketch size (minhash).")
    sketch_compare.add_argument("--precision", type=int, default=10, help="Precision for HLL.")
    sketch_compare.add_argument("--json", type=Path, help="Optional output path.")
    sketch_compare.set_defaults(func=command_sketch_compare)

    motif_cmd = subparsers.add_parser("motif", help="Motif discovery helpers.")
    motif_sub = motif_cmd.add_subparsers(dest="motif_command", required=True)
    motif_find = motif_sub.add_parser("find", help="Discover motifs via EM/STEME/online.")
    motif_find.add_argument("--fasta", type=Path, required=True, help="FASTA file with sequences.")
    motif_find.add_argument("--width", type=int, required=True, help="Motif width (k).")
    motif_find.add_argument("--solver", choices=["em", "steme", "online"], default="em", help="Solver to use (default: em).")
    motif_find.add_argument("--iterations", type=int, default=50, help="Iterations (EM/STEME).")
    motif_find.add_argument("--restarts", type=int, default=5, help="STEME random restarts (default: 5).")
    motif_find.add_argument("--learning-rate", type=float, default=0.3, help="Online learning rate (default: 0.3).")
    motif_find.add_argument("--passes", type=int, default=3, help="Online passes over the data (default: 3).")
    motif_find.add_argument("--json", type=Path, help="Optional JSON output path.")
    motif_find.add_argument("--plot", type=Path, help="Optional PWM plot path (requires matplotlib).")
    motif_find.add_argument("--plot-viz-spec", type=Path, help="Optional viz-spec JSON path for --plot.")
    motif_find.set_defaults(func=command_motif_find)

    veribiota_cmd = subparsers.add_parser("veribiota", help="Lean/VeriBiota export helpers.")
    veribiota_sub = veribiota_cmd.add_subparsers(dest="veribiota_command", required=True)

    veribiota_export = veribiota_sub.add_parser("export", help="Export an edit DAG payload to a Lean module.")
    veribiota_export.add_argument("--input", type=Path, required=True, help="Path to a Helix edit DAG JSON.")
    veribiota_export.add_argument("--out", type=Path, help="Lean output file (defaults to stdout).")
    veribiota_export.add_argument("--dag-name", default="example_dag", help="Lean identifier for the DAG value.")
    veribiota_export.add_argument("--module-name", default="VeriBiotaBridge", help="Lean namespace for the output.")
    veribiota_export.add_argument(
        "--lean-import",
        default="VeriBiota",
        help="Lean module to import for EditDAG definitions (default: VeriBiota).",
    )
    veribiota_export.add_argument(
        "--skip-eval",
        action="store_true",
        help="Skip emitting `#eval VeriBiota.check <dag>`.",
    )
    veribiota_export.add_argument(
        "--skip-theorem",
        action="store_true",
        help="Skip inserting the placeholder theorem/`admit` block.",
    )
    veribiota_export.set_defaults(func=command_veribiota_export)

    veribiota_check = veribiota_sub.add_parser(
        "lean-check", help="Generate lean-check summaries for downstream Lean validation."
    )
    veribiota_check.add_argument("--input", type=Path, required=True, help="Edit DAG JSON input.")
    veribiota_check.add_argument("--dag-name", help="Override DAG name stored in the check file.")
    veribiota_check.add_argument("--out", type=Path, help="Lean-check JSON output (defaults to stdout).")
    veribiota_check.add_argument(
        "--skip-validate",
        action="store_true",
        help="Skip verifying the generated lean-check payload.",
    )
    veribiota_check.add_argument(
        "--prob-tol",
        type=float,
        default=5e-3,
        help="Probability tolerance for validation (default: 5e-3).",
    )
    veribiota_check.set_defaults(func=command_veribiota_lean_check)

    veribiota_preflight = veribiota_sub.add_parser(
        "preflight", help="Validate lean-check JSON files (and optional source payloads)."
    )
    veribiota_preflight.add_argument("--checks", type=Path, nargs="+", required=True, help="Lean-check JSON paths.")
    veribiota_preflight.add_argument(
        "--payloads",
        type=Path,
        nargs="+",
        help="Optional edit DAG JSONs to cross-check seq hashes (must align with --checks).",
    )
    veribiota_preflight.add_argument(
        "--prob-tol",
        type=float,
        default=5e-3,
        help="Probability tolerance for validation (default: 5e-3).",
    )
    veribiota_preflight.set_defaults(func=command_veribiota_preflight)

    veribiota_multi = veribiota_sub.add_parser(
        "export-dags", help="Export multiple edit DAG payloads into a single Lean module."
    )
    veribiota_multi.add_argument("--inputs", type=Path, nargs="+", required=True, help="List of edit DAG JSON files.")
    veribiota_multi.add_argument("--out", type=Path, help="Lean output file (defaults to stdout).")
    veribiota_multi.add_argument("--dag-names", nargs="+", help="Optional Lean identifiers (defaults to file stems).")
    veribiota_multi.add_argument("--module-name", default="VeriBiotaBridge", help="Lean namespace for the output.")
    veribiota_multi.add_argument(
        "--lean-import",
        default="VeriBiota",
        help="Lean module to import for EditDAG definitions (default: VeriBiota).",
    )
    veribiota_multi.add_argument("--list-name", default="exampleDags", help="Lean identifier for the DAG list.")
    veribiota_multi.add_argument(
        "--theorem-name",
        help="Identifier for the aggregate theorem (default: <list_name>_all_checked).",
    )
    veribiota_multi.add_argument(
        "--eval",
        action="store_true",
        help="Emit `#eval VeriBiota.check <dag>` for each DAG (disabled by default).",
    )
    veribiota_multi.add_argument(
        "--skip-theorem",
        action="store_true",
        help="Skip inserting the aggregate theorem stub.",
    )
    veribiota_multi.set_defaults(func=command_veribiota_export_dags)

    veribiota_suite = veribiota_sub.add_parser(
        "export-suite", help="Export DAG JSONs directly into a VeriBiota checkout."
    )
    veribiota_suite.add_argument("--inputs", type=Path, nargs="+", required=True, help="DAG JSON inputs.")
    veribiota_suite.add_argument("--dag-names", nargs="+", help="Optional Lean identifiers for each DAG.")
    veribiota_suite.add_argument("--veribiota-root", type=Path, required=True, help="Path to VeriBiota repo root.")
    veribiota_suite.add_argument(
        "--module-path",
        type=Path,
        required=True,
        help="Relative path inside VeriBiota repo for the Lean module (e.g., Biosim/VeriBiota/Helix/CrisprMicro.lean).",
    )
    veribiota_suite.add_argument("--module-name", help="Override Lean module name (defaults to module path).")
    veribiota_suite.add_argument(
        "--lean-import",
        default="Biosim.VeriBiota.EditDAG",
        help="Lean module to import for EditDAG definitions (default: Biosim.VeriBiota.EditDAG).",
    )
    veribiota_suite.add_argument("--list-name", default="allDags", help="Identifier for the suite's DAG list.")
    veribiota_suite.add_argument(
        "--theorem-name",
        help="Identifier for the aggregate theorem (default: <list_name>_all_checked).",
    )
    veribiota_suite.add_argument(
        "--eval",
        action="store_true",
        help="Emit `#eval VeriBiota.check <dag>` for each DAG (disabled by default).",
    )
    veribiota_suite.add_argument(
        "--skip-theorem",
        action="store_true",
        help="Skip inserting the aggregate theorem stub.",
    )
    veribiota_suite.set_defaults(func=command_veribiota_export_suite)

    return parser


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
        args.func(args)
    except ValueError as exc:
        parser.error(str(exc))


if __name__ == "__main__":  # pragma: no cover - manual invocation path
    main()
