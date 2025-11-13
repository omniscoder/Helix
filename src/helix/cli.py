"""Unified Helix CLI that wraps DNA, peptide, RNA, protein, and triage helpers."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import re
import sys
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, TextIO
from types import SimpleNamespace

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
from .prime.simulator import simulate_prime_edit
from .genome.digital import DigitalGenome as CoreDigitalGenome
from .crispr.dag_api import build_crispr_edit_dag
from .prime.dag_api import build_prime_edit_dag
from .pcr.model import Primer, PrimerPair, PCRConfig
from .pcr.dag_api import pcr_edit_dag
from .edit.dag import EditDAG, EditDAGFrame, EditEdge, EditNode, dag_from_payload
from .edit.report import render_html_report
from .experiment.spec import (
    CrisprExperimentSpec,
    ExperimentSpecError,
    PrimeExperimentSpec,
    load_experiment_spec,
)

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


def _current_command_str() -> str:
    return "helix " + " ".join(sys.argv[1:])


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


def _ensure_viz_available(feature: str = "visualization") -> None:
    if not VIZ_AVAILABLE:
        raise SystemExit(
            f"{feature} requires the 'viz' extra. Install with `pip install \"veri-helix[viz]\"`."
        )


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
            command=_current_command_str(),
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
            "command": _current_command_str(),
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
    try:
        outcomes = simulate_prime_edit(
            genome,
            editor,
            peg,
            max_outcomes=args.max_outcomes,
        ) or []
    except NotImplementedError as exc:
        raise SystemExit(str(exc))

    payload = {
        "schema": {"kind": "prime.edit_sim", "spec_version": SPEC_VERSION},
        "meta": {
            "helix_version": HELIX_VERSION,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "command": _current_command_str(),
        },
        "params": {"max_outcomes": args.max_outcomes},
        "editor": _serialize_prime_editor(editor),
        "peg": _serialize_peg(peg),
        "genome": _digital_genome_summary(genome, source=args.genome),
        "outcomes": [_serialize_prime_outcome(outcome) for outcome in outcomes],
    }
    validated = validate_viz_payload("prime.edit_sim", payload)
    _write_json_output(validated, args.json)
    if not outcomes:
        print("No prime editing outcomes were produced with the current parameters.")


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
            command=_current_command_str(),
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
            command=_current_command_str(),
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
            command=_current_command_str(),
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
            command=_current_command_str(),
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
            command=_current_command_str(),
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
            command=_current_command_str(),
            viz_spec_path=spec_path,
        )


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
            command=_current_command_str(),
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
                command=_current_command_str(),
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
                command=_current_command_str(),
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
                command=_current_command_str(),
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
                command=_current_command_str(),
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
                command=_current_command_str(),
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
                command=_current_command_str(),
                viz_spec_path=viz_path,
            )
    print(f"Demo visualizations written to {output_dir}")


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
    prime_sim.add_argument("--max-outcomes", type=int, default=16, help="Maximum number of outcomes to emit (default: 16).")
    prime_sim.add_argument("--json", type=Path, help="Optional output JSON path (defaults to stdout).")
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
