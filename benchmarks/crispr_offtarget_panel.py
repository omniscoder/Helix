"""Benchmark CRISPR off-target scoring against synthetic assay counts.

This harness treats each entry as a guide + off-target assay over a shared
digital genome window. For each guide it:

  * Enumerates off-target candidates via helix.crispr.score.enumerate_off_targets.
  * Scores hits with helix.crispr.score.score_off_targets.
  * Compares scores to synthetic read counts from a TSV.

All sequences are purely computational DNA strings; no wet-lab conditions
are modeled here.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import yaml

from helix import bioinformatics
from helix.crispr import score as crispr_score

PANEL_KIND = "helix.crispr.offtarget_panel.v1"


@dataclass(frozen=True)
class GenomeSpec:
    fasta: Path


@dataclass(frozen=True)
class OffTargetAssayColumns:
    chrom: str
    start: str
    strand: str
    read_count: str


@dataclass(frozen=True)
class OffTargetAssaySpec:
    file: Path
    columns: OffTargetAssayColumns
    assay_type: Optional[str]


@dataclass(frozen=True)
class OffTargetSearchSpec:
    max_mm: int
    max_gap: int


@dataclass(frozen=True)
class OffTargetGuideSpec:
    id: str
    sequence: str
    pam: str
    search: OffTargetSearchSpec
    assay: Optional[OffTargetAssaySpec]


@dataclass(frozen=True)
class OffTargetPanelSpec:
    name: str
    description: Optional[str]
    kind: str
    genome: GenomeSpec
    pam_default: str
    max_mm_default: int
    max_gap_default: int
    guides: List[OffTargetGuideSpec]


def _resolve_path(base: Path, value: str) -> Path:
    path = Path(str(value).strip())
    if not path.is_absolute():
        path = (base / path).resolve()
    return path


def _parse_genome(base_dir: Path, raw: Mapping[str, Any]) -> GenomeSpec:
    fasta_raw = raw.get("fasta")
    if not fasta_raw:
        raise ValueError("Off-target panel requires genome.fasta.")
    fasta_path = _resolve_path(base_dir, str(fasta_raw))
    if not fasta_path.exists():
        raise FileNotFoundError(f"Genome FASTA '{fasta_path}' not found.")
    return GenomeSpec(fasta=fasta_path)


def _parse_assay(base_dir: Path, raw: Mapping[str, Any]) -> OffTargetAssaySpec:
    file_raw = raw.get("file")
    if not file_raw:
        raise ValueError("Off-target assay requires a 'file' path.")
    file_path = _resolve_path(base_dir, str(file_raw))
    cols_raw = raw.get("columns") or {}
    if not isinstance(cols_raw, Mapping):
        cols_raw = {}
    chrom_col = str(cols_raw.get("chrom") or "chrom")
    start_col = str(cols_raw.get("start") or "start")
    strand_col = str(cols_raw.get("strand") or "strand")
    read_col = str(cols_raw.get("read_count") or cols_raw.get("reads") or "reads")
    columns = OffTargetAssayColumns(
        chrom=chrom_col,
        start=start_col,
        strand=strand_col,
        read_count=read_col,
    )
    assay_type_raw = raw.get("type")
    assay_type = str(assay_type_raw).strip() or None if assay_type_raw is not None else None
    return OffTargetAssaySpec(file=file_path, columns=columns, assay_type=assay_type)


def _parse_search(defaults: Mapping[str, Any], raw: Mapping[str, Any]) -> OffTargetSearchSpec:
    max_mm = int(raw.get("max_mm", defaults.get("max_mm", 1)))
    max_gap = int(raw.get("max_gap", defaults.get("max_gap", 0)))
    return OffTargetSearchSpec(max_mm=max_mm, max_gap=max_gap)


def _parse_guide(base_dir: Path, defaults: Mapping[str, Any], raw: Mapping[str, Any]) -> OffTargetGuideSpec:
    guide_id = str(raw.get("id") or "").strip() or "guide"
    seq = bioinformatics.normalize_sequence(str(raw.get("sequence") or ""))
    if not seq:
        raise ValueError(f"Guide '{guide_id}' requires a non-empty 'sequence'.")
    pam = str(raw.get("pam") or defaults.get("pam") or "SpCas9-NGG").strip() or "SpCas9-NGG"

    search_raw = raw.get("search") or {}
    if not isinstance(search_raw, Mapping):
        search_raw = {}
    search = _parse_search(defaults, search_raw)

    assay_raw = raw.get("assay")
    assay: Optional[OffTargetAssaySpec]
    if isinstance(assay_raw, Mapping):
        assay = _parse_assay(base_dir, assay_raw)
    else:
        assay = None

    return OffTargetGuideSpec(
        id=guide_id,
        sequence=seq,
        pam=pam,
        search=search,
        assay=assay,
    )


def load_panel_spec(path: Path) -> OffTargetPanelSpec:
    cfg_path = Path(path)
    if not cfg_path.exists():
        raise FileNotFoundError(f"Off-target panel config '{cfg_path}' not found.")
    data = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    if not isinstance(data, Mapping):
        raise ValueError("Off-target panel config must be a YAML mapping.")

    kind = str(data.get("kind") or "").strip()
    if kind and kind != PANEL_KIND:
        raise ValueError(f"Unexpected panel kind '{kind}'. Expected '{PANEL_KIND}'.")

    name = str(data.get("name") or "unnamed_offtarget_panel")
    desc_raw = data.get("description")
    description = str(desc_raw).strip() or None if desc_raw is not None else None

    base_dir = cfg_path.parent
    genome_raw = data.get("genome") or {}
    if not isinstance(genome_raw, Mapping):
        genome_raw = {}
    genome = _parse_genome(base_dir, genome_raw)

    defaults_raw = data.get("defaults") or {}
    if not isinstance(defaults_raw, Mapping):
        defaults_raw = {}
    pam_default = str(defaults_raw.get("pam") or "SpCas9-NGG").strip() or "SpCas9-NGG"
    max_mm_default = int(defaults_raw.get("max_mm", 1))
    max_gap_default = int(defaults_raw.get("max_gap", 0))

    guides_raw = data.get("guides") or []
    if not isinstance(guides_raw, list) or not guides_raw:
        raise ValueError("Off-target panel config requires a non-empty 'guides' list.")
    guides = [
        _parse_guide(base_dir, defaults_raw, entry)
        for entry in guides_raw
        if isinstance(entry, Mapping)
    ]

    return OffTargetPanelSpec(
        name=name,
        description=description,
        kind=PANEL_KIND,
        genome=genome,
        pam_default=pam_default,
        max_mm_default=max_mm_default,
        max_gap_default=max_gap_default,
        guides=guides,
    )


def _load_genome_sequence(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    seq_lines = [line for line in lines if not line.startswith(">")]
    return bioinformatics.normalize_sequence("".join(seq_lines))


def _load_assay_counts(spec: OffTargetAssaySpec) -> Dict[Tuple[str, int], int]:
    if not spec.file.exists():
        return {}
    counts: Dict[Tuple[str, int], int] = {}
    with spec.file.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            chrom = str(row.get(spec.columns.chrom) or "").strip()
            start_raw = row.get(spec.columns.start)
            strand = str(row.get(spec.columns.strand) or "").strip() or "+"
            reads_raw = row.get(spec.columns.read_count)
            if start_raw is None or reads_raw is None:
                continue
            try:
                start = int(start_raw)
                reads = int(reads_raw)
            except ValueError:
                continue
            key = (strand, start)
            counts[key] = counts.get(key, 0) + max(0, reads)
    return counts


def _pearson(x: Sequence[float], y: Sequence[float]) -> Optional[float]:
    if len(x) != len(y) or len(x) < 2:
        return None
    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)
    num = 0.0
    den_x = 0.0
    den_y = 0.0
    for xi, yi in zip(x, y):
        dx = xi - mean_x
        dy = yi - mean_y
        num += dx * dy
        den_x += dx * dx
        den_y += dy * dy
    if den_x <= 0.0 or den_y <= 0.0:
        return None
    return num / math.sqrt(den_x * den_y)


def _run_for_guide(
    genome_seq: str,
    guide: OffTargetGuideSpec,
) -> Dict[str, Any]:
    guide_dict = {"id": guide.id, "sequence": guide.sequence}
    hits = crispr_score.enumerate_off_targets(
        genome_seq,
        guide_dict,
        guide.pam,
        max_mm=guide.search.max_mm,
        max_gap=guide.search.max_gap,
    )
    weights = crispr_score.load_weights(None)
    off_params = weights.get("off_target", {})
    scored_hits = crispr_score.score_off_targets(hits, off_params)

    # Map predicted hits by (strand, start) for comparison.
    pred_map: Dict[Tuple[str, int], Dict[str, Any]] = {}
    for hit in scored_hits:
        strand = str(hit.get("strand") or "").strip() or "+"
        start = int(hit.get("start", 0))
        key = (strand, start)
        pred_map[key] = {
            "strand": strand,
            "start": start,
            "end": int(hit.get("end", start)),
            "distance": int(hit.get("distance", 0)),
            "pam_ok": bool(hit.get("pam_ok", True)),
            "score": float(hit.get("score", 0.0)),
        }

    entry: Dict[str, Any] = {
        "guide_id": guide.id,
        "sequence": guide.sequence,
        "pam": guide.pam,
        "max_mm": guide.search.max_mm,
        "max_gap": guide.search.max_gap,
        "predicted_hits": sorted(
            pred_map.values(),
            key=lambda h: float(h.get("score", 0.0)),
            reverse=True,
        ),
    }

    if guide.assay is None:
        entry["status"] = "missing_assay"
        return entry

    obs_counts = _load_assay_counts(guide.assay)
    if not obs_counts:
        entry["status"] = "missing_counts"
        return entry

    # Build overlapping pairs (score vs log1p(reads)) for correlation.
    scores: List[float] = []
    logs: List[float] = []
    overlap = 0
    for key, count in obs_counts.items():
        hit = pred_map.get(key)
        if hit is None:
            continue
        overlap += 1
        scores.append(float(hit["score"]))
        logs.append(math.log1p(max(0, count)))

    metrics: Dict[str, Any] = {
        "observed_sites": len(obs_counts),
        "predicted_sites": len(pred_map),
        "overlap_sites": overlap,
    }
    if scores and logs:
        r = _pearson(scores, logs)
        if r is not None:
            metrics["pearson_r_score_vs_log_reads"] = r

    entry["status"] = "ok"
    entry["metrics"] = metrics
    # Also expose observed sites for downstream plots.
    entry["observed_sites"] = [
        {
            "strand": key[0],
            "start": key[1],
            "reads": count,
        }
        for key, count in sorted(obs_counts.items(), key=lambda kv: kv[1], reverse=True)
    ]
    return entry


def run_panel(panel: OffTargetPanelSpec) -> Dict[str, Any]:
    genome_seq = _load_genome_sequence(panel.genome.fasta)
    results: List[Dict[str, Any]] = []
    pearsons: List[float] = []

    for guide in panel.guides:
        entry = _run_for_guide(genome_seq, guide)
        results.append(entry)
        metrics = entry.get("metrics") or {}
        r = metrics.get("pearson_r_score_vs_log_reads")
        if isinstance(r, (int, float)):
            pearsons.append(float(r))

    summary: Dict[str, Any] = {
        "panel_name": panel.name,
        "panel_description": panel.description,
        "guides": len(panel.guides),
        "evaluated": sum(1 for e in results if e.get("status") == "ok"),
    }
    if pearsons:
        summary["mean_pearson_r"] = sum(pearsons) / len(pearsons)

    payload: Dict[str, Any] = {
        "schema": {"kind": "crispr_offtarget_panel", "spec_version": "1.0"},
        "panel": {
            "name": panel.name,
            "description": panel.description,
            "kind": panel.kind,
        },
        "summary": summary,
        "results": results,
    }
    return payload


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare Helix CRISPR off-target scores to assay-style hit counts."
    )
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Off-target panel YAML config (see templates/crispr_offtarget_panel.helix.yml).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output JSON path for benchmark results.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_arg_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    panel = load_panel_spec(args.config)
    payload = run_panel(panel)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"[OK] CRISPR off-target panel results written to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

