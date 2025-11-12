"""Schema validation utilities for Helix visualization payloads."""
from __future__ import annotations

import json
import logging
from importlib import resources
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

_SPEC_MANIFEST = json.loads(
    resources.files(__package__).joinpath("spec_manifest.json").read_text(encoding="utf-8")
)
SPEC_VERSION: str = _SPEC_MANIFEST["spec_version"]

try:  # pragma: no cover - optional dependency
    from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator

    PYDANTIC_AVAILABLE = True
except ImportError:  # pragma: no cover - fallback path
    BaseModel = object  # type: ignore
    PYDANTIC_AVAILABLE = False

LOGGER = logging.getLogger(__name__)
if not PYDANTIC_AVAILABLE:  # pragma: no cover - informational log
    LOGGER.debug("pydantic not installed; viz payload validation is disabled.")


class SchemaError(RuntimeError):
    """Raised when a visualization payload fails schema validation."""


def _parse_version(value: str) -> List[int]:
    try:
        return [int(part) for part in value.split(".")]
    except ValueError:
        return []


def _check_spec_version(payload: Dict[str, Any], kind: str) -> None:
    version: Optional[str] = None
    if isinstance(payload, dict):
        if "spec_version" in payload:
            version = str(payload["spec_version"])
        elif isinstance(payload.get("meta"), dict) and "spec_version" in payload["meta"]:
            version = str(payload["meta"]["spec_version"])
    if version:
        provided = _parse_version(version)
        required = _parse_version(SPEC_VERSION)
        if not provided:
            raise SchemaError(f"{kind}: schema mismatch (spec_version '{version}' is invalid).")
        if provided < required:
            raise SchemaError(
                f"{kind}: schema mismatch (provided spec_version {version}, expected >= {SPEC_VERSION})."
            )


if PYDANTIC_AVAILABLE:

    class _BaseSchemaModel(BaseModel):
        model_config = ConfigDict(extra="allow", protected_namespaces=(), populate_by_name=True)

    class MinimizersPayload(_BaseSchemaModel):
        sequence_length: int
        minimizers: List[Any]


    class SeedAnchor(_BaseSchemaModel):
        ref_start: int
        qry_start: int
        ref_end: Optional[int] = None
        qry_end: Optional[int] = None
        len: Optional[int] = None

        @model_validator(mode="after")
        def ensure_endpoints(cls, values: "SeedAnchor") -> "SeedAnchor":
            if values.ref_end is None and values.qry_end is None and values.len is None:
                raise ValueError("Provide ref_end/qry_end or len for each anchor.")
            return values


    class SeedChainPayload(_BaseSchemaModel):
        ref_length: int
        qry_length: int
        chains: List[List[SeedAnchor]]
        meta: Optional[Dict[str, Any]] = None


    class RNADotplotPayload(_BaseSchemaModel):
        posterior: List[List[float]]
        meta: Optional[Dict[str, Any]] = None

        @model_validator(mode="after")
        def ensure_square(cls, values: "RNADotplotPayload") -> "RNADotplotPayload":
            size = len(values.posterior)
            if any(len(row) != size for row in values.posterior):
                raise ValueError("posterior matrix must be square.")
            return values


    class AlignmentResult(_BaseSchemaModel):
        read_id: Optional[str] = None
        read_length: Optional[int] = None
        seed_hits: Optional[int] = None
        alignments: List[Dict[str, Any]]


    class AlignmentRibbonPayload(_BaseSchemaModel):
        ref_length: Optional[int] = None
        qry_length: Optional[int] = None
        cigar: Optional[str] = None
        ref_start: Optional[int] = None
        qry_start: Optional[int] = None
        metadata: Optional[Dict[str, Any]] = None
        results: Optional[List[AlignmentResult]] = None
        meta: Optional[Dict[str, Any]] = None

        @model_validator(mode="after")
        def ensure_alignment(cls, values: "AlignmentRibbonPayload") -> "AlignmentRibbonPayload":
            if values.results:
                return values
            missing = [
                field
                for field in ("ref_length", "qry_length", "cigar", "ref_start", "qry_start")
                if getattr(values, field) is None
            ]
            if missing:
                raise ValueError(f"alignment payload missing fields: {', '.join(missing)}")
            return values


    class DistanceHeatmapPayload(_BaseSchemaModel):
        labels: List[str]
        matrix: List[List[float]]
        method: Optional[str] = None
        meta: Optional[Dict[str, Any]] = None

        @model_validator(mode="after")
        def ensure_shape(cls, values: "DistanceHeatmapPayload") -> "DistanceHeatmapPayload":
            size = len(values.matrix)
            if size != len(values.labels):
                raise ValueError("labels length must match matrix size.")
            if any(len(row) != size for row in values.matrix):
                raise ValueError("matrix must be square.")
            return values


    PWMColumn = Union[Dict[str, float], Sequence[float]]


    class MotifLogoPayload(_BaseSchemaModel):
        alphabet: Optional[List[str]] = None
        pwm: List[PWMColumn]
        background: Optional[Sequence[float]] = None
        consensus: Optional[str] = None
        meta: Optional[Dict[str, Any]] = None
        spec_version: Optional[str] = None

        @model_validator(mode="after")
        def ensure_columns(cls, values: "MotifLogoPayload") -> "MotifLogoPayload":
            alphabet = values.alphabet or ["A", "C", "G", "T"]
            width = len(alphabet)
            for column in values.pwm:
                if isinstance(column, dict):
                    continue
                if len(list(column)) != width:
                    raise ValueError("PWM columns must match alphabet length.")
            if values.background is not None and len(list(values.background)) != width:
                raise ValueError("Background distribution must match alphabet length.")
            return values


    class CRISPRPamSite(_BaseSchemaModel):
        start: int
        end: int


    class CRISPRGuideEntry(_BaseSchemaModel):
        id: Optional[str] = None
        start: int
        end: int
        strand: str
        pam_site: CRISPRPamSite
        gc_content: float
        sequence: Optional[str] = None
        metadata: Optional[Dict[str, Any]] = None


    class CRISPRGuidesPayload(_BaseSchemaModel):
        schema_alias: Dict[str, Any] = Field(alias="schema")
        input_sha256: str
        sequence_length: int
        pam: Dict[str, Any]
        params: Dict[str, Any]
        guides: List[CRISPRGuideEntry]
        meta: Optional[Dict[str, Any]] = None


    class CRISPROffTargetHit(_BaseSchemaModel):
        start: int
        end: int
        strand: str
        pam_ok: bool
        mismatches: List[Dict[str, Any]] = []
        distance: Optional[int] = None
        score: Optional[float] = None


    class CRISPROffTargetGuideHits(_BaseSchemaModel):
        guide_id: Optional[str] = None
        hits: List[CRISPROffTargetHit]
        on_target_score: Optional[float] = None


    class CRISPROffTargetsPayload(_BaseSchemaModel):
        schema_alias: Dict[str, Any] = Field(alias="schema")
        input_sha256: str
        pam: Dict[str, Any]
        params: Dict[str, Any]
        guides: List[CRISPROffTargetGuideHits]
        meta: Optional[Dict[str, Any]] = None


    class CRISPRSimOutcome(_BaseSchemaModel):
        label: str
        count: int
        probability: float
        diff: Optional[Dict[str, Any]] = None


    class CRISPRSimPayload(_BaseSchemaModel):
        schema_alias: Dict[str, Any] = Field(alias="schema")
        site: Dict[str, Any]
        guide: Dict[str, Any]
        priors: Dict[str, Any]
        draws: int
        outcomes: List[CRISPRSimOutcome]
        meta: Optional[Dict[str, Any]] = None


    _SCHEMA_MODELS: Dict[str, Any] = {
        "viz_minimizers": MinimizersPayload,
        "viz_seed_chain": SeedChainPayload,
        "viz_rna_dotplot": RNADotplotPayload,
        "viz_alignment_ribbon": AlignmentRibbonPayload,
        "viz_distance_heatmap": DistanceHeatmapPayload,
        "viz_motif_logo": MotifLogoPayload,
        "crispr.guides": CRISPRGuidesPayload,
        "crispr.offtargets": CRISPROffTargetsPayload,
        "crispr.sim": CRISPRSimPayload,
    }
else:
    _SCHEMA_MODELS = {}


def validate_viz_payload(kind: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    """Validate a viz payload, raising SchemaError with helpful context."""
    _check_spec_version(payload, kind)
    if not PYDANTIC_AVAILABLE:
        return payload
    model = _SCHEMA_MODELS.get(kind)
    if not model:
        return payload
    try:
        validated = model.model_validate(payload)
    except ValidationError as exc:  # pragma: no cover - exercised via tests
        raise SchemaError(f"{kind} payload invalid: {exc}") from exc
    return validated.model_dump(exclude_none=True, by_alias=True)


def describe_schema(kind: str | None = None) -> str:
    """Return a JSON description of available schemas or a specific schema."""
    if not PYDANTIC_AVAILABLE:
        return "pydantic is not installed; install the 'schema' extra to enable schema descriptions."
    if kind is None:
        options = ", ".join(sorted(_SCHEMA_MODELS))
        return f"Available viz schemas (spec_version={SPEC_VERSION}): {options}"
    model = _SCHEMA_MODELS.get(kind)
    if not model:
        raise SchemaError(f"Unknown schema '{kind}'.")
    schema = model.model_json_schema()
    return json.dumps(schema, indent=2)


def manifest() -> Dict[str, Any]:
    """Return a manifest of all registered schemas."""
    data: Dict[str, Any] = {
        "spec_version": SPEC_VERSION,
        "schemas": {},
        "validator_enabled": PYDANTIC_AVAILABLE,
    }
    for kind, model in _SCHEMA_MODELS.items():
        entry: Dict[str, Any] = {
            "spec_version": SPEC_VERSION,
            "validator": PYDANTIC_AVAILABLE,
        }
        if PYDANTIC_AVAILABLE:
            entry["schema"] = model.model_json_schema()
        data["schemas"][kind] = entry
    return data


def load_manifest(path: Union[str, Path]) -> Dict[str, Any]:
    """Load a manifest JSON from disk."""
    with Path(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _prop_keys(entry: Dict[str, Any]) -> set:
    schema = entry.get("schema") or {}
    props = schema.get("properties")
    if isinstance(props, dict):
        return set(props.keys())
    return set()


def diff_manifests(base: Dict[str, Any], target: Dict[str, Any]) -> Dict[str, Any]:
    """Compute a simple diff between two manifest dictionaries."""
    diff: Dict[str, Any] = {
        "spec_versions": {"from": base.get("spec_version"), "to": target.get("spec_version")},
        "schemas": {"added": [], "removed": [], "modified": {}},
    }
    base_schemas = base.get("schemas", {})
    target_schemas = target.get("schemas", {})
    base_keys = set(base_schemas.keys())
    target_keys = set(target_schemas.keys())
    diff["schemas"]["added"] = sorted(target_keys - base_keys)
    diff["schemas"]["removed"] = sorted(base_keys - target_keys)
    for key in sorted(base_keys & target_keys):
        base_entry = base_schemas[key]
        target_entry = target_schemas[key]
        entry_diff: Dict[str, Any] = {}
        base_version = base_entry.get("spec_version")
        target_version = target_entry.get("spec_version")
        if base_version != target_version:
            entry_diff["spec_version"] = {"from": base_version, "to": target_version}
        base_props = _prop_keys(base_entry)
        target_props = _prop_keys(target_entry)
        added_props = sorted(target_props - base_props)
        removed_props = sorted(base_props - target_props)
        if added_props:
            entry_diff["properties_added"] = added_props
        if removed_props:
            entry_diff["properties_removed"] = removed_props
        if entry_diff:
            diff["schemas"]["modified"][key] = entry_diff
    return diff


def format_manifest_diff(diff: Dict[str, Any], fmt: str = "table") -> str:
    if fmt == "json":
        return json.dumps(diff, indent=2)
    lines: List[str] = []
    spec_from = diff["spec_versions"].get("from")
    spec_to = diff["spec_versions"].get("to")
    lines.append(f"Schema diff (from {spec_from} to {spec_to})")
    added = diff["schemas"]["added"]
    removed = diff["schemas"]["removed"]
    modified = diff["schemas"]["modified"]
    if not added and not removed and not modified:
        lines.append("No schema changes detected.")
        return "\n".join(lines)
    if added:
        lines.append("Added schemas:")
        for key in added:
            lines.append(f"  + {key}")
    if removed:
        lines.append("Removed schemas:")
        for key in removed:
            lines.append(f"  - {key}")
    for key, entry in modified.items():
        lines.append(f"Modified schema: {key}")
        spec_change = entry.get("spec_version")
        if spec_change:
            lines.append(f"    spec_version: {spec_change['from']} -> {spec_change['to']}")
        if entry.get("properties_added"):
            lines.append(f"    properties added: {', '.join(entry['properties_added'])}")
        if entry.get("properties_removed"):
            lines.append(f"    properties removed: {', '.join(entry['properties_removed'])}")
    return "\n".join(lines)
