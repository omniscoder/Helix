"""Convert Helix edit DAG payloads into Lean/VeriBiota declarations."""
from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Sequence

_LEAN_QUOTE_PATTERN = re.compile(r'(["\\])')


def _lean_string(value: str) -> str:
    escaped = _LEAN_QUOTE_PATTERN.sub(r"\\\1", value)
    escaped = escaped.replace("\n", "\\n")
    return f"\"{escaped}\""


def _lean_float(value: float) -> str:
    if math.isnan(value):
        return "Float.nan"
    if math.isinf(value):
        return "Float.inf" if value > 0 else "Float.negInf"
    text = f"{value:.12g}"
    if "." not in text and "e" not in text.lower():
        text += ".0"
    return text


def _format_assoc_list(mapping: Mapping[str, Any]) -> str:
    if not mapping:
        return "[]"
    entries = []
    for key in sorted(mapping):
        value = mapping[key]
        if isinstance(value, (dict, list)):
            value_str = json.dumps(value, sort_keys=True)
        else:
            value_str = str(value)
        entries.append(f"({_lean_string(str(key))}, {_lean_string(value_str)})")
    return "[ " + ", ".join(entries) + " ]"


def _format_string_list(items: Sequence[str]) -> str:
    if not items:
        return "[]"
    return "[ " + ", ".join(_lean_string(str(item)) for item in items) + " ]"


def _format_event(event: Mapping[str, Any]) -> str:
    metadata = event.get("metadata") or {}
    fields = [
        ("chrom", _lean_string(str(event.get("chrom", "")))),
        ("start", str(int(event.get("start", 0)))),
        ("end", str(int(event.get("end", 0)))),
        ("replacement", _lean_string(str(event.get("replacement", "")))),
        ("metadata", _format_assoc_list(metadata)),
    ]
    return _format_record(fields)


def _format_record(fields: Iterable[tuple[str, str]]) -> str:
    lines = [f"{name} := {value}" for name, value in fields]
    if len(lines) == 1:
        return "{ " + lines[0] + " }"
    indented = ",\n".join("    " + line for line in lines)
    return "{\n" + indented + "\n}"


def _indent(text: str, spaces: int) -> str:
    prefix = " " * spaces
    return "\n".join(prefix + line if line else "" for line in text.splitlines())


def _format_nodes(nodes: Mapping[str, Any]) -> str:
    if not nodes:
        return "[]"
    formatted = []
    for node_id in sorted(nodes):
        node = nodes[node_id]
        metadata = node.get("metadata") or {}
        parent_ids = node.get("parent_ids") or node.get("parents") or []
        seq_hashes = node.get("seq_hashes") or {}
        diffs = [_format_event(entry) for entry in node.get("diffs", [])]
        sequences = node.get("sequences") or {}

        fields = [
            ("id", _lean_string(str(node_id))),
            ("logProb", _lean_float(float(node.get("log_prob", 0.0)))),
            ("metadata", _format_assoc_list(metadata)),
            ("parentIds", _format_string_list(parent_ids)),
            ("seqHashes", _format_assoc_list(seq_hashes)),
            (
                "diffs",
                "[]"
                if not diffs
                else "[\n"
                + ",\n".join(_indent(entry, 8) for entry in diffs)
                + "\n    ]",
            ),
            ("sequences", _format_assoc_list(sequences)),
        ]
        formatted.append(_indent(_format_record(fields), 4))
    return "[\n" + ",\n".join(formatted) + "\n  ]"


def _format_edges(edges: Sequence[Mapping[str, Any]]) -> str:
    if not edges:
        return "[]"
    formatted = []
    for edge in edges:
        metadata = edge.get("metadata") or {}
        event = edge.get("event") or {}
        fields = [
            ("source", _lean_string(str(edge.get("source", "")))),
            ("target", _lean_string(str(edge.get("target", "")))),
            ("rule", _lean_string(str(edge.get("rule", edge.get("rule_name", ""))))),
            ("event", _format_event(event)),
            ("metadata", _format_assoc_list(metadata)),
        ]
        formatted.append(_indent(_format_record(fields), 4))
    return "[\n" + ",\n".join(formatted) + "\n  ]"


def _safe_identifier(value: str, fallback: str) -> str:
    if not value:
        return fallback
    clean = re.sub(r"[^0-9A-Za-z_]", "_", value)
    if clean and clean[0].isdigit():
        clean = f"_{clean}"
    return clean or fallback


def _sanitize_namespace(namespace: str) -> list[str]:
    parts = [part for part in re.split(r"[.:]", namespace) if part]
    if not parts:
        parts = ["VeriBiotaBridge"]
    sanitized = [_safe_identifier(part, f"n{idx}") for idx, part in enumerate(parts)]
    return sanitized


def _build_dag_record(payload: Mapping[str, Any]) -> str:
    nodes = payload.get("nodes", {})
    edges = payload.get("edges", [])
    root_id = payload.get("root_id") or (nodes and sorted(nodes)[0]) or ""
    return _format_record(
        [
            ("nodes", _format_nodes(nodes)),
            ("edges", _format_edges(edges)),
            ("rootId", _lean_string(str(root_id))),
        ]
    )


def dag_payload_to_lean(
    payload: Mapping[str, Any],
    *,
    dag_name: str = "example_dag",
    module_name: str = "VeriBiotaBridge",
    import_module: str = "VeriBiota",
    include_eval: bool = True,
    include_theorem: bool = True,
) -> str:
    """Render a Lean module that recreates the provided DAG payload."""

    namespace_parts = _sanitize_namespace(module_name)
    safe_dag = _safe_identifier(dag_name, "example_dag")

    dag_record = _build_dag_record(payload)

    lines = [
        f"import {import_module}",
        "",
        f"open {import_module}",
        "",
    ]
    for part in namespace_parts:
        lines.extend([f"namespace {part}", ""])
    lines.extend(
        [
        f"/-- Auto-generated from a Helix edit DAG payload. -/",
        f"def {safe_dag} : EditDAG :=",
        _indent(dag_record, 2),
        "",
        ]
    )
    if include_eval:
        lines.extend(
            [
                f"#eval VeriBiota.check {safe_dag}",
                "",
            ]
        )
    if include_theorem:
        lines.extend(
            [
                "/-!",
                f"Replace the stub below with a proof once VeriBiota automation is available.",
                "",
                f"theorem {safe_dag}_verified : VeriBiota.check {safe_dag} := by",
                "  -- TODO: fill in proof.",
                "  admit",
                "-/",
                "",
            ]
        )
    for part in reversed(namespace_parts):
        lines.extend([f"end {part}", ""])
    return "\n".join(lines)


def dag_payloads_to_lean(
    payloads: Sequence[Mapping[str, Any]],
    dag_names: Sequence[str],
    *,
    module_name: str = "VeriBiotaBridge",
    import_module: str = "VeriBiota",
    include_eval: bool = False,
    include_theorem: bool = False,
    list_name: str = "exampleDags",
    theorem_name: str | None = None,
) -> str:
    """Render a Lean module defining several DAGs and an aggregate list."""

    if len(payloads) != len(dag_names):
        raise ValueError("payloads and dag_names must be the same length.")

    namespace_parts = _sanitize_namespace(module_name)
    safe_names = [_safe_identifier(name, f"dag{idx}") for idx, name in enumerate(dag_names)]
    dag_records = [_build_dag_record(payload) for payload in payloads]

    list_identifier = _safe_identifier(list_name, "exampleDags")
    theorem_identifier = _safe_identifier(theorem_name or f"{list_identifier}_all_checked", "allExamplesChecked")

    lines = [
        f"import {import_module}",
        "",
        f"open {import_module}",
        "",
    ]
    for part in namespace_parts:
        lines.extend([f"namespace {part}", ""])

    for dag_name, dag_record in zip(safe_names, dag_records):
        lines.extend(
            [
                f"/-- Auto-generated from a Helix edit DAG payload. -/",
                f"def {dag_name} : EditDAG :=",
                _indent(dag_record, 2),
                "",
            ]
        )
        if include_eval:
            lines.extend([f"#eval VeriBiota.check {dag_name}", ""])

    dag_list = (
        "[\n"
        + ",\n".join(_indent(name, 4) for name in safe_names)
        + ("\n  ]" if safe_names else "\n  ]")
    )
    lines.extend(
        [
            f"/-- Collection of auto-generated DAG examples. -/",
            f"def {list_identifier} : List EditDAG :=",
            _indent(dag_list, 2),
            "",
        ]
    )

    if include_theorem:
        lines.extend(
            [
                "/-!",
                f"Replace the stub below with a proof once VeriBiota automation is available.",
                "",
                f"theorem {theorem_identifier} :",
                f"  ∀ dag ∈ {list_identifier}, VeriBiota.check dag := by",
                "  intro dag hDag",
                "  -- TODO: provide proof for each DAG in the list.",
                "  admit",
                "-/",
                "",
            ]
        )

    for part in reversed(namespace_parts):
        lines.extend([f"end {part}", ""])
    return "\n".join(lines)


def module_name_from_path(module_path: Path) -> str:
    """Derive a Lean module name from a filesystem path."""

    stem_path = module_path.with_suffix("")
    parts = [part for part in stem_path.parts if part not in ("", ".")]
    return ".".join(parts)
