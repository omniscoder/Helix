"""Minimal workflow runner used by the Helix CLI tests.

The original project split workflows into a separate package. To keep this
repository self-contained (and ensure `helix workflows` remains usable),
we provide a lightweight implementation that is sufficient for the tests:

- Supports simple YAML configs with `workflows:` and `steps:`.
- Each step has a `command` (string or list of CLI tokens), an `args` dict
  mapped directly to `--key value` flags, and an optional `stdout` filename.
- Optional `schema` entries:
    schema:
      kind: viz_alignment_ribbon
      output: map.json

  If present, we compute the schema kind/version and SHA-256 hash from the
  referenced JSON artifact and surface those back to the CLI for reporting.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

import yaml


@dataclass
class WorkflowStepResult:
    command: str
    output_path: Optional[Path]
    schema_kind: Optional[str]
    schema_version: Optional[str]
    schema_hash: Optional[str]


@dataclass
class WorkflowRunResult:
    name: str
    output_dir: Path
    steps: List[WorkflowStepResult]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _run_cli(command: Iterable[str], *, cwd: Path, stdout_path: Optional[Path]) -> None:
    env = os.environ.copy()
    src_path = _repo_root() / "src"
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(src_path) + (os.pathsep + existing if existing else "")

    if stdout_path is not None:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
        with stdout_path.open("w", encoding="utf-8") as handle:
            proc = subprocess.run(
                [sys.executable, "-m", "helix.cli", *command],
                cwd=cwd,
                env=env,
                text=True,
                stdout=handle,
                stderr=subprocess.PIPE,
            )
    else:
        proc = subprocess.run(
            [sys.executable, "-m", "helix.cli", *command],
            cwd=cwd,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    if proc.returncode != 0:
        raise RuntimeError(f"Workflow step failed ({' '.join(command)}): {proc.stderr}")


def _extract_schema_info(step_cfg: Mapping[str, Any], wf_dir: Path) -> tuple[Optional[str], Optional[str], Optional[str]]:
    schema_cfg = step_cfg.get("schema") or {}
    if not schema_cfg:
        return None, None, None
    artifact_rel = schema_cfg.get("output")
    if not artifact_rel:
        return schema_cfg.get("kind"), None, None
    artifact_path = wf_dir / str(artifact_rel)
    if not artifact_path.exists():
        return schema_cfg.get("kind"), None, None
    payload = json.loads(artifact_path.read_text(encoding="utf-8"))
    meta = payload.get("meta", {})
    kind = schema_cfg.get("kind") or meta.get("schema", {}).get("kind")
    version = meta.get("spec_version")
    sha256 = hashlib.sha256(artifact_path.read_bytes()).hexdigest()
    return kind, version, sha256


def run_workflow_config(
    config: Path,
    output_dir: Path,
    selected: Optional[str] = None,
) -> List[WorkflowRunResult]:
    config_path = Path(config)
    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    workflows_cfg = cfg.get("workflows") or []
    results: List[WorkflowRunResult] = []

    for workflow_cfg in workflows_cfg:
        name = str(workflow_cfg.get("name") or "workflow")
        if selected and selected != name:
            continue
        wf_dir = Path(output_dir) / name
        wf_dir.mkdir(parents=True, exist_ok=True)
        step_results: List[WorkflowStepResult] = []

        for step_cfg in workflow_cfg.get("steps") or []:
            command_entry = step_cfg.get("command")
            if isinstance(command_entry, str):
                cli_tokens: List[str] = [command_entry]
            else:
                cli_tokens = [str(token) for token in (command_entry or [])]

            args_cfg = step_cfg.get("args") or {}
            for key, value in args_cfg.items():
                flag = f"--{str(key).replace('_', '-')}"
                cli_tokens.append(flag)
                cli_tokens.append(str(value))

            stdout_rel = step_cfg.get("stdout")
            stdout_path = wf_dir / stdout_rel if stdout_rel else None

            _run_cli(cli_tokens, cwd=wf_dir, stdout_path=stdout_path)

            schema_kind, schema_version, schema_hash = _extract_schema_info(step_cfg, wf_dir)
            step_results.append(
                WorkflowStepResult(
                    command=" ".join(cli_tokens),
                    output_path=stdout_path,
                    schema_kind=schema_kind,
                    schema_version=schema_version,
                    schema_hash=schema_hash,
                )
            )

        results.append(WorkflowRunResult(name=name, output_dir=wf_dir, steps=step_results))

    return results

