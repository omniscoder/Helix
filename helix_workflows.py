"""YAML-driven Helix workflow runner."""
from __future__ import annotations

import argparse
import contextlib
import io
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence

import helix_cli

try:  # pragma: no cover - optional dependency
    import yaml
except ImportError:  # pragma: no cover - handled at runtime
    yaml = None


INPUT_PATH_KEYS = {"input", "spectrum_file"}
OUTPUT_PATH_KEYS = {
    "dot_output",
    "json",
    "output",
    "frameshift_csv",
    "clusters_csv",
    "orfs_csv",
    "orfs_fasta",
    "orf_fasta",
    "orf_csv",
}


@dataclass
class StepResult:
    command: str
    argv: List[str]
    stdout: str
    output_path: Path | None


@dataclass
class WorkflowResult:
    name: str
    output_dir: Path
    steps: List[StepResult]


def _ensure_yaml_available() -> None:
    if yaml is None:
        raise ImportError("PyYAML is required for workflow configs. Install it with 'pip install pyyaml'.")


def _load_config(path: Path) -> dict:
    _ensure_yaml_available()
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _normalize_path(value: str, base_dir: Path) -> str:
    path = Path(value)
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    return str(path)


def _dict_to_cli(args_dict: Dict[str, object], input_base: Path, output_base: Path) -> List[str]:
    argv: List[str] = []
    for key, value in args_dict.items():
        flag = f"--{key.replace('_', '-')}"
        if isinstance(value, bool):
            if value:
                argv.append(flag)
            continue
        if value is None:
            continue
        if isinstance(value, (list, tuple)):
            for item in value:
                argv.extend([flag, str(item)])
            continue
        if key in INPUT_PATH_KEYS and isinstance(value, str):
            value = _normalize_path(value, input_base)
        elif key in OUTPUT_PATH_KEYS and isinstance(value, str):
            value = _normalize_path(value, output_base)
        argv.extend([flag, str(value)])
    return argv


def _run_step(parser: argparse.ArgumentParser, command: str, argv: Sequence[str]) -> str:
    args = parser.parse_args([command, *argv])
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(buffer):
        args.func(args)
    return buffer.getvalue()


def run_workflow_config(
    config_path: Path,
    *,
    output_dir: Path,
    selected: str | None = None,
) -> List[WorkflowResult]:
    config = _load_config(config_path)
    workflows = config.get("workflows", [])
    if selected:
        workflows = [wf for wf in workflows if wf.get("name") == selected]
        if not workflows:
            raise ValueError(f"No workflow named '{selected}' in {config_path}.")

    parser = helix_cli.build_parser()
    results: List[WorkflowResult] = []
    base_dir = config_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    for wf in workflows:
        name = wf.get("name")
        if not name:
            raise ValueError("Workflow missing 'name'.")
        wf_dir = output_dir / name
        wf_dir.mkdir(parents=True, exist_ok=True)

        steps_config = wf.get("steps", [])
        step_results: List[StepResult] = []
        for idx, step in enumerate(steps_config, start=1):
            command = step.get("command")
            if not command:
                raise ValueError(f"Workflow '{name}' step {idx} missing 'command'.")
            args_dict = step.get("args", {})
            if not isinstance(args_dict, dict):
                raise ValueError(f"Workflow '{name}' step {idx} args must be a mapping.")
            argv = _dict_to_cli(args_dict, base_dir, wf_dir)
            stdout_text = _run_step(parser, command, argv)

            stdout_file = step.get("stdout")
            output_path: Path | None = None
            if stdout_file:
                output_path = wf_dir / stdout_file
                output_path.parent.mkdir(parents=True, exist_ok=True)
                output_path.write_text(stdout_text, encoding="utf-8")

            log_path = wf_dir / f"step_{idx}_{command}.log"
            log_path.write_text(stdout_text, encoding="utf-8")
            step_results.append(StepResult(command=command, argv=list(argv), stdout=stdout_text, output_path=output_path))

        results.append(WorkflowResult(name=name, output_dir=wf_dir, steps=step_results))

    return results
