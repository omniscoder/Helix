#!/usr/bin/env python3
"""Generate the Prime Editing micro-universe DAG and lean-check payload."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _run_cli(*args: str) -> None:
    cmd = [sys.executable, "-m", "helix.cli", *args]
    subprocess.run(cmd, check=True, cwd=_repo_root())


def main() -> None:
    root = _repo_root()
    workdir = root / "veribiota_work"
    workdir.mkdir(parents=True, exist_ok=True)

    genome = root / "tests/veribiota/prime_micro.fna"
    peg_cfg = root / "tests/veribiota/prime_peg.json"
    editor_cfg = root / "tests/veribiota/prime_editor.json"

    dag_path = workdir / "prime_micro.dag.json"
    check_path = workdir / "prime_micro.lean-check.json"

    _run_cli(
        "prime",
        "dag",
        "--genome",
        str(genome),
        "--peg-config",
        str(peg_cfg),
        "--editor-config",
        str(editor_cfg),
        "--max-depth",
        "2",
        "--json",
        str(dag_path),
    )
    _run_cli(
        "veribiota",
        "lean-check",
        "--input",
        str(dag_path),
        "--out",
        str(check_path),
    )
    print(f"[veribiota] Generated Prime micro DAG → {dag_path}")
    print(f"[veribiota] Wrote lean-check summary → {check_path}")


if __name__ == "__main__":
    main()
