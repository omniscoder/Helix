#!/usr/bin/env python3
"""Generate the PCR micro-universe DAG and lean-check payload."""

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

    genome = root / "tests/veribiota/pcr_micro.fna"
    primer_cfg = root / "tests/veribiota/pcr_primers.json"
    pcr_cfg = root / "tests/veribiota/pcr_config.json"

    dag_path = workdir / "pcr_micro.dag.json"
    check_path = workdir / "pcr_micro.lean-check.json"

    _run_cli(
        "pcr",
        "dag",
        "--genome",
        str(genome),
        "--primer-config",
        str(primer_cfg),
        "--pcr-config",
        str(pcr_cfg),
        "--out",
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
    print(f"[veribiota] Generated PCR micro DAG → {dag_path}")
    print(f"[veribiota] Wrote lean-check summary → {check_path}")


if __name__ == "__main__":
    main()
