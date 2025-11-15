from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

from helix.cli import command_veribiota_export_suite


def _write_dummy_dag(tmp_path: Path, name: str) -> Path:
    path = tmp_path / f"{name}.json"
    payload = {
        "root_id": "root",
        "nodes": {
            "root": {
                "log_prob": 0.0,
                "metadata": {"stage": "root"},
                "parent_ids": [],
                "seq_hashes": {},
            }
        },
        "edges": [],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_export_suite_writes_into_target(tmp_path: Path) -> None:
    dag_path = _write_dummy_dag(tmp_path, "example")
    veri_root = tmp_path / "VeriBiota"
    module_rel = Path("Biosim/VeriBiota/Helix/CrisprMicro.lean")
    args = SimpleNamespace(
        inputs=[dag_path],
        dag_names=None,
        veribiota_root=veri_root,
        module_path=module_rel,
        module_name=None,
        lean_import="Biosim.VeriBiota.EditDAG",
        list_name="allDags",
        theorem_name=None,
        eval=False,
        skip_theorem=True,
    )
    command_veribiota_export_suite(args)
    target = veri_root / module_rel
    assert target.exists()
    text = target.read_text()
    assert "namespace Biosim" in text
    assert "def allDags : List EditDAG" in text
