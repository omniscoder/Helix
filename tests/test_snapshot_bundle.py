from __future__ import annotations

import json
from pathlib import Path

import pytest

from helix.snapshot_bundle import SnapshotBundle, SnapshotBundleError, pack_snapshot_bundle, select_run_entry


def _write_json(path: Path, payload: dict) -> Path:
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_pack_snapshot_bundle_and_read(tmp_path: Path) -> None:
    session = {"version": 2, "state": {"run_kind": "CRISPR"}}
    session_path = _write_json(tmp_path / "session.json", session)
    run_payload = {
        "state": {
            "run_kind": "CRISPR",
            "run_id": 1,
            "config": {"guide": {"id": "G1"}},
            "outcomes": [],
        }
    }
    run_path = _write_json(tmp_path / "run1.json", run_payload)
    bundle_path = tmp_path / "bundle.hxs"

    manifest = pack_snapshot_bundle(session_path=session_path, run_paths=[run_path], out_path=bundle_path)

    assert bundle_path.exists()
    assert manifest["runs"][0]["id"] == "1"

    bundle = SnapshotBundle(bundle_path)
    manifest_from_disk = bundle.load_manifest()
    entry = select_run_entry(manifest_from_disk, run_id="1")
    snapshot = bundle.read_json(entry["params"]["path"])
    assert snapshot["state"]["run_id"] == 1


def test_select_run_entry_by_kind(tmp_path: Path) -> None:
    session_path = _write_json(tmp_path / "session.json", {"version": 2, "state": {}})
    run_a = _write_json(
        tmp_path / "a.json",
        {"state": {"run_kind": "CRISPR", "run_id": 10, "outcomes": []}},
    )
    run_b = _write_json(
        tmp_path / "b.json",
        {"state": {"run_kind": "PRIME", "run_id": 11, "outcomes": []}},
    )
    bundle_path = tmp_path / "bundle.hxs"
    manifest = pack_snapshot_bundle(session_path=session_path, run_paths=[run_a, run_b], out_path=bundle_path)

    entry = select_run_entry(manifest, kind="prime")
    assert entry["id"] == "11"

    with pytest.raises(SnapshotBundleError):
        select_run_entry(manifest)
