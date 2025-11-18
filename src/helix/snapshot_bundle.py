from __future__ import annotations

import json
import uuid
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from . import __version__ as HELIX_VERSION
from .schema import SPEC_VERSION


class SnapshotBundleError(Exception):
    """Raised when snapshot bundle operations fail."""


def _canonical_bytes(payload: Mapping[str, Any]) -> bytes:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def _sha256_bytes(data: bytes) -> str:
    import hashlib

    return hashlib.sha256(data).hexdigest()


def _load_json(path: Path) -> Mapping[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:  # pragma: no cover - json errors bubble up
        raise SnapshotBundleError(f"Failed to read JSON from {path}: {exc}") from exc


def _snapshot_id() -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"hxs_{stamp}_{uuid.uuid4().hex[:6]}"


def _normalize_run_id(raw_id: object, fallback: str) -> str:
    token = str(raw_id or fallback).strip() or fallback
    return token.replace(" ", "_")


@dataclass
class SnapshotBundle:
    path: Path

    def load_manifest(self) -> Mapping[str, Any]:
        with zipfile.ZipFile(self.path, "r") as bundle:
            try:
                with bundle.open("manifest.json") as handle:
                    return json.load(handle)
            except KeyError as exc:  # pragma: no cover - invalid bundle
                raise SnapshotBundleError("manifest.json missing from bundle") from exc

    def read_bytes(self, relative_path: str) -> bytes:
        with zipfile.ZipFile(self.path, "r") as bundle:
            try:
                return bundle.read(relative_path)
            except KeyError as exc:  # pragma: no cover - invalid path
                raise SnapshotBundleError(f"{relative_path} missing from bundle") from exc

    def read_json(self, relative_path: str) -> Mapping[str, Any]:
        data = self.read_bytes(relative_path)
        try:
            return json.loads(data.decode("utf-8"))
        except Exception as exc:  # pragma: no cover - json errors bubble up
            raise SnapshotBundleError(f"Invalid JSON payload at {relative_path}: {exc}") from exc


def select_run_entry(
    manifest: Mapping[str, Any], *, run_id: str | None = None, kind: str | None = None
) -> Mapping[str, Any]:
    runs = manifest.get("runs") or []
    if not runs:
        raise SnapshotBundleError("Snapshot bundle contains no runs.")
    normalized_kind = kind.upper() if isinstance(kind, str) else None
    if run_id:
        for entry in runs:
            if str(entry.get("id")) == run_id:
                if normalized_kind and str(entry.get("kind")).upper() != normalized_kind:
                    continue
                return entry
        raise SnapshotBundleError(f"Run '{run_id}' not found in manifest.")
    if normalized_kind:
        for entry in runs:
            if str(entry.get("kind")).upper() == normalized_kind:
                return entry
        raise SnapshotBundleError(f"No run matching kind '{kind}' in manifest.")
    if len(runs) == 1:
        return runs[0]
    raise SnapshotBundleError("Multiple runs available; specify --run-id or --kind.")


def pack_snapshot_bundle(
    *,
    session_path: Path,
    run_paths: Sequence[Path],
    out_path: Path,
    extra_assets: Sequence[Path] | None = None,
) -> Mapping[str, Any]:
    if not run_paths:
        raise SnapshotBundleError("At least one run snapshot is required to build a bundle.")

    files_to_write: list[tuple[str, bytes]] = []

    def add_file(arcname: str, data: bytes) -> dict[str, Any]:
        files_to_write.append((arcname, data))
        return {"path": arcname, "sha256": _sha256_bytes(data), "size": len(data)}

    manifest: dict[str, Any] = {
        "snapshot_spec": "1.0.0",
        "snapshot_id": _snapshot_id(),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "studio_version": HELIX_VERSION,
        "cli_version": HELIX_VERSION,
        "schema_version": f"helix.schema/{SPEC_VERSION}",
        "runs": [],
        "compare": [],
        "reports": [],
        "assets": [],
    }

    session_bytes = Path(session_path).read_bytes()
    session_entry = add_file("session/session.json", session_bytes)
    manifest["assets"].append({"path": session_entry["path"], "sha256": session_entry["sha256"], "kind": "session", "role": "session"})

    for run_idx, run_path in enumerate(run_paths, start=1):
        normalized_path = run_path
        if normalized_path.is_dir():
            candidate = normalized_path / "snapshot.json"
            if candidate.exists():
                normalized_path = candidate
        payload = _load_json(normalized_path)
        state = payload.get("state") if isinstance(payload, Mapping) else {}
        run_id = _normalize_run_id(state.get("run_id"), f"run_{run_idx}")
        kind = str((state or {}).get("run_kind", "UNKNOWN")).upper()
        rel_path = f"runs/{run_id}/snapshot.json"
        entry = add_file(rel_path, json.dumps(payload, indent=2).encode("utf-8"))
        manifest["runs"].append(
            {
                "id": run_id,
                "kind": kind,
                "engine": {"name": "helix.session", "version": HELIX_VERSION},
                "params": entry,
                "artifacts": [],
                "status": {"outcome": "success"},
            }
        )

    for asset_path in extra_assets or []:
        data = Path(asset_path).read_bytes()
        sha = _sha256_bytes(data)
        ext = asset_path.suffix
        rel = f"assets/{sha}{ext}"
        ref = add_file(rel, data)
        manifest["assets"].append(
            {
                "path": ref["path"],
                "sha256": ref["sha256"],
                "kind": "asset",
                "role": asset_path.stem,
            }
        )

    manifest["runs"].sort(key=lambda item: str(item.get("id")))
    manifest["assets"].sort(key=lambda item: item.get("path", ""))

    manifest_copy = dict(manifest)
    manifest_copy.pop("manifest_sha256", None)
    manifest_sha = _sha256_bytes(_canonical_bytes(manifest_copy))
    manifest["manifest_sha256"] = manifest_sha
    files_to_write.append(("manifest.json", _canonical_bytes(manifest)))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(out_path, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
        for arcname, data in files_to_write:
            bundle.writestr(arcname, data)
    return manifest


__all__ = [
    "SnapshotBundle",
    "SnapshotBundleError",
    "pack_snapshot_bundle",
    "select_run_entry",
]
