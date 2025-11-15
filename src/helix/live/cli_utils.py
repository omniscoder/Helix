"""Helpers shared between CLI commands and LiveLab."""

from __future__ import annotations

import argparse
import json
import queue
import re
import shlex
import sys
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence

from .. import __version__ as HELIX_VERSION
from ..live import StateReducer
from ..live.realtime import RealtimeHook, RealtimeServer


def current_command_str() -> str:
    return " ".join(shlex.quote(arg) for arg in sys.argv)


def sanitize_run_name(name: str) -> str:
    cleaned = re.sub(r"[^a-zA-Z0-9_-]+", "-", name.strip().lower())
    cleaned = cleaned.strip("-")
    return cleaned or "livegraph"


def compose_run_id(model_name: str, variant: Optional[str], seed: Optional[int]) -> str:
    variant_token = (variant or "default").replace(" ", "_").lower()
    seed_token = f"seed{seed}" if seed is not None else "seed0"
    return f"{sanitize_run_name(model_name)}_{variant_token}_{seed_token}"


def prepare_live_bundle_dir(bundle: Optional[Path], model_name: str) -> Path:
    if bundle:
        target = Path(bundle)
    else:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        target = Path("live_runs") / f"{sanitize_run_name(model_name)}-{stamp}"
    target.mkdir(parents=True, exist_ok=True)
    return target


def maybe_create_realtime_queue(enabled: bool, endpoint: Optional[str] = None) -> Optional[RealtimeHook]:
    if not enabled:
        return None
    q: queue.Queue = queue.Queue(maxsize=512)
    server: Optional[RealtimeServer] = None
    try:
        server = RealtimeServer(endpoint or "tcp://127.0.0.1:8765")
        print(f"[realtime] serving feed on {server.label}")
    except Exception as exc:  # pragma: no cover - transport failures
        print(f"[realtime] failed to start server: {exc}", file=sys.stderr)
        server = None

    def _consumer():
        while True:
            payload = q.get()
            if payload is None:
                break
            runtime = payload.get("runtime") or {}
            slice_name = runtime.get("slice") or runtime.get("model") or "-"
            variant = runtime.get("variant") or "-"
            t = payload.get("time")
            try:
                t_display = f"{float(t):.2f}"
            except (TypeError, ValueError):
                t_display = "?"
            delta = payload.get("delta") or {}
            added = len(delta.get("added", {}))
            updated = len(delta.get("updated", {}))
            removed = len(delta.get("removed", {}))
            print(f"[realtime:{slice_name}/{variant}] t={t_display} +{added} ~{updated} -{removed}")
            if server:
                server.broadcast(dict(payload))

    threading.Thread(target=_consumer, daemon=True).start()
    return RealtimeHook(queue=q, server=server)


def dump_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _write_live_bundle(
    bundle_dir: Path,
    *,
    hgx_path: Path,
    model_name: str,
    scheduler,
    reducer: StateReducer,
    events: Sequence[Any],
    meta: Mapping[str, Any],
) -> None:
    model_contents = hgx_path.read_text(encoding="utf-8") if hgx_path.exists() else ""
    (bundle_dir / "model.hgx").write_text(model_contents, encoding="utf-8")
    dump_json(bundle_dir / "snapshots.json", scheduler.snapshots)
    dump_json(bundle_dir / "deltas.json", reducer.history)
    dump_json(bundle_dir / "events.json", list(events))
    enriched_meta = dict(meta)
    enriched_meta.update(
        {
            "model": model_name,
            "snapshots": len(scheduler.snapshots),
            "helix_version": HELIX_VERSION,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "command": current_command_str(),
        }
    )
    dump_json(bundle_dir / "meta.json", enriched_meta)


def finalize_live_run(
    args: argparse.Namespace,
    *,
    scheduler,
    spec_path: Path,
    model_name: str,
    duration: float,
    sync_dt: float,
    dt_overrides: Mapping[str, float],
    const_inputs: Mapping[str, Mapping[str, float]],
    input_series_file: Optional[str],
    hz: float,
    reducer: Optional[StateReducer],
    event_bus,
    realtime_hook: Optional[RealtimeHook] = None,
    meta_extra: Optional[Mapping[str, Any]] = None,
) -> Path:
    wall_start = datetime.now()
    scheduler.run_until(duration, wall_budget_ms=args.wall_ms)
    wall_time_sec = (datetime.now() - wall_start).total_seconds()
    islands_meta: List[Dict[str, Any]] = []
    for island in getattr(scheduler, "islands", []):
        nodes = [getattr(node, "name", None) for node in getattr(island, "nodes", []) if getattr(node, "name", None)]
        islands_meta.append({"name": getattr(island, "name", "island"), "nodes": nodes, "dt": getattr(island, "dt", None)})

    bundle_dir = prepare_live_bundle_dir(args.bundle, model_name)
    meta: Dict[str, Any] = {
        "kind": "helix.live.run.v1",
        "seed": args.seed,
        "duration": duration,
        "sync_dt": sync_dt,
        "default_dt": args.default_dt,
        "dt_overrides": dt_overrides,
        "input_constants": const_inputs,
        "input_series_file": input_series_file,
        "hz": hz,
        "wall_budget_ms": args.wall_ms,
        "wall_time_sec": wall_time_sec,
        "islands": islands_meta,
        "model": model_name,
        "spec_path": str(spec_path),
        "run_id": scheduler.runtime_meta.get("run_id") if scheduler.runtime_meta else None,
    }
    if meta_extra:
        meta.update(meta_extra)
    reducer_obj = reducer or StateReducer(hz=hz)
    events = event_bus.drain() if event_bus else []
    _write_live_bundle(
        bundle_dir,
        hgx_path=spec_path,
        model_name=model_name,
        scheduler=scheduler,
        reducer=reducer_obj,
        events=events,
        meta=meta,
    )
    print(f"Live run bundle written to {bundle_dir}")
    print(
        f"  snapshots={len(scheduler.snapshots)} islands={len(getattr(scheduler, 'islands', []))} "
        f"dt=[{', '.join(f'{isl.name}:{isl.dt}' for isl in getattr(scheduler, 'islands', []))}]"
    )
    if realtime_hook:
        realtime_hook.queue.put(None)
        if realtime_hook.server:
            realtime_hook.server.close()
    return bundle_dir


def print_plan_summary(
    *,
    title: str,
    model_name: str,
    sync_dt: float,
    default_dt: float,
    overrides: Mapping[str, float],
    islands: Sequence[Mapping[str, Any]],
    node_count: int,
    edge_count: int,
) -> None:
    print(title)
    print(f"  model: {model_name}")
    print(f"  sync_dt: {sync_dt}")
    print(f"  default_dt: {default_dt}")
    if overrides:
        print("  dt overrides:")
        for name, value in sorted(overrides.items()):
            print(f"    - {name} = {value}")
    else:
        print("  dt overrides: (none)")
    print(f"  graph: {node_count} nodes / {edge_count} edges")
    if not islands:
        print("  islands: (none)")
        return
    print("  islands:")
    for idx, island in enumerate(islands, 1):
        name = island.get("name") or f"island_{idx-1}"
        dt = island.get("dt")
        nodes = ", ".join(island.get("nodes") or [])
        print(f"    {idx}. {name} (dt={dt})")
        if nodes:
            print(f"       nodes: {nodes}")
