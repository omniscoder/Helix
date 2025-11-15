"""Interactive LiveLab helpers (LiveSession + CLI shell)."""

from __future__ import annotations

import shlex
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence

from ..core.graph import Edge, GraphIR
from ..core.scheduler import Island, LiveScheduler
from ..dsl import HGXModel, load_hgx, dump_hgx, build_graph_from_hgx
from ..nodes import ABMNode, CouplerNode, FieldNode, ObserverNode
from ..live import StateReducer
from ..live.cli_utils import (
    finalize_live_run,
    maybe_create_realtime_queue,
    prepare_live_bundle_dir,
    print_plan_summary,
)

_NODE_FACTORIES = {
    "field": FieldNode,
    "abm": ABMNode,
    "observer": ObserverNode,
    "coupler": CouplerNode,
}


def _coerce_value(text: str) -> Any:
    lowered = text.lower()
    if lowered in {"true", "false"}:
        return lowered == "true"
    try:
        if "." in text:
            return float(text)
        return int(text)
    except ValueError:
        return text


class LiveSession:
    """Mutable graph + scheduler builder."""

    def __init__(self, name: str = "LiveSession", sync_dt: float = 0.25, default_dt: float = 0.1, seed: int = 0):
        self.name = name
        self.sync_dt = sync_dt
        self.default_dt = default_dt
        self.seed = seed
        self.graph = GraphIR()
        self.dt_overrides: Dict[str, float] = {}
        self._node_specs: Dict[str, Dict[str, Any]] = {}
        self._scheduler: Optional[LiveScheduler] = None
        self._dirty = True

    @classmethod
    def empty(cls, **kwargs) -> "LiveSession":
        return cls(**kwargs)

    @classmethod
    def from_hgx(cls, path: Path, *, sync_dt: Optional[float] = None, default_dt: Optional[float] = None) -> "LiveSession":
        model = load_hgx(path)
        session = cls(
            name=model.name,
            sync_dt=sync_dt if sync_dt is not None else 0.25,
            default_dt=default_dt if default_dt is not None else 0.1,
        )
        graph = build_graph_from_hgx(model)
        session.graph = graph
        for node in model.nodes:
            session._node_specs[node["name"]] = {"type": node["type"], "params": dict(node.get("params", {}))}
        session._dirty = True
        return session

    def list_nodes(self) -> List[str]:
        return sorted(node.name for node in self.graph.iter_nodes())

    def list_edges(self) -> List[Edge]:
        return sorted(self.graph.edges, key=lambda e: (e.source[0], e.target[0]))

    def add_node(self, name: str, kind: str, params: Optional[Mapping[str, Any]] = None) -> None:
        if name in self._node_specs:
            raise ValueError(f"Node {name} already exists.")
        cls = _NODE_FACTORIES.get(kind.lower())
        if not cls:
            raise ValueError(f"Unknown node kind '{kind}'. Available: {', '.join(sorted(_NODE_FACTORIES))}")
        params = dict(params or {})
        node = cls(name=name, **params)
        self.graph.add_node(node)
        self._node_specs[name] = {"type": type(node).__name__, "params": dict(params)}
        self._dirty = True

    def remove_node(self, name: str) -> None:
        if name not in self._node_specs:
            raise ValueError(f"Node {name} not found.")
        self._node_specs.pop(name, None)
        self.dt_overrides.pop(name, None)
        nodes = [node for node in self.graph.iter_nodes() if node.name != name]
        edges = [
            edge
            for edge in self.graph.edges
            if edge.source[0] != name and edge.target[0] != name
        ]
        new_graph = GraphIR()
        for node in nodes:
            new_graph.add_node(node)
        for edge in edges:
            new_graph.connect(edge.source[0], edge.source[1], edge.target[0], edge.target[1])
        self.graph = new_graph
        self._dirty = True

    def connect(self, src: str, src_port: str, dst: str, dst_port: str) -> None:
        self.graph.connect(src, src_port, dst, dst_port)
        self._dirty = True

    def disconnect(self, src: str, dst: str) -> None:
        edges = [
            edge for edge in self.graph.edges if not (edge.source[0] == src and edge.target[0] == dst)
        ]
        nodes = list(self.graph.iter_nodes())
        new_graph = GraphIR()
        for node in nodes:
            new_graph.add_node(node)
        for edge in edges:
            new_graph.connect(edge.source[0], edge.source[1], edge.target[0], edge.target[1])
        self.graph = new_graph
        self._dirty = True

    def set_param(self, node_name: str, param: str, value: Any) -> None:
        try:
            node = self.graph.node(node_name)
        except KeyError as exc:
            raise ValueError(f"Node {node_name} not found.") from exc
        setattr(node, param, value)
        spec = self._node_specs.setdefault(node_name, {"type": type(node).__name__, "params": {}})
        spec.setdefault("params", {})[param] = value
        self._dirty = True

    def set_dt_override(self, node_name: str, dt: float) -> None:
        if node_name not in self._node_specs:
            raise ValueError(f"Node {node_name} does not exist.")
        self.dt_overrides[node_name] = float(dt)
        self._dirty = True

    def plan(self) -> None:
        graph = self.graph
        components = graph.strongly_connected_components()
        if not components:
            components = [(node.name,) for node in graph.iter_nodes()]
        islands = []
        for idx, comp in enumerate(components):
            comp_dt = min(self.dt_overrides.get(name, self.default_dt) for name in comp)
            islands.append({"name": f"island_{idx}", "dt": comp_dt, "nodes": list(comp)})
        print_plan_summary(
            title=f"Live plan for {self.name}",
            model_name=self.name,
            sync_dt=self.sync_dt,
            default_dt=self.default_dt,
            overrides=self.dt_overrides,
            islands=islands,
            node_count=len(self._node_specs),
            edge_count=len(self.graph.edges),
        )

    def run(
        self,
        duration: float,
        *,
        realtime: bool = False,
        realtime_endpoint: Optional[str] = None,
        hz: float = 60.0,
        bundle: Optional[Path] = None,
    ) -> Path:
        scheduler = self._ensure_scheduler()
        realtime_hook = maybe_create_realtime_queue(realtime, realtime_endpoint)
        realtime_queue = realtime_hook.queue if realtime_hook else None
        reducer = StateReducer(hz=hz, realtime_queue=realtime_queue)
        scheduler.state_reducer = reducer
        scheduler.runtime_meta = {"model": self.name, "run_mode": "livelab"}
        if realtime_hook and realtime_hook.server:
            realtime_hook.server.set_command_handler(lambda payload: None)
        args = SimpleNamespace(
            bundle=bundle,
            seed=self.seed,
            default_dt=self.default_dt,
            wall_ms=5.0,
        )
        path = finalize_live_run(
            args,
            scheduler=scheduler,
            spec_path=Path(f"{self.name}.hgx"),
            model_name=self.name,
            duration=duration,
            sync_dt=self.sync_dt,
            dt_overrides=self.dt_overrides,
            const_inputs={},
            input_series_file=None,
            hz=hz,
            reducer=reducer,
            event_bus=None,
            realtime_hook=realtime_hook,
            meta_extra={"target_type": "livelab"},
        )
        return path

    def export_hgx(self, path: Path) -> None:
        nodes = []
        for node_name, spec in sorted(self._node_specs.items()):
            nodes.append(
                {
                    "name": node_name,
                    "type": spec["type"],
                    "params": dict(spec.get("params", {})),
                }
            )
        edges = []
        for edge in sorted(self.graph.edges, key=lambda e: (e.source[0], e.target[0])):
            edges.append(
                {
                    "source": {"node": edge.source[0], "port": edge.source[1]},
                    "target": {"node": edge.target[0], "port": edge.target[1]},
                }
            )
        dump_hgx(HGXModel(name=self.name, nodes=nodes, edges=edges), path)

    def _ensure_scheduler(self) -> LiveScheduler:
        if self._scheduler is None or self._dirty:
            self._scheduler = self._build_scheduler()
            self._dirty = False
        return self._scheduler

    def _build_scheduler(self) -> LiveScheduler:
        graph = self.graph
        components = graph.strongly_connected_components()
        if not components:
            components = [(node.name,) for node in graph.iter_nodes()]
        islands: List[Island] = []
        for idx, comp in enumerate(components):
            nodes = [graph.node(name) for name in comp]
            node_dts = [self.dt_overrides.get(name, self.default_dt) for name in comp]
            dt = min(node_dts) if node_dts else self.default_dt
            islands.append(Island(name=f"island_{idx}", nodes=nodes, dt=dt))
        return LiveScheduler(graph=graph, islands=islands, sync_dt=self.sync_dt, state_reducer=None)


@dataclass
class LiveCommand:
    raw: str


class LiveDevShell:
    """Simple REPL for LiveSession."""

    def __init__(self, session: LiveSession):
        self.session = session
        self.history: List[LiveCommand] = []

    def loop(self) -> None:
        print("Helix LiveLab – type :help for commands. Prefix commands with ':'.")
        while True:
            try:
                line = input("live> ").strip()
            except EOFError:
                print()
                break
            if not line:
                continue
            if not line.startswith(":"):
                print("Commands must start with ':'. See :help.")
                continue
            self.history.append(LiveCommand(raw=line))
            try:
                if not self._dispatch(line[1:].strip()):
                    break
            except Exception as exc:  # pragma: no cover - interactive shell
                print(f"error: {exc}")

    def _dispatch(self, payload: str) -> bool:
        if not payload:
            return True
        parts = shlex.split(payload)
        cmd = parts[0].lower()
        args = parts[1:]
        if cmd in {"quit", "exit"}:
            return False
        if cmd == "help":
            self._print_help()
        elif cmd == "list":
            self._handle_list(args)
        elif cmd == "plan":
            self.session.plan()
        elif cmd == "history":
            for entry in self.history:
                print(entry.raw)
        elif cmd == "add":
            self._handle_add(args)
        elif cmd in {"del", "remove"}:
            self._handle_remove(args)
        elif cmd == "connect":
            self._handle_connect(args)
        elif cmd in {"disconnect", "unlink"}:
            self._handle_disconnect(args)
        elif cmd == "set":
            self._handle_set(args)
        elif cmd == "dt":
            self._handle_dt(args)
        elif cmd == "run":
            self._handle_run(args)
        elif cmd == "export":
            self._handle_export(args)
        else:
            print(f"Unknown command '{cmd}'. Type :help for options.")
        return True

    def _print_help(self) -> None:
        print(
            """
Commands:
  :help                       Show this help.
  :list nodes|edges           List graph nodes or edges.
  :add node <name> kind=<kind> [param=value ...]
  :del node <name>            Remove a node.
  :connect src.port dst.port  Connect nodes (e.g., egf.tile cells.field_signal).
  :disconnect src dst         Remove edges between nodes.
  :set node.param value       Update a node parameter.
  :dt node=<name> value=<dt>  Override per-node dt.
  :plan                       Show SCC → island breakdown.
  :run [t=5] [realtime]       Run for duration (seconds); realtime streams deltas.
  :export hgx <path>          Write current graph to HGX.
  :export bundle <dir>        Run + export bundle into directory.
  :history                    Show recent commands.
  :quit                       Exit the shell.
"""
        )

    def _handle_list(self, args: Sequence[str]) -> None:
        if not args:
            print("Usage: :list nodes|edges")
            return
        what = args[0].lower()
        if what == "nodes":
            for node in self.session.list_nodes():
                spec = self.session._node_specs.get(node, {})
                print(f"- {node} ({spec.get('type')})")
        elif what == "edges":
            for edge in self.session.list_edges():
                print(f"- {edge.source[0]}.{edge.source[1]} -> {edge.target[0]}.{edge.target[1]}")
        else:
            print("Usage: :list nodes|edges")

    def _handle_add(self, args: Sequence[str]) -> None:
        if not args or args[0] != "node" or len(args) < 2:
            print("Usage: :add node <name> kind=<kind> [param=value ...]")
            return
        name = args[1]
        params = {}
        kind = None
        for token in args[2:]:
            if "=" not in token:
                continue
            key, value = token.split("=", 1)
            if key == "kind":
                kind = value
            else:
                params[key] = _coerce_value(value)
        if not kind:
            print("Provide kind=<Field|ABM|Observer|Coupler>.")
            return
        self.session.add_node(name, kind, params)
        print(f"Added node {name} ({kind}).")

    def _handle_remove(self, args: Sequence[str]) -> None:
        if len(args) < 2 or args[0] != "node":
            print("Usage: :del node <name>")
            return
        self.session.remove_node(args[1])
        print(f"Removed node {args[1]}.")

    def _handle_connect(self, args: Sequence[str]) -> None:
        if len(args) < 2:
            print("Usage: :connect src.port dst.port")
            return
        src = args[0]
        dst = args[1]
        if "." not in src or "." not in dst:
            print("Specify ports as node.port.")
            return
        src_node, src_port = src.split(".", 1)
        dst_node, dst_port = dst.split(".", 1)
        self.session.connect(src_node, src_port, dst_node, dst_port)
        print(f"Connected {src} -> {dst}.")

    def _handle_disconnect(self, args: Sequence[str]) -> None:
        if len(args) < 2:
            print("Usage: :disconnect <src> <dst>")
            return
        self.session.disconnect(args[0], args[1])
        print(f"Disconnected edges from {args[0]} to {args[1]}.")

    def _handle_set(self, args: Sequence[str]) -> None:
        if len(args) < 2 or "." not in args[0]:
            print("Usage: :set node.param value")
            return
        node, param = args[0].split(".", 1)
        value = _coerce_value(args[1])
        self.session.set_param(node, param, value)
        print(f"Set {node}.{param} = {value}")

    def _handle_dt(self, args: Sequence[str]) -> None:
        if not args:
            print("Usage: :dt node=<name> value=<dt>")
            return
        kv = dict(arg.split("=", 1) for arg in args if "=" in arg)
        node = kv.get("node")
        value = kv.get("value")
        if not node or value is None:
            print("Usage: :dt node=<name> value=<dt>")
            return
        self.session.set_dt_override(node, float(value))
        print(f"Set dt override for {node} = {value}")

    def _handle_run(self, args: Sequence[str]) -> None:
        duration = 5.0
        realtime = False
        bundle = None
        for token in args:
            if token == "realtime":
                realtime = True
            elif token.startswith("t="):
                duration = float(token.split("=", 1)[1])
            elif token.startswith("bundle="):
                bundle = Path(token.split("=", 1)[1])
        path = self.session.run(duration, realtime=realtime, bundle=bundle)
        print(f"Run complete. Bundle at {path}")

    def _handle_export(self, args: Sequence[str]) -> None:
        if len(args) < 2:
            print("Usage: :export hgx <path> | :export bundle <dir>")
            return
        kind, target = args[0], Path(args[1])
        if kind.lower() == "hgx":
            self.session.export_hgx(target)
            print(f"Saved HGX to {target}")
        elif kind.lower() == "bundle":
            path = self.session.run(duration=0.0, bundle=target)
            print(f"Bundle written to {path}")
        else:
            print("Usage: :export hgx <path> | :export bundle <dir>")
