"""Multi-rate scheduler for solver islands."""

from __future__ import annotations

from time import perf_counter, sleep
from threading import Event
from typing import Callable, Dict, Iterable, List, Mapping, MutableMapping, Optional

from .cache import ContentAddressedCache
from .graph import GraphIR, GraphNode
from .invariants import assert_invariants

InputResolver = Callable[[str, float], Mapping[str, float]]
ExternalInputProvider = Callable[[str, float], Mapping[str, float]]


class Island:
    """Group of strongly-coupled nodes that share a solver."""

    def __init__(self, name: str, nodes: Iterable[GraphNode], dt: float, cache: Optional[ContentAddressedCache] = None):
        self.name = name
        self.nodes = list(nodes)
        self.dt = dt
        self.cache = cache or ContentAddressedCache()

    def step(self, t: float, input_resolver: InputResolver) -> Dict[str, Mapping[str, float]]:
        outputs: Dict[str, Mapping[str, float]] = {}
        for node in self.nodes:
            node_inputs = input_resolver(node.name, t)
            cached = self.cache.get(node.signature(), node_inputs, self.dt)
            if cached is not None:
                outputs[node.name] = cached
                continue
            result = node.step(t, self.dt, node_inputs)
            self.cache.put(node.signature(), node_inputs, self.dt, result)
            outputs[node.name] = result
        return outputs


class LiveScheduler:
    """Drives all islands while keeping a steady sync cadence."""

    def __init__(
        self,
        graph: GraphIR,
        islands: Iterable[Island],
        sync_dt: float,
        state_reducer,
        external_inputs: Optional[ExternalInputProvider] = None,
        event_bus=None,
        hot_swap_manager=None,
    ) -> None:
        self.graph = graph
        self.islands = list(islands)
        self.sync_dt = sync_dt
        self.state_reducer = state_reducer
        self.event_bus = event_bus
        self.hot_swap_manager = hot_swap_manager
        self.external_inputs = external_inputs or (lambda _node, _t: {})
        self.node_buffers: Dict[str, Mapping[str, float]] = {}
        self.island_buffers: Dict[str, Mapping[str, Mapping[str, float]]] = {}
        self._island_time: Dict[str, float] = {island.name: 0.0 for island in self.islands}
        self._snapshot: Dict[str, Mapping[str, float]] = {}
        self.snapshots: List[Dict[str, Mapping[str, float]]] = []
        self.runtime_meta: Dict[str, object] = {}
        self._pause_event: Event = Event()
        self._pause_event.set()

    def _gather_inputs(self, node_name: str, t: float) -> Dict[str, float]:
        values: Dict[str, float] = dict(self.external_inputs(node_name, t))
        for edge in self.graph.incoming_edges(node_name):
            src_node, src_port = edge.source
            dst_port = edge.target[1]
            upstream = self.node_buffers.get(src_node, {})
            if src_port in upstream:
                values[dst_port] = upstream[src_port]
        return values

    def run_until(self, T: float, wall_budget_ms: float = 5.0) -> None:
        """Advance the scheduler until simulation time reaches T."""

        t_sync = 0.0
        wall_anchor = perf_counter()
        while t_sync < T:
            self._pause_event.wait()
            sync_target = min(T, t_sync + self.sync_dt)
            for island in self.islands:
                island_time = self._island_time[island.name]
                while island_time < sync_target:
                    t_current = island_time
                    outputs = island.step(t_current, lambda node_name, tt=t_current: self._gather_inputs(node_name, tt))
                    self.node_buffers.update(outputs)
                    self.island_buffers[island.name] = outputs
                    island_time += island.dt
                self._island_time[island.name] = island_time
            t_sync = sync_target
            snapshot = self.collect_snapshot()
            self.snapshots.append(snapshot)
            if self.state_reducer:
                meta = {"time": t_sync, "runtime": self.runtime_meta}
                self.state_reducer.push(snapshot, meta=meta)
            if self.event_bus:
                self.event_bus.publish({"t": t_sync, "snapshot": snapshot})
            assert_invariants(self.graph, snapshot)
            elapsed_ms = (perf_counter() - wall_anchor) * 1e3
            if elapsed_ms < wall_budget_ms:
                sleep((wall_budget_ms - elapsed_ms) * 1e-3)
            wall_anchor = perf_counter()

    def collect_snapshot(self) -> Dict[str, Mapping[str, float]]:
        """Return the latest node outputs keyed by node name."""

        snapshot: Dict[str, Mapping[str, float]] = {}
        for node in self.graph.iter_nodes():
            if node.name in self.node_buffers:
                snapshot[node.name] = dict(self.node_buffers[node.name])
        self._snapshot = snapshot
        return snapshot.copy()

    def apply_hot_swap(self, scope: str, payload) -> None:
        if not self.hot_swap_manager:
            return
        self.hot_swap_manager.apply(scope, payload)

    def update_external(self, node_name: str, values: Mapping[str, float]) -> None:
        """Override the external input provider for a single node."""

        base_provider = self.external_inputs

        def provider(name: str, t: float) -> Mapping[str, float]:
            if name == node_name:
                return values
            return base_provider(name, t)

        self.external_inputs = provider

    def pause(self) -> None:
        self._pause_event.clear()

    def resume(self) -> None:
        self._pause_event.set()

    def toggle_pause(self) -> None:
        if self._pause_event.is_set():
            self._pause_event.clear()
        else:
            self._pause_event.set()

    def paused(self) -> bool:
        return not self._pause_event.is_set()
