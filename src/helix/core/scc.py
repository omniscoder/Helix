"""Strongly connected component utilities."""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Set, Tuple


class Graph:
    """Lightweight directed graph wrapper used for SCC computations."""

    def __init__(self) -> None:
        self.adj: Dict[str, Set[str]] = defaultdict(set)

    def add_edge(self, u: str, v: str) -> None:
        self.adj[u].add(v)

    def vertices(self) -> Set[str]:
        verts = set(self.adj.keys())
        for targets in self.adj.values():
            verts.update(targets)
        return verts


def tarjan_scc(g: Graph) -> List[Tuple[str, ...]]:
    """Tarjan's SCC algorithm."""

    index: Dict[str, int] = {}
    lowlink: Dict[str, int] = {}
    onstack: Set[str] = set()
    stack: List[str] = []
    idx = 0
    sccs: List[Tuple[str, ...]] = []

    def strongconnect(v: str) -> None:
        nonlocal idx
        index[v] = idx
        lowlink[v] = idx
        idx += 1
        stack.append(v)
        onstack.add(v)

        for w in g.adj.get(v, set()):
            if w not in index:
                strongconnect(w)
                lowlink[v] = min(lowlink[v], lowlink[w])
            elif w in onstack:
                lowlink[v] = min(lowlink[v], index[w])

        if lowlink[v] == index[v]:
            comp: List[str] = []
            while True:
                w = stack.pop()
                onstack.remove(w)
                comp.append(w)
                if w == v:
                    break
            sccs.append(tuple(comp))

    for vertex in g.vertices():
        if vertex not in index:
            strongconnect(vertex)
    return sccs


def condensation_dag(g: Graph) -> Tuple[List[Tuple[str, ...]], Graph]:
    """Return SCCs and the condensation DAG."""

    comps = tarjan_scc(g)
    comp_index = {v: i for i, comp in enumerate(comps) for v in comp}
    dag = Graph()
    for u, targets in g.adj.items():
        for v in targets:
            cu, cv = comp_index[u], comp_index[v]
            if cu != cv:
                dag.add_edge(cu, cv)
    return comps, dag
