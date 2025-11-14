"""Helix LiveGraph orchestrator for the EGFR/GRB2 MVP slice."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, Tuple

import yaml

from helix.core.graph import GraphIR
from helix.core.scheduler import Island, LiveScheduler
from helix.live import StateReducer
from helix.live.metadata import describe_graph_nodes
from helix.live.hot_swap import HotSwapManager
from helix.nodes import ABMNode, CouplerNode, FieldNode, ObserverNode, SBMLNode

CONFIG_PATH = Path("models/egfr_grb2/config.yaml")


def _load_yaml(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _load_params_file(path: Path) -> Dict[str, float]:
    if not path:
        return {}
    if path.exists():
        return _load_yaml(path) or {}
    return {}


def build_from_config(config_path: Path, variant: str = "wt") -> Tuple[GraphIR, LiveScheduler, Dict]:
    """Construct the Helix LiveGraph and scheduler for a given variant config."""

    config_path = Path(config_path)
    cfg = _load_yaml(config_path)
    if variant not in cfg["variants"]:
        raise ValueError(f"Unknown variant '{variant}'. Available: {list(cfg['variants'])}")
    variant_cfg = cfg["variants"][variant]
    params = _load_params_file(Path(variant_cfg.get("params", "")))

    sim_cfg = cfg.get("sim", {})
    rd_cfg = sim_cfg.get("rd", {})
    abm_cfg = sim_cfg.get("abm", {})
    size = abm_cfg.get("size", [32, 32])
    field = FieldNode(
        name="egf_field",
        diffusion=float(rd_cfg.get("D", 0.05)),
        shape=size,
        baseline=float(rd_cfg.get("baseline", 0.0)),
        stripe_value=float(rd_cfg.get("stripe_value", 1.0)),
        stripe_span=rd_cfg.get("stripe_span", 0.25),
    )

    output_name = sim_cfg.get("sbml_output", "pERK")
    sbml = SBMLNode(
        name="egfr_rules",
        sbml_path=variant_cfg["sbml"],
        outputs={output_name: "out"},
        amici_module=variant_cfg.get("amici_module"),
        parameters=params,
        initial_pERK=0.05,
    )

    coupler = CouplerNode(name="perk_to_behavior")
    if isinstance(size, (list, tuple)) and len(size) == 2:
        inferred_agents = int(size[0] * size[1])
    else:
        inferred_agents = 100
    agents = int(abm_cfg.get("agents", inferred_agents))
    cells = ABMNode(name="epithelium", agents=agents, layout=size)
    perk_report = ObserverNode(name="perk_report")
    cell_report = ObserverNode(name="cell_report")

    graph = GraphIR()
    for node in (field, sbml, coupler, cells, perk_report, cell_report):
        graph.add_node(node)

    output_port = output_name
    graph.connect("egf_field", "tile", "egfr_rules", "ligand")
    graph.connect("egfr_rules", output_port, "perk_to_behavior", "input")
    graph.connect("perk_to_behavior", "output", "epithelium", "field_signal")
    graph.connect("egf_field", "map", "epithelium", "field_map")
    graph.connect("epithelium", "agent_delta", "egf_field", "signal")
    graph.connect("egfr_rules", output_port, "perk_report", "signal")
    graph.connect("epithelium", "agents", "cell_report", "signal")

    node_meta = describe_graph_nodes(graph)
    reducer = StateReducer(hz=30, node_meta=node_meta)
    sync_dt = float(sim_cfg.get("sync_dt", 1.0))
    islands = [
        Island("rd", [field], dt=0.5),
        Island("pathway", [sbml], dt=0.2),
        Island("coupler", [coupler], dt=0.2),
        Island("cells", [cells], dt=1.0),
        Island("observers", [perk_report, cell_report], dt=1.0),
    ]
    hot_swap_manager = HotSwapManager()
    scheduler = LiveScheduler(
        graph=graph,
        islands=islands,
        sync_dt=sync_dt,
        state_reducer=reducer,
        hot_swap_manager=hot_swap_manager,
    )

    def _apply_variant(payload: Dict[str, object]) -> None:
        variant_name = payload.get("variant")
        if not isinstance(variant_name, str):
            return
        if variant_name not in cfg["variants"]:
            return
        target = cfg["variants"][variant_name]
        params_override = _load_params_file(Path(target.get("params", "")))
        sbml.params.update({k: float(v) for k, v in params_override.items()})
        scheduler.runtime_meta["variant"] = variant_name

    hot_swap_manager.register("variant", _apply_variant)

    stripe_value = float(rd_cfg.get("stripe_value", 1.0))
    scheduler.update_external("egf_field", {"control": stripe_value})
    return graph, scheduler, cfg


def build_egfr_grb2(config_path: Path = CONFIG_PATH, variant: str = "wt") -> Tuple[GraphIR, LiveScheduler]:
    graph, scheduler, _cfg = build_from_config(config_path, variant)
    return graph, scheduler


def run_variant(variant: str, duration: float) -> LiveScheduler:
    """Run a variant and return the populated scheduler (for inspection/tests)."""

    _, scheduler, cfg = build_from_config(CONFIG_PATH, variant=variant)
    scheduler.run_until(duration)
    return scheduler


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the EGFR/GRB2 MVP orchestrator.")
    parser.add_argument("--variant", choices=["wt", "grb2_ko"], default="wt")
    parser.add_argument("--duration", type=float, default=None, help="Simulation horizon in minutes (default from config).")
    args = parser.parse_args()

    cfg = _load_yaml(CONFIG_PATH)
    duration = args.duration or float(cfg["simulation"].get("duration", 60.0))
    scheduler = run_variant(args.variant, duration)
    if scheduler.snapshots:
        perk = scheduler.snapshots[-1].get("egfr_rules", {}).get(cfg["simulation"].get("sbml_output", "pERK"), None)
        print(f"Completed {args.variant} run (T={duration} min). Final pERK={perk}")


if __name__ == "__main__":
    main()
