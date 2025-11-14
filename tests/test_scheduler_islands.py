from helix.core.graph import GraphIR
from helix.core.scheduler import Island, LiveScheduler
from helix.live import StateReducer
from helix.nodes import FieldNode, ObserverNode, RuleNetNode


def test_scheduler_runs_multi_rate_islands():
    graph = GraphIR()
    rule = RuleNetNode("rule", decay=0.1, initial=0.0)
    field = FieldNode("field", diffusion=0.4)
    observer = ObserverNode("obs")

    for node in [rule, field, observer]:
        graph.add_node(node)

    graph.connect("rule", "rate", "field", "signal")
    graph.connect("field", "tile", "obs", "signal")

    islands = [
        Island("ode", [rule], dt=0.1),
        Island("rd", [field, observer], dt=0.2),
    ]
    reducer = StateReducer(hz=60)

    scheduler = LiveScheduler(
        graph=graph,
        islands=islands,
        sync_dt=0.5,
        state_reducer=reducer,
    )

    scheduler.update_external("rule", {"drive": 1.0})
    scheduler.run_until(1.0)

    assert rule.state["value"] > 0.0
    assert field.state["tile"] > 0.0
    assert reducer.history  # ensures snapshots emitted
