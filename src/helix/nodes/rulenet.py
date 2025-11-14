"""RuleNet nodes wrap ODE/SSA style reactions."""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

from .base import BaseLiveNode, port

try:  # pragma: no cover - optional dependency
    import amici  # type: ignore
except ImportError:  # pragma: no cover
    amici = None  # type: ignore


class RuleNetNode(BaseLiveNode):
    def __init__(self, name: str, decay: float = 0.25, initial: float = 0.0) -> None:
        ports = {
            "drive": port("drive", "in"),
            "rate": port("rate", "out"),
        }
        super().__init__(name=name, kind="rulenet", ports=ports, state={"value": initial})
        self.decay = decay
        self.metadata["hash"] = f"rulenet:{name}:{decay}:{initial}"

    def step(self, t: float, dt: float, inputs: Mapping[str, float]):  # type: ignore[override]
        value = self.get_state("value")
        drive = inputs.get("drive", 0.0)
        dv = drive - self.decay * value
        value = value + dv * dt
        self.set_state("value", value)
        return {"rate": value}


class SBMLNode(BaseLiveNode):
    """Lightweight SBML-backed node with an optional AMICI backend."""

    DEFAULT_PARAMS: Dict[str, float] = {
        "basal_activation": 0.02,
        "dephos_rate": 0.3,
        "grb2_scale": 1.0,
        "response_scale": 1.0,
        "km": 0.5,
    }

    def __init__(
        self,
        name: str,
        sbml_path: str,
        outputs: Optional[Mapping[str, str]] = None,
        amici_module: Optional[str] = None,
        parameters: Optional[Mapping[str, float]] = None,
        initial_pERK: float = 0.05,
    ) -> None:
        port_defs = {"ligand": port("ligand", "in")}
        self.output_names = list(outputs.keys()) if outputs else ["pERK"]
        for out_name in self.output_names:
            port_defs[out_name] = port(out_name, "out")
        super().__init__(name=name, kind="sbml", ports=port_defs, state={"pERK": initial_pERK})
        self.sbml_path = str(sbml_path)
        self.params = dict(self.DEFAULT_PARAMS)
        if parameters:
            self.params.update({k: float(v) for k, v in parameters.items()})
        self.metadata["hash"] = f"sbml:{Path(sbml_path).name}"
        self._backend = "simple"
        self._amici_model = None
        self._amici_solver = None
        if amici_module and amici:
            try:  # pragma: no cover - heavy optional dependency
                module = importlib.import_module(amici_module)
                model: Any = module.getModel()  # type: ignore[attr-defined]
                solver = model.getSolver()
                self._amici_model = model
                self._amici_solver = solver
                self._backend = "amici"
            except Exception:  # pragma: no cover - fallback on failure
                self._backend = "simple"

    def step(self, t: float, dt: float, inputs: Mapping[str, float]):  # type: ignore[override]
        ligand = max(inputs.get("ligand", 0.0), 0.0)
        if self._backend == "amici":
            return self._amici_step(t, dt, ligand)
        return self._simple_step(dt, ligand)

    def _simple_step(self, dt: float, ligand: float) -> Dict[str, float]:
        perk = self.get_state("pERK")
        basal = self.params["basal_activation"]
        km = self.params["km"]
        scale = self.params["grb2_scale"]
        response = self.params["response_scale"]
        activation = basal + response * (scale * ligand / (km + ligand + 1e-9))
        perk += dt * (activation - self.params["dephos_rate"] * perk)
        perk = max(perk, 0.0)
        self.set_state("pERK", perk)
        outputs = {name: perk for name in self.output_names}
        return outputs

    def _amici_step(self, t: float, dt: float, ligand: float) -> Dict[str, float]:
        model = self._amici_model
        solver = self._amici_solver
        if model is None or solver is None:
            return self._simple_step(dt, ligand)
        # Update model parameters for ligand-dependent input.
        if hasattr(model, "setParameterById"):  # pragma: no cover
            try:
                model.setParameterById("ligand_input", ligand)
            except Exception:
                pass
        solver.setMaxSteps(1000)
        model.setTimepoints([t + dt])
        r = amici.runAmiciSimulation(model, solver)  # pragma: no cover
        perk_idx = None
        for idx, name in enumerate(model.getStateIds()):
            if name == "pERK":
                perk_idx = idx
                break
        if perk_idx is None:
            return self._simple_step(dt, ligand)
        perk = float(r.x[0][perk_idx])
        self.set_state("pERK", perk)
        return {name: perk for name in self.output_names}
