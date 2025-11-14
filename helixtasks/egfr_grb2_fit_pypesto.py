"""Calibration scaffold for the EGFR/GRB2 slice using AMICI + pyPESTO."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import numpy as np
import yaml

CFG_PATH = Path("models/egfr_grb2/config.yaml")


def _load_csv(csv_path: Path) -> np.ndarray:
    return np.genfromtxt(csv_path, delimiter=",", names=True, dtype=None, encoding="utf-8")


def _split_by_condition(data: np.ndarray) -> Dict[str, np.ndarray]:
    if "condition" not in data.dtype.names:
        return {"default": data}
    groups: Dict[str, np.ndarray] = {}
    for cond in np.unique(data["condition"]):
        groups[str(cond)] = data[data["condition"] == cond]
    return groups


def fit_variant(variant: str, variant_cfg: Dict, data: np.ndarray) -> Dict[str, float]:
    try:
        import amici  # type: ignore
        import pypesto
        import pypesto.optimize
    except ImportError as exc:  # pragma: no cover
        raise SystemExit("Install amici + pypesto to run calibration (`pip install amici pypesto`).") from exc

    sbml_path = Path(variant_cfg["sbml"])
    importer = amici.SbmlImporter(str(sbml_path))
    model_name = Path(sbml_path).stem
    out_dir = Path("build/amici") / model_name
    model = importer.sbml2amici(model_name, str(out_dir), verbose=False)
    solver = model.getSolver()
    solver.setMaxSteps(10000)

    edatas = []
    grouped = _split_by_condition(data)
    for cond, subset in grouped.items():
        edata = amici.ExpData(model.get())
        edata.setTimepoints(subset["time_min"])
        edata.setObservedData(subset["pERK"])
        edatas.append(edata)

    objective = pypesto.AmiciObjective(model=model, solver=solver, edatas=edatas)
    problem = pypesto.Problem(
        objective=objective,
        lb=[-3, -3, -3, -3],
        ub=[3, 3, 3, 3],
        x_names=["basal_activation", "dephos_rate", "grb2_scale", "response_scale"],
    )
    optimizer = pypesto.ScipyOptimizer(method="L-BFGS-B")
    result = pypesto.minimize(problem, optimizer=optimizer, n_starts=5)
    best = result.optimize_result.list[0]
    return dict(zip(problem.x_names, 10 ** best.x))


def main() -> None:
    parser = argparse.ArgumentParser(description="Fit EGFR/GRB2 pathway parameters using pyPESTO.")
    parser.add_argument("--config", type=Path, default=CFG_PATH, help="Config YAML for EGFR/GRB2 slice.")
    args = parser.parse_args()
    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    data = _load_csv(Path(cfg["data"]["synthetic"]))
    for variant, spec in cfg["variants"].items():
        params = fit_variant(variant, spec, data)
        out_path = Path(spec["params"])
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(yaml.safe_dump(params), encoding="utf-8")
        print(f"[{variant}] wrote calibrated parameters to {out_path}")


if __name__ == "__main__":
    main()
