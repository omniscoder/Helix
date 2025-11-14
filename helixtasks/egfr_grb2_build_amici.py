"""Utility script to generate AMICI models for the EGFR/GRB2 SBML files."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

CFG_PATH = Path("models/egfr_grb2/config.yaml")


def _module_to_dir(module_path: str) -> Path:
    parts = module_path.split(".")
    if len(parts) < 2:
        raise ValueError(f"Invalid AMICI module path '{module_path}'. Expected dotted module.")
    return Path(*parts[:-1])


def build_variant(name: str, sbml_path: Path, module_path: str) -> None:
    try:
        import amici  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise SystemExit("Install AMICI (`pip install amici`) to build models.") from exc

    importer = amici.SbmlImporter(str(sbml_path))
    out_dir = _module_to_dir(module_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    model_name = module_path.split(".")[-1]
    importer.sbml2amici(model_name, str(out_dir), verbose=False)
    print(f"[{name}] Generated AMICI module '{module_path}' from {sbml_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate AMICI modules for all EGFR/GRB2 variants.")
    parser.add_argument(
        "--config",
        type=Path,
        default=CFG_PATH,
        help="Path to the EGFR/GRB2 config.yaml",
    )
    args = parser.parse_args()
    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    for variant, spec in cfg.get("variants", {}).items():
        build_variant(
            variant,
            Path(spec["sbml"]),
            spec["amici_module"],
        )


if __name__ == "__main__":
    main()
