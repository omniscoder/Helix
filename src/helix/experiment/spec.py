"""
Experiment specification helpers (YAML → structured config).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Union

try:  # pragma: no cover - optional dependency
    import yaml
except Exception:  # pragma: no cover
    yaml = None


class ExperimentSpecError(ValueError):
    """Raised when an experiment config is invalid."""


def _ensure_yaml_available() -> None:
    if yaml is None:
        raise ExperimentSpecError(
            "PyYAML is required for experiment configs. Install it with 'pip install pyyaml'."
        )


def _resolve_path(base: Path, value: str) -> Path:
    path = Path(value.strip())
    if not path.is_absolute():
        path = (base / path).resolve()
    return path


@dataclass(frozen=True)
class GenomeSpec:
    fasta: Path
    region: Optional[str]


@dataclass(frozen=True)
class CrisprSimulationSpec:
    max_depth: int = 2
    min_prob: float = 1e-4
    max_sites: int = 20
    seed: int = 0
    use_gpu: bool = False


@dataclass(frozen=True)
class PrimeSimulationSpec:
    max_depth: int = 2
    min_prob: float = 1e-4
    seed: int = 0


@dataclass(frozen=True)
class CrisprExperimentSpec:
    kind: str
    name: str
    description: Optional[str]
    genome: GenomeSpec
    cas_config: Path
    guide_sequence: str
    guide_name: Optional[str]
    guide_pam: Optional[str]
    simulation: CrisprSimulationSpec


@dataclass(frozen=True)
class PrimeExperimentSpec:
    kind: str
    name: str
    description: Optional[str]
    genome: GenomeSpec
    editor_config: Path
    peg_config: Optional[Path]
    peg_inline: Optional[Dict[str, str]]
    simulation: PrimeSimulationSpec


ExperimentSpec = Union[CrisprExperimentSpec, PrimeExperimentSpec]


def load_experiment_spec(path: Path) -> ExperimentSpec:
    _ensure_yaml_available()
    cfg_path = Path(path)
    if not cfg_path.exists():
        raise ExperimentSpecError(f"Experiment config '{cfg_path}' not found.")
    data = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ExperimentSpecError("Experiment config must be a YAML mapping.")
    kind = str(data.get("kind", "")).strip()
    if kind == "helix.crispr.experiment.v1":
        return _parse_crispr_spec(cfg_path, data)
    if kind == "helix.prime.experiment.v1":
        return _parse_prime_spec(cfg_path, data)
    raise ExperimentSpecError(
        "Unknown experiment kind. Supported kinds: 'helix.crispr.experiment.v1', 'helix.prime.experiment.v1'."
    )


def _parse_genome_spec(cfg_path: Path, data: Dict[str, Any]) -> GenomeSpec:
    genome = data.get("genome")
    if not isinstance(genome, dict):
        raise ExperimentSpecError("Experiment config requires a 'genome' section.")
    fasta = genome.get("fasta")
    if not fasta:
        raise ExperimentSpecError("Genome section requires a 'fasta' path.")
    fasta_path = _resolve_path(cfg_path.parent, str(fasta))
    if not fasta_path.exists():
        raise ExperimentSpecError(f"Genome FASTA '{fasta_path}' not found.")
    region = genome.get("region")
    region_value = str(region).strip() if region else None
    return GenomeSpec(fasta=fasta_path, region=region_value or None)


def _parse_crispr_spec(cfg_path: Path, data: Dict[str, Any]) -> CrisprExperimentSpec:
    genome = _parse_genome_spec(cfg_path, data)
    cas = data.get("cas")
    if not isinstance(cas, dict) or not cas.get("config"):
        raise ExperimentSpecError("CRISPR experiments require 'cas.config'.")
    cas_config = _resolve_path(cfg_path.parent, str(cas["config"]))
    if not cas_config.exists():
        raise ExperimentSpecError(f"Cas config '{cas_config}' not found.")

    guide = data.get("guide")
    if not isinstance(guide, dict) or not guide.get("sequence"):
        raise ExperimentSpecError("CRISPR experiments require a 'guide.sequence'.")
    guide_sequence = str(guide["sequence"]).strip().upper()
    if not guide_sequence:
        raise ExperimentSpecError("Guide sequence must be non-empty.")
    guide_name = str(guide.get("name")).strip() if guide.get("name") else None
    guide_pam = str(guide.get("pam")).strip() if guide.get("pam") else None

    simulation = _parse_crispr_simulation(data.get("simulation", {}))
    return CrisprExperimentSpec(
        kind="helix.crispr.experiment.v1",
        name=str(data.get("name") or "unnamed_crispr_experiment"),
        description=str(data.get("description") or "").strip() or None,
        genome=genome,
        cas_config=cas_config,
        guide_sequence=guide_sequence,
        guide_name=guide_name,
        guide_pam=guide_pam,
        simulation=simulation,
    )


def _parse_crispr_simulation(section: Dict[str, Any]) -> CrisprSimulationSpec:
    if not isinstance(section, dict):
        section = {}
    return CrisprSimulationSpec(
        max_depth=int(section.get("max_depth", 2)),
        min_prob=float(section.get("min_prob", 1e-4)),
        max_sites=int(section.get("max_sites", 20)),
        seed=int(section.get("seed", 0)),
        use_gpu=bool(section.get("use_gpu", False)),
    )


def _parse_prime_spec(cfg_path: Path, data: Dict[str, Any]) -> PrimeExperimentSpec:
    genome = _parse_genome_spec(cfg_path, data)
    editor = data.get("editor")
    if not isinstance(editor, dict) or not editor.get("config"):
        raise ExperimentSpecError("Prime experiments require 'editor.config'.")
    editor_config = _resolve_path(cfg_path.parent, str(editor["config"]))
    if not editor_config.exists():
        raise ExperimentSpecError(f"Prime editor config '{editor_config}' not found.")

    peg_section = data.get("peg")
    peg_config: Optional[Path] = None
    peg_inline: Optional[Dict[str, str]] = None
    if isinstance(peg_section, dict):
        if peg_section.get("config"):
            peg_config = _resolve_path(cfg_path.parent, str(peg_section["config"]))
            if not peg_config.exists():
                raise ExperimentSpecError(f"peg config '{peg_config}' not found.")
        elif all(key in peg_section for key in ("spacer", "pbs", "rtt")):
            peg_inline = {
                "spacer": str(peg_section["spacer"]).strip().upper(),
                "pbs": str(peg_section["pbs"]).strip().upper(),
                "rtt": str(peg_section["rtt"]).strip().upper(),
            }
            if not all(peg_inline.values()):
                raise ExperimentSpecError("Inline peg definitions require non-empty spacer/pbs/rtt.")
            if peg_section.get("name"):
                peg_inline["name"] = str(peg_section["name"]).strip()
        else:
            raise ExperimentSpecError("Prime experiments require either 'peg.config' or inline spacer/pbs/rtt.")
    else:
        raise ExperimentSpecError("Prime experiments require a 'peg' section.")

    simulation = _parse_prime_simulation(data.get("simulation", {}))
    return PrimeExperimentSpec(
        kind="helix.prime.experiment.v1",
        name=str(data.get("name") or "unnamed_prime_experiment"),
        description=str(data.get("description") or "").strip() or None,
        genome=genome,
        editor_config=editor_config,
        peg_config=peg_config,
        peg_inline=peg_inline,
        simulation=simulation,
    )


def _parse_prime_simulation(section: Dict[str, Any]) -> PrimeSimulationSpec:
    if not isinstance(section, dict):
        section = {}
    return PrimeSimulationSpec(
        max_depth=int(section.get("max_depth", 2)),
        min_prob=float(section.get("min_prob", 1e-4)),
        seed=int(section.get("seed", 0)),
    )
