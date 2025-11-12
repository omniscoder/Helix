"""PAM definitions and simple validators."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Mapping, Sequence

_IUPAC: Mapping[str, Sequence[str]] = {
    "A": ("A",),
    "C": ("C",),
    "G": ("G",),
    "T": ("T", "U"),
    "U": ("T", "U"),
    "R": ("A", "G"),
    "Y": ("C", "T", "U"),
    "S": ("G", "C"),
    "W": ("A", "T", "U"),
    "K": ("G", "T", "U"),
    "M": ("A", "C"),
    "B": ("C", "G", "T", "U"),
    "D": ("A", "G", "T", "U"),
    "H": ("A", "C", "T", "U"),
    "V": ("A", "C", "G"),
    "N": ("A", "C", "G", "T", "U"),
}


@dataclass(frozen=True)
class PAM:
    """Container describing a PAM pattern."""

    name: str
    pattern: str
    orientation: str = "3prime"  # guide-3' adjacency (SpCas9 default)
    notes: str | None = None

    def as_dict(self) -> Dict[str, str]:
        return {"name": self.name, "pattern": self.pattern, "orientation": self.orientation, "notes": self.notes}


_DEFAULT_PAMS: Dict[str, PAM] = {
    "SpCas9-NGG": PAM(name="SpCas9-NGG", pattern="NGG", orientation="3prime", notes="canonical SpCas9"),
    "SaCas9-NNGRRT": PAM(name="SaCas9-NNGRRT", pattern="NNGRRT", orientation="3prime"),
    "SpRY-NNN": PAM(name="SpRY-NNN", pattern="NNN", orientation="3prime", notes="near-PAMless"),
}


def list_pams() -> Sequence[str]:
    """Return the registered PAM names."""

    return tuple(sorted(_DEFAULT_PAMS))


def get_pam(name: str) -> Dict[str, str]:
    """Return a PAM definition by name."""

    try:
        pam = _DEFAULT_PAMS[name]
    except KeyError as exc:  # pragma: no cover - defensive
        raise ValueError(f"Unknown PAM '{name}'. Known PAMs: {', '.join(list_pams())}.") from exc
    return pam.as_dict()


def register_pam(pam: PAM) -> None:
    """Register a PAM at runtime (useful for notebooks/tests)."""

    _DEFAULT_PAMS[pam.name] = pam


def match_pam(seq: str, pam: Mapping[str, str] | PAM, pos: int = 0) -> bool:
    """Return True if `seq[pos:pos+len(pattern)]` satisfies the PAM rule."""

    if isinstance(pam, PAM):
        pattern = pam.pattern
    else:
        pattern = str(pam.get("pattern", ""))
    if not pattern:
        return False
    region = seq[pos : pos + len(pattern)].upper()
    if len(region) != len(pattern):
        return False
    for base, symbol in zip(region, pattern.upper()):
        allowed = _IUPAC.get(symbol)
        if not allowed or base.upper() not in allowed:
            return False
    return True


_IUPAC_COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "U": "A",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "D": "H",
    "H": "D",
    "V": "B",
    "N": "N",
}


def reverse_complement_pattern(pattern: str) -> str:
    """Return the reverse complement for an IUPAC pattern string."""

    letters = []
    for char in reversed(pattern.upper()):
        letters.append(_IUPAC_COMPLEMENT.get(char, "N"))
    return "".join(letters)
