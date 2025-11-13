"""Hash utilities for deterministic edit DAGs."""
from __future__ import annotations

import hashlib
from typing import Iterable


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def hash_event(chrom: str, start: int, end: int, replacement: str, stage: str) -> str:
    payload = f"{chrom}|{start}|{end}|{replacement}|{stage}".encode("utf-8")
    return _sha256(payload)[:16]


def hash_node_id(parent_ids: Iterable[str], event_hash: str, stage: str, time_step: int | None = None) -> str:
    parents = "|".join(sorted(parent_ids)).encode("utf-8")
    time_bytes = str(time_step or 0).encode("utf-8")
    payload = b"N|" + parents + b"|" + event_hash.encode("utf-8") + b"|" + stage.encode("utf-8") + b"|" + time_bytes
    return "n_" + _sha256(payload)[:16]
