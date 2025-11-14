"""Content-addressed cache for solver islands."""

from __future__ import annotations

import hashlib
import json
import threading
from typing import Any, Dict, Mapping, Optional, Tuple


def _stable_payload(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


class ContentAddressedCache:
    """Simple thread-safe content addressed cache."""

    def __init__(self) -> None:
        self._store: Dict[str, Any] = {}
        self._lock = threading.Lock()

    def _key(self, signature: str, inputs: Mapping[str, Any], dt: float) -> str:
        digest = hashlib.sha256()
        digest.update(signature.encode("utf-8"))
        digest.update(str(dt).encode("utf-8"))
        digest.update(_stable_payload(inputs).encode("utf-8"))
        return digest.hexdigest()

    def get(self, signature: str, inputs: Mapping[str, Any], dt: float) -> Optional[Any]:
        key = self._key(signature, inputs, dt)
        with self._lock:
            return self._store.get(key)

    def put(self, signature: str, inputs: Mapping[str, Any], dt: float, value: Any) -> None:
        key = self._key(signature, inputs, dt)
        with self._lock:
            self._store[key] = value

    def clear(self) -> None:
        with self._lock:
            self._store.clear()
