"""Hot swap manager for live edits."""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Callable, DefaultDict, List


class HotSwapManager:
    def __init__(self) -> None:
        self._callbacks: DefaultDict[str, List[Callable[[Any], None]]] = defaultdict(list)

    def register(self, scope: str, callback: Callable[[Any], None]) -> None:
        self._callbacks[scope].append(callback)

    def apply(self, scope: str, payload: Any) -> None:
        for callback in self._callbacks.get(scope, []):
            callback(payload)
