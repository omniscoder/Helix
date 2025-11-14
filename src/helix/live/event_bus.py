"""Lock-free style event bus using a bounded deque."""

from __future__ import annotations

from collections import deque
from typing import Any, Callable, Deque, List


class EventBus:
    def __init__(self, capacity: int = 1024):
        self.capacity = capacity
        self._events: Deque[Any] = deque(maxlen=capacity)
        self._subscribers: List[Callable[[Any], None]] = []

    def publish(self, event: Any) -> None:
        self._events.append(event)
        for callback in list(self._subscribers):
            callback(event)

    def subscribe(self, callback: Callable[[Any], None]) -> None:
        self._subscribers.append(callback)

    def drain(self) -> List[Any]:
        events = list(self._events)
        self._events.clear()
        return events

    def __len__(self) -> int:  # pragma: no cover - trivial
        return len(self._events)
