"""Transport helpers for realtime viz."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

from .feed import BundleFeed, LiveFrame, RealtimeFeed


class VizChannel:
    """Thin wrapper that normalizes realtime frames/control traffic."""

    def __init__(self, *, endpoint: Optional[str], bundle: Optional[Path]):
        if bundle:
            self._feed = BundleFeed(Path(bundle))
            self.endpoint = None
            self.interactive = False
        else:
            normalized = endpoint or "tcp://127.0.0.1:8765"
            self._feed = RealtimeFeed(normalized)
            self.endpoint = normalized
            self.interactive = True
        self.slice: str = "-"
        self.variant: str = "-"
        self.run_id: Optional[str] = None
        self.variant_options: List[str] = []
        self._latest: Optional[LiveFrame] = None

    def poll(self) -> Optional[LiveFrame]:
        frame = self._feed.poll()
        if frame:
            runtime = frame.runtime or {}
            self.slice = runtime.get("slice") or self.slice
            self.variant = runtime.get("variant") or self.variant
            run_id = runtime.get("run_id") or runtime.get("model")
            if run_id:
                self.run_id = str(run_id)
            variants = runtime.get("variants")
            if isinstance(variants, (list, tuple)):
                self.variant_options = [str(v) for v in variants]
            self._latest = frame
        return frame

    def send_control(self, payload) -> None:
        if not self.interactive:
            return
        self._feed.send_command(payload)

    def close(self) -> None:
        self._feed.close()
