"""Realtime feed transport utilities for helix live runs."""

from __future__ import annotations

import queue
import threading
from dataclasses import dataclass
from multiprocessing.connection import Client, Listener
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

AUTH_KEY = b"helix-realtime"


def parse_endpoint(endpoint: Optional[str]) -> Tuple[Union[str, Tuple[str, int]], str]:
    """
    Normalize a human-friendly endpoint string into a multiprocessing address.

    Supports:
      - ``tcp://host:port`` (e.g., tcp://127.0.0.1:8765)
      - ``ipc://path`` for UNIX domain sockets.
    """

    if not endpoint:
        return ("127.0.0.1", 8765), "tcp://127.0.0.1:8765"
    endpoint = endpoint.strip()
    if endpoint.startswith("ipc://"):
        path = endpoint[len("ipc://") :]
        path = path or "/tmp/helix-realtime.sock"
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        return path, f"ipc://{path}"
    if endpoint.startswith("tcp://"):
        host_port = endpoint[len("tcp://") :]
        if ":" not in host_port:
            host_port += ":8765"
        host, port_text = host_port.rsplit(":", 1)
        port = int(port_text)
        return (host or "127.0.0.1", port), f"tcp://{host or '127.0.0.1'}:{port}"
    if ":" in endpoint:
        host, port_text = endpoint.rsplit(":", 1)
        port = int(port_text)
        return (host or "127.0.0.1", port), f"tcp://{host or '127.0.0.1'}:{port}"
    raise ValueError(f"Unsupported realtime endpoint '{endpoint}'.")


@dataclass
class RealtimeHook:
    queue: "queue.Queue"
    server: Optional["RealtimeServer"] = None


class RealtimeServer:
    """Fan-out server that streams reducer payloads to external clients."""

    def __init__(self, endpoint: str):
        self.address, self.label = parse_endpoint(endpoint)
        self._listener = Listener(self.address, authkey=AUTH_KEY)
        self._clients: List = []
        self._lock = threading.Lock()
        self._running = threading.Event()
        self._running.set()
        self._on_command: Optional[Callable[[Dict[str, Any]], None]] = None
        threading.Thread(target=self._accept_loop, daemon=True).start()

    def _accept_loop(self) -> None:
        while self._running.is_set():
            try:
                conn = self._listener.accept()
            except (OSError, EOFError):
                break
            with self._lock:
                self._clients.append(conn)
            threading.Thread(target=self._command_loop, args=(conn,), daemon=True).start()

    def _command_loop(self, conn) -> None:
        while self._running.is_set():
            try:
                if not conn.poll(0.1):
                    continue
                payload = conn.recv()
            except (EOFError, OSError):
                break
            if isinstance(payload, dict) and self._on_command:
                self._on_command(payload)
        self._drop(conn)

    def broadcast(self, payload: Dict[str, Any]) -> None:
        with self._lock:
            alive = []
            for conn in self._clients:
                try:
                    conn.send(payload)
                    alive.append(conn)
                except (EOFError, OSError):
                    pass
            self._clients = alive

    def set_command_handler(self, handler: Optional[Callable[[Dict[str, Any]], None]]) -> None:
        self._on_command = handler

    def close(self) -> None:
        self._running.clear()
        try:
            self._listener.close()
        except Exception:
            pass
        with self._lock:
            for conn in self._clients:
                try:
                    conn.close()
                except Exception:
                    pass
            self._clients = []

    def _drop(self, conn) -> None:
        with self._lock:
            if conn in self._clients:
                self._clients.remove(conn)
        try:
            conn.close()
        except Exception:
            pass


class RealtimeClient:
    """Client helper that consumes broadcast payloads and exposes a thread-safe queue."""

    def __init__(self, endpoint: str, *, max_queue: int = 256):
        address, label = parse_endpoint(endpoint)
        self.label = label
        self._conn = Client(address, authkey=AUTH_KEY)
        self._queue: "queue.Queue" = queue.Queue(maxsize=max_queue)
        self._closed = threading.Event()
        threading.Thread(target=self._recv_loop, daemon=True).start()

    def _recv_loop(self) -> None:
        while not self._closed.is_set():
            try:
                payload = self._conn.recv()
            except (EOFError, OSError):
                break
            try:
                self._queue.put(payload, timeout=0.1)
            except queue.Full:
                pass
        self.close()

    def poll(self) -> Optional[Dict[str, Any]]:
        try:
            return self._queue.get_nowait()
        except queue.Empty:
            return None

    def send_command(self, payload: Dict[str, Any]) -> None:
        try:
            self._conn.send(payload)
        except (EOFError, OSError):
            self.close()

    def close(self) -> None:
        if self._closed.is_set():
            return
        self._closed.set()
        try:
            self._conn.close()
        except Exception:
            pass
