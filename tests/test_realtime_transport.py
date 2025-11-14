"""Realtime server/client coordination tests."""

from __future__ import annotations

import time

from helix.live.realtime import RealtimeClient, RealtimeServer, parse_endpoint


def test_parse_endpoint_tcp():
    address, label = parse_endpoint("tcp://127.0.0.1:9001")
    assert address == ("127.0.0.1", 9001)
    assert label == "tcp://127.0.0.1:9001"


def test_realtime_server_client_roundtrip(tmp_path):
    endpoint = f"ipc://{tmp_path / 'helix.sock'}"
    server = RealtimeServer(endpoint)
    received = []

    def _handler(command):
        received.append(command)

    server.set_command_handler(_handler)
    client = RealtimeClient(endpoint)

    sample = {"time": 0.0, "runtime": {"variant": "wt"}, "delta": {"added": {}, "updated": {}, "removed": {}}}
    server.broadcast(sample)
    deadline = time.time() + 1.0
    message = None
    while time.time() < deadline and message is None:
        message = client.poll()
    assert message == sample

    client.send_command({"command": "pause"})
    deadline = time.time() + 1.0
    while time.time() < deadline and not received:
        time.sleep(0.01)
    assert received[0]["command"] == "pause"
    client.close()
    server.close()
