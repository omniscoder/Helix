from helix.live import StateReducer, diff_snapshots


def test_diff_snapshots_tracks_changes():
    prev = {"a": {"x": 1.0}}
    current = {"a": {"x": 2.0}, "b": {"y": 3.0}}
    delta = diff_snapshots(current, prev)
    assert "b" in delta["added"]
    assert "a" in delta["updated"]
    assert not delta["removed"]


def test_state_reducer_records_history():
    reducer = StateReducer(hz=30)
    reducer.push({"n1": {"out": 1.0}})
    reducer.push({"n1": {"out": 1.0}, "n2": {"out": 2.0}})
    assert len(reducer.history) == 2
    assert "n2" in reducer.history[-1]["delta"]["added"]
