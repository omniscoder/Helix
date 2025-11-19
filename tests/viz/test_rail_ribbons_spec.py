from __future__ import annotations

import numpy as np

from helix.gui.modern.qt import HelixModernWidget
from helix.gui.modern.spec import EditEvent, EditEventType, EditVisualizationSpec


def _make_simple_spec() -> EditVisualizationSpec:
    sequence = "A" * 20
    cut_idx = 10
    events = [
        EditEvent(t=0.2, type=EditEventType.RECOGNITION, index=cut_idx),
        EditEvent(t=0.35, type=EditEventType.NICK_PRIMARY, index=cut_idx),
        EditEvent(t=0.9, type=EditEventType.REPAIR_COMPLETE, index=cut_idx),
    ]
    return EditVisualizationSpec(
        sequence=sequence,
        pam_index=0,
        guide_range=(0, len(sequence)),
        edit_type="test",
        events=events,
    )


def test_param_ribbons_cut_and_repair_shapes() -> None:
    spec = _make_simple_spec()
    interleaved, ranges, head_interleaved = HelixModernWidget._build_param_ribbons_from_spec(spec)

    # Expect two ribbons (cut + repair), each with 40 samples.
    assert len(ranges) == 2
    total_vertices = sum(count for _, count in ranges)
    assert interleaved.shape == (total_vertices, 5)

    pos = interleaved[:, :3]
    evt_time = interleaved[:, 3]
    kinds = interleaved[:, 4]

    first0, count0 = ranges[0]
    first1, count1 = ranges[1]

    cut_arc = pos[first0:first0 + count0]
    repair_arc = pos[first1:first1 + count1]
    cut_t = float(evt_time[first0])
    repair_t = float(evt_time[first1])
    cut_kind = float(kinds[first0])
    repair_kind = float(kinds[first1])

    rail_spacing = 0.3

    # Cut ribbon: symmetric arch above the top rail, never dips below y_top.
    assert cut_kind == 0.0
    assert cut_arc[:, 1].min() >= rail_spacing - 1e-4

    # Repair ribbon: runs from top rail down to bottom rail.
    assert repair_kind == 2.0
    assert abs(repair_arc[0, 1] - rail_spacing) < 1e-4
    assert abs(repair_arc[-1, 1] + rail_spacing) < 1e-4
    assert repair_arc[:, 1].min() < rail_spacing - 0.05

    # Event times propagate into the weight channel used for time gating.
    assert abs(cut_t - 0.35) < 1e-6
    assert abs(repair_t - 0.9) < 1e-6

    # Arrowheads are present and their tips coincide with ribbon endpoints.
    assert head_interleaved is not None
    heads = head_interleaved.reshape(-1, 5)
    assert heads.shape[0] == 6  # 2 ribbons * 3 vertices per triangle
    head_pos = heads[:, :3]
    tips = head_pos[0::3]
    assert any(
        np.allclose(tip, cut_arc[-1], atol=1e-5)
        or np.allclose(tip, repair_arc[-1], atol=1e-5)
        for tip in tips
    )

