"""§118 — do stage-1 near-misses cause stage 2's excess falls? Reads the stored exact-CME run."""

import json
import pathlib

STORE = pathlib.Path("results/near_misses.json")


def _d():
    return json.loads(STORE.read_text())


def _t16(om):
    return next(r for r in _d()["cells"][str(om)] if r["t"] == 16.0)


def test_tagging_is_wired():
    assert _d()["p1_rel"] < 1e-5


def test_never_dipped_set_is_the_q_process():
    """P2: no plateau near 0.84 — it approaches 1 from above like §117's absorbing chain."""
    for om in (30, 55):
        seq = [r["r_never"] for r in _d()["cells"][str(om)]]
        assert seq == sorted(seq, reverse=True)
        assert 1.0 < seq[-1] < 1.03


def test_the_excess_sits_in_the_dipped_set():
    for om in (30, 55):
        r = _t16(om)
        assert r["dipped"] > 5 * r["never"]
        assert r["r_never"] - r["r_free"] > 0.15


def test_clipping_at_the_saddle_does_nothing():
    """P4: the null that turned out to be mis-designed (the saddle is below stage 2's edge)."""
    for om in (30, 55):
        r = _t16(om)
        assert abs(r["r_clipped"] - r["r_free"]) < 0.01


def test_selection_by_approach_refuted():
    """P6: stage 2 is mostly still HIGH when stage 1 first crosses, except at Omega = 14."""
    a = {x["omega"]: x["p_low_at_crossing"] for x in _d()["approach"]}
    assert a[30] < 0.3 and a[55] < 0.3
    assert a[14] > 0.6


def _closed(om, ratio):
    r = _t16(om)
    return (ratio - r["r_free"]) / (r["r_never"] - r["r_free"])


def test_hiding_the_sub_saddle_phase_closes_little():
    """P7 failed: the sub-saddle phase is a minor part of the cause."""
    for om in (30, 55):
        x = next(v for v in _d()["p7"][str(om)]["as_healthy"] if v["t"] == 16.0)
        assert 0.05 < _closed(om, x["ratio"]) < 0.25


def test_hiding_post_dip_input_below_the_edge_closes_most():
    """P8: ~70% closed, the same at both volumes; its 0.03 criterion still failed."""
    fr = []
    for om in (30, 55):
        x = next(v for v in _d()["p8"][str(om)] if v["t"] == 16.0)
        assert _t16(om)["r_never"] - x["ratio"] > 0.03
        fr.append(_closed(om, x["ratio"]))
    assert all(0.6 < f < 0.8 for f in fr)
    assert abs(fr[0] - fr[1]) < 0.05


def test_dose_response_saturates_below_the_prediction():
    """P9b failed: monotone in the threshold, saturating near 0.75 rather than reaching 0.9."""
    for om in (30, 55):
        o = _d()["p9"][str(om)]
        assert abs(o["self_ratio"] / o["r_free"] - 1) < 1e-9, "P9a: the tagged chain is wired"
        fr = [_closed(om, o[str(lv)]) for lv in (1.3, 1.5795, 2.0, 2.5)]
        assert fr == sorted(fr)
        assert 0.6 < fr[-1] < 0.85
        assert fr[-1] - fr[-2] < 0.03, "saturated"


def test_already_low_split_is_not_the_line():
    """P10b failed: the still-high-at-dip subgroup keeps falling well above never-dipped."""
    for om in (30, 55):
        o = _d()["p10"][str(om)]
        assert o["p_hi_at_dip"] - o["p_never"] > 0.03


def test_state_at_dip_plus_post_dip_input_is_sufficient_and_splits():
    """P11 / P11b: the pair closes the gap (over-correcting at Om = 30); redraw alone is ~1/4."""
    for om in (30, 55):
        assert _d()["p11"][str(om)]["closed"] > 0.95
        assert 0.2 < _d()["p11b"][str(om)]["closed"] < 0.3
    assert _d()["p11"]["30"]["closed"] > 1.05, "over-correction, which is why P11b exists"
