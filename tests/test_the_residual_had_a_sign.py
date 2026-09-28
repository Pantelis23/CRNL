"""§116 — §107's residual drift with stage 2's rate replaced by the exact gap."""

import json
import pathlib

import pytest

STORE = pathlib.Path("results/the_residual_had_a_sign.json")


def _d():
    return json.loads(STORE.read_text())


def test_the_rebuilt_closure_reproduces_section_107_at_the_solver_floor():
    """P1 printed FAILS against a 1e-9 gate; the agreement is at the eigensolver's floor
    (T-NUM-a), which is the most the instrument supports."""
    d = _d()
    assert 1e-9 < d["p1_worst"] < 1e-5, d["p1_worst"]


def test_rebuild_is_live_not_just_stored():
    from experiments.chain_without_a_joint_solve import chain_operating_points
    from experiments.the_corrected_closure import geometric_rate, pi_low
    from experiments.the_residual_had_a_sign import closure_with

    meas = {r["omega"]: r["meas_ratio"] for r in json.loads(
        pathlib.Path("results/where_the_expansion_frays.json").read_text())["cells"]}
    mus, _ = chain_operating_points(14, 2)
    r = closure_with(14, geometric_rate(14, mus[0]), pi_low(14, mus[0])) / meas[14]
    assert r == pytest.approx(0.7902761288779301, rel=1e-6)


def test_the_exact_rate_has_the_sign_section_107_said_was_impossible():
    rows = _d()["rows"]
    ratios = [r["k_exact"] / r["k_geom"] for r in rows]
    assert ratios == sorted(ratios), "k_exact / k_geom rises with Omega"
    assert ratios[-1] > 1.5


def test_two_thirds_of_the_log_drift_is_stage_twos_rate():
    import math
    d = _d()
    assert d["span_geom"] == pytest.approx(1.833, abs=0.005)
    assert math.log(d["span_exact"]) / math.log(d["span_geom"]) < 0.4


def test_the_remainder_is_u_shaped_not_gone():
    """§116.2: the span criterion passed; the cells say the residual changed shape."""
    ys = [r["ratio_exact"] for r in _d()["rows"]]
    i = ys.index(min(ys))
    assert 0 < i < len(ys) - 1, "minimum is interior"
    assert ys[-1] / min(ys) > 1.15, "and it rises again at large Omega"
    assert all(y < 1.0 for y in ys), "the model is low at every Omega"


def test_p2_failed_at_omega_14_and_stays_recorded():
    rows = _d()["rows"]
    assert rows[0]["omega"] == 14 and rows[0]["k_exact"] < rows[0]["k_geom"]
    assert _d()["p2_all_down"] is False


def test_published_conclusions_audit():
    a = _d()["audit"]
    assert a["s110_refutation"] is True
    assert a["s109_exits_bracket"] is False
    assert a["s102_fast_end_Om14"] is True and a["s102_fast_end_Om30"] is False
    assert a["s107_toward_fast_end"] is False


def test_matched_seed_removes_the_low_omega_arm_but_not_the_large_omega_rise():
    """§116.5: the seed was half the story -- the stored verdict's 'exonerated' is wrong."""
    rows = _d()["seed_kill"]["rows"]
    lo = [r["ratio_exact_matched"] for r in rows if r["omega"] <= 40]
    hi = [r["ratio_exact_matched"] for r in rows if r["omega"] >= 40]
    assert max(lo) / min(lo) < 1.03, "flat to ~2% over Omega = 14-40"
    assert hi == sorted(hi) and hi[-1] / hi[0] > 1.15, "and a rise beyond 40 survives"
    assert 1.15 < min(lo) < 1.25, "on a flat ~21% offset"
