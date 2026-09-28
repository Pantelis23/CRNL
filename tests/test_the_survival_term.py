"""§115 — the survival term, the window, and the exact Q-process rate."""

import json
import math
import pathlib

import pytest

from experiments.the_survival_term import conditioned_gap, invert, stage1_death
from experiments.two_axes_for_the_position import PUBLISHED_RAILS, Elem

STORE = pathlib.Path("results/the_survival_term.json")


def _d():
    return json.loads(STORE.read_text())


# ---------------------------------------------------------------- live checks (cheap cells)

def test_ignoring_survival_reproduces_section_109():
    """P1: the new return path leaves v unchanged, so the old inversion is recovered exactly."""
    pub = {r["omega"]: r for r in
           json.loads(pathlib.Path("results/it_was_the_seed.json").read_text())["p3"]}
    e = Elem("published", *PUBLISHED_RAILS, cap_mult=1.25)
    mus, _ = e.operating_points(14, 2)
    pl = e.pi_low(14, mus[0])
    v, s1 = e.joint_absorbing(14, 2.0, return_survival=True)
    # NOT `==`: two identical expm_multiply calls differ in the last digit (T-NUM-a). This line
    # was first written as `==` -- the error §114.8c and §114.10 had already recorded, repeated
    # in the section written immediately after them -- and failed on 0.1355656168063782 vs
    # 0.13556561680637907.
    assert v == pytest.approx(e.joint_absorbing(14, 2.0), rel=1e-12)
    km, _, _ = e.candidate_averages(14)
    assert invert(v, 1.0, pl, 2.0) / km == pytest.approx(pub[14]["qsd"], rel=1e-6)
    assert 0.7 < s1 < 0.85, "stage 1 dies appreciably at Omega = 14 -- the whole point"


def test_lambda0_is_stage_one_death_exactly():
    """Stage 2 does not feed back into stage 1, so the conditioned process's leading decay
    must be stage 1's own QSD decay rate."""
    e = Elem("published", *PUBLISHED_RAILS, cap_mult=2.0)
    l0, l1, _ = conditioned_gap(e, 14)
    assert l0 == pytest.approx(stage1_death(e, 14), rel=1e-9)
    assert l1 > l0


# ---------------------------------------------------------------- the stored measurement

def test_survival_correction_shrinks_the_windowed_drift_by_a_quarter():
    p3 = _d()["p3"]
    assert p3["span_published"] == pytest.approx(1.1756, abs=1e-3)
    assert 0.70 < p3["span_corrected"] / p3["span_published"] < 0.80


def test_the_position_depends_on_the_unswept_window():
    """§115.2: the published t = 2 is one point on a converging sequence."""
    pos = _d()["p8"]["positions"]
    for om in ("55", "70"):
        seq = pos[om]
        assert seq == sorted(seq), "rises with the window"
        steps = [seq[i + 1] - seq[i] for i in range(len(seq) - 1)]
        assert steps == sorted(steps, reverse=True), "and converges -- increments shrink"
    assert pos["70"][-1] - pos["70"][0] > 0.5


def test_exact_gap_passes_its_three_checks():
    d = _d()
    assert all(m["corr"] > 0.9 for m in d["p11"]["mode"]), "lambda1 must be stage 2's inter-well mode"
    assert all(abs(b["cap3"] - b["cap2"]) < 1e-4 for b in d["p11"]["box"]), "box-converged"
    for r in d["p10"]["v3"]:
        if r["ratio"] > 1.1:        # where the occupancy bias is appreciable
            assert abs(r["windowed_cond"] - r["exact"]) < abs(r["windowed_pinned"] - r["exact"])
            assert abs(r["windowed_cond"] - r["exact"]) < 0.1


def test_exact_drift_is_real_half_size_and_stays_in_the_bracket():
    v1 = _d()["p10"]["v1"]
    pos = [r["pos_exact_cap2"] for r in v1]
    assert pos == sorted(pos), "the drift is real: still monotone in Omega"
    assert 0.45 < (pos[-1] - pos[0]) / 1.1756 < 0.55
    assert all(0.0 < p < 1.0 for p in pos), "never exits the bracket (§109 said it did)"


def test_section_114_3_outlier_was_the_instrument():
    """The 3.2x outlier sits below its bracketing pair under the exact rate, and it had the
    largest pinned-vs-conditioned occupancy ratio in the grid."""
    d = _d()
    rail = {r["cell"]: r["pos_exact"] for r in d["p9"]["rail"]}
    pub30 = next(r for r in d["p10"]["v1"] if r["omega"] == 30)["pos_exact_cap2"]
    assert rail["E-084 Om=14"] < min(rail["E-070 Om=30"], pub30)
    assert max(rail["E-070 Om=30"], pub30) / rail["E-084 Om=14"] < 1.5
    ratios = {r["cell"]: r["ratio"] for r in d["p10"]["v3"]}
    assert max(ratios, key=ratios.get) == "E-084 Om=14"


def test_rail_axis_no_longer_shallower_and_a_omega_now_leads():
    c = _d()["p12"]["candidates"]
    assert 0.9 < c["A*Omega"]["slope_ratio"] < 1.3
    assert c["A*Omega"]["rms_over_span"] < c["(margin/sigma)^2"]["rms_over_span"]


def test_failed_predictions_stay_recorded():
    """Rule 3: P4b and P5 failed and the stored verdicts must say so."""
    d = _d()
    assert d["p4b"]["holds"] is False
    p5 = d["p5"]
    assert all(abs(v["slope_corr"] - 1) > abs(v["slope_uncorr"] - 1) for v in p5.values())
