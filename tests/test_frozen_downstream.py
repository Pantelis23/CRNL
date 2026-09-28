"""§120 -- predicting the near-miss constant from single-element data: frozen fails, memory nearly."""

import json
import pathlib

STORE = pathlib.Path("results/frozen_downstream.json")
POSITIVE_ZONE = "SXbZac"


def _d():
    return json.loads(STORE.read_text())


def _cells(arms):
    return {k: v for k, v in _d().items() if "/" not in k and k.split(":")[0] in arms}


def test_frozen_hazard_overpredicts_and_is_not_a_constant():
    """R2 FAILED: memoryless stage 2 over-predicts p_d everywhere, by 1.1x to over 100x."""
    r = [v["ratio_d"] for v in _cells("SXbZacdN").values()]
    assert min(r) > 1.0
    assert max(r) / min(r) > 50


def test_frozen_hazard_swaps_a_and_c():
    """R1 FAILED: a and c in the wrong order at both Omega."""
    d = _d()
    for om in (40, 55):
        assert d[f"a:{om}"]["p_d"] < d[f"c:{om}"]["p_d"]
        assert d[f"a:{om}"]["p1_red"] > d[f"c:{om}"]["p1_red"]


def test_memory_reduction_runs_low_by_a_near_constant_factor():
    """R4 failed by its letter (a at 0.34-0.39) but every positive-zone cell sits in [0.3, 0.7]
    while the exact p_d spans ~90x."""
    c = _cells(POSITIVE_ZONE)
    r = [v["ratio_ode"] for v in c.values()]
    assert all(0.3 < x < 0.7 for x in r)
    pd = [v["p_d"] for v in c.values()]
    assert max(pd) / min(pd) > 50
    assert min(v["ratio_ode"] for k, v in c.items() if k.startswith("a:")) < 0.5, "the letter failed"


def test_memory_reduction_orders_the_arms_at_40():
    """R5 held at Omega = 40, including a < Z (the equal-width pair zone width could not split)."""
    d = _d()
    order = "NdacZbXS"
    ode = [d[f"{k}:40"]["p_d_ode"] for k in order]
    ex = [d[f"{k}:40"]["p_d"] for k in order]
    assert ex == sorted(ex) and ode == sorted(ode)


def test_memory_reduction_step_converged():
    d = _d()
    assert abs(d["X:40/tau_half"]["p_d_ode"] / d["X:40"]["p_d_ode"] - 1) < 0.02


def test_stage_two_noise_closes_most_of_the_gap_and_never_overshoots():
    """R6: stochastic stage 2 lifts every cell into [0.66, 1.06] of the exact p_d, across an
    ~1800x range; X misses R6's 0.7 floor, and no positive-zone cell overshoots 1."""
    c = _cells("SXbZacdN")
    r = {k: v["ratio_noisy"] for k, v in c.items()}
    assert all(0.6 < x < 1.1 for x in r.values())
    assert all(r[k] < 1.0 for k in r if k[0] in POSITIVE_ZONE)
    assert min(r["X:40"], r["X:55"]) < 0.7, "the letter failed in X"
    pd = [v["p_d"] for v in c.values()]
    assert max(pd) / min(pd) > 1000
