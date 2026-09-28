"""§117 — the closure's remainder, factor by factor, and the near-miss channel."""

import json
import pathlib

STORE = pathlib.Path("results/where_the_remainder_lives.json")


def _d():
    return json.loads(STORE.read_text())


def test_factor_split_is_an_identity():
    assert _d()["p1_worst"] < 1e-6


def test_p_transmit_was_innocent():
    """T-CASC-aj's suspect: one Omega=30 value used everywhere. Right to <= 2% at every Omega."""
    assert all(0.99 < r["r_pt"] < 1.03 for r in _d()["rows"])


def test_the_large_omega_rise_was_the_pinned_occupancy():
    rows = _d()["rows"]
    pin = [r["r_los1_pinned"] for r in rows]
    cond = [r["r_los1_cond"] for r in rows]
    assert pin[-1] < pin[0] - 0.1, "pinned pi falls away at large Omega"
    assert cond[-1] > cond[0], "the consistent pair does not"


def test_stage_one_residual_is_a_window_transient():
    for r in _d()["p5"]:
        assert r["r_los0_t8"] < r["r_los0_t2"]
        assert abs(r["r_los0_t8"] - 1) < 0.02


def test_stage_two_residual_is_not():
    """The kill test split: stage 2's conditioned factor gets WORSE with the window."""
    for r in _d()["p5"]:
        if r["omega"] >= 20:
            assert r["r_los1_cond_t8"] < r["r_los1_cond_t2"]
            assert r["r_los1_cond_t8"] < 0.9


def test_the_q_process_is_exact_for_the_absorbing_chain():
    seq = [x["ratio"] for x in _d()["p7"]["seq"]]
    assert seq == sorted(seq, reverse=True), "monotone approach"
    steps = [seq[i] - seq[i + 1] for i in range(len(seq) - 1)]
    assert all(steps[i + 1] < steps[i] for i in range(len(steps) - 1)), "increments shrink"
    assert abs(seq[-1] - 1) < 0.01


def test_the_free_chain_plateaus_below_it_the_near_miss_channel():
    free30 = next(x for x in _d()["p6"] if x["omega"] == 30)["seq"]
    by_t = {x["t"]: x["ratio"] for x in free30}
    assert abs(by_t[16.0] - by_t[32.0]) < 0.01, "a plateau, not a transient"
    assert 0.80 < by_t[32.0] < 0.87
    absorbing = {x["t"]: x["ratio"] for x in _d()["p7"]["seq"]}
    assert absorbing[32.0] - by_t[32.0] > 0.15, "same seed, same window, different chains"
