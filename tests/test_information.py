"""Cost per bit delivered: the comparator-free measure, and its one trap."""

import numpy as np
import pytest

from crnl.information import (
    cost_per_bit, mutual_information_bits, shannon_bits,
)


# -- the information primitives --------------------------------------------

def test_shannon_bits_of_a_fair_coin_is_one():
    assert shannon_bits([0.5, 0.5]) == pytest.approx(1.0)
    assert shannon_bits([1.0, 0.0]) == pytest.approx(0.0)
    assert shannon_bits([0.25] * 4) == pytest.approx(2.0)


def test_mutual_information_is_one_bit_for_a_perfect_channel():
    """Disjoint conditionals carry the whole input bit."""
    assert mutual_information_bits([1.0, 0.0], [0.0, 1.0]) == pytest.approx(1.0)


def test_mutual_information_is_zero_when_the_output_ignores_the_input():
    assert mutual_information_bits([0.3, 0.7], [0.3, 0.7]) == pytest.approx(0.0)


def test_mutual_information_never_exceeds_one_bit_for_a_binary_input():
    rng = np.random.default_rng(0)
    for _ in range(50):
        a = rng.dirichlet(np.ones(5))
        b = rng.dirichlet(np.ones(5))
        assert -1e-12 <= mutual_information_bits(a, b) <= 1.0 + 1e-12


# -- the cascade measure ----------------------------------------------------

def test_information_decays_and_cost_accumulates_with_depth():
    rows = cost_per_bit(0.15, 30, 8.0, 12)
    info = [r["I_bits"] for r in rows]
    cost = [r["ds"] for r in rows]
    assert info == sorted(info, reverse=True)          # monotone decay
    assert cost == sorted(cost)                        # monotone accumulation
    assert info[0] < 1.0                               # the channel acts first


def test_cost_per_bit_grows_superlinearly_with_depth():
    """Cost accumulates linearly while information decays exponentially, so
    kT/bit must grow faster than depth. Measured at gamma=0.05, Omega=60,
    t=32: 42 / 534 / 2199 kT per bit at depth 1 / 10 / 30."""
    rows = cost_per_bit(0.05, 60, 32.0, 30)
    d1, d10, d30 = rows[0], rows[9], rows[29]
    assert d10["kT_per_bit"] / d1["kT_per_bit"] > 10.0
    assert d30["kT_per_bit"] / d10["kT_per_bit"] > 3.0


def test_a_bit_is_cheapest_at_strong_drive_when_depth_is_demanding():
    """The naive reading -- weak drive is cheap -- is exactly backwards once the
    bit has to survive: weak drive delivers no bits, so its cost per bit
    diverges. Measured at Omega=60, depth 30, t=32: 2199 kT/bit at gamma=0.05
    versus 15953 at gamma=0.30."""
    strong = cost_per_bit(0.05, 60, 32.0, 30)[-1]["kT_per_bit"]
    weak = cost_per_bit(0.30, 60, 32.0, 30)[-1]["kT_per_bit"]
    assert strong < weak / 3.0


def test_reliability_is_bought_superlinearly_in_population():
    """Cost per bit RISES with Omega: cost is extensive while information
    saturates. Measured (gamma=0.05, t=32, depth 30): 1240 / 2199 / 4198 kT per
    bit at Omega = 30 / 60 / 120, for 0.52 / 0.60 / 0.63 bits."""
    c = [cost_per_bit(0.05, om, 32.0, 30)[-1] for om in (30, 60, 120)]
    assert [r["kT_per_bit"] for r in c] == sorted(r["kT_per_bit"] for r in c)
    assert [r["I_bits"] for r in c] == sorted(r["I_bits"] for r in c)


# -- THE trap ---------------------------------------------------------------

def test_depth_one_is_degenerate_and_depth_thirty_is_not():
    """The regression guard on this module's central caveat.

    At depth 1 a stage that barely runs scores BEST, because one channel
    application from a rail hardly damages the bit: t_stage=0.05 costs 0.89
    kT/bit against 20.2 at t_stage=16. That is the do-nothing degeneracy this
    project has now met three times, and it is why cost_per_bit must never be
    quoted at depth 1.

    At depth 30 the ordering reverses -- doing nothing becomes catastrophic
    (5493 kT/bit at t=0.05) and an interior optimum exists. If this test ever
    stops discriminating, the measure has lost the property that makes it
    well-posed.
    """
    lazy_d1 = cost_per_bit(0.15, 30, 0.05, 1)[-1]["kT_per_bit"]
    busy_d1 = cost_per_bit(0.15, 30, 16.0, 1)[-1]["kT_per_bit"]
    assert lazy_d1 < busy_d1, "depth 1 should reward doing nothing"

    lazy_d30 = cost_per_bit(0.15, 30, 0.05, 30)[-1]["kT_per_bit"]
    busy_d30 = cost_per_bit(0.15, 30, 16.0, 30)[-1]["kT_per_bit"]
    assert busy_d30 < lazy_d30 / 2.0, "depth 30 must punish doing nothing"


def test_stage_time_has_an_interior_optimum_at_demanding_depth():
    """Too short does not restore; too long pays for idle cycling."""
    ts = [1.0, 4.0, 16.0, 64.0, 256.0]
    costs = [cost_per_bit(0.05, 30, t, 30)[-1]["kT_per_bit"] for t in ts]
    best = int(np.argmin(costs))
    assert 0 < best < len(ts) - 1, dict(zip(ts, costs))


def test_depth_must_be_positive():
    with pytest.raises(ValueError, match="depth"):
        cost_per_bit(0.15, 30, 8.0, 0)


def test_mutual_information_is_clamped_at_zero_not_slightly_negative():
    """When the bit is destroyed, I is a difference of two ~equal entropies and
    lands at float noise (-8.88e-16 measured). It must read 0, not a negative
    number that turns one cell's cost into `inf` while its neighbour reports a
    large finite value."""
    p = np.array([0.1, 0.2, 0.3, 0.4])
    assert mutual_information_bits(p, p.copy()) == 0.0
    assert mutual_information_bits(p, p.copy()) >= 0.0


def test_meaningfully_negative_mutual_information_raises():
    """The guard must not silently swallow a real misalignment."""
    import crnl.information as info_mod
    real_H = info_mod.shannon_bits
    # call order inside mutual_information_bits is (mixture, p_plus, p_minus);
    # make the mixture look LESS entropic than its parts, which is impossible
    # for real distributions and so must be reported rather than clamped
    calls = iter([0.0, 5.0, 5.0])
    try:
        info_mod.shannon_bits = lambda p: next(calls)
        with pytest.raises(ValueError, match="not float noise"):
            info_mod.mutual_information_bits([0.5, 0.5], [0.5, 0.5])
    finally:
        info_mod.shannon_bits = real_H
