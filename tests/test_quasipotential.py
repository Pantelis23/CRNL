"""Tests for the exact quasipotential.

The one that matters most is `test_bottleneck_on_a_known_graph`: the flood is
the only nontrivial algorithm here, and on the real simplex there is nothing to
check it against, so it is pinned on a five-node path where the min-max barrier
can be read off by eye.
"""
from __future__ import annotations

import numpy as np
import pytest

from crnl.networks.am_asymmetric import am_asymmetric, basin_boundary
from crnl.networks.am_reversible import am_reversible, delta_star, lambda_antisym
from crnl.quasipotential import (
    P_FLOOR,
    barriers,
    bottleneck,
    landscape,
    ridge_profile,
)


def kappa_theory(gamma: float) -> float:
    """lambda / D0 with D0 = (2/9)(1+gamma): FINDINGS 15's corrected coefficient."""
    return lambda_antisym(gamma) / ((2.0 / 9.0) * (1.0 + gamma))


def test_bottleneck_on_a_known_graph():
    """Min-max over a path graph, where the answer is visible: 7."""
    W = np.array([0.0, 5.0, 2.0, 7.0, 1.0])
    nb = [np.array([1]), np.array([0, 2]), np.array([1, 3]),
          np.array([2, 4]), np.array([3])]
    land = {"W": W, "states": np.arange(5).reshape(-1, 1)}
    assert bottleneck(land, nb, 0, 4)["level"] == pytest.approx(7.0)
    # and the barrier between adjacent low nodes is the lower pass
    assert bottleneck(land, nb, 0, 2)["level"] == pytest.approx(5.0)


def test_bottleneck_is_symmetric_in_its_endpoints():
    W = np.array([0.0, 5.0, 2.0, 7.0, 1.0])
    nb = [np.array([1]), np.array([0, 2]), np.array([1, 3]),
          np.array([2, 4]), np.array([3])]
    land = {"W": W, "states": np.arange(5).reshape(-1, 1)}
    assert (bottleneck(land, nb, 0, 4)["level"]
            == bottleneck(land, nb, 4, 0)["level"])


def test_landscape_is_a_shifted_log_of_a_distribution():
    land = landscape(am_reversible(0.4), 60, 60.0)
    assert land["p"].sum() == pytest.approx(1.0)
    assert land["W"][np.isfinite(land["W"])].min() == pytest.approx(0.0)
    assert land["dW_max"] == pytest.approx(-np.log(P_FLOOR) / 60.0)
    # W is exactly -(1/Omega) ln p wherever it is resolved
    ok = land["resolved"]
    ref = -np.log(land["p"][ok]) / 60.0
    assert np.allclose(land["W"][ok] - land["W"][ok].min(), ref - ref.min())


def test_symmetric_network_has_equal_barriers():
    """beta = 0 must give dW_x == dW_y to solver precision, not merely close."""
    for gamma in [0.35, 0.45]:
        br = barriers(am_reversible(gamma), 120, 120.0)
        assert br["dW_x"] == pytest.approx(br["dW_y"], rel=1e-9)
        assert abs(br["saddle_eps"]) < 2.0 / 120


def test_barrier_level_dominates_both_minima():
    br = barriers(am_reversible(0.4), 120, 120.0)
    assert br["W_saddle"] >= max(br["W_x"], br["W_y"])
    assert br["dW_x"] > 0 and br["dW_y"] > 0


def test_out_of_range_barrier_raises_rather_than_returning_inf():
    """A barrier past the probability floor must be an error, not a number."""
    with pytest.raises(ValueError, match="exceeds this instrument's range"):
        barriers(am_reversible(0.05), 400, 400.0)


def test_ridge_profile_is_a_transverse_minimum():
    land = landscape(am_reversible(0.4), 80, 80.0)
    eps, Wm = ridge_profile(land)
    states, W = land["states"], land["W"]
    d = states[:, 0] - states[:, 1]
    for e, w in list(zip(eps, Wm))[::7]:
        sel = d == int(round(e * 80))
        assert w == pytest.approx(W[sel].min())


@pytest.mark.parametrize("gamma,total", [(0.35, 200), (0.40, 240), (0.45, 400)])
def test_ridge_curvature_matches_lambda_over_D0(gamma, total):
    """The measured coefficient must track lambda/D0, not (9/2)*lambda.

    Checked at ONE finite Omega, so the tolerance is 5% rather than the 0.1%
    that FINDINGS 15 reports after extrapolating Omega -> infinity. Each Omega
    here is inside that section's instrument window (lattice-limited from below,
    probability-floor-limited from above); Omega = 200 at gamma = 0.45 is NOT,
    and reading a 9% error there as a failure was this test's first bug.

    The second assertion is the one that actually discriminates, and it is not
    tolerance-sensitive: the two candidate formulas differ by 1/(1+gamma), which
    is 25-30% here, and the measurement sits on top of one of them.
    """
    land = landscape(am_reversible(gamma), total, float(total))
    eps, W = ridge_profile(land)
    ok = np.isfinite(W) & (np.abs(eps) < 0.25 * delta_star(gamma))
    k = -float(np.polyfit(eps[ok], W[ok], 2)[0])
    assert k == pytest.approx(kappa_theory(gamma), rel=0.05)
    assert abs(k - kappa_theory(gamma)) < abs(k - 1.5 * (1 - 2 * gamma))


def test_tilt_deepens_the_favoured_basin():
    """beta > 0 must make X's barrier out larger and Y's smaller."""
    gamma, N = 0.42, 150
    sym = barriers(am_reversible(gamma), N, float(N))
    for beta in [0.05, 0.10]:
        net = am_asymmetric(gamma, beta)
        br = barriers(net, N, float(N), split=basin_boundary(gamma, beta))
        assert br["dW_x"] > sym["dW_x"], (beta, br["dW_x"], sym["dW_x"])
        assert br["dW_y"] < sym["dW_y"], (beta, br["dW_y"], sym["dW_y"])


def test_tilt_is_antisymmetric_in_the_barriers():
    """Relabelling X<->Y is beta -> -beta, so the two barriers must swap."""
    gamma, N, beta = 0.42, 150, 0.08
    pos = barriers(am_asymmetric(gamma, beta), N, float(N),
                   split=basin_boundary(gamma, beta))
    neg = barriers(am_asymmetric(gamma, -beta), N, float(N),
                   split=basin_boundary(gamma, -beta))
    assert pos["dW_x"] == pytest.approx(neg["dW_y"], rel=1e-6)
    assert pos["dW_y"] == pytest.approx(neg["dW_x"], rel=1e-6)
