"""Tests for the biased-source design rule (FINDINGS 17).

The load-bearing one is `test_symmetric_source_is_symmetric`: FINDINGS 16's first
run reported P(ok|X) = 0.638 against P(ok|Y) = 0.667 at beta = 0, where symmetry
forces equality, because integer division gave the two inputs different
magnitudes. That artifact was 20% of the effect under test, so it gets a test
rather than a comment.
"""
from __future__ import annotations

import numpy as np
import pytest

from experiments.biased_source import errors, information, optimal_beta


@pytest.mark.parametrize("omega", [100, 150, 200])
@pytest.mark.parametrize("gamma", [0.35, 0.42])
def test_symmetric_source_is_symmetric(gamma, omega):
    """beta = 0 must give e+ == e- exactly, not approximately."""
    e_plus, e_minus = errors(gamma, 0.0, omega, 0.25)
    assert e_plus == pytest.approx(e_minus, rel=1e-9), (e_plus, e_minus)


@pytest.mark.parametrize("omega", [100, 200])
def test_tilt_trades_one_error_for_the_other(omega):
    """beta > 0 must lower the X error and raise the Y error, monotonically."""
    prev = None
    for beta in (0.0, 0.01, 0.02, 0.04):
        ep, em = errors(0.35, beta, omega, 0.25)
        if prev is not None:
            assert ep < prev[0], (beta, ep, prev[0])
            assert em > prev[1], (beta, em, prev[1])
        prev = (ep, em)


def test_information_is_bounded_by_the_source_entropy():
    for p in (0.5, 0.7, 0.9):
        h = -p * np.log2(p) - (1 - p) * np.log2(1 - p)
        for e in (1e-6, 0.01, 0.2):
            I = information(p, e, e)
            assert 0.0 <= I <= h + 1e-12, (p, e, I, h)
    # a useless channel carries nothing
    assert information(0.7, 0.5, 0.5) == pytest.approx(0.0, abs=1e-12)
    # a perfect channel carries the whole source
    assert information(0.7, 0.0, 0.0) == pytest.approx(
        -0.7 * np.log2(0.7) - 0.3 * np.log2(0.3), rel=1e-9)


def test_symmetric_source_wants_no_tilt():
    """FINDINGS 16: beta = 0 is the maximum, not merely stationary."""
    r = optimal_beta(0.35, 150, 0.25, 0.50)
    assert r is not None
    assert r["beta_star"] < 2e-3, r["beta_star"]


@pytest.mark.parametrize("p", [0.70, 0.90])
def test_biased_source_wants_a_tilt_toward_the_likely_symbol(p):
    """FINDINGS 17's headline: beta* > 0, and it grows with the prior."""
    r = optimal_beta(0.35, 150, 0.25, p)
    assert r is not None
    assert r["beta_star"] > 5e-3, r
    # and the tilt stays far from the fold -- the rule is a gentle one
    assert r["beta_star_over_bc"] < 0.25, r["beta_star_over_bc"]


def test_the_optimum_matches_error_log_ratio_to_prior_log_odds():
    """P2's SHAPE: ln(e-/e+) proportional to ln(p/(1-p)) with zero intercept.

    The coefficient is NOT asserted to be 1 -- it is 0.76 at Omega = 200 and its
    Omega -> infinity limit is undetermined (FINDINGS 17.2). What is pinned here
    is the proportionality and the vanishing intercept, which is the part that
    holds to R^2 = 0.9999.
    """
    ps = [0.60, 0.75, 0.90]
    xs, ys = [], []
    for p in ps:
        r = optimal_beta(0.35, 150, 0.25, p)
        assert r is not None
        xs.append(np.log(p / (1 - p)))
        ys.append(r["P2_measured"])
    slope, intercept = np.polyfit(xs, ys, 1)
    resid = np.array(ys) - (slope * np.array(xs) + intercept)
    r2 = 1 - resid.var() / np.array(ys).var()
    assert r2 > 0.999, (r2, ys)
    assert abs(intercept) < 0.05, intercept
    assert 0.5 < slope < 1.0, slope
