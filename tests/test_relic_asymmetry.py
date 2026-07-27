"""Tests for the relic-asymmetry result (FINDINGS 18) and the inverted tilt rule.

`test_symmetric_start_is_a_coin_flip` is the one that guards the whole section:
with beta = 0 and an exactly symmetric start, P(X) must be 0.5 by symmetry, so
any parity or rounding bug in building the initial state shows up here rather
than as a plausible-looking asymmetry (which is how FINDINGS 16's harness bug
first presented).
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import norm

from experiments.relic_asymmetry import (
    beta_for_u,
    crossover_beta_root_omega,
    p_survive,
    seed_ratio,
)
from experiments.tilt_rule_limit import matched_prior


@pytest.mark.parametrize("omega", [60, 120, 240])
def test_symmetric_start_is_a_coin_flip(omega):
    assert p_survive(0.05, 0.0, omega) == pytest.approx(0.5, abs=1e-9)


def test_seed_ratio_is_linear_in_beta_and_root_omega():
    """u = (g/lambda)/sigma must scale as beta*sqrt(Omega), the collapse variable."""
    base = seed_ratio(0.05, 0.01, 100.0)
    assert seed_ratio(0.05, 0.02, 100.0) == pytest.approx(2 * base, rel=1e-12)
    assert seed_ratio(0.05, 0.01, 400.0) == pytest.approx(2 * base, rel=1e-12)
    # and beta_for_u inverts it
    for u in (0.5, 2.0):
        assert seed_ratio(0.05, beta_for_u(0.05, u, 240.0), 240.0) == pytest.approx(
            u, rel=1e-9)


def test_crossover_uses_the_corrected_diffusion():
    """The threshold rests on FINDINGS 15's D_0 = (1+gamma)/9, not (1/9).

    With the uncorrected diffusion the crossover would lose its 1/(1-gamma), so
    this pins the gamma-dependence rather than just the number.
    """
    assert crossover_beta_root_omega(0.05) == pytest.approx(0.8204, rel=1e-3)
    for g in (0.0, 0.05, 0.2, 0.35):
        assert crossover_beta_root_omega(g) == pytest.approx(
            np.sqrt(3) / 2 * (1 - 2 * g) / (1 - g), rel=1e-12)
    # at gamma = 0 the two conventions agree, so the test must not pass trivially
    assert crossover_beta_root_omega(0.0) == pytest.approx(np.sqrt(3) / 2, rel=1e-12)


@pytest.mark.parametrize("u", [0.5, 1.0, 2.0])
def test_p_survive_collapses_across_omega_onto_phi(u):
    """FINDINGS 18.1: P(X) depends on (beta, Omega) only through u, and equals Phi(u)."""
    ps = [p_survive(0.05, beta_for_u(0.05, u, om), om) for om in (60, 120, 240)]
    assert max(ps) - min(ps) < 0.006, ps            # measured worst case 0.0028
    assert np.mean(ps) == pytest.approx(norm.cdf(u), abs=0.012)  # measured 0.0065


def test_p_survive_is_monotone_in_the_tilt():
    ps = [p_survive(0.05, beta_for_u(0.05, u, 120), 120)
          for u in (0.0, 0.5, 1.0, 2.0)]
    assert ps == sorted(ps)
    assert ps[-1] > 0.95


def test_matched_prior_inverts_the_optimality_condition():
    """A tilt that helps X must be matched to a prior favouring X, and vice versa.

    dln(e+)/dbeta < 0 and dln(e-)/dbeta > 0 always (the tilt trades one error for
    the other), so the sign convention is what this pins.
    """
    p = matched_prior(1e-2, 5e-2, -30.0, 30.0)
    assert p is not None and 0.5 < p < 1.0
    # a larger gap between the two errors demands a more extreme prior
    p2 = matched_prior(5e-3, 8e-2, -30.0, 30.0)
    assert p2 is not None and p2 > p
    # mirroring the two errors mirrors the prior about 1/2
    pm = matched_prior(5e-2, 1e-2, -30.0, 30.0)
    assert pm is not None and pm == pytest.approx(1.0 - p, rel=1e-9)
    # an error gap too wide to be optimal for ANY prior has no root, and the
    # function says so instead of clamping to p = 1
    assert matched_prior(1e-3, 1e-1, -30.0, 30.0) is None


def test_matched_prior_is_symmetric_at_equal_errors():
    """Equal errors and equal-and-opposite sensitivities -> no prior preference."""
    p = matched_prior(2e-2, 2e-2, -25.0, 25.0)
    assert p is not None
    assert p == pytest.approx(0.5, abs=1e-6)
