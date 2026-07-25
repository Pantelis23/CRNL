"""The crossover between the restoration wall and the channel floor.

FINDINGS 12. These tests pin the two limits of the saddle-point formula and the
fact that they are limits of ONE expression, plus the regime where the formula
fits for the wrong reason.
"""

import numpy as np
import pytest

from crnl.information import (
    crossover_omega, flip_probability, predicted_exponent, wall_coefficient,
)
from crnl.networks.am_reversible import GAMMA_C, delta_star, lambda_antisym


# -- the coefficient --------------------------------------------------------

def test_wall_coefficient_reduces_to_the_irreversible_value():
    """design.md 9 derives kappa = 3/2 for irreversible AM."""
    assert wall_coefficient(0.0) == pytest.approx(1.5)


def test_wall_coefficient_tracks_the_restoring_gain():
    """kappa = (9/2) * lambda_antisym, and both vanish together at gamma_c --
    which is what makes restoration impossible there rather than merely slow."""
    for gamma in (0.0, 0.1, 0.25, 0.4, 0.49):
        assert wall_coefficient(gamma) == pytest.approx(
            4.5 * lambda_antisym(gamma), rel=1e-12)
    assert wall_coefficient(GAMMA_C) == pytest.approx(0.0, abs=1e-15)


# -- the two limits of one formula -----------------------------------------

def test_small_omega_limit_is_the_restoration_wall():
    """2 kappa Omega sigma^2 << 1: the exponent is linear in Omega, i.e.
    P(error) ~ exp(-kappa e^2 Omega), FINDINGS 1-2."""
    gamma, nf = 0.05, 0.02          # tiny noise -> wall regime out to large Omega
    kappa, d = wall_coefficient(gamma), delta_star(gamma)
    for omega in (4, 8, 16):
        assert predicted_exponent(gamma, omega, nf) == pytest.approx(
            kappa * omega * d ** 2, rel=0.02)


def test_large_omega_limit_is_the_omega_independent_channel_floor():
    """2 kappa Omega sigma^2 >> 1: the exponent saturates at delta*^2/(2 sigma^2)
    and more molecules buy nothing."""
    gamma, nf = 0.05, 0.45
    d = delta_star(gamma)
    ceiling = d ** 2 / (2.0 * (nf * d) ** 2)
    big = [predicted_exponent(gamma, om, nf) for om in (4000, 40000)]
    assert big[1] > big[0]
    assert big[1] == pytest.approx(ceiling, rel=0.01)
    assert all(b < ceiling for b in big)          # approached from below


def test_crossover_is_where_the_two_terms_balance():
    for gamma in (0.05, 0.3):
        for nf in (0.15, 0.35):
            om_x = crossover_omega(gamma, nf)
            # at Omega_x the denominator is exactly 2, i.e. the exponent is half
            # its unsaturated value
            kappa, d = wall_coefficient(gamma), delta_star(gamma)
            assert predicted_exponent(gamma, om_x, nf) == pytest.approx(
                0.5 * kappa * om_x * d ** 2, rel=1e-9)


def test_crossover_moves_to_larger_population_as_noise_falls():
    xs = [crossover_omega(0.05, nf) for nf in (0.45, 0.35, 0.20, 0.10)]
    assert xs == sorted(xs)


# -- measured behaviour -----------------------------------------------------

def test_flip_probability_falls_exponentially_in_the_wall_regime():
    """At low channel noise the population-limited side is alive: p falls by
    orders of magnitude with Omega. Measured at gamma=0.05, sigma/delta*=0.10:
    2.0e-3 at Omega=4 down to 1.8e-14 at Omega=96."""
    ps = [flip_probability(0.05, om, 16.0, 0.10, depth=40) for om in (4, 12, 24)]
    assert all(np.isfinite(ps))
    assert ps[0] / ps[-1] > 1e3                    # measured ~8500x
    assert ps == sorted(ps, reverse=True)


def test_flip_probability_saturates_in_the_channel_regime():
    """At high channel noise more molecules buy almost nothing. Measured at
    gamma=0.05, sigma/delta*=0.45: 3.4e-2 -> 1.5e-2 over a 30x population."""
    ps = [flip_probability(0.05, om, 16.0, 0.45, depth=40) for om in (4, 32, 96)]
    assert all(np.isfinite(ps))
    assert ps[0] / ps[-1] < 4.0


def test_the_two_regimes_differ_by_orders_of_magnitude_at_the_same_population():
    """The whole point: at fixed Omega, which side you are on decides whether a
    bigger population helps at all."""
    quiet = flip_probability(0.05, 64, 16.0, 0.10, depth=40)
    loud = flip_probability(0.05, 64, 16.0, 0.45, depth=40)
    assert loud / quiet > 1e6


def test_collapse_holds_at_strong_drive():
    """-ln p against the parameter-free prediction, one gamma."""
    xs, ys = [], []
    for nf in (0.15, 0.28, 0.45):
        for om in (8, 16, 32, 64):
            p = flip_probability(0.05, om, 16.0, nf, depth=40)
            if np.isfinite(p) and p > 1e-15:
                xs.append(predicted_exponent(0.05, om, nf))
                ys.append(-np.log(p))
    xs, ys = np.array(xs), np.array(ys)
    slope, intercept = np.polyfit(xs, ys, 1)
    r2 = 1 - np.sum((ys - (slope * xs + intercept)) ** 2) / np.sum(
        (ys - ys.mean()) ** 2)
    assert r2 > 0.95, r2
    assert 0.5 < slope < 1.5, slope


def test_near_gamma_c_the_flip_is_not_channel_driven():
    """The caveat that keeps FINDINGS 12 honest.

    At gamma=0.45 the formula still fits (R^2=0.97) but for the WRONG reason:
    p is nearly independent of the channel noise (measured 0.089 -> 0.121 across
    a 4.5x change in sigma at Omega=96), so the channel term is not doing the
    work. There the stage simply fails to hold the state on its own -- the
    landscape is too shallow and t_stage is well under one relaxation time.
    A good fit in that block is not evidence for the mechanism.
    """
    ps = [flip_probability(0.45, 96, 16.0, nf, depth=40)
          for nf in (0.10, 0.45)]
    assert all(np.isfinite(ps))
    assert ps[1] / ps[0] < 2.0            # 4.5x more noise, <2x more flipping

    # contrast: at strong drive the same noise change moves p by >1e6
    strong = [flip_probability(0.05, 96, 16.0, nf, depth=40)
              for nf in (0.10, 0.45)]
    assert strong[1] / strong[0] > 1e6
