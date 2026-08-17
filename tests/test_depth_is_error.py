"""FINDINGS §76: the depth ceiling is the per-stage error rate, exactly."""
from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import norm

from experiments.cascade_saturated import d_max_closed
from experiments.depth_is_error import b_star, c_star


def test_c_star_is_the_information_half_life_constant():
    from scipy.optimize import brentq

    def h(p):
        return -(p * np.log2(p) + (1 - p) * np.log2(1 - p))

    assert h((1 + b_star()) / 2) == pytest.approx(0.5, abs=1e-12)
    assert c_star() == pytest.approx(0.124266404564, abs=1e-11)


def test_depth_times_error_converges_to_c_star():
    """§76: D_max * eps -> c*, monotonically, over seven decades."""
    cs = c_star()
    dev = [abs(d_max_closed(e, e) * e / cs - 1) for e in (1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8)]
    assert all(np.diff(dev) < 0)
    assert dev[-1] < 1e-6


def test_the_constant_is_the_same_across_substrates():
    """§76: Schloegl and a bare step function give the same c*, to under 1%."""
    from experiments.cascade_saturated import eps_pair
    from experiments.ceiling_is_it_the_element import eps_from, grid, pc_step

    cs = c_star()
    vals = []
    for lam in (1.0, 4.0):
        e = eps_pair(3600, 0.1 * lam, 1.0 * lam, 1.9 * lam, 0.35 * 0.9 * lam)
        vals.append(d_max_closed(*e) * 0.5 * (e[0] + e[1]) / cs)
    x, _ = grid(0.1, 1.9, 1.6, 3600)
    for f in (0.45, 0.28):
        e = eps_from(pc_step(x, 0.1, 1.0, 1.9, 3600), x, 0.1, 1.9, f * 0.9)
        vals.append(d_max_closed(*e) * 0.5 * (e[0] + e[1]) / cs)
    assert max(abs(v - 1) for v in vals) < 0.01, vals


def test_arithmetic_mean_is_the_right_epsilon_for_an_asymmetric_channel():
    """The decay rate is 1 - e_hi - e_lo, so the SUM is what enters."""
    cs = c_star()
    e_hi, e_lo = 1e-3, 3e-3
    d = d_max_closed(e_hi, e_lo)
    ar = d * 0.5 * (e_hi + e_lo) / cs
    ge = d * np.sqrt(e_hi * e_lo) / cs
    ha = d * 2 / (1 / e_hi + 1 / e_lo) / cs
    assert abs(ar - 1) < abs(ge - 1) and abs(ar - 1) < abs(ha - 1)


def test_section_12s_factor_is_an_algebraic_prefactor_that_rises():
    """§76.1: the ratio to exp(z^2/2)/4 grows linearly in z -- it is NOT a constant 3."""
    cs = c_star()
    ratios = []
    for f in (0.45, 0.35, 0.28):
        z = 1.0 / f
        ratios.append((cs / float(norm.cdf(-z))) / (np.exp(z ** 2 / 2) / 4))
    assert all(np.diff(ratios) > 0), ratios
    # and it tracks the predicted 4 c* sqrt(2 pi) z to within 20%
    for f, r in zip((0.45, 0.35, 0.28), ratios):
        pred = 4 * cs * np.sqrt(2 * np.pi) / f
        assert abs(r - pred) / pred < 0.20, (f, r, pred)
