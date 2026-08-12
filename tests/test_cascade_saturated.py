"""FINDINGS §72: the saturated-regime depth ceiling, and the instruments behind it."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.cascade_saturated import (d_max_saturated, eps_pair,
                                           mutual_info_depth, p_cross)


def test_p_cross_is_a_proper_splitting_probability():
    pc, lo, hi = p_cross(600, 0.1, 1.0, 1.9)
    assert pc[lo] == pytest.approx(1.0, abs=1e-12)
    assert pc[hi] == pytest.approx(0.0, abs=1e-12)
    assert (pc >= -1e-12).all() and (pc <= 1 + 1e-12).all()
    assert np.all(np.diff(pc[lo:hi + 1]) <= 1e-12)


def test_p_cross_vectorised_matches_the_quadratic_reference():
    """§72.1: the O(cap) suffix sum must equal the O(cap^2) per-state log-sum-exp."""
    from experiments.cascade_schlogl import rates, schlogl_consts

    om, r1, r2, r3 = 200, 0.1, 1.0, 1.9
    c = schlogl_consts(r1, r2, r3)
    cap = int(np.ceil(1.6 * r3 * om))
    lam, mu = rates(om, c, cap)
    lo, hi = int(round(r1 * om)), int(round(r3 * om))
    lp = np.full(cap + 1, -np.inf)
    acc = 0.0
    lp[lo] = 0.0
    for k in range(lo + 1, hi + 1):
        acc += np.log(mu[k]) - np.log(lam[k])
        lp[k] = acc

    def lse(v):
        v = v[np.isfinite(v)]
        m = v.max()
        return m + np.log(np.exp(v - m).sum())

    den = lse(lp[lo:hi])
    ref = np.array([1.0 if n <= lo else (0.0 if n >= hi
                                         else float(np.exp(lse(lp[n:hi]) - den)))
                    for n in range(cap + 1)])
    got, _, _ = p_cross(om, r1, r2, r3)
    assert np.abs(ref - got).max() < 1e-12


def test_asymmetric_channel_is_not_a_symmetric_one():
    """The element is asymmetric, so eps_hi != eps_lo and (1-2eps)^D would be wrong."""
    e_hi, e_lo = eps_pair(900, 0.1, 1.0, 1.9, 0.35 * 0.9)
    assert e_hi > 0 and e_lo > 0
    assert abs(e_hi - e_lo) / max(e_hi, e_lo) > 1e-3


def test_depth_ceiling_saturates_and_matches_am_within_tens_of_percent():
    """§72: ratio to exp(Delta^2/2 sigma^2)/4 lands on AM's, on a substrate sharing nothing."""
    D = 0.9
    for f, am_ratio in ((0.35, 50.0 / 14.8), (0.28, 489.0 / 147.0)):
        ds = [d_max_saturated(*eps_pair(om, 0.1, 1.0, 1.9, f * D))
              for om in (7200, 14400, 28800)]
        assert (max(ds) - min(ds)) / np.mean(ds) < 0.05          # saturated
        ratio = ds[-1] / (np.exp(1.0 / (2 * f ** 2)) / 4.0)
        assert abs(ratio - am_ratio) / am_ratio < 0.30, (f, ratio, am_ratio)


def test_mutual_information_is_monotone_and_bounded():
    I = mutual_info_depth(0.02, 0.03, 50)
    assert (I <= 1 + 1e-12).all() and (I >= -1e-12).all()
    assert np.all(np.diff(I) <= 1e-12)
