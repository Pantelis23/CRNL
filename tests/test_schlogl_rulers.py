"""§67's verdict rule and instrument, tested before the experiment ran (§66's convention)."""
from __future__ import annotations

import numpy as np
import pytest

from crnl.networks.am_reversible import cycle_affinity, reverse_pairing
from experiments.schlogl_rulers import (AM_FLOOR, birth_death, cell, schlogl,
                                        sigma_local, verdict_transfer)


def test_verdict_reaches_both_branches():
    assert verdict_transfer(2.05, 2.0, 0.25, "G")[0] == "transfers"
    assert verdict_transfer(7.10, 2.0, 0.25, "G")[0] == "does-not"
    # exactly at the tolerance edge counts as transferring, and does not crash
    assert verdict_transfer(2.5, 2.0, 0.25, "G")[0] == "transfers"


def test_schlogl_fixed_points_are_where_they_were_placed():
    for x0, m in ((1.0, 0.5), (2.0, 0.4), (0.7, 0.3)):
        net = schlogl(x0, m)
        for x in (x0 - m, x0, x0 + m):
            assert float(net.rhs(np.array([x]))[0]) == pytest.approx(0.0, abs=1e-12)


def test_affinity_floor_is_two_ln_three_and_independent_of_x0():
    """FINDINGS §67: Schloegl's affinity floor at the death of bistability is 2 ln 3.

    Checked against the engine's `cycle_affinity`, which takes the null space of the
    per-pair stoichiometry and knows nothing of the closed form.
    """
    for x0 in (0.4, 1.0, 2.5):
        net = schlogl(x0, 1e-7)
        A = cycle_affinity(net, reverse_pairing(net))
        assert A == pytest.approx(2.0 * np.log(3.0), abs=1e-10), x0
    assert 2.0 * np.log(3.0) != pytest.approx(AM_FLOOR, abs=1e-3)   # NOT the same as AM's


def test_sigma_local_is_nonnegative_and_vanishes_only_at_equilibrium():
    """Entropy production is nonnegative pointwise; a sign error here would be invisible."""
    s = sigma_local(400, 1.0, 0.5, cap=1200)
    assert (s[3:] >= -1e-12).all(), s.min()
    assert s[3:].max() > 0


def test_chain_drift_differs_from_the_field_at_exactly_first_order():
    """The chain uses falling factorials, mass action uses powers: they differ by O(1/Omega).

    That gap IS the discreteness this project measures. A test demanding equality would be
    demanding something false -- which the first version of §67's P1 gate did.
    """
    x0, m = 1.0, 0.5
    net = schlogl(x0, m)
    scaled = []
    for om in (200, 400, 800, 1600):
        lam, mu, _ = birth_death(om, x0, m, cap=int(3 * (x0 + m) * om))
        n = int(1.2 * x0 * om)
        gap = abs((lam[n] - mu[n]) / om - float(net.rhs(np.array([n / om]))[0]))
        assert gap > 0                                   # it must NOT be zero
        scaled.append(om * gap)
    assert np.ptp(scaled) / np.mean(scaled) < 0.05, scaled


def test_cap_is_not_a_reflecting_wall():
    """RULE 10: a cap would push probability back and manufacture restoration."""
    a = cell(200, 1.0, 0.5, 0.35, 0.80, cap_mult=3.0)
    b = cell(200, 1.0, 0.5, 0.35, 0.80, cap_mult=6.0)
    for k in ("p_down", "mean_T", "Sigma", "Q"):
        assert a[k] == pytest.approx(b[k], rel=1e-9), k
