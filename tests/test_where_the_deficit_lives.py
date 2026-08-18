"""FINDINGS §85: locating §84's deficit. Two suspects refuted, the third confirmed."""
from __future__ import annotations

import numpy as np
import pytest

from crnl.networks.am_reversible import delta_star
from experiments.nu_is_two import A_nf
from experiments.where_the_deficit_lives import (
    MEASURED, A_fp, A_wkb, am_fast, rates_u, slow_manifold, u_star,
)


class TestTheReductionIsFaithful:
    """P1: any deficit must be the reduction, not a bug in it."""

    @pytest.mark.parametrize("g", [0.30, 0.40, 0.46])
    def test_slow_manifold_hits_b_star_at_the_attractor(self, g):
        assert slow_manifold(float(delta_star(g)), g) == pytest.approx(g / (1 + g), abs=1e-14)

    @pytest.mark.parametrize("g", [0.30, 0.40, 0.46])
    def test_reduced_drift_equals_the_exact_factored_drift_everywhere(self, g):
        """Not just at the fixed points -- identically in u."""
        for u in np.linspace(0.01, float(delta_star(g)) * 0.99, 10):
            lam, mu = rates_u(u, g)
            exact = (1 + g) * u * (slow_manifold(u, g) - g / (1 + g))
            assert (lam - mu) == pytest.approx(exact, abs=1e-14)


class TestSuspectsOneAndTwoAreRefuted:
    def test_the_exact_manifold_is_worse_than_84s_truncation(self):
        """S1 refuted, and backwards: truncation is not the deficit."""
        for g in (0.30, 0.40, 0.46):
            m = MEASURED[g]
            assert abs(A_nf(g) / m - 1) < abs(A_fp(g) / m - 1)

    def test_wkb_and_fokker_planck_agree_far_inside_the_deficit(self):
        """S2 refuted: lam_u/mu_u stays near 1, so the cubic term never gets going."""
        for g in (0.30, 0.40, 0.46):
            assert abs(A_wkb(g) / A_fp(g) - 1) < 0.01
            assert abs(A_wkb(g) / A_fp(g) - 1) < abs(1 - A_fp(g) / MEASURED[g]) / 5

    def test_wkb_exceeds_fp_as_the_cubic_term_requires(self):
        """ln r - 2(r-1)/(r+1) ~ (r-1)^3/12 > 0: the sign is right even though the size is not."""
        for g in (0.30, 0.40, 0.46):
            assert A_wkb(g) > A_fp(g)

    def test_84s_closer_agreement_is_a_cancellation_not_an_improvement(self):
        gs = sorted(MEASURED)
        r = [A_nf(g) / A_fp(g) for g in gs]
        assert not (all(np.diff(r) > 0) or all(np.diff(r) < 0))   # non-monotone


class TestTheSpeedUpIsNeutralWhereItMustBe:
    """P4's test is only decisive if M really does leave the lead alone."""

    def test_the_scaled_pair_does_not_change_the_leads_rates_at_fixed_b(self):
        """f1/r1 have Delta u = 0. rates_u depends on M ONLY through the slow manifold."""
        u, g = 0.3, 0.40
        b = slow_manifold(u, g, 1.0)
        s = 1 - b
        x, y = (s + u) / 2, (s - u) / 2
        expect = (b * x + g * y * y, g * x * x + b * y)
        assert rates_u(u, g, 1.0) == pytest.approx(expect, abs=1e-14)

    def test_M_leaves_the_cycle_affinity_exactly_unchanged(self):
        from crnl.networks.am_reversible import cycle_affinity, reverse_pairing
        base = am_fast(0.40, 1.0)
        a0 = cycle_affinity(base, reverse_pairing(base))
        for M in (2.0, 8.0, 16.0):
            net = am_fast(0.40, M)
            assert cycle_affinity(net, reverse_pairing(net)) == pytest.approx(a0, abs=1e-13)

    def test_M_does_move_the_landscape_so_only_the_ratio_is_comparable(self):
        """Stated because it would be an error to read A across M as one system."""
        assert u_star(0.40, 8.0) > u_star(0.40, 1.0) * 1.05
