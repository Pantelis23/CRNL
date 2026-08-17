"""FINDINGS §83: identical ODEs, different reliability.

The neutrality gate is the section, so it is tested hardest: if a pair perturbs the drift at
all, every number in §83 is measuring a landscape change instead.
"""
from __future__ import annotations

import numpy as np
import pytest

from experiments.ode_does_not_determine_it import (
    A_mfpt, A_quad, am_A, am_neutral, drift, eta_lna, lam_mu,
)


class TestTheNeutralPairIsActuallyNeutral:
    def test_schlogl_drift_is_untouched_to_machine_precision(self):
        xs = np.linspace(0.02, 2.5, 61)
        f0 = drift(xs, 0.1, 1.0, 1.9, 0.0)
        for c in (0.5, 2.0, 20.0, 1000.0):
            assert np.abs(drift(xs, 0.1, 1.0, 1.9, c) - f0).max() < 1e-11

    def test_am_drift_is_untouched_on_the_whole_simplex(self):
        rng = np.random.default_rng(1)
        base = am_neutral(0.35, 0.0)
        S0, f = base.stoichiometry_matrix(), base.fluxes
        for _ in range(100):
            x = rng.dirichlet([1.0, 1.0, 1.0])
            f0 = S0 @ f(x)
            for c in (0.25, 1.0, 10.0):
                net = am_neutral(0.35, c)
                assert np.abs(net.stoichiometry_matrix() @ net.fluxes(x) - f0).max() < 1e-12

    def test_am_neutral_pair_preserves_the_conservation_law(self):
        net = am_neutral(0.35, 0.7)
        assert np.abs(np.ones(3) @ net.stoichiometry_matrix()).max() < 1e-12

    def test_but_the_propensities_do_change(self):
        """The point of the section: the drift is invariant, lam and mu are not."""
        l0, m0 = lam_mu(1.9, 0.1, 1.0, 1.9, 0.0)
        l1, m1 = lam_mu(1.9, 0.1, 1.0, 1.9, 5.0)
        assert (l1 - l0) == pytest.approx(m1 - m0, rel=1e-12)   # same shift in both
        assert l1 / l0 > 1.5 and m1 / m0 > 1.5                  # and it is a large one


class TestTheMeasurement:
    @pytest.mark.parametrize("c", [0.0, 2.0, 20.0])
    def test_quadrature_matches_the_exact_first_passage_at_every_c(self, c):
        """Rule 16, absolute: a quadrature drifting with c would fabricate the whole result."""
        Aq, err = A_quad(0.1, 1.0, 1.9, c)
        assert err < 1e-9
        assert Aq == pytest.approx(A_mfpt(0.1, 1.0, 1.9, c), rel=5e-3)

    def test_A_falls_monotonically_at_a_fixed_ODE(self):
        As = [A_quad(0.1, 1.0, 1.9, c)[0] for c in (0.0, 0.5, 1.0, 2.0, 5.0, 20.0)]
        assert all(np.diff(As) < 0)
        assert As[0] / As[-1] > 5

    def test_eta_falls_in_the_same_direction(self):
        """Two independent routes -- quadrature over the barrier, Lyapunov at the rail. A
        disagreement in DIRECTION would mean one of them is computed wrongly."""
        etas = [eta_lna(0.1, 1.0, 1.9, c)[0] for c in (0.0, 1.0, 5.0, 20.0)]
        assert all(np.diff(etas) < 0)

    def test_am_collapses_the_same_way_and_agrees_with_82(self):
        a0, _ = am_A(0.35, 0.0)
        a1, _ = am_A(0.35, 0.5)
        assert a0 == pytest.approx(0.07158, abs=2e-4)    # §82's stored value, same instrument
        assert a1 < a0 and a0 / a1 > 1.5
