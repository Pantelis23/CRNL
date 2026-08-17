"""FINDINGS §84: nu = 2 exactly, from the pitchfork normal form with no fitted parameter."""
from __future__ import annotations

import numpy as np
import pytest

from crnl.networks.am_reversible import GAMMA_C, delta_star, lambda_antisym
from experiments.nu_is_two import A_nf, b2, effective_nu, eps, nf_delta_star, omegas_for


class TestTheEliminationIsRight:
    """P1. No measurement can rescue a wrong elimination, so this is checked against closed
    forms the module computes independently."""

    @pytest.mark.parametrize("g", [0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.49, 0.499])
    def test_the_normal_form_drift_is_exactly_lambda_antisym(self, g):
        assert (1.0 + g) * eps(g) == pytest.approx(lambda_antisym(g), abs=1e-15)

    def test_nf_delta_star_converges_to_the_exact_one(self):
        """Rule 20: it is a LEADING-ORDER form. It must converge, not sit within a tolerance --
        at gamma = 0 it is 18% off and that is not a failure."""
        resid = [abs(nf_delta_star(g) / float(delta_star(g)) - 1.0)
                 for g in (0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.49, 0.499)]
        assert all(np.diff(resid) < 0)
        assert resid[0] > 0.15 and resid[-1] < 1e-3

    def test_the_saddle_is_at_one_third_for_every_gamma(self):
        """b0 = 1/3 exactly for all gamma is what makes eps closed-form; it is the reason the
        prefactor cancels and A has no free parameter."""
        for g in (0.05, 0.25, 0.45, 0.499):
            h = -(1 + g) / 2 + (2 + g) / 3 + (3 / 2) * (g - 1) / 9
            assert h == pytest.approx(0.0, abs=1e-15)


class TestNuIsExactlyTwo:
    def test_A_nf_vanishes_quadratically_at_gamma_c(self):
        """With gamma = 1/2 - d, (1-2g) = 2d, so A -> 2*(2d)^2/(2*0.5*1.5^2) = (8/2.25) d^2."""
        for d in (1e-2, 1e-3, 1e-4):
            assert A_nf(GAMMA_C - d) / d ** 2 == pytest.approx(8.0 / 2.25, rel=0.02)

    def test_the_local_exponent_converges_to_2_from_below(self):
        nus = [effective_nu(np.linspace(GAMMA_C - 2 * w, GAMMA_C - w, 9),
                            [A_nf(g) for g in np.linspace(GAMMA_C - 2 * w, GAMMA_C - w, 9)])
               for w in (0.10, 0.03, 0.01, 0.003)]
        assert all(v < 2.0 for v in nus)
        assert all(np.diff(nus) > 0)
        assert nus[-1] > 1.99

    def test_every_finite_window_reads_below_2(self):
        """§63.2 excluded 2 over [0.20, 0.45]. An effective exponent under 2 is exactly what
        nu = 2 with a correction to scaling produces, so that window cannot exclude it."""
        gs = np.linspace(0.20, 0.45, 9)
        assert effective_nu(gs, [A_nf(g) for g in gs]) < 1.9


class TestTheInstrumentWindow:
    def test_omega_band_keeps_the_barrier_inside_the_usable_range(self):
        """P5: too small and WKB has not started; too large and ln T nears the ~35 ceiling
        that turned one local slope NEGATIVE at Omega = 650 in scouting."""
        for g in (0.30, 0.38, 0.44):
            prod = [A_nf(g) * om for om in omegas_for(g)]
            assert min(prod) > 4.0 and max(prod) < 22.0

    def test_b2_is_negative_so_the_pitchfork_is_supercritical(self):
        for g in (0.0, 0.25, 0.49):
            assert b2(g) < 0
