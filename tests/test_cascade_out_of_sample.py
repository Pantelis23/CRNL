"""FINDINGS §97: out of sample -- the derivations transfer, the fitted slope does not."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.cascade_out_of_sample import (
    A_new, F_new, HILL_K_NEW, HILL_N_NEW, RAILS_NEW, SLOPE_91, down_roots, hill_new,
    lna_width, stage1_law,
)


class TestTheElementIsGenuinelyDifferent:
    def test_it_is_not_the_calibration_element(self):
        assert RAILS_NEW != (0.15, 1.0, 3.1827)
        assert (HILL_N_NEW, HILL_K_NEW) != (4.0, 1.0)
        assert A_new() / 0.190241 > 2                 # a much deeper barrier

    def test_neutral_at_the_rail_and_it_transmits(self):
        """§91 P1(b): a neutrality gate alone would not catch a coupling that carries no
        signal, and this section's P1 checks both halves."""
        r1, r2, r3 = RAILS_NEW
        assert F_new(r3) == pytest.approx(r3, abs=1e-9)
        assert len(down_roots(r1)) == 1               # loses its high rail

    def test_the_margin_is_outside_the_calibrated_range(self):
        """Rule 19: §91's slope was traced over 1.81-4.70 sigma. Applying it here is an
        extrapolation, and that is part of why P4 failed."""
        r1, r2, r3 = RAILS_NEW
        up, pi1, cap = stage1_law(30)
        x1 = up / 30
        mu1 = float((pi1 * x1).sum())
        sd1 = float(np.sqrt((pi1 * (x1 - mu1) ** 2).sum()))
        xs = np.linspace(r1, r3, 2001)
        xc = next(float(x) for x in xs[::-1] if len(down_roots(x)) < 3)
        assert (r3 - xc) / sd1 > 4.70


class TestWhatTransfers:
    def test_the_static_transfer_mean_is_essentially_exact(self):
        """§97 P2: 0.00% on an element it was never calibrated on."""
        r1, r2, r3 = RAILS_NEW
        up, pi1, cap = stage1_law(30)
        x1 = up / 30
        mu1 = float((pi1 * x1).sum())
        pred = float(sum(w * F_new(x) for w, x in zip(pi1, x1)
                         if np.isfinite(F_new(x)))) + (mu1 - r3)
        assert pred == pytest.approx(4.26781, rel=2e-4)

    def test_the_lna_width_lands_inside_its_own_accuracy(self):
        assert lna_width(4.26779, 30) == pytest.approx(0.50236, rel=0.05)


class TestWhatDoesNot:
    def test_the_fitted_slope_is_far_wrong_here(self):
        """§97 P4. The slope is NOT re-fitted; it is applied as published and it misses."""
        r1, r2, r3 = RAILS_NEW
        up, pi1, cap = stage1_law(30)
        x1 = up / 30
        mu1 = float((pi1 * x1).sum())
        sd1 = float(np.sqrt((pi1 * (x1 - mu1) ** 2).sum()))
        xs = np.linspace(r1, r3, 2001)
        xc = next(float(x) for x in xs[::-1] if len(down_roots(x)) < 3)
        sd2 = lna_width(4.26779, 30)
        pen = float(np.exp(SLOPE_91 * ((r3 - xc) / sd2 - (r3 - xc) / sd1)))
        assert pen < 0.5 * 4.4742                     # 74% low

    def test_the_barrier_depth_differs_which_the_margin_law_cannot_see(self):
        """The axis §91 never swept: every one of its 14 points had the same A*Omega."""
        assert abs(A_new() * 30 - 13.2) < 0.5
        assert abs(0.190241 * 30 - 5.7) < 0.5
