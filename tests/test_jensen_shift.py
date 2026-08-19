"""FINDINGS §95: the width reduces to the mean; curvature supplies 42% of the mean shift."""
from __future__ import annotations

import numpy as np
import pytest

import experiments.chemical_cascade as cc
from experiments.jensen_shift import F, MU1, MU2, SD1, SD2, d2F, dF, lna_width
from experiments.margin_law import R3


class TestTheMap:
    def test_neutral_at_the_rail(self):
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        assert F(R3) == pytest.approx(R3, abs=1e-12)

    def test_the_derivatives_are_converged_in_the_step(self):
        """Rule 13: the finite-difference step is a second axis, not a detail."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        vals = [d2F(MU1, h) for h in (1e-2, 1e-3, 5e-4)]
        assert abs(vals[-1] - vals[-2]) < 0.02 * abs(vals[-1])

    def test_the_map_is_concave_at_the_operating_point(self):
        """Without concavity Jensen pushes the other way and the mechanism is dead."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        assert d2F(MU1, 1e-3) < 0
        assert 0 < dF(MU1, 1e-3) < 1          # attenuating, as a restoring element must


class TestTheMeanShift:
    def test_jensen_has_the_right_sign_and_closes_part_of_the_gap(self):
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        d_intr = MU1 - R3
        no_j = F(MU1) + d_intr
        with_j = no_j + 0.5 * d2F(MU1, 1e-3) * SD1 ** 2
        assert with_j < no_j                                   # lowers the operating point
        assert abs(MU2 - with_j) < abs(MU2 - no_j)             # and closes some of the gap
        frac = 1 - abs(MU2 - with_j) / abs(MU2 - no_j)
        assert 0.3 < frac < 0.6                                # partial, ~42%, not all of it

    def test_it_is_only_partial_and_that_is_the_claim(self):
        """§95 reports 42%, not 100%. If a later change makes it look complete, that is a
        signal something was tuned."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        with_j = F(MU1) + (MU1 - R3) + 0.5 * d2F(MU1, 1e-3) * SD1 ** 2
        assert abs(MU2 - with_j) > 0.01                        # a real residual remains


class TestTheWidthFollowsTheMean:
    def test_the_LNA_at_the_measured_mean_gives_the_measured_width(self):
        """§95's useful half: the width is not a separate problem."""
        assert lna_width(MU2, 30) == pytest.approx(SD2, rel=0.03)

    def test_the_LNA_error_on_stage_one_is_comparable(self):
        """The 1.5% residual at stage 2 is inside the LNA's own accuracy, measured on a stage
        whose input is a noiseless chemostat."""
        assert abs(lna_width(MU1, 30) / SD1 - 1) < 0.05

    def test_moving_the_operating_point_down_widens_the_stage(self):
        assert lna_width(2.90, 30) > lna_width(3.05, 30)
