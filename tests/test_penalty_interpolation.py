"""FINDINGS §98: the position between the two limits is a function of the barrier depth."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.penalty_interpolation import (
    A_CAL, NEW_AOM, NEW_FAST, NEW_FROZEN, NEW_MEAS, measure, position,
)

_MEMO: dict = {}


def _m(om):
    if om not in _MEMO:
        _MEMO[om] = measure(om)
    return _MEMO[om]


class TestThePosition:
    def test_it_is_zero_at_the_fast_limit_and_one_at_the_frozen(self):
        assert position(5.0, 1.0, 1.0) == pytest.approx(0.0)
        assert position(5.0, 1.0, 5.0) == pytest.approx(1.0)

    def test_it_falls_with_the_barrier_depth(self):
        """§98 P2: deeper barriers sit further from the frozen limit."""
        ps = [position(*_m(om)) for om in (20, 40, 70)]
        assert all(np.diff(ps) < 0)

    def test_the_bracket_is_not_guaranteed_at_a_shallow_barrier(self):
        """§98 P2's aside: at A*Omega = 2.7 the position exceeds 1 -- the two limits are
        asymptotic, and reporting them as a bracket everywhere would be wrong."""
        assert position(*_m(14)) > 1.0


class TestTheSecondElementLandsOnTheCurve:
    def test_out_of_sample_within_eight_percent(self):
        """§98 P3, and nothing is fitted: the curve comes from ONE element's Omega sweep and
        the other element's point is checked against it."""
        oms = [20, 30, 40, 55, 70, 85]
        aoms = np.array([A_CAL * o for o in oms])
        ps = np.array([position(*_m(o)) for o in oms])
        assert np.all(np.diff(aoms) > 0)                  # np.interp needs increasing x
        pred = float(np.interp(NEW_AOM, aoms, ps))
        assert abs(pred - position(NEW_FROZEN, NEW_FAST, NEW_MEAS)) < 0.08

    def test_reversed_x_would_have_returned_the_first_element(self):
        """§98.1: the bug that printed FAILS. np.interp with decreasing x returns fp[0], and
        the tell is that the 'prediction' is exactly a row of the table."""
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([10.0, 20.0, 30.0])
        assert np.interp(2.5, x[::-1], y[::-1]) == pytest.approx(y[0])
        assert np.interp(2.5, x, y) == pytest.approx(25.0)
