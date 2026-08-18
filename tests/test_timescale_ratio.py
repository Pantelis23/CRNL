"""FINDINGS §93: the penalty depends on the timescale ratio; the collapse is algebra."""
from __future__ import annotations

import numpy as np
import pytest

import experiments.chemical_cascade as cc
from experiments.margin_law import R2, R3, stage1_stationary
from experiments.timescale_ratio import penalty, pinned_reference


class TestTheProperTimeGate:
    def test_the_diagonal_is_constant_to_machine_precision(self):
        """P1: with the window in the downstream's own clock, s_up = s_dn = c is the SAME
        physics relabelled. Anything else is a bug, not an effect."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        _, pi1 = stage1_stationary(30)
        vals = [penalty(30, c, c, 2.0, pi1)[0] for c in (0.25, 1.0, 4.0)]
        assert max(vals) / min(vals) < 1.0001

    def test_wall_time_normalisation_is_NOT_flat(self):
        """§93 P1's failure, kept: at fixed WALL time the diagonal runs 26.3 -> 2.9, because
        scaling the downstream changes how much of its own clock has elapsed."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        _, pi1 = stage1_stationary(30)
        # penalty() divides t0 by s_dn, so to get the old behaviour pass t0 = 2*s_dn
        a = penalty(30, 0.25, 0.25, 2.0 * 0.25, pi1)[0]
        b = penalty(30, 4.0, 4.0, 2.0 * 4.0, pi1)[0]
        assert a / b > 5


class TestTheCollapseIsAnIdentity:
    def test_equal_ratios_agree_exactly(self):
        """Q t = t0[(s_up/s_dn) Q1 + Q2]. A wiring check: it would catch the coupling being
        placed in the wrong block."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        _, pi1 = stage1_stationary(30)
        a = penalty(30, 1.0, 4.0, 2.0, pi1)[0]
        b = penalty(30, 4.0, 16.0, 2.0, pi1)[0]
        assert a == pytest.approx(b, rel=1e-9)


class TestTheShapeIsThePhysics:
    def test_plateau_below_ratio_one_and_fall_above(self):
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        _, pi1 = stage1_stationary(30)
        slow = [penalty(30, r, 1.0, 2.0, pi1)[0] for r in (0.0625, 0.25, 1.0)]
        fast = penalty(30, 16.0, 1.0, 2.0, pi1)[0]
        assert max(slow) / min(slow) < 1.10        # plateau
        assert fast < 0.5 * max(slow)              # motional narrowing
        assert slow[-1] == max(slow)               # turns at ratio 1

    def test_the_pinned_reference_does_not_depend_on_the_downstream_clock(self):
        """With the window in the downstream's proper time the reference is invariant --
        which is what makes the penalty a clean ratio."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        a = pinned_reference(30, 1.0, 2.0)
        b = pinned_reference(30, 4.0, 2.0 / 4.0)
        assert a == pytest.approx(b, rel=1e-9)
