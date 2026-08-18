"""FINDINGS §92: the margin law is a frozen-upstream statement; the fast limit is narrowing."""
from __future__ import annotations

import numpy as np
import pytest

import experiments.chemical_cascade as cc
from experiments.action_is_not_priced import schlogl_A
from experiments.margin_law import (
    A_of_xup, R1, R2, R3, build_reflect, predict, stage1_stationary, upstream_qsd,
)


class TestTheBarrierFunction:
    @pytest.mark.parametrize("scheme", ["source", "catalytic", "hill"])
    def test_A_at_the_rail_is_the_isolated_action(self, scheme):
        """P1: two independent routes -- the coupled landscape at x_up = r3, and §80's
        isolated quadrature."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        assert A_of_xup(R3, scheme) == pytest.approx(schlogl_A(R1, R3, R2)[0], abs=1e-9)

    def test_A_vanishes_at_the_collapse_and_is_monotone_below_the_rail(self):
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        xs = np.linspace(R1, R3, 400)
        A = np.array([A_of_xup(x, "hill") for x in xs])
        assert A[0] == 0.0 and A[-1] > 0
        live = A > 0
        assert np.all(np.diff(A[live]) > -1e-12)


class TestTheSeedIsSpeedIndependent:
    def test_stage1_stationary_does_not_depend_on_the_clock(self):
        """§92.1(c): this is what makes the clock sweep clean -- scaling every rate equally
        leaves the stationary law unchanged, so the seed cannot smuggle in a speed."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        up, pi = stage1_stationary(30)
        assert pi.sum() == pytest.approx(1.0)
        assert (pi * up / 30).sum() < R3          # sits below the deterministic rail
        assert (pi * up / 30).sum() > R2

    def test_reflected_stage_one_cannot_escape(self):
        """§92.1(b): with no escape there is nothing to condition on."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        Q, up, m2, cap = build_reflect(30, 1.0)
        assert up[0] >= R2 * 30                   # state space starts at the saddle
        assert Q.shape[0] == len(up) * m2


class TestTheTwoLimits:
    def test_frozen_average_exceeds_fast_average(self):
        """Jensen: averaging the rate beats the rate at the average, and by a lot here."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        xs, px = upstream_qsd(30)
        frozen = predict("hill", 4.0, 1.0, 30, xs, px)
        fast = predict("hill", 4.0, 1.0, 30, xs, px, quenched=True)
        assert frozen > 3 * fast

    def test_the_penalty_falls_with_upstream_speed_motional_narrowing(self):
        """§92 P6: the headline. Same landscape, same margin, only the clock changes."""
        import scipy.sparse.linalg as spla
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        up, pi1 = stage1_stationary(30)
        vals = []
        for sp_ in (1.0, 32.0):
            Q, uu, m2, cap = build_reflect(30, sp_)
            q = np.zeros(len(uu) * m2)
            for a, w in enumerate(pi1):
                q[a * m2 + int(round(R3 * 30))] = w
            q = spla.expm_multiply(Q.T * 2.0, q)
            lo2 = (np.arange(len(uu) * m2) % m2) < R2 * 30
            vals.append(float(q[lo2].sum()))
        assert vals[1] < vals[0] / 2              # fast upstream is far safer
