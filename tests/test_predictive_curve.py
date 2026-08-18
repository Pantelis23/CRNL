"""FINDINGS §87: the escape curve from the tilted generator, with no stationary solve."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.predictive_curve import H, action, curve, fast_fixed, _grad
from experiments.where_the_deficit_lives import MEASURED, slow_manifold, u_star


class TestTheZeroMomentumSheet:
    def test_H_vanishes_identically_at_zero_momentum(self):
        """P1, and it is why the joint solve collapsed: the whole p=0 sheet is a root
        manifold, for EVERY s, not just the deterministic one."""
        for g in (0.40, 0.44):
            for u in np.linspace(0.05, 0.45, 5):
                for s in (0.5, 0.62, 0.7, 0.8):
                    assert H(u, s, 0.0, 0.0, g) == pytest.approx(0.0, abs=1e-14)

    def test_zero_momentum_flow_recovers_the_deterministic_slow_manifold(self):
        for g in (0.40, 0.44):
            for u in np.linspace(0.05, 0.40, 5):
                sdet = 1.0 - slow_manifold(u, g)
                assert _grad(u, sdet, 0.0, 0.0, g)[0] == pytest.approx(0.0, abs=1e-8)


class TestThePrediction:
    @pytest.mark.parametrize("g", [0.40, 0.44])
    def test_the_tilted_action_beats_the_deterministic_manifold(self, g):
        A, rows = action(g, n=60)
        m = MEASURED[g]
        assert len(rows) > 40
        assert abs(1 - A / m) < 0.03
        det = {0.40: 0.9027, 0.44: 0.9348}[g]
        assert abs(1 - A / m) < abs(1 - det)

    def test_the_action_is_grid_converged(self):
        """Rule 20: the residual must be shown to converge, not gated at one resolution."""
        vals = [action(0.40, n)[0] for n in (30, 60, 120)]
        assert abs(vals[-1] - vals[-2]) < abs(vals[1] - vals[0])
        assert abs(vals[-1] / vals[-2] - 1) < 1e-3

    def test_the_escape_momentum_is_negative(self):
        """§87.1: the first scan swept pu > 0 and returned 'no solution' -- a search that
        cannot reach the answer's sign reports cleanly and looks like a result."""
        rows = curve(0.40, n=20)
        assert rows and all(r["pu"] < 0 for r in rows)
        mid = rows[len(rows) // 2]
        assert 0.02 < abs(mid["pu"]) < 0.5

    def test_the_instanton_curve_is_displaced_OPPOSITE_to_86s_ridge(self):
        """§87 P2: §86's ridge is displaced UP (+0.00254 at gamma=0.40); the instanton's
        curve is displaced DOWN. The measurement in §86 stands; its mechanism does not."""
        for g, ridge_up in ((0.40, 0.00254), (0.44, 0.00104)):
            rows = curve(g, n=40)
            disp = float(np.mean([r["b"] - r["b_det"] for r in rows]))
            assert ridge_up > 0
            assert disp < 0
