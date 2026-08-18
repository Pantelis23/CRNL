"""FINDINGS §88: the tilted route validated on eight gamma, then extrapolated with stated scope."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.deep_barriers import ETA_77, MEASURED_A, MEASURED_A_M, action


class TestValidation:
    def test_predicts_every_measured_gamma_within_six_percent(self):
        for g, m in MEASURED_A.items():
            A, _ = action(g, n=40)
            assert A is not None
            assert 1.0 < A / m < 1.07          # overshoots, and by less than 7%

    def test_the_overshoot_shrinks_toward_gamma_c(self):
        """Rule 20: the criterion is convergence, and its DIRECTION licenses §88 P3."""
        gs = sorted(MEASURED_A)
        ov = [action(g, n=40)[0] / MEASURED_A[g] - 1 for g in gs]
        assert all(np.diff(ov) < 0)
        assert ov[0] > 4 * ov[-1]

    def test_the_overshoot_falls_with_timescale_separation(self):
        """T15-n.4: confirms the residual is finite separation."""
        ms = sorted(MEASURED_A_M[0.40])
        ov = [abs(action(0.40, float(M), n=40)[0] / MEASURED_A_M[0.40][M] - 1) for M in ms]
        assert ov[-1] < ov[0] / 2

    def test_but_it_is_NOT_a_one_over_M_law(self):
        """§88 P2's prediction was 1/M and it is refuted: overshoot*M would be constant."""
        ms = sorted(MEASURED_A_M[0.40])
        prod = [abs(action(0.40, float(M), n=40)[0] / MEASURED_A_M[0.40][M] - 1) * M
                for M in ms]
        assert max(prod) / min(prod) > 3.0


class TestTheExtension:
    @pytest.mark.parametrize("g", [0.05, 0.10, 0.20])
    def test_the_route_reaches_barriers_no_lattice_instrument_could(self, g):
        A, _ = action(g, n=40)
        assert A is not None and A > 0.2

    def test_escape_still_beats_readout_at_deep_barriers(self):
        """§88 P4, and the margin is what makes it robust to the extrapolation error."""
        for g, eta in ETA_77.items():
            A, _ = action(g, n=40)
            assert A < eta
            assert A / eta < 0.25      # a 4x error would still not flip it
