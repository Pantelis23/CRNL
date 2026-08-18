"""FINDINGS §90: the measured A's own ansatz spread exceeds §89's residual."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.deep_barriers import MEASURED_A
from experiments.richardson_A import LADDER, ln_T_ladder, richardson, slopes
from experiments.the_dropped_term import dropped_term


# The ladder costs several exact CME solves at Omega up to 1700; without this memo the
# parametrised cases re-run it a dozen times and the module times out.
_LADDER_MEMO: dict = {}
_DROP_MEMO: dict = {}


def _pts(g):
    if g not in _LADDER_MEMO:
        _LADDER_MEMO[g] = ln_T_ladder(g, LADDER[g])
    return _LADDER_MEMO[g]


def _drop(g):
    if g not in _DROP_MEMO:
        _DROP_MEMO[g] = dropped_term(g, n=60)
    return _DROP_MEMO[g]


def _cands(g):
    sl = slopes(_pts(g))
    x = np.array([1.0 / m for _, m in sl])
    y = np.array([s for s, _ in sl])
    return {"last": float(y[-1]),
            "lin_all": float(np.polyfit(x, y, 1)[-1]),
            "lin3": float(np.polyfit(x[-3:], y[-3:], 1)[-1]),
            "quad": float(np.polyfit(x, y, 2)[-1])}


class TestTheBiasIsReal:
    @pytest.mark.parametrize("g", [0.35, 0.40, 0.44])
    def test_local_slopes_are_linear_in_one_over_omega(self, g):
        """P2: there IS a 1/Omega bias to extrapolate, not scatter."""
        sl = slopes(_pts(g))
        x = np.array([1.0 / m for _, m in sl])
        y = np.array([s for s, _ in sl])
        assert abs(np.corrcoef(x, y)[0, 1]) > 0.9

    @pytest.mark.parametrize("g", [0.35, 0.40, 0.44])
    def test_the_two_point_slope_is_biased_low(self, g):
        A_rich, _ = richardson(_pts(g))
        assert A_rich > MEASURED_A[g]
        assert 0.001 < A_rich / MEASURED_A[g] - 1 < 0.006


class TestButItIsNotResolvable:
    @pytest.mark.parametrize("g", [0.35, 0.40, 0.44])
    def test_the_ansaetze_straddle_one(self, g):
        """Rule 15: the candidates disagree, so the quantity is unresolved and says so."""
        d = _drop(g)
        corrected = d["A"] - d["C"]
        ratios = [corrected / c for c in _cands(g).values()]
        assert min(ratios) < 1.0 < max(ratios)

    @pytest.mark.parametrize("g", [0.35, 0.40, 0.44])
    def test_the_ansatz_spread_exceeds_89s_residual(self, g):
        """The instrument cannot resolve the thing it was brought in to judge."""
        v = list(_cands(g).values())
        spread = max(v) / min(v) - 1
        d = _drop(g)
        residual_89 = abs((d["A"] - d["C"]) / MEASURED_A[g] - 1)
        assert spread > residual_89
