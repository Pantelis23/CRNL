"""FINDINGS §94: the width sequence does not follow the LNA recursion.

The §94 numbers are at Omega = 30, where a D = 3 solve is ~10 minutes -- far too slow for the
suite. These tests protect the machinery and the QUALITATIVE structure at a small Omega, plus the
deterministic-gain facts, which are instant. The Omega = 30 quantities live in FINDINGS and are
reproduced by `experiments/depth_compounding.py`.
"""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse.linalg as spla

import experiments.chemical_cascade as cc
from experiments.cascade_schlogl import schlogl_consts
from experiments.depth_compounding import build_chain, seed, stage_stats
from experiments.margin_law import R1, R2, R3, stage1_stationary

C = schlogl_consts(R1, R2, R3)
OM = 14
_MEMO: dict = {}


def _run(D, all_reflected=False, om=OM, t=2.0):
    key = (D, all_reflected, om, t)
    if key not in _MEMO:
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        _, pi1 = stage1_stationary(om)
        Q, ref, dims, strides, cap = build_chain(om, D, all_reflected=all_reflected)
        _MEMO[key] = (spla.expm_multiply(Q.T * t, seed(om, ref, dims, strides, pi1)),
                      ref, dims, strides)
    return _MEMO[key]


def _hi(x):
    cc.HILL_N, cc.HILL_K = 4.0, 1.0
    r = cc.downstream_roots(x, C, R3, "hill")
    return r[-1] if len(r) else np.nan


def _g(x, h=1e-4):
    return abs((_hi(x + h) - _hi(x - h)) / (2 * h))


class TestTheCouplingIsOneWay:
    def test_a_downstream_stage_changes_nothing_upstream(self):
        """P1: without this, a depth comparison is meaningless."""
        p2, r2, d2, s2 = _run(2)
        p3, r3_, d3, s3 = _run(3)
        a = stage_stats(p2, OM, r2, d2, s2, 0)
        b = stage_stats(p3, OM, r3_, d3, s3, 0)
        assert a[0] == pytest.approx(b[0], abs=1e-9)
        assert a[1] == pytest.approx(b[1], abs=1e-9)


class TestTheWidthSequence:
    def test_widths_rise_monotonically_and_the_mean_drifts_down(self):
        p, ref, dims, strides = _run(3, all_reflected=True)
        st = [stage_stats(p, OM, ref, dims, strides, k, all_reflected=True)
              for k in range(3)]
        sd = [x[1] for x in st]
        mu = [x[0] for x in st]
        assert sd[1] > sd[0] and sd[2] > sd[1]
        assert mu[1] < mu[0] and mu[0] < R3          # drifting down from the rail

    def test_the_LNA_recursion_predicts_sigma3(self):
        """§96.1(a) corrected this. The original test asserted that sigma_3 OVERSHOOTS the
        fixed point, which was the seed bug: the last stage was being seeded at its saddle and
        relaxing upward, inflating its width. With every stage seeded at its rail the recursion
        predicts sigma_3 to 1.3% here and 4.3% at Omega = 30, and sigma_3 sits essentially AT
        the fixed point rather than far past it."""
        p, ref, dims, strides = _run(3, all_reflected=True)
        sd = [stage_stats(p, OM, ref, dims, strides, k, all_reflected=True)[1]
              for k in range(3)]
        g2 = (sd[1] / sd[0]) ** 2 - 1.0
        assert 0 < g2 < 1.0
        pred = sd[0] * np.sqrt(1.0 + g2 + g2 ** 2)
        assert pred == pytest.approx(sd[2], rel=0.05)
        assert sd[2] < 1.05 * sd[0] / np.sqrt(1 - g2)      # at the fixed point, not past it
