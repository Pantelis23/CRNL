"""FINDINGS §96: the mean shift is the static-transfer average; two §94/§95 bugs corrected."""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse.linalg as spla

import experiments.chemical_cascade as cc
from experiments.depth_compounding import build_chain, seed, stage_stats
from experiments.jensen_shift import F, MU1, MU2, SD1, d2F
from experiments.margin_law import R2, R3, stage1_stationary
from experiments.static_transfer_limit import G_of_h, build_two_reflected, stats

OM = 30


class TestTheTwoBugsStayFixed:
    def test_every_stage_seeds_at_its_rail_when_all_reflected(self):
        """§96.1(a): the last stage was seeded at a mean of 1.1333 against a rail at 3.1827 --
        essentially at its saddle -- and relaxed upward through the window."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        _, pi1 = stage1_stationary(OM)
        Q, ref, dims, strides, cap = build_chain(OM, 3, all_reflected=True)
        p0 = seed(OM, ref, dims, strides, pi1, all_reflected=True)
        for k in range(1, 3):
            mu, _ = stage_stats(p0, OM, ref, dims, strides, k, all_reflected=True)
            assert abs(mu - R3) < 0.05

    def test_a_reflected_stage_is_not_high_side_filtered(self):
        """§96.1(b): a reflected stage cannot escape, so filtering counts > R2*Omega merely
        drops its boundary site. Stage 1 must sit exactly at its stationary mean."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        ref1, pi1 = stage1_stationary(OM)
        exact = float((pi1 * ref1).sum()) / OM
        Q, ref, dims, strides, cap = build_chain(OM, 2, all_reflected=True)
        p = spla.expm_multiply(
            Q.T * 2.0, seed(OM, ref, dims, strides, pi1, all_reflected=True))
        mu, _ = stage_stats(p, OM, ref, dims, strides, 0, all_reflected=True)
        assert mu == pytest.approx(exact, abs=1e-6)

    def test_two_independent_builds_of_the_same_system_agree(self):
        """The P1 gate that caught (b): they disagreed, including on stage 1, which is
        provably invariant."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        _, pi1 = stage1_stationary(OM)
        Q, ref, dims, strides, cap = build_chain(OM, 2, all_reflected=True)
        pa = spla.expm_multiply(
            Q.T * 2.0, seed(OM, ref, dims, strides, pi1, all_reflected=True))
        Qb, refb, n = build_two_reflected(OM, 1.0)
        pb = np.zeros(n * n)
        pos = list(refb).index(int(round(R3 * OM)))
        for i, w in enumerate(pi1):
            pb[i * n + pos] = w
        pb = spla.expm_multiply(Qb.T * 2.0, pb)
        for k in range(2):
            a = stage_stats(pa, OM, ref, dims, strides, k, all_reflected=True)
            b = stats(pb, refb, n, k)
            assert a[0] == pytest.approx(b[0] / OM, abs=1e-6)
            assert a[1] == pytest.approx(b[1] / OM, abs=1e-6)


class TestTheStaticTransferAverage:
    def test_the_exact_frozen_average_beats_the_second_order_truncation(self):
        """§96.2: the truncation was 89% of §95's residual."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        ref1, pi1 = stage1_stationary(OM)
        xs = ref1 / OM
        d_intr = MU1 - R3
        exact = float(sum(w * F(x) for w, x in zip(pi1, xs) if np.isfinite(F(x)))) + d_intr
        trunc = F(MU1) + 0.5 * d2F(MU1, 1e-3) * SD1 ** 2 + d_intr
        assert abs(MU2 - exact) < 0.2 * abs(MU2 - trunc)
        assert abs(MU2 - exact) / MU2 < 0.005            # 0.12% of the mean

    def test_the_two_limits_are_different_and_bracket_the_measurement(self):
        """Frozen averages the OUTPUT, fast averages the INPUT -- Jensen on different
        functions, so they must differ, and the measurement must lie between."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        ref1, pi1 = stage1_stationary(OM)
        xs = ref1 / OM
        d_intr = MU1 - R3
        frozen = float(sum(w * F(x) for w, x in zip(pi1, xs) if np.isfinite(F(x)))) + d_intr
        fast = G_of_h(float(sum(w * cc.hill(x, R3) for w, x in zip(pi1, xs)))) + d_intr
        assert fast > frozen
        assert frozen < MU2 < fast

    def test_the_upstream_clock_moves_the_operating_point(self):
        """P2: if it did not, the correlation time would not be what the residual is made of."""
        cc.HILL_N, cc.HILL_K = 4.0, 1.0
        _, pi1 = stage1_stationary(OM)
        mus = []
        for s in (0.25, 32.0):
            Q, ref, n = build_two_reflected(OM, s)
            p = np.zeros(n * n)
            pos = list(ref).index(int(round(R3 * OM)))
            for i, w in enumerate(pi1):
                p[i * n + pos] = w
            p = spla.expm_multiply(Q.T * 2.0, p)
            mus.append(stats(p, ref, n, 1)[0] / OM)
        assert abs(mus[0] - mus[1]) > 0.002
