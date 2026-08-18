"""FINDINGS §91: a chemically-coupled cascade. Transmission is half the gate."""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse.linalg as spla

import experiments.chemical_cascade as cc
from experiments.cascade_schlogl import schlogl_consts

R1, R2, R3 = cc.RAILS
C = schlogl_consts(R1, R2, R3)


class TestBothHalvesOfTheGate:
    @pytest.mark.parametrize("scheme", ["source", "catalytic", "hill"])
    def test_neutral_at_the_rail(self, scheme):
        """All three reproduce the isolated element when the upstream is correct."""
        for nn in (5, 30, 60, 90):
            a, b = cc.rates_stage(float(nn), R3 * 30, 30, C, R3, False, scheme)
            a0, b0 = cc.rates_stage(float(nn), 0.0, 30, C, R3, True, scheme)
            assert a == pytest.approx(a0, abs=1e-12)
            assert b == pytest.approx(b0, abs=1e-12)

    def test_source_coupling_does_not_transmit(self):
        """§91.1: it passes neutrality and still carries no signal -- the downstream keeps
        its high rail when the upstream is at r1."""
        assert len(cc.downstream_roots(R1, C, R3, "source")) == 3

    @pytest.mark.parametrize("scheme", ["catalytic", "hill"])
    def test_the_others_transmit(self, scheme):
        lo = cc.downstream_roots(R1, C, R3, scheme)
        assert len(lo) == 1 and lo[0] < R2


class TestTheArtifactIsReproducible:
    def test_source_coupling_shows_fake_filtering(self):
        """The false result §91.1 caught: sublinear accumulation from disconnection.
        P(down low | up low) is what distinguishes it from real filtering."""
        Q, cap, m, strides, _ = cc.build(2, 30, "source")
        p = spla.expm_multiply(Q.T * 2.0, cc.seed_high(2, 30, m, strides, R3))
        mk = cc.masks(2, m, strides, 30, R2)
        P1, PD = float(p[mk[0]].sum()), float(p[mk[1]].sum())
        cond = float(p[mk[1] & mk[0]].sum()) / P1
        assert PD / (1 - (1 - P1) ** 2) < 0.7      # looks like filtering
        assert cond < 0.05                          # but nothing propagated


class TestTheMarginControlsIt:
    def _run(self, n, K):
        cc.HILL_N, cc.HILL_K = n, K
        xs = np.linspace(R1, R3, 2001)
        xc = next(x for x in xs[::-1] if len(cc.downstream_roots(x, C, R3, "hill")) < 3)
        Q, cap, m, strides, _ = cc.build(2, 30, "hill")
        p = spla.expm_multiply(Q.T * 2.0, cc.seed_high(2, 30, m, strides, R3))
        mk = cc.masks(2, m, strides, 30, R2)
        P1 = float(p[mk[0]].sum())
        own = float(p[mk[1] & ~mk[0]].sum()) / P1
        return float(xc), own

    def test_penalty_falls_with_margin_over_an_order_of_magnitude(self):
        try:
            xc_lo, own_lo = self._run(4.0, 1.6)     # small margin
            xc_hi, own_hi = self._run(4.0, 0.6)     # large margin
        finally:
            cc.HILL_N, cc.HILL_K = 4.0, 1.0
        assert xc_lo > xc_hi                        # smaller margin
        assert own_lo / own_hi > 8                  # far worse composition

    def test_the_two_knobs_agree_at_a_matched_margin(self):
        """The collapse: moving the margin by exponent or by half-max must give the same
        penalty. Compared at nearly equal x_crit, not across different ranges (§91.3)."""
        try:
            xc_n, own_n = self._run(8.0, 1.0)
            xc_k, own_k = self._run(4.0, 0.8)
        finally:
            cc.HILL_N, cc.HILL_K = 4.0, 1.0
        assert abs(xc_n - xc_k) < 0.05              # same margin, different knob
        assert own_n == pytest.approx(own_k, rel=0.15)
