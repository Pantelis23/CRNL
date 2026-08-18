"""FINDINGS §89: the reduction's error is the fast-pair term it drops. Computed, not fitted."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.deep_barriers import MEASURED_A, MEASURED_A_M
from experiments.the_dropped_term import dropped_term


class TestTheDroppedTermExplainsTheOvershoot:
    def test_C_has_the_right_sign_and_size_at_every_gamma(self):
        for g, m in MEASURED_A.items():
            d = dropped_term(g, n=60)
            over = d["A"] - m
            assert over > 0                      # the reduction overshoots
            assert d["C"] > 0                    # and the dropped term subtracts
            assert 0.8 < d["C"] / over < 1.1     # same size, no fit

    def test_subtracting_C_corrects_the_action_to_a_fifth_of_a_percent(self):
        for g, m in MEASURED_A.items():
            d = dropped_term(g, n=60)
            raw = abs(d["A"] / m - 1)
            cor = abs((d["A"] - d["C"]) / m - 1)
            assert cor < raw / 4
            assert cor < 0.003

    def test_it_reproduces_the_non_one_over_M_shape_without_fitting_a_power(self):
        """§88's overshoot falls but not like 1/M and RISES at the first step. The correction
        must follow that shape, including the rise."""
        ms = sorted(MEASURED_A_M[0.40])
        over, cor = [], []
        for M in ms:
            d = dropped_term(0.40, float(M), n=60)
            m = MEASURED_A_M[0.40][M]
            over.append(d["A"] - m)
            cor.append(abs((d["A"] - d["C"]) / m - 1))
        assert over[1] > over[0]                 # the raw overshoot RISES first
        assert all(c < 0.003 for c in cor)       # the corrected value does not care

    def test_the_residual_is_flat_which_is_why_it_is_only_a_suspect(self):
        """§89.1: the correction removes the gamma-structure and leaves a constant, which
        points at the measurement side. Recorded as a suspect, not explained."""
        r = [(dropped_term(g, n=60)["A"] - dropped_term(g, n=60)["C"]) / m
             for g, m in sorted(MEASURED_A.items())]
        assert all(v > 1.0 for v in r)
        assert max(r) - min(r) < 0.002           # flat across a 7x spread in raw overshoot
