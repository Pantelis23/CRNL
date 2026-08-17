"""FINDINGS §81: the escape action redone -- integral vs exact, transfer, and the underflow."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.action_redo import A_am, A_integral, A_mfpt
from experiments.derive_eta import schlogl_V


def test_the_integral_predicts_A_on_landscapes_it_never_saw():
    """§81: -int ln(mu/lam) dx from the rate functions, against the exact MFPT."""
    for rails in ((0.4, 1.0, 2.2), (0.3, 1.0, 1.7), (0.05, 1.0, 2.5)):
        Ai, err = A_integral(*rails)
        assert err < 1e-9
        assert Ai == pytest.approx(A_mfpt(*rails), rel=2e-3), rails


def test_A_does_not_transfer():
    As = [A_integral(*r)[0] for r in ((0.1, 1.0, 1.9), (0.3, 1.0, 1.7), (0.05, 1.0, 2.5))]
    assert max(As) / min(As) > 5


def test_A_over_eta_is_not_constant_so_80_is_not_a_relabelling():
    """§81: the two coefficients respond differently to landscape shape."""
    rr = []
    for r1, r2, r3 in ((0.1, 1.0, 1.9), (0.4, 1.0, 2.2), (0.3, 1.0, 1.7), (0.05, 1.0, 2.5)):
        A = A_integral(r1, r2, r3)[0]
        eta = ((r3 - r1) / 2) ** 2 / (2 * schlogl_V(r1, r2, r3))
        rr.append(A / eta)
    assert max(rr) / min(rr) > 1.5


def test_underflowed_quasipotential_cells_are_excluded():
    """§81.1: at gamma=0.20 pi(saddle) hits the 1e-300 floor and 'A' becomes 691/Omega."""
    assert len(A_am(0.20)) == 0                       # excluded, not fitted
    ser = A_am(0.35)
    assert len(ser) >= 2
    vals = [v for _, v in ser]
    assert all(np.diff(vals) < 0)                     # converging downward, not 1/Omega
    prods = [om * v for om, v in ser]
    assert max(prods) / min(prods) > 1.2              # barrier GROWS with Omega, unlike the floor
