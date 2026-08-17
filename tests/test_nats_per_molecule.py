"""FINDINGS §77: nats of reliability per molecule -- the last free number."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.nats_per_molecule import am_eta, ln_inv_eps, schlogl_eta


def _eta(fn, oms):
    ys = [fn(o) for o in oms]
    assert all(y is not None for y in ys)
    return [(ys[i + 1][0] - ys[i][0]) / (oms[i + 1] - oms[i]) for i in range(len(oms) - 1)]


def test_eta_settles_so_ln_inv_eps_is_linear_in_omega():
    """§77 P1: gated on convergence, not a tolerance (rule 20)."""
    sl = _eta(lambda o: schlogl_eta(1.0, o), [1600, 6400, 25600])
    assert abs(sl[-1] - sl[-2]) / abs(sl[-1]) < 0.05
    sl_am = _eta(lambda o: am_eta(0.05, o), [120, 240, 400])
    assert abs(sl_am[-1] - sl_am[-2]) / abs(sl_am[-1]) < 0.01


def test_schlogl_eta_scales_exactly_with_lambda_as_75s_collapse_requires():
    """§77 P3: lambda = 4 must give 4x lambda = 1, or §75 and §77 disagree."""
    a = _eta(lambda o: schlogl_eta(1.0, o), [6400, 25600])[0]
    b = _eta(lambda o: schlogl_eta(4.0, o), [6400, 25600])[0]
    assert b / a == pytest.approx(4.0, rel=0.01)


def test_eta_does_not_transfer_and_moves_with_the_landscape():
    """§77.1: a factor of ~600 across elements, and it moves inside each substrate too."""
    am_lo = _eta(lambda o: am_eta(0.20, o), [240, 400])[0]
    am_hi = _eta(lambda o: am_eta(0.05, o), [240, 400])[0]
    sch = _eta(lambda o: schlogl_eta(1.0, o), [6400, 25600])[0]
    assert am_hi / am_lo > 4          # landscape dependence within AM
    assert am_hi / sch > 100          # and no transfer between substrates


def test_am_eta_reports_rather_than_swallows_an_untrustworthy_solve():
    """§77.2: the engine guards the stationary solve; a bare `except` hid the exclusion."""
    assert am_eta(0.30, 240) is None          # the cell the first version dropped silently
    assert am_eta(0.20, 120) is not None


def test_ln_inv_eps_is_computed_in_logs():
    """A tail of e^-4000 must not underflow to inf/nan."""
    v = ln_inv_eps(1.0, 0.01)
    assert np.isfinite(v) and v > 4000
