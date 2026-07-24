"""Tests for the scaling/theory experiment helpers.

These cover the pure analysis functions (grid construction, power-law fitting,
collapse residual, crossings) rather than the Monte Carlo itself, which is
covered by the engine tests.
"""

from __future__ import annotations

import numpy as np
import pytest

from experiments.quasipotential import (
    KAPPA_THEORY,
    omega_grid_for,
    fit_power_law,
)
from experiments.freezeout_scaling import collapse_residual, fit_collapse, crossing
from experiments.radix_scaling import locate_band, omega_required


# --------------------------------------------------------------------------- #
# quasipotential                                                              #
# --------------------------------------------------------------------------- #

def test_omega_grid_scales_inversely_with_eps_squared():
    # c ~ eps^2, so the visible-wall window sits at Omega ~ 1/eps^2:
    # halving eps must push the grid ~4x higher.
    small = omega_grid_for(0.05)
    big = omega_grid_for(0.10)
    assert min(small) > 3 * min(big)
    assert max(small) > 3 * max(big)


def test_fit_power_law_recovers_known_exponent():
    eps = np.array([0.04, 0.06, 0.08, 0.10, 0.14, 0.20])
    c = 1.5 * eps ** 2                      # the exact theory curve
    p, kappa = fit_power_law(eps, c)
    assert p == pytest.approx(2.0, abs=1e-6)
    assert kappa == pytest.approx(KAPPA_THEORY, rel=1e-6)


def test_fit_power_law_ignores_bad_points():
    eps = np.array([0.05, 0.1, 0.2])
    c = np.array([1.5 * 0.05 ** 2, np.nan, 1.5 * 0.2 ** 2])
    p, kappa = fit_power_law(eps, c)
    assert p == pytest.approx(2.0, abs=1e-6)  # nan dropped, not propagated


# --------------------------------------------------------------------------- #
# freeze-out finite-size scaling                                              #
# --------------------------------------------------------------------------- #

def _synthetic_fss(Hc=0.05, a=0.4, omegas=(40, 160, 640)):
    """Data built to obey D = F((H-Hc)*Omega^a) exactly, F a logistic."""
    by = {}
    for w in omegas:
        rows = []
        for H in np.geomspace(0.02, 0.4, 18):
            u = (H - Hc) * w ** a
            rows.append({"hubble": float(H), "order": float(1.0 / (1.0 + np.exp(3 * u))),
                         "p_frozen": 0.0, "relic": 0.0})
        by[w] = rows
    return by


def test_collapse_residual_is_minimal_at_true_parameters():
    by = _synthetic_fss(Hc=0.05, a=0.4)
    true = collapse_residual(by, 0.05, 0.4)
    # any sizeable displacement in either parameter must collapse worse
    assert true < collapse_residual(by, 0.09, 0.4)
    assert true < collapse_residual(by, 0.05, 0.1)
    assert true < collapse_residual(by, 0.02, 0.7)


def test_fit_collapse_recovers_synthetic_parameters():
    by = _synthetic_fss(Hc=0.05, a=0.4)
    fit = fit_collapse(by)
    assert fit["Hc"] == pytest.approx(0.05, abs=0.015)
    assert fit["a"] == pytest.approx(0.4, abs=0.06)


def test_crossing_interpolates_and_reports_nan_when_uncrossed():
    rows = [{"hubble": 0.05, "order": 0.9},
            {"hubble": 0.10, "order": 0.6},
            {"hubble": 0.20, "order": 0.3}]
    h = crossing(rows, level=0.5)
    assert 0.10 < h < 0.20
    never = [{"hubble": 0.05, "order": 0.9}, {"hubble": 0.1, "order": 0.8}]
    assert np.isnan(crossing(never, level=0.5))


# --------------------------------------------------------------------------- #
# radix scaling                                                               #
# --------------------------------------------------------------------------- #

def _pt(omega, loss, trials=1000):
    return {"omega": omega, "trials": trials,
            "champion_wins": int(round(trials * (1 - loss)))}


def test_locate_band_picks_the_error_rich_window():
    rows = [_pt(20, 0.60), _pt(50, 0.30), _pt(100, 0.15),
            _pt(200, 0.05), _pt(400, 0.001)]
    lo, hi = locate_band(rows)
    assert lo == 50 and hi == 200     # drops the saturated and the vanished ends


def test_locate_band_falls_back_when_window_is_empty():
    # all points saturated: must still return a usable range, not crash
    rows = [_pt(20, 0.60), _pt(50, 0.55), _pt(100, 0.52)]
    lo, hi = locate_band(rows)
    assert lo <= hi


def test_omega_required_inverts_the_fitted_wall():
    # loss = exp(intercept - c*Omega); at Omega* loss == target
    fit = {"c": 0.02, "intercept": 0.0}
    om = omega_required(fit, target=0.05)
    assert np.exp(fit["intercept"] - fit["c"] * om) == pytest.approx(0.05)
    assert np.isnan(omega_required(None))
