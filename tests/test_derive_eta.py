"""FINDINGS §78: eta = Delta^2/(2V) from the LNA, against §77's measured values."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.chemical_channel_noise import am_rail_width, rail_width
from experiments.derive_eta import AM_ETA, SCH_ETA, am_V, schlogl_V


def test_lna_variance_is_the_exact_rail_width_in_the_large_omega_limit():
    """§78 P1: the content of the section. Convergence, not a tolerance (rule 20)."""
    V, _ = am_V(0.20)
    dev = []
    for om in (60, 120, 240, 400):
        sd, _ = am_rail_width(0.20, om)
        dev.append(abs(sd / np.sqrt(V / om) - 1))
    assert all(np.diff(dev) < 0), dev
    assert dev[-1] < 2e-3

    Vs = schlogl_V(0.1, 1.0, 1.9)
    devs = [abs(rail_width(o, 0.1, 1.0, 1.9)["sd_exact"] / np.sqrt(Vs / o) - 1)
            for o in (1600, 6400, 25600)]
    assert all(np.diff(devs) < 0), devs
    assert devs[-1] < 1e-3


def test_eta_predicted_matches_section_77s_measurements():
    """§78 P2: absolute, against stored numbers, nothing fitted."""
    for g, meas in AM_ETA.items():
        V, ds = am_V(g)
        assert ds ** 2 / (2 * V) == pytest.approx(meas, rel=0.01), g
    for spread, meas in SCH_ETA.items():
        V = schlogl_V(1.0 - spread, 1.0, 1.0 + spread)
        assert spread ** 2 / (2 * V) == pytest.approx(meas, rel=0.01), spread


def test_the_lna_error_is_worst_where_the_rail_is_shallowest():
    """§78 P4: the residual must grow toward gamma_c, not away from it."""
    err = {}
    for g, meas in AM_ETA.items():
        V, ds = am_V(g)
        err[g] = abs(ds ** 2 / (2 * V) / meas - 1)
    assert err[0.30] > err[0.20] > err[0.05]


def test_eta_scales_exactly_with_lambda():
    """§78 P5: §75's collapse demands eta ~ lambda, and the derivation must reproduce it."""
    a = 0.9 ** 2 / (2 * schlogl_V(0.1, 1.0, 1.9))
    b = (4 * 0.9) ** 2 / (2 * schlogl_V(0.4, 4.0, 7.6))
    assert b / a == pytest.approx(4.0, rel=1e-6)
