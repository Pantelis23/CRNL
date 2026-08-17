"""FINDINGS §80: escape beats readout in a chemically-coupled cascade."""
from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import norm

from experiments.cascade_schlogl import schlogl_consts
from experiments.chemical_channel_noise import rail_width
from experiments.derive_eta import schlogl_V
from experiments.two_exponentials import ln_mfpt

RAILS = (0.1, 1.0, 1.9)


def test_mfpt_is_positive_and_grows_exponentially_with_omega():
    """§80 P1: the first instrument returned a NEGATIVE time. Guard it."""
    lts = [ln_mfpt(o, *RAILS) for o in (400, 800, 1600, 3200)]
    assert all(np.isfinite(t) and t > 0 for t in lts), lts
    loc = [(lts[i + 1] - lts[i]) / o for i, o in enumerate((400, 800, 1600))]
    assert all(0.02 < v < 0.03 for v in loc), loc
    assert abs(loc[-1] - loc[-2]) / loc[-1] < 0.01     # converging


def test_escape_exponent_is_smaller_than_the_readout_exponent():
    """§80 P2: A < eta, so escape dominates and §75-§79 priced the subdominant mode."""
    D = (RAILS[2] - RAILS[0]) / 2
    eta = D ** 2 / (2 * schlogl_V(*RAILS))
    A = (ln_mfpt(3200, *RAILS) - ln_mfpt(1600, *RAILS)) / 1600
    assert A < eta
    assert eta / A == pytest.approx(2.29, rel=0.05)


def test_escape_dominates_by_hundreds_of_nats_at_realistic_omega():
    for om in (1600, 6400):
        sd = rail_width(om, *RAILS)["sd_exact"]
        ln_read = float(norm.logcdf(-((RAILS[2] - RAILS[0]) / 2) / sd))
        ln_esc = np.log(2.0) - ln_mfpt(om, *RAILS)
        assert ln_esc - ln_read > 50, (om, ln_esc - ln_read)


def test_the_escape_action_is_a_deterministic_side_integral():
    """§80.1: A = -int ln(mu/lam) dx from the RATE FUNCTIONS -- no chain, no CME."""
    k1a, k1r, k2b, k2r = schlogl_consts(*RAILS)
    val, err = quad(lambda x: np.log((k1r * x ** 3 + k2r * x) / (k1a * x ** 2 + k2b)),
                    RAILS[1], RAILS[2], limit=200)
    A_exact = (ln_mfpt(3200, *RAILS) - ln_mfpt(1600, *RAILS)) / 1600
    assert -val == pytest.approx(A_exact, rel=1e-3)


def test_an_external_channel_still_binds_on_readout_so_71_72_stand():
    """§80.1: the two regimes are different cascades; §71/§72 used the wide external channel."""
    D = (RAILS[2] - RAILS[0]) / 2
    ln_read_ext = float(norm.logcdf(-D / (0.35 * D)))
    ln_esc = np.log(2.0) - ln_mfpt(1600, *RAILS)
    assert ln_read_ext > ln_esc + 20
