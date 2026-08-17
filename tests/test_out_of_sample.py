"""FINDINGS §79: the ODE+LNA prediction on systems the derivation never saw."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.chemical_channel_noise import am_rail_width, rail_width
from experiments.derive_eta import am_V, schlogl_V
from experiments.out_of_sample import (quartic_consts, quartic_roots, quartic_V,
                                       quartic_width)


def _dev(pred_V, exact_fn, oms):
    return [abs(exact_fn(o) / np.sqrt(pred_V / o) - 1) for o in oms]


def test_quartic_autocatalysis_is_predicted_from_the_ode_alone():
    """§79: a DIFFERENT reaction order (3X<->4X), never used in the derivation."""
    c = quartic_consts(m=0.8)
    r = quartic_roots(c)
    assert len(r) >= 3
    V = quartic_V(c, r)
    dev = _dev(V, lambda o: quartic_width(o, c, r), (1600, 6400, 25600))
    assert all(np.diff(dev) < 0), dev
    assert dev[-1] < 0.01


def test_asymmetric_rails_are_predicted_from_the_ode_alone():
    V = schlogl_V(0.4, 1.0, 2.2)
    dev = _dev(V, lambda o: rail_width(o, 0.4, 1.0, 2.2)["sd_exact"], (1600, 6400, 25600))
    assert all(np.diff(dev) < 0), dev
    assert dev[-1] < 0.01


def test_the_prediction_fails_where_it_was_predicted_to_fail():
    """§79 P2: near gamma_c the rail is shallow and the LNA must break -- chosen in advance."""
    V, _ = am_V(0.45)
    dev = _dev(V, lambda o: am_rail_width(0.45, o)[0], (60, 120, 240, 400))
    assert not all(np.diff(dev) < 0)      # does NOT converge
    assert dev[-1] > 0.10                 # and is badly wrong

    Vg, _ = am_V(0.10)                    # while a deep rail is fine
    devg = _dev(Vg, lambda o: am_rail_width(0.10, o)[0], (60, 120, 240))
    assert devg[-1] < 0.01
    assert dev[-1] > 10 * devg[-1]


def test_depth_is_exponentially_sensitive_so_79_1s_caveat_is_real():
    """§79.1: a sub-percent sigma error is a large-factor D_max error."""
    from scipy.stats import norm
    from experiments.depth_is_error import c_star

    V, ds = am_V(0.35)
    om = 400
    sd_exact, _ = am_rail_width(0.35, om)
    ln_pred = np.log(c_star()) - norm.logcdf(-ds / np.sqrt(V / om))
    ln_exact = np.log(c_star()) - norm.logcdf(-ds / sd_exact)
    sigma_err = abs(sd_exact / np.sqrt(V / om) - 1)
    assert sigma_err < 0.02                       # sigma is good
    assert np.exp(ln_pred - ln_exact) > 5         # D_max is not
