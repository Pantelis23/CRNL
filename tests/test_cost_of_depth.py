"""FINDINGS §74: what rail separation costs, and the closed form that made it computable."""
from __future__ import annotations

import numpy as np
import pytest

from crnl.networks.am_reversible import delta_star
from experiments.cascade_saturated import (d_max_closed, d_max_saturated,
                                           mutual_info_at, mutual_info_depth)
from experiments.cost_of_depth import depth_at, schlogl_affinity


def test_closed_form_matches_the_iterative_chain_72_published():
    """§74.3: the O(log) closed form must not move §72's numbers."""
    for e_hi, e_lo in ((0.02, 0.03), (0.001, 0.0015), (1e-4, 3e-4)):
        for D in (5, 37, 210):
            assert mutual_info_at(e_hi, e_lo, D) == pytest.approx(
                mutual_info_depth(e_hi, e_lo, D)[-1], abs=1e-12)
        assert np.ceil(d_max_closed(e_hi, e_lo)) == d_max_saturated(e_hi, e_lo)


def test_am_rail_separation_is_bounded_by_conservation():
    """§74.1: delta* rises with drive but is capped at 1, so depth saturates."""
    ds = [float(delta_star(g)) for g in (0.45, 0.20, 0.05, 0.002, 1e-4)]
    assert all(np.diff(ds) > 0)
    assert ds[-1] < 1.0
    assert ds[-1] > 0.999
    # and the depth ceiling it implies is finite
    assert depth_at(1.0, 0.15, -1.0, 1.0, 0.0) < 1e11


def test_schlogl_affinity_is_exactly_scale_invariant():
    """§74.2: r -> lambda r leaves ln[e1 e2/e3] untouched, so Delta is free in affinity."""
    base = (0.5, 1.0, 1.5)
    a0 = schlogl_affinity(*base)
    for lam in (0.25, 4.0, 64.0, 1000.0):
        assert schlogl_affinity(*(lam * v for v in base)) == pytest.approx(a0, abs=1e-13)


def test_open_element_passes_the_conservative_ceiling_at_matched_affinity():
    """§74.2: the headline -- same affinity, unbounded depth, paid in material."""
    sig = 0.15
    am_ceiling = depth_at(1.0, sig, -1.0, 1.0, 0.0)
    deep = depth_at(2.0, sig, 2.0, 6.0, 4.0)          # lambda = 4, same affinity
    assert deep is None or deep > 10 * am_ceiling     # None means the 1e18 cap was passed


def test_depth_grid_is_not_clamped_at_zero():
    """§74.3: AM's coordinate is signed; a positivity clamp deleted its whole low rail."""
    d = depth_at(0.4, 0.15, -0.4, 0.4, 0.0)
    assert d is not None and d > 1
