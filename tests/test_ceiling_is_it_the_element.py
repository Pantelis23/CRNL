"""FINDINGS §73: the depth ceiling does not see the element."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.cascade_saturated import d_max_saturated
from experiments.ceiling_is_it_the_element import (eps_from, grid, pc_langevin,
                                                   pc_schlogl, pc_sigmoid, pc_step)

R1, R2, R3, OM = 0.1, 1.0, 1.9, 3600
DELTA = (R3 - R1) / 2.0


def _ratio(pc, f):
    x, _ = grid(R1, R3, 1.6, OM)
    d = d_max_saturated(*eps_from(pc, x, R1, R3, f * DELTA))
    return d / (np.exp(1.0 / (2 * f ** 2)) / 4.0)


def test_a_step_function_reproduces_the_ceiling():
    """The headline: no chemistry, no dynamics, same ratio -- so §72's reading is deflated."""
    x, _ = grid(R1, R3, 1.6, OM)
    for f in (0.45, 0.35, 0.28):
        chem = _ratio(pc_schlogl(x, R1, R2, R3, OM), f)
        step = _ratio(pc_step(x, R1, R2, R3, OM), f)
        assert abs(step - chem) / chem < 0.30, (f, chem, step)


def test_a_slope_matched_sigmoid_matches_the_chemistry_closely():
    x, _ = grid(R1, R3, 1.6, OM)
    pc = pc_schlogl(x, R1, R2, R3, OM)
    sad = int(np.argmin(np.abs(x - R2)))
    slope = float(-(pc[sad + 1] - pc[sad - 1]) / (x[sad + 1] - x[sad - 1]))
    for f in (0.45, 0.35):
        assert _ratio(pc_sigmoid(x, R1, R2, R3, OM, 4.0 * slope), f) == pytest.approx(
            _ratio(pc, f), rel=0.05)


def test_langevin_scale_density_is_exp_U_over_D_not_2U_over_D():
    """§73.2: the first version had the factor of 2 and a cumsum that produced 0/0.

    The bug pointed toward the conclusion just published, so it is pinned.
    """
    x, _ = grid(R1, R3, 1.6, OM)
    pc = pc_langevin(x, R1, R2, R3, OM)
    assert pc is not None and np.all(np.isfinite(pc))
    lo = int(np.argmin(np.abs(x - R1)))
    hi = int(np.argmin(np.abs(x - R3)))
    assert pc[lo] == pytest.approx(1.0, abs=1e-9)
    assert pc[hi] == pytest.approx(0.0, abs=1e-9)
    assert np.all(np.diff(pc[lo:hi + 1]) <= 1e-12)
    assert _ratio(pc, 0.35) == pytest.approx(_ratio(pc_schlogl(x, R1, R2, R3, OM), 0.35),
                                             rel=0.15)
