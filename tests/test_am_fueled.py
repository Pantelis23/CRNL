"""Tests for the depleting drive (FINDINGS 20).

`test_reduces_to_am_reversible_at_fixed_fuel` is the anchor, in the style of
FINDINGS 19's 0/300 check: hold the tank levels fixed and the (X, Y, B) drift must
be `am_reversible(gamma_inf * w/f)` with time rescaled by f, exactly. Everything
else about the fueled network is a statement about how that ratio MOVES, so if the
frozen-tank case drifts, nothing downstream means anything.
"""
from __future__ import annotations

import numpy as np
import pytest

from crnl.expanding import common_order
from crnl.networks.am_fueled import (
    GAMMA_C, am_fueled, death_waste_fraction, gamma_effective, initial_counts,
)
from crnl.networks.am_reversible import am_reversible, reverse_pairing
from crnl.vectorized import compile_network


@pytest.mark.parametrize("f,w", [(0.9, 0.1), (0.8, 0.2), (0.7, 0.3), (0.6, 0.4)])
def test_reduces_to_am_reversible_at_fixed_fuel(f, w):
    """The anchor: frozen tank == am_reversible(gamma_eff), time-rescaled by f."""
    g_eff = w / f
    nf, nr = am_fueled(1.0), am_reversible(min(g_eff, 0.999))
    Sf, Sr = nf.stoichiometry_matrix(), nr.stoichiometry_matrix()
    rng = np.random.default_rng(0)
    for _ in range(120):
        x, y = rng.uniform(0.05, 0.6, 2)
        if x + y > 0.95:
            continue
        b = 1.0 - x - y
        dfl = (Sf @ nf.fluxes(np.array([x, y, b, f, w])))[:3]
        drv = (Sr @ nr.fluxes(np.array([x, y, b]))) * f
        assert np.allclose(dfl, drv, atol=1e-13), (f, w, x, y, dfl, drv)


def test_death_point_is_parameter_free():
    """gamma_eff = 1/2 at w/Phi = 1/(1+2*gamma_inf); 1/3 at gamma_inf = 1."""
    assert death_waste_fraction(1.0) == pytest.approx(1 / 3)
    for gi in (0.5, 1.0, 2.0, 4.0):
        wf = death_waste_fraction(gi)
        # at that waste fraction the effective drive is exactly gamma_c
        assert gamma_effective(1 - wf, wf, gi) == pytest.approx(GAMMA_C, rel=1e-12)


def test_gamma_effective_rises_as_the_tank_empties():
    """The mirror of FINDINGS 19's cooling, which drove gamma DOWN."""
    gs = [float(gamma_effective(100 - w, w, 1.0)) for w in (1, 10, 25, 33, 50)]
    assert gs == sorted(gs)
    assert gs[0] < GAMMA_C < gs[-1]
    # an empty tank is inf, not 1 -- a different statement, not to be clipped
    assert np.isinf(float(gamma_effective(0, 100, 1.0)))


def test_every_reaction_conserves_both_totals():
    net = am_fueled(1.0)
    S = net.stoichiometry_matrix()
    for j in range(net.n_reactions):
        assert S[:3, j].sum() == 0, net.reactions[j].name    # X + Y + B
        assert S[3:, j].sum() == 0, net.reactions[j].name    # F + W


def test_fuel_is_an_independent_coordinate():
    """A full cycle returns (X,Y,B) to the start while burning three fuel.

    This is why the fixed-gamma model is a projection rather than an equivalent
    description: n_F cannot be recovered from (X, Y, B).
    """
    net = am_fueled(1.0)
    S = net.stoichiometry_matrix()
    cycle = S[:, 0] + S[:, 1] + S[:, 2]        # f1 + f2 + f3
    assert list(cycle[:3]) == [0, 0, 0]
    assert cycle[3] == -3 and cycle[4] == 3


def test_network_is_uniform_order_three():
    """It leaves the 2->2 class but stays uniform, so the FINDINGS 19 machinery
    still accepts it -- with lambda = 2H rather than H."""
    assert common_order(compile_network(am_fueled(1.0), 100.0)) == 3


def test_reverse_pairing_still_resolves():
    """The forward/reverse pairs must be derivable, not assumed."""
    net = am_fueled(1.0)
    pair = reverse_pairing(net)
    assert list(pair) == [3, 4, 5, 0, 1, 2]


@pytest.mark.parametrize("bad", [0.0, -1.0, np.nan, np.inf])
def test_rejects_bad_gamma_inf(bad):
    with pytest.raises(ValueError):
        am_fueled(bad)


def test_initial_counts_are_valid_and_committed():
    n0 = initial_counts(120, 1200, waste0=277, gamma_inf=1.0)
    assert n0.sum() == 120 + 1200
    assert n0[:3].sum() == 120 and n0[3] + n0[4] == 1200
    assert (n0 >= 0).all()
    assert n0[0] > n0[1]                       # committed to X
