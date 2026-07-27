"""Tests for the tilted AM network.

The load-bearing ones are `test_drift_matches_engine` (the hand-derived ODE in
am_asymmetric is the basis of every fixed point, so it must agree with the
generic mass-action engine reaction-for-reaction) and
`test_affinity_is_beta_independent` (the whole reason for scaling each reverse
by its own forward; if it failed, beta would be confounded with the drive and
every beta sweep would be a gamma sweep in disguise).
"""
from __future__ import annotations

import numpy as np
import pytest

from crnl.networks.am_asymmetric import (
    am_asymmetric,
    basin_boundary,
    beta_critical,
    drift,
    interior_fixed_points,
)
from crnl.networks.am_reversible import (
    GAMMA_C,
    am_reversible,
    cycle_affinity,
    delta_star,
    reverse_pairing,
)

GAMMAS = [0.05, 0.2, 0.35, 0.45]


@pytest.mark.parametrize("gamma", GAMMAS)
@pytest.mark.parametrize("beta", [-0.4, -0.1, 0.0, 0.1, 0.4])
def test_drift_matches_engine(gamma, beta):
    """The hand-written drift must equal S @ fluxes from the generic engine."""
    net = am_asymmetric(gamma, beta)
    S = net.stoichiometry_matrix()
    rng = np.random.default_rng(0)
    for _ in range(20):
        x, y = rng.uniform(0.02, 0.9, 2)
        if x + y > 0.98:
            continue
        state = np.array([x, y, 1.0 - x - y])
        engine = S @ net.fluxes(state)
        assert np.allclose(drift([x, y], gamma, beta), engine[:2], atol=1e-12)


@pytest.mark.parametrize("gamma", GAMMAS)
@pytest.mark.parametrize("beta", [0.0, 0.15, 0.5, 0.9])
def test_affinity_is_beta_independent(gamma, beta):
    """A = -3 ln gamma for every tilt: the tilt costs no thermodynamic force."""
    net = am_asymmetric(gamma, beta)
    a = cycle_affinity(net, reverse_pairing(net))
    assert a == pytest.approx(-3.0 * np.log(gamma), rel=1e-12)


@pytest.mark.parametrize("gamma", GAMMAS)
def test_reduces_to_symmetric(gamma):
    """beta = 0 must reproduce am_reversible, rate for rate and root for root."""
    a, b = am_asymmetric(gamma, 0.0), am_reversible(gamma)
    assert len(a.reactions) == len(b.reactions)
    for ra, rb in zip(a.reactions, b.reactions):
        assert ra.reactants == rb.reactants
        assert ra.products == rb.products
        assert ra.k == pytest.approx(rb.k, rel=1e-15)
    pts = interior_fixed_points(gamma, 0.0)
    assert len(pts) == 3
    assert pts[-1]["x"] - pts[-1]["y"] == pytest.approx(delta_star(gamma), rel=1e-7)
    mid = pts[1]
    assert mid["kind"] == "saddle"
    assert (mid["x"], mid["y"]) == pytest.approx((1 / 3, 1 / 3), abs=1e-9)


@pytest.mark.parametrize("gamma", [0.05, 0.2, 0.35])
def test_tilt_moves_the_boundary_not_the_attractors(gamma):
    """The bias lives in the saddle. Positive beta must favour X (boundary < 0)."""
    assert basin_boundary(gamma, 0.0) == pytest.approx(0.0, abs=1e-9)
    prev = 0.0
    for beta in [0.02, 0.05, 0.1]:
        bnd = basin_boundary(gamma, beta)
        assert bnd < prev, f"boundary must move toward Y as beta grows: {bnd} !< {prev}"
        prev = bnd


@pytest.mark.parametrize("gamma", [0.05, 0.2, 0.35])
def test_tilt_is_antisymmetric_in_beta(gamma):
    """Relabelling X<->Y is exactly beta -> -beta; nothing may break that."""
    for beta in [0.05, 0.2]:
        pos = interior_fixed_points(gamma, beta)
        neg = interior_fixed_points(gamma, -beta)
        assert len(pos) == len(neg)
        # mirroring (x,y) -> (y,x) maps one set onto the other
        for p in pos:
            assert any(abs(p["x"] - q["y"]) < 1e-7 and abs(p["y"] - q["x"]) < 1e-7
                       for q in neg)
        assert basin_boundary(gamma, -beta) == pytest.approx(
            -basin_boundary(gamma, beta), abs=1e-7)


@pytest.mark.parametrize("gamma", [0.2, 0.35, 0.45])
def test_beta_critical_brackets_the_fold(gamma):
    """3 interior fixed points just below beta_c, 1 just above."""
    bc = beta_critical(gamma)
    assert 0.0 < bc < 1.0
    assert len(interior_fixed_points(gamma, bc - 1e-3)) == 3
    assert len(interior_fixed_points(gamma, bc + 1e-3)) == 1
    assert np.isnan(basin_boundary(gamma, bc + 1e-3))


def test_beta_critical_falls_as_the_landscape_weakens():
    """A shallower landscape tolerates less tilt. Monotone in gamma."""
    bcs = [beta_critical(g) for g in [0.05, 0.2, 0.3, 0.35, 0.4, 0.45]]
    assert all(np.diff(bcs) < 0), bcs


def test_beta_critical_undefined_above_gamma_c():
    """No fold to find when there was never a landscape: nan, not 0.0."""
    for gamma in [GAMMA_C, 0.6, 0.9]:
        assert np.isnan(beta_critical(gamma))


def test_fixed_point_count_is_grid_independent():
    """The multi-start grid must not be what sets the number of roots."""
    import crnl.networks.am_asymmetric as mod

    for gamma, beta in [(0.2, 0.0), (0.2, 0.3), (0.35, 0.1), (0.05, 0.5)]:
        coarse = len(interior_fixed_points(gamma, beta))
        # rerun the same root find on a seed grid ~3x denser
        seeds = [(gx, gy)
                 for gx in np.linspace(0.01, 0.97, 40)
                 for gy in np.linspace(0.01, 0.97, 40)
                 if gx + gy < 0.995]
        from scipy.optimize import fsolve
        found = []
        for s in seeds:
            sol, _, ier, _ = fsolve(mod.drift, s, args=(gamma, beta),
                                    fprime=mod._jac, full_output=True)
            if ier != 1:
                continue
            x, y = float(sol[0]), float(sol[1])
            if min(x, y, 1 - x - y) <= 1e-7:
                continue
            if np.max(np.abs(mod.drift(sol, gamma, beta))) > 1e-9:
                continue
            if any(np.hypot(x - f[0], y - f[1]) < 1e-6 for f in found):
                continue
            found.append((x, y))
        assert len(found) == coarse, (gamma, beta, len(found), coarse)


@pytest.mark.parametrize("beta", [-1.0, 1.0, 1.5, np.nan, np.inf])
def test_beta_domain_is_enforced(beta):
    with pytest.raises(ValueError):
        am_asymmetric(0.2, beta)


def test_equilibrium_has_no_landscape_at_any_tilt():
    """gamma = 1 is detailed balance for every beta; one fixed point, stable."""
    for beta in [0.0, 0.3, 0.7]:
        pts = interior_fixed_points(1.0, beta)
        assert len(pts) == 1
        assert pts[0]["stable"]
