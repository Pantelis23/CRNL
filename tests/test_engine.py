"""Verification suite -- the physics checks are the tests (design.md §6 build order).

Run with:  pytest -q     (or  python -m pytest tests/ -q)
"""

from __future__ import annotations

import numpy as np
import pytest

from crnl import approximate_majority, integrate, gillespie, seed_for
from crnl import classify
from crnl.reactions import Reaction, ReactionNetwork


# --------------------------------------------------------------------------- #
# reactions.py -- stoichiometry, fluxes, units convention (§2.1, §3.2, §3.3)   #
# --------------------------------------------------------------------------- #

def test_am_stoichiometry_is_conservative():
    net = approximate_majority()
    S = net.stoichiometry_matrix()
    # every AM reaction is 2->2, so each column sums to zero (count conserved)
    assert np.allclose(S.sum(axis=0), 0.0)


def test_am_rhs_matches_reduced_odes():
    # design.md §2.2: dx/dt = x(1 - x - 2y), dy/dt = y(1 - 2x - y)
    net = approximate_majority()
    rng = np.random.default_rng(0)
    for _ in range(20):
        x, y = rng.random(2) * 0.4
        b = 1 - x - y
        r = net.rhs([x, y, b])
        assert r[0] == pytest.approx(x * (1 - x - 2 * y), abs=1e-9)
        assert r[1] == pytest.approx(y * (1 - 2 * x - y), abs=1e-9)


def test_homodimer_units_convention():
    # design.md §3.2/§3.3: for 2A -> B with k, the stochastic constant is 2k/Omega
    # and the propensity is c * n(n-1)/2. The two facts together must reproduce
    # the deterministic drift d[A]/dt = -2k[A]^2.
    net = ReactionNetwork(
        species=["A", "B"],
        reactions=[Reaction({"A": 2}, {"B": 1}, k=1.0, name="2A->B")],
    )
    omega = 100.0
    c = net.stochastic_constants(omega)[0]
    assert c == pytest.approx(2 * 1.0 / omega)  # 2k/Omega, NOT k/Omega

    nA = 40
    a = net.propensities(np.array([nA, 0]), np.array([c]))[0]
    assert a == pytest.approx(c * nA * (nA - 1) / 2)  # c n(n-1)/2, NOT c n^2

    # mean stochastic drift of [A] matches the ODE -2k[A]^2 to O(1/n)
    A = nA / omega
    ode_drift = -2 * 1.0 * A ** 2
    ssa_mean_drift = -2 * a / omega  # each firing changes n_A by -2
    assert ssa_mean_drift == pytest.approx(ode_drift, rel=1 / nA + 1e-9)


def test_heterobimolecular_units_convention():
    net = ReactionNetwork(
        species=["A", "B", "C"],
        reactions=[Reaction({"A": 1, "B": 1}, {"C": 2}, k=1.0)],
    )
    omega = 50.0
    c = net.stochastic_constants(omega)[0]
    assert c == pytest.approx(1.0 / omega)  # k/Omega
    a = net.propensities(np.array([10, 7, 0]), np.array([c]))[0]
    assert a == pytest.approx(c * 10 * 7)


# --------------------------------------------------------------------------- #
# deterministic.py -- fixed points and eigenvalues (§2.3)                      #
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize(
    "x, expected_kind, expected_eigs",
    [
        ([1.0, 0.0, 0.0], "stable", [-1.0, -1.0]),
        ([0.0, 1.0, 0.0], "stable", [-1.0, -1.0]),
        ([1 / 3, 1 / 3, 1 / 3], "saddle", [-1.0, 1 / 3]),
        ([0.0, 0.0, 1.0], "unstable", [1.0, 1.0]),
    ],
)
def test_am_fixed_point_classification(x, expected_kind, expected_eigs):
    net = approximate_majority()
    fp = classify.classify_fixed_point(net, x)
    assert fp.kind == expected_kind
    got = sorted(np.round(fp.eigenvalues.real, 4))
    exp = sorted(round(e, 4) for e in expected_eigs)
    assert np.allclose(got, exp, atol=1e-3)


def test_conservation_drift_is_tiny():
    net = approximate_majority()
    traj = integrate(net, [0.51, 0.49, 0.0], t_span=(0, 200))
    # total must stay 1 to numerical tolerance; drift is watched, not clamped
    assert traj.conserved_drift < 1e-6


def test_ode_glides_to_x_rail_from_slight_bias():
    # design.md §4: the deterministic curve is "the lie" -- 51/49 always -> X.
    net = approximate_majority()
    final = classify.find_stable_endpoint(net, [0.51, 0.49, 0.0])
    assert final[0] == pytest.approx(1.0, abs=1e-3)
    assert final[1] == pytest.approx(0.0, abs=1e-3)


# --------------------------------------------------------------------------- #
# stochastic.py -- SSA converges to the ODE as Omega grows (§6)               #
# --------------------------------------------------------------------------- #

def test_ssa_converges_to_ode_as_omega_grows():
    # design.md §6: the single best test that the §3.2 units convention is right.
    # A wrong stochastic constant leaves the ODE untouched and shows up ONLY as
    # the two engines failing to converge. Strongly-biased start so the mean SSA
    # outcome is unambiguous; error probability -> 0 as Omega grows.
    net = approximate_majority()
    x0 = np.array([0.60, 0.40, 0.0])
    ode_final = classify.find_stable_endpoint(net, x0)
    assert ode_final[0] == pytest.approx(1.0, abs=1e-3)  # ODE says X wins

    prev_err = 1.0
    for omega in (50, 200):
        n0 = np.round(x0 * omega).astype(int)
        n0[2] = int(omega) - n0[0] - n0[1]
        wrong = 0
        trials = 300
        for t in range(trials):
            res = gillespie(net, n0, omega, seed_for(omega, t))
            if classify.classify_am_outcome(res) == "Y":
                wrong += 1
        err = wrong / trials
        # error fraction must fall as Omega grows (restoration wall)
        assert err <= prev_err + 0.02
        prev_err = err
    assert prev_err < 0.05  # by Omega=200 a 60/40 start almost never errs


def test_ssa_absorbs_into_three_bins():
    # design.md §3.4: X-wins, Y-wins, and the all-blank bin the ODE denies.
    # All-blank is a genuine finite-Omega outcome but only non-negligible at the
    # low-Omega end (empirically ~0.5% at Omega=10, vanishing by Omega=20), so we
    # probe it where it actually lives.
    net = approximate_majority()
    omega = 10
    from collections import Counter

    bins = Counter()
    for t in range(3000):
        # start balanced so all three outcomes are reachable, incl. all-blank
        n0 = np.array([omega // 2, omega - omega // 2, 0])
        res = gillespie(net, n0, omega, seed_for(omega, t, base=7))
        assert res.absorbed
        bins[classify.classify_am_outcome(res)] += 1
    assert bins["X"] > 0 and bins["Y"] > 0  # decisions always reachable
    assert bins["B"] > 0  # the all-blank bin the deterministic repeller denies
    assert bins["undecided"] == 0  # every AM trajectory absorbs


def test_seed_is_replayable():
    net = approximate_majority()
    omega = 40
    n0 = np.array([21, 19, 0])
    a = gillespie(net, n0, omega, seed_for(omega, 123))
    b = gillespie(net, n0, omega, seed_for(omega, 123))
    assert np.array_equal(a.n_final, b.n_final)
    assert a.steps == b.steps
