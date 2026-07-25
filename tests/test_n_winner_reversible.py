"""Reversible n-winner AM and the affinity floor (FINDINGS 13)."""

import numpy as np
import pytest

from crnl.networks.n_winner_reversible import (
    affinity_critical, breaking_mode, gamma_critical, lambda_breaking,
    n_winner_reversible, symmetric_state,
)
from crnl.deterministic import jacobian


# -- the n=2 case must reproduce FINDINGS 9.1 exactly ----------------------

@pytest.mark.parametrize("gamma", [0.0, 0.1, 0.25, 0.4, 0.5, 0.6, 0.9])
def test_n2_reproduces_the_known_closed_form(gamma):
    """lambda_antisym = (1-2 gamma)/3 from am_reversible. If this drifts, the
    generalisation is not a generalisation."""
    assert lambda_breaking(2, gamma) == pytest.approx((1 - 2 * gamma) / 3, abs=1e-12)


def test_n2_recovers_gamma_c_one_half_and_three_ln_two():
    assert gamma_critical(2) == pytest.approx(0.5, abs=1e-12)
    assert affinity_critical(2) == pytest.approx(3 * np.log(2), abs=1e-11)


def test_n2_symmetric_point_is_gamma_independent_at_one_third():
    """The degeneracy that is SPECIAL to n=2 and does not survive n>=3."""
    for gamma in (0.02, 0.3, 0.8):
        x, b = symmetric_state(2, gamma)
        assert x == pytest.approx(1 / 3, abs=1e-12)
        assert b == pytest.approx(1 / 3, abs=1e-12)


# -- structure -------------------------------------------------------------

def test_network_shape_and_reversibility():
    n = 4
    net = n_winner_reversible(n, 0.3)
    assert net.n_species == n + 1
    assert net.n_reactions == 2 * (n * (n - 1) // 2 + n)   # forwards + reverses
    fwd = net.n_reactions // 2
    for j in range(fwd):
        f, r = net.reactions[j], net.reactions[j + fwd]
        assert f.reactants == r.products and f.products == r.reactants
        assert r.k == pytest.approx(0.3 * f.k)


@pytest.mark.parametrize("n", [2, 3, 4, 6])
def test_cycle_dimension_is_n_choose_2(n):
    """Counting each reversible PAIR as one edge -- the counting under which
    FINDINGS 9.1 called AM's cycle space one-dimensional. C(n,2)+n edges against
    rank(S)=n leaves exactly C(n,2): 1, 3, 6, 15.

    Mixing this with the count that treats forward and reverse as separate edges
    (4, 9, 16) is how an earlier draft of this module claimed n=3 had a
    9-dimensional cycle space while comparing it to 9.1's 1.
    """
    net = n_winner_reversible(n, 0.3)
    forward = net.n_reactions // 2
    S_forward = net.stoichiometry_matrix()[:, :forward]
    rank = np.linalg.matrix_rank(S_forward)
    assert rank == n
    assert forward - rank == n * (n - 1) // 2


def test_symmetric_state_is_a_fixed_point():
    for n in (2, 3, 5, 9):
        for gamma in (0.01, 0.2, 0.7):
            x, b = symmetric_state(n, gamma)
            net = n_winner_reversible(n, gamma)
            state = np.concatenate([np.full(n, x), [b]])
            assert np.abs(net.rhs(state)).max() < 1e-12
            assert x * n + b == pytest.approx(1.0, abs=1e-12)


def test_symmetric_state_moves_with_gamma_for_n_at_least_three():
    """Unlike n=2. Measured at n=3: x = 0.2050 / 0.2269 / 0.2431."""
    xs = [symmetric_state(3, g)[0] for g in (0.02, 0.2, 0.6)]
    assert xs == sorted(xs)
    assert xs[-1] - xs[0] > 0.03
    assert not any(x == pytest.approx(1 / 4, abs=1e-3) for x in xs)  # not 1/(n+1)


def test_breaking_mode_is_a_true_eigenvector():
    """Otherwise the Rayleigh quotient would be a bound, not the eigenvalue."""
    for n in (3, 8, 16):
        gamma = 0.05
        x, b = symmetric_state(n, gamma)
        J = jacobian(n_winner_reversible(n, gamma),
                     np.concatenate([np.full(n, x), [b]]))
        v = breaking_mode(n)
        Jv = J @ v
        lam = v @ Jv
        assert np.linalg.norm(Jv - lam * v) / np.linalg.norm(Jv) < 1e-10


# -- the measured law ------------------------------------------------------

def test_gamma_c_falls_much_faster_than_one_over_n():
    """Kills the pre-run prediction gamma_c = 1/n (which would give A_c=3 ln n).
    Measured at n=32: 3.58e-5 against 1/n = 0.031, ~870x below."""
    assert gamma_critical(32) < 0.031 / 100


def test_affinity_floor_approaches_nine_ln_n():
    """A_c/ln n climbs 3.0 -> ~9. Measured 8.9894 at n=256."""
    ratios = [affinity_critical(n) / np.log(n) for n in (2, 8, 32, 128)]
    assert ratios == sorted(ratios)
    assert ratios[0] == pytest.approx(3.0, abs=1e-6)     # exact at n=2
    assert 8.9 < ratios[-1] < 9.0


def test_the_landscape_really_dies_above_gamma_c():
    """The eigenvalue crossing is a bifurcation, not just a sign change:
    below gamma_c a broken state survives, above it symmetry is restored."""
    from scipy.integrate import solve_ivp

    n = 4
    gc = gamma_critical(n)
    for factor, expect_broken in ((0.7, True), (1.4, False)):
        gamma = gc * factor
        x, b = symmetric_state(n, gamma)
        net = n_winner_reversible(n, gamma)
        y0 = np.concatenate([np.full(n, x), [b]])
        y0[0] += 1e-6
        y0[1] -= 1e-6
        sol = solve_ivp(lambda t, y: net.rhs(y), (0, 3000.0), y0,
                        method="LSODA", rtol=1e-11, atol=1e-14)
        spread = sol.y[:n, -1].max() - sol.y[:n, -1].min()
        # bool(): numpy comparisons return np.bool_, and np.bool_(True) is not
        # the Python singleton, so `is expect_broken` fails on a passing case
        assert bool(spread > 0.5) is expect_broken, (gamma, spread)


def test_rejects_degenerate_inputs():
    with pytest.raises(ValueError, match="n >= 2"):
        n_winner_reversible(1, 0.3)
    with pytest.raises(ValueError, match="gamma"):
        n_winner_reversible(3, -0.1)
