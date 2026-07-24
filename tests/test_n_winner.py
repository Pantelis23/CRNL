from __future__ import annotations

import numpy as np
import pytest

from crnl import approximate_majority
from crnl.networks.n_winner import n_winner
from crnl import classify


def test_n_winner_2_reduces_to_am():
    # n_winner(2) has species [X1,X2,B] in the same order as AM's [X,Y,B],
    # so its RHS must equal AM's on every state.
    nw = n_winner(2)
    am = approximate_majority()
    assert nw.species == ["X1", "X2", "B"]
    rng = np.random.default_rng(0)
    for _ in range(20):
        x, y = rng.random(2) * 0.4
        state = [x, y, 1 - x - y]
        assert np.allclose(nw.rhs(state), am.rhs(state), atol=1e-12)


def test_n_winner_reaction_counts():
    # disagreement = C(n,2) pairs, recruitment = n, total = C(n,2)+n
    for n in (2, 3, 5):
        nw = n_winner(n)
        assert nw.n_species == n + 1
        assert nw.n_reactions == n * (n - 1) // 2 + n


def test_n_winner_all_heterobimolecular_no_homodimer():
    # spec FRAGILE-1: every reaction has two DISTINCT reactant species, coeff 1.
    nw = n_winner(4)
    for r in nw.reactions:
        assert sum(r.reactants.values()) == 2
        assert len(r.reactants) == 2  # two distinct species, not one with coeff 2
        assert all(c == 1 for c in r.reactants.values())


def test_n_winner_conserves_count():
    nw = n_winner(5)
    S = nw.stoichiometry_matrix()
    assert np.allclose(S.sum(axis=0), 0.0)  # every reaction is 2->2


def test_n_winner_rails_are_stable():
    # each all-Xi corner classifies stable within the stoichiometric subspace
    n = 3
    nw = n_winner(n)
    for i in range(n):
        x = np.zeros(n + 1)
        x[i] = 1.0
        fp = classify.classify_fixed_point(nw, x)
        assert fp.kind == "stable", f"rail X{i+1} was {fp.kind}"


def test_n_winner_symmetric_point_is_unstable():
    # Interior symmetric fixed point: x_i = 1/(2n-1), b = (n-1)/(2n-1).
    # Derivation: at x_i=x for all i, dXi/dt = -(n-1)x^2 + b*x = 0 => b=(n-1)x,
    # and n*x + b = 1 => x = 1/(2n-1). At n=2 this is (1/3,1/3,1/3) (AM saddle).
    for n in (2, 3, 4):
        nw = n_winner(n)
        x = 1.0 / (2 * n - 1)
        b = (n - 1) / (2 * n - 1)
        state = np.array([x] * n + [b])
        # confirm it really is a fixed point
        assert np.allclose(nw.rhs(state), 0.0, atol=1e-9)
        fp = classify.classify_fixed_point(nw, state)
        assert fp.kind in ("unstable", "saddle"), f"n={n} symmetric was {fp.kind}"
