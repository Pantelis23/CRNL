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


from crnl.reactions import Reaction, ReactionNetwork
from crnl.networks.am import approximate_majority as _am  # noqa: F401


def _reference_props(net, n, omega):
    return net.propensities(np.asarray(n), net.stochastic_constants(omega))


@pytest.mark.parametrize("net", [
    ReactionNetwork(["A"], [Reaction({"A": 1}, {}, 1.0)]),                  # unimolecular
    ReactionNetwork(["A", "B", "C"], [Reaction({"A": 1, "B": 1}, {"C": 2}, 1.0)]),  # hetero
    ReactionNetwork(["A", "B"], [Reaction({"A": 2}, {"B": 1}, 1.0)]),       # homodimer
    approximate_majority(),
    n_winner(3),
])
def test_vectorized_matches_reference_on_small_states(net):
    # spec FRAGILE-2: cover the boundary n_i < coeff regime EXHAUSTIVELY for
    # small states (each count in 0..3), not just random large counts -- this is
    # where the homodimer/falling-factorial zeroing lives.
    from itertools import product
    from crnl.vectorized import compile_network, propensities_fast

    omega = 10.0
    compiled = compile_network(net, omega)
    ns = net.n_species
    for combo in product(range(4), repeat=ns):
        state = np.array(combo, dtype=np.int64)
        got = propensities_fast(compiled, state)
        exp = _reference_props(net, state, omega)
        # tight rtol; bit-exact == is not portable across multiply order (spec §6.2)
        assert np.allclose(got, exp, rtol=1e-12, atol=1e-15), (combo, got, exp)


def test_vectorized_matches_reference_on_random_large_states():
    from crnl.vectorized import compile_network, propensities_fast
    net = n_winner(6)
    omega = 200.0
    compiled = compile_network(net, omega)
    rng = np.random.default_rng(1)
    for _ in range(200):
        state = rng.integers(0, 60, size=net.n_species)
        got = propensities_fast(compiled, state)
        exp = _reference_props(net, state, omega)
        assert np.allclose(got, exp, rtol=1e-12, atol=1e-15)


def test_gillespie_fast_matches_reference_distribution():
    # The fast SSA is the same Markov chain as the reference (propensities agree
    # to 1e-12); assert their AM outcome DISTRIBUTIONS agree within statistics.
    # Bit-identical trajectories are not asserted (float multiply-order can flip
    # a measure-zero selection tie), but the chains are statistically identical.
    from crnl.vectorized import compile_network, gillespie_fast
    from crnl import gillespie, seed_for
    from crnl import classify

    net = approximate_majority()
    omega = 50
    n0 = np.array([26, 24, 0])
    trials = 800

    ref_x = 0
    for t in range(trials):
        r = gillespie(net, n0, omega, seed_for(omega, t))
        ref_x += classify.classify_am_outcome(r) == "X"

    compiled = compile_network(net, omega)
    fast_x = 0
    for t in range(trials):
        r = gillespie_fast(compiled, n0, seed_for(omega, t), species=list(net.species))
        fast_x += classify.classify_am_outcome(r) == "X"

    p_ref = ref_x / trials
    p_fast = fast_x / trials
    se = (p_ref * (1 - p_ref) / trials) ** 0.5
    assert abs(p_fast - p_ref) <= 4 * se + 1e-9, (p_ref, p_fast)


def test_gillespie_fast_absorbs_am():
    from crnl.vectorized import compile_network, gillespie_fast
    from crnl import seed_for
    net = approximate_majority()
    omega = 40
    compiled = compile_network(net, omega)
    r = gillespie_fast(compiled, np.array([21, 19, 0]), seed_for(omega, 3))
    assert r.absorbed
    assert r.n_final.sum() == omega  # count conserved


from crnl.stochastic import SSAResult


def _result(species, counts, absorbed=True):
    return SSAResult(t_final=1.0, n_final=np.array(counts), steps=1,
                     absorbed=absorbed, species=list(species))


def test_classify_winner_bins():
    from crnl.classify import classify_winner
    sp = ["X1", "X2", "X3", "B"]
    assert classify_winner(_result(sp, [7, 0, 0, 3])) == "X1"
    assert classify_winner(_result(sp, [0, 0, 0, 10])) == "blank"
    assert classify_winner(_result(sp, [4, 3, 0, 3])) == "coexist"
    assert classify_winner(_result(sp, [4, 0, 0, 6], absorbed=False)) == "undecided"


def test_classify_am_outcome_unchanged():
    # existing behavior preserved: AM single winners, blank -> "B", undecided
    from crnl.classify import classify_am_outcome
    sp = ["X", "Y", "B"]
    assert classify_am_outcome(_result(sp, [10, 0, 0])) == "X"
    assert classify_am_outcome(_result(sp, [0, 8, 0])) == "Y"
    assert classify_am_outcome(_result(sp, [0, 0, 12])) == "B"
    assert classify_am_outcome(_result(sp, [5, 5, 2], absorbed=False)) == "undecided"
