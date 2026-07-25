from __future__ import annotations

import numpy as np
import pytest

from crnl.networks.am_reversible import am_reversible, reverse_pairing, GAMMA_C


def test_structure_and_order():
    net = am_reversible(0.3)
    assert net.species == ["X", "Y", "B"]
    assert net.n_reactions == 6
    # documented order: f1, f2, f3, r1, r2, r3
    assert net.reactions[0].reactants == {"X": 1, "Y": 1}
    assert net.reactions[0].products == {"B": 2}
    assert net.reactions[1].reactants == {"B": 1, "X": 1}
    assert net.reactions[1].products == {"X": 2}
    assert net.reactions[2].reactants == {"B": 1, "Y": 1}
    assert net.reactions[2].products == {"Y": 2}
    assert net.reactions[3].reactants == {"B": 2}
    assert net.reactions[3].products == {"X": 1, "Y": 1}
    assert net.reactions[4].reactants == {"X": 2}
    assert net.reactions[4].products == {"B": 1, "X": 1}
    assert net.reactions[5].reactants == {"Y": 2}
    assert net.reactions[5].products == {"B": 1, "Y": 1}


def test_rates_scale_with_gamma():
    net = am_reversible(0.25, k=2.0)
    assert [r.k for r in net.reactions[:3]] == [2.0, 2.0, 2.0]
    assert [r.k for r in net.reactions[3:]] == [0.5, 0.5, 0.5]


def test_all_reactions_are_two_to_two():
    net = am_reversible(0.3)
    S = net.stoichiometry_matrix()
    assert np.allclose(S.sum(axis=0), 0.0)          # count conserved
    for r in net.reactions:
        assert sum(r.reactants.values()) == 2
        assert sum(r.products.values()) == 2


def test_reverses_are_homodimers():
    # the first live use of the engine's homodimer path
    net = am_reversible(0.3)
    for r in net.reactions[3:]:
        assert len(r.reactants) == 1                 # one species...
        assert list(r.reactants.values()) == [2]     # ...with coefficient 2


def test_reverse_pairing_is_derived():
    net = am_reversible(0.3)
    assert list(reverse_pairing(net)) == [3, 4, 5, 0, 1, 2]


def test_reverse_pairing_reports_unpaired():
    # an irreversible network has no pairs at all -> all -1
    from crnl.networks import approximate_majority
    assert list(reverse_pairing(approximate_majority())) == [-1, -1, -1]


def test_gamma_c_constant():
    assert GAMMA_C == 0.5


def test_rejects_negative_gamma():
    with pytest.raises(ValueError):
        am_reversible(-0.1)


def test_gamma_zero_reproduces_irreversible_am_ode():
    """At gamma = 0 the reversible ODE must BE the irreversible AM ODE.

    A deterministic control, exact and cheap: it catches a sign error or a
    mis-built reverse reaction before any stochastic run. Parallels
    test_n_winner_2_reduces_to_am.
    """
    from crnl.networks import approximate_majority
    rev = am_reversible(0.0)
    am = approximate_majority()
    rng = np.random.default_rng(0)
    for _ in range(20):
        x, y = rng.random(2) * 0.4
        state = [x, y, 1 - x - y]
        assert np.allclose(rev.rhs(state), am.rhs(state), atol=1e-12)
