from __future__ import annotations

import math

import numpy as np
import pytest

from crnl.networks.am_reversible import am_reversible, reverse_pairing, cycle_affinity
from crnl.thermo import ln_multinomial, entropy_step, decompose


def test_ln_multinomial_matches_factorials():
    for n in ([2, 3, 4], [0, 0, 5], [10, 0, 1]):
        N = sum(n)
        expect = math.lgamma(N + 1) - sum(math.lgamma(c + 1) for c in n)
        assert ln_multinomial(np.array(n)) == pytest.approx(expect)


def test_entropy_step_matches_exact_identity_at_every_gamma():
    """dS = lnW(n') - lnW(n) + s*ln(1/gamma), s = +1 forward / -1 reverse.

    Verified in the spec to 9.3e-15 over 1662 jumps. This is a stronger oracle
    than a gamma=1-only check and catches both homodimer factor-of-2 failure
    modes at every gamma.
    """
    for N in (7, 12):
        for gamma in (0.05, 0.3, 0.77, 1.0):
            net = am_reversible(gamma)
            pairing = reverse_pairing(net)
            S = net.stoichiometry_matrix()
            cs = net.stochastic_constants(float(N))
            checked = 0
            for nX in range(N + 1):
                for nY in range(N + 1 - nX):
                    n = np.array([nX, nY, N - nX - nY], dtype=np.int64)
                    a = net.propensities(n, cs)
                    for j in range(6):
                        if a[j] <= 0:
                            continue
                        n2 = n + S[:, j].astype(np.int64)
                        if n2.min() < 0:
                            continue
                        if net.propensities(n2, cs)[pairing[j]] <= 0:
                            continue
                        got = entropy_step(net, pairing, j, n, n2, cs)
                        s = 1.0 if j < 3 else -1.0
                        want = (ln_multinomial(n2) - ln_multinomial(n)
                                + s * math.log(1.0 / gamma))
                        assert got == pytest.approx(want, abs=1e-12)
                        checked += 1
            assert checked > 100          # the sweep really ran


def test_entropy_step_vanishes_on_a_cycle_at_equilibrium():
    # gamma = 1 is detailed balance: f1->f2->f3 returns to the start state and
    # the summed entropy production must be exactly 0.
    net = am_reversible(1.0)
    pairing = reverse_pairing(net)
    S = net.stoichiometry_matrix()
    cs = net.stochastic_constants(30.0)
    n = np.array([10, 10, 10], dtype=np.int64)
    total, state = 0.0, n
    for j in (0, 1, 2):
        nxt = state + S[:, j].astype(np.int64)
        total += entropy_step(net, pairing, j, state, nxt, cs)
        state = nxt
    assert np.array_equal(state, n)                # closed the cycle
    assert total == pytest.approx(0.0, abs=1e-12)


def test_cycle_sum_is_minus_three_ln_gamma():
    # away from equilibrium the same cycle yields the affinity, at ANY state
    gamma = 0.3
    net = am_reversible(gamma)
    pairing = reverse_pairing(net)
    S = net.stoichiometry_matrix()
    cs = net.stochastic_constants(60.0)
    A = cycle_affinity(net, pairing)
    for start in ([20, 20, 20], [30, 20, 10], [5, 5, 50]):
        state = np.array(start, dtype=np.int64)
        total = 0.0
        for j in (0, 1, 2):
            nxt = state + S[:, j].astype(np.int64)
            total += entropy_step(net, pairing, j, state, nxt, cs)
            state = nxt
        assert total == pytest.approx(A, abs=1e-12)
        assert total == pytest.approx(-3 * np.log(gamma), abs=1e-12)


def test_entropy_step_rejects_unpaired_reaction():
    from crnl.networks import approximate_majority
    net = approximate_majority()                   # irreversible: no pairs
    pairing = reverse_pairing(net)
    cs = net.stochastic_constants(30.0)
    n = np.array([10, 10, 10], dtype=np.int64)
    with pytest.raises(ValueError):
        entropy_step(net, pairing, 0, n, n, cs)


def test_decompose_sums_to_total():
    n0 = np.array([46, 34, 40], dtype=np.int64)
    n_stop = np.array([90, 10, 20], dtype=np.int64)
    A = -3 * np.log(0.3)
    d = decompose(n0, n_stop, net_reaction_firings=375.9, affinity=A)
    assert d["boundary"] == pytest.approx(ln_multinomial(n_stop) - ln_multinomial(n0))
    assert d["cycle"] == pytest.approx((A / 3.0) * 375.9)
    assert d["total"] == pytest.approx(d["boundary"] + d["cycle"])


def test_decompose_applies_the_factor_three_exactly_once():
    # regression: the /3 belongs to decompose alone. If a caller also divides,
    # results are 3x too small. Pin the arithmetic so the convention is explicit.
    n0 = np.array([10, 10, 10], dtype=np.int64)
    A = -3 * np.log(0.5)                      # = 3*ln2
    d = decompose(n0, None, net_reaction_firings=30.0, affinity=A, boundary=0.0)
    assert d["cycle"] == pytest.approx(np.log(2) * 30.0)   # (A/3)*30 = ln2*30
    assert d["total"] == pytest.approx(np.log(2) * 30.0)


def test_decompose_prefers_an_explicit_boundary():
    n0 = np.array([10, 10, 10], dtype=np.int64)
    d = decompose(n0, None, 0.0, 1.0, boundary=7.5)
    assert d["boundary"] == pytest.approx(7.5)
    with pytest.raises(ValueError):
        decompose(n0, None, 0.0, 1.0)         # neither boundary nor n_stop
