from __future__ import annotations

import numpy as np
import pytest

from crnl.networks.am_reversible import am_reversible, reverse_pairing, cycle_affinity
from crnl.cme import enumerate_states, generator, stationary, ep_rate


def test_enumerate_states_covers_the_simplex():
    states, index = enumerate_states(3, 5)
    assert len(states) == (5 + 1) * (5 + 2) // 2        # 21 for N=5
    assert all(s.sum() == 5 for s in states)
    assert len(index) == len(states)
    assert index[tuple(states[0])] == 0


def test_generator_rows_sum_to_zero():
    net = am_reversible(0.3)
    Q = generator(net, 6, 6.0)
    assert np.allclose(np.asarray(Q.sum(axis=1)).ravel(), 0.0, atol=1e-12)


def test_stationary_is_a_probability_vector_in_the_kernel():
    net = am_reversible(0.3)
    N = 8
    Q = generator(net, N, float(N))
    p = stationary(net, N, float(N))
    assert p.min() >= -1e-12
    assert p.sum() == pytest.approx(1.0)
    assert np.max(np.abs(p @ Q.toarray())) < 1e-9      # p Q = 0


def test_equilibrium_at_gamma_one_is_uniform_multinomial():
    # exact result: at gamma = 1 the stationary distribution is the uniform
    # multinomial, and detailed balance holds edge by edge.
    from crnl.thermo import ln_multinomial
    N = 9
    net = am_reversible(1.0)
    p = stationary(net, N, float(N))
    states, index = enumerate_states(3, N)
    ln_expect = np.array([ln_multinomial(s) for s in states])
    expect = np.exp(ln_expect - ln_expect.max())
    expect /= expect.sum()
    assert np.allclose(p, expect, atol=1e-10)


def test_ep_rate_equals_affinity_times_flux():
    """Schnakenberg: sigma == A * J exactly, with J the single cycle flux.

    This asserts the IDENTITY, not just that sigma is positive and rising -- the
    identity is what makes the EP machinery trustworthy, and it is free to check
    because the net fluxes come straight from the stationary distribution.
    """
    from crnl.cme import enumerate_states, stationary
    from crnl.networks.am_reversible import cycle_affinity

    N = 12
    net1 = am_reversible(1.0)
    assert ep_rate(net1, N, float(N), reverse_pairing(net1)) == pytest.approx(0.0, abs=1e-9)

    prev = -1.0
    for gamma in (0.1, 0.2, 0.3):
        net = am_reversible(gamma)
        pairing = reverse_pairing(net)
        sigma = ep_rate(net, N, float(N), pairing)
        assert sigma > 0.0
        assert sigma > prev                     # rises with gamma in this range
        prev = sigma

        # independent net cycle fluxes J_j = <a_forward - a_reverse>
        states, _ = enumerate_states(3, N)
        p = stationary(net, N, float(N))
        cs = net.stochastic_constants(float(N))
        S = net.stoichiometry_matrix().astype(int)
        J = np.zeros(3)
        for i, n in enumerate(states):
            a = net.propensities(n, cs)
            for j in range(3):
                if a[j] > 0 and (n + S[:, j]).min() >= 0:
                    J[j] += p[i] * a[j]
                rev = int(pairing[j])
                if a[rev] > 0 and (n + S[:, rev]).min() >= 0:
                    J[j] -= p[i] * a[rev]
        # a single cycle forces all three net fluxes equal
        assert J[0] == pytest.approx(J[1], rel=1e-8)
        assert J[1] == pytest.approx(J[2], rel=1e-8)
        A = cycle_affinity(net, pairing)
        assert sigma == pytest.approx(A * J[0], rel=1e-8)
