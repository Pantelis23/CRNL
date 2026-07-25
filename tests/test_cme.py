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


def test_stationary_raises_for_a_reducible_chain():
    net = am_reversible(0.3)
    for total in (0, 1, 2):
        with pytest.raises(ValueError):
            stationary(net, total, float(max(total, 1)))
    # total = 3 is the smallest irreducible case and must succeed
    p = stationary(net, 3, 3.0)
    assert p.sum() == pytest.approx(1.0)


def test_stationary_guard_prevents_a_silent_nan_return():
    """Without the total < 3 guard this would emit a bare MatrixRankWarning and
    return an all-NaN vector instead of raising."""
    import warnings
    net = am_reversible(0.3)
    with pytest.raises(ValueError):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            stationary(net, 1, 1.0)


def test_first_passage_matches_a_hand_solvable_case():
    # Absorb as soon as X is gone or Y is gone, from a symmetric start on a tiny
    # simplex: mean time must be positive and finite, split ~ 1/2 by symmetry.
    from crnl.cme import first_passage
    net = am_reversible(0.3)
    N = 9
    res = first_passage(
        net, N, float(N), start=np.array([3, 3, 3]),
        is_absorbing=lambda n: n[0] == 0 or n[1] == 0,
    )
    assert res["valid"]
    assert res["mean_time"] > 0.0
    assert np.isfinite(res["mean_time"])
    assert res["split"] == pytest.approx(0.5, abs=0.05)


def test_first_passage_split_follows_the_bias():
    from crnl.cme import first_passage
    net = am_reversible(0.3)
    N = 12
    biased = first_passage(
        net, N, float(N), start=np.array([8, 2, 2]),
        is_absorbing=lambda n: n[0] == 0 or n[1] == 0,
    )
    assert biased["split"] > 0.8              # X strongly favoured


def test_first_passage_guard_rejects_a_broken_solve():
    """The guard must FAIL a case that is actually broken.

    Not an if/else: this asserts valid is False. The geometry matters -- the
    pathology needs the Part-B setup (start at the attractor, absorb on the far
    side), where the solve returns a NEGATIVE mean time (~-6e16) at Omega=30,
    gamma=0.02. With the easier "X or Y extinct" absorbing set the same
    parameters solve perfectly happily, which is why an earlier version of this
    test had zero coverage of a mandatory guard.
    """
    from crnl.cme import first_passage
    from crnl.networks.am_reversible import fixed_points
    gamma, N = 0.02, 30
    net = am_reversible(gamma)
    att = max((f for f in fixed_points(gamma) if f["kind"] == "attractor"),
              key=lambda f: f["x"])
    n_x = int(round(att["x"] * N))
    n_b = int(round(att["b"] * N))
    start = np.array([n_x, N - n_x - n_b, n_b])
    res = first_passage(net, N, float(N), start,
                        lambda n: int(n[0]) - int(n[1]) <= -0.3 * N)
    assert res["valid"] is False
    assert (res["mean_time"] < 0
            or res["mean_time"] > 1e10
            or res["residual"] > 1e-8)


def test_first_passage_accepts_a_healthy_solve():
    # the complement: a well-conditioned problem must pass, with a tiny residual
    from crnl.cme import first_passage
    net = am_reversible(0.3)
    res = first_passage(net, 12, 12.0, np.array([7, 5, 0]),
                        lambda n: abs(int(n[0]) - int(n[1])) >= 6)
    assert res["valid"] is True
    assert res["residual"] < 1e-10
    assert 0.0 < res["mean_time"] < 1e6


def test_decomposition_composes_with_first_passage():
    """THE test that pins the two halves together.

    `decompose` applies affinity/3; `first_passage` returns raw per-reaction
    firings. If either side also divided by 3, every dissipation number would be
    exactly 3x wrong -- and no test of either function alone can see it. This
    checks the composed result against an independent solve of Qtt V = -u with
    u_i = sum_j a_j * dS_j, i.e. entropy production integrated directly.
    """
    import scipy.sparse.linalg as spla
    from crnl.cme import first_passage, generator, enumerate_states
    from crnl.thermo import decompose, entropy_step

    gamma, N = 0.3, 24
    net = am_reversible(gamma)
    pairing = reverse_pairing(net)
    A = cycle_affinity(net, pairing)
    start = np.array([14, 10, 0], dtype=np.int64)
    theta_counts = 12

    def absorbing(n):
        return abs(int(n[0]) - int(n[1])) >= theta_counts

    fp = first_passage(net, N, float(N), start, absorbing, pairing)
    assert fp["valid"]
    composed = decompose(start, None, fp["net_reaction_firings"], A,
                         boundary=fp["boundary"])["total"]

    # independent reference: integrate entropy production directly
    states, index = enumerate_states(3, N)
    absorb_mask = np.array([absorbing(s) for s in states])
    Q = generator(net, N, float(N))
    trans = np.where(~absorb_mask)[0]
    tmap = {int(i): k for k, i in enumerate(trans)}
    Qtt = Q[trans][:, trans].tocsr()
    cs = net.stochastic_constants(float(N))
    S = net.stoichiometry_matrix().astype(np.int64)
    u = np.zeros(len(trans))
    for i in trans:
        a = net.propensities(states[i], cs)
        acc = 0.0
        for j in range(6):
            if a[j] <= 0:
                continue
            n2 = states[i] + S[:, j]
            if n2.min() < 0 or net.propensities(n2, cs)[pairing[j]] <= 0:
                continue
            acc += a[j] * entropy_step(net, pairing, j, states[i], n2, cs)
        u[tmap[int(i)]] = acc
    V = spla.spsolve(Qtt, -u)
    reference = float(V[tmap[int(index[tuple(start)])]])

    assert composed == pytest.approx(reference, rel=1e-6)


def test_first_passage_requires_reachable_absorbing_set():
    from crnl.cme import first_passage
    net = am_reversible(0.3)
    with pytest.raises(ValueError):
        first_passage(net, 9, 9.0, start=np.array([3, 3, 3]),
                      is_absorbing=lambda n: False)      # nothing absorbs
