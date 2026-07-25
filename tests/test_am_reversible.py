from __future__ import annotations

import warnings

import numpy as np
import pytest

from crnl.networks.am_reversible import am_reversible, reverse_pairing, GAMMA_C
from crnl.reactions import Reaction, ReactionNetwork


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


def test_rejects_nan_gamma():
    with pytest.raises(ValueError):
        am_reversible(float("nan"))


def test_reverse_pairing_raises_on_ambiguity():
    # two reverses of the same forward: A->B has two candidate reverses
    # (B->A twice, different rates), so the pairing cannot be resolved.
    net = ReactionNetwork(
        species=["A", "B"],
        reactions=[
            Reaction({"A": 1}, {"B": 1}, 1.0, name="f:A->B"),
            Reaction({"B": 1}, {"A": 1}, 0.5, name="r1:B->A"),
            Reaction({"B": 1}, {"A": 1}, 0.7, name="r2:B->A"),
        ],
        name="ambiguous",
    )
    with pytest.raises(ValueError):
        reverse_pairing(net)


def test_gamma_zero_reproduces_irreversible_am_ode():
    """At gamma = 0 the reversible ODE must BE the irreversible AM ODE.

    This is a deterministic, exact and cheap control: it establishes that the
    forward triple matches irreversible AM and that the reverse contributions
    vanish at gamma = 0 (every reverse rate is exactly 0.0, so the reverse
    stoichiometry never enters `rhs`). It does NOT check that the reverse
    stoichiometry itself is correctly built -- that is pinned separately by
    test_reverses_are_homodimers and test_structure_and_order. Parallels
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


def test_homodimer_units_on_reverse_reactions():
    # FRAGILE: first live use of the engine's homodimer path. The reverse
    # constants must be 2*gamma*k/Omega, not gamma*k/Omega.
    gamma, k, omega = 0.3, 1.0, 100.0
    net = am_reversible(gamma, k)
    cs = net.stochastic_constants(omega)
    assert cs[0] == pytest.approx(k / omega)              # hetero forward
    assert cs[1] == pytest.approx(k / omega)
    assert cs[2] == pytest.approx(k / omega)
    for j in (3, 4, 5):
        assert cs[j] == pytest.approx(2 * gamma * k / omega)   # homodimer


def test_homodimer_propensity_is_combinatorial():
    # a = c*n(n-1)/2, NOT c*n^2. n=1 is the discriminating case: the correct
    # form is exactly 0 there, a c*n^2 bug is not.
    gamma, k, omega = 0.3, 1.0, 100.0
    net = am_reversible(gamma, k)
    cs = net.stochastic_constants(omega)
    c_r2 = cs[4]
    for nX in (0, 1, 2, 50):
        a = net.propensities(np.array([nX, 0, 0]), cs)[4]   # r2: 2X -> B+X
        assert a == pytest.approx(c_r2 * nX * (nX - 1) / 2)
    # explicitly pin the discriminating case
    assert net.propensities(np.array([1, 0, 0]), cs)[4] == 0.0
    # (a c*n^2 bug would give c_r2 * 1**2 = c_r2 != 0 here, unlike the correct
    # c*n(n-1)/2 form pinned above)


def test_cycle_affinity_from_network_k():
    from crnl.networks.am_reversible import cycle_affinity
    for gamma in (0.05, 0.3, 0.5, 0.9):
        net = am_reversible(gamma)
        A = cycle_affinity(net, reverse_pairing(net))
        assert A == pytest.approx(-3.0 * np.log(gamma))
    # equilibrium: no drive
    net1 = am_reversible(1.0)
    assert cycle_affinity(net1, reverse_pairing(net1)) == pytest.approx(0.0)


def test_cycle_affinity_is_not_the_stochastic_constant_value():
    # the trap: computing the affinity from stochastic_constants instead of k
    # is off by -3*ln(2) because of the homodimer factor.
    from crnl.networks.am_reversible import cycle_affinity
    gamma, omega = 0.3, 100.0
    net = am_reversible(gamma)
    A = cycle_affinity(net, reverse_pairing(net))
    cs = net.stochastic_constants(omega)
    A_wrong = float(np.log(cs[0] * cs[1] * cs[2]) - np.log(cs[3] * cs[4] * cs[5]))
    assert A == pytest.approx(-3 * np.log(gamma))
    assert A_wrong == pytest.approx(A - 3 * np.log(2))
    assert abs(A - A_wrong) > 2.0        # the error is ~2.079, not subtle


def test_affinity_at_gamma_c_is_three_ln_two():
    from crnl.networks.am_reversible import cycle_affinity
    net = am_reversible(GAMMA_C)
    A = cycle_affinity(net, reverse_pairing(net))
    assert A == pytest.approx(3 * np.log(2))
    assert A == pytest.approx(2.0794415417, abs=1e-9)


def test_cycle_affinity_at_gamma_zero_is_deliberately_inf():
    # gamma = 0 makes every reverse rate 0, so the cycle is irreversible and
    # the drive is unbounded. This must be a deliberate `inf` return, not a
    # stumble into divide-by-zero in np.log (which would fail under
    # filterwarnings = error).
    from crnl.networks.am_reversible import cycle_affinity
    net = am_reversible(0.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        A = cycle_affinity(net, reverse_pairing(net))
    assert A == float("inf")


def test_cycle_affinity_raises_with_no_reversible_pairs():
    from crnl.networks.am_reversible import cycle_affinity
    from crnl.networks import approximate_majority
    net = approximate_majority()
    with pytest.raises(ValueError):
        cycle_affinity(net, reverse_pairing(net))


def test_cycle_affinity_raises_on_disjoint_cycles():
    # two independent reversible pairs (A<->B, C<->D) sum to zero just like a
    # real cycle does, so a naive "does the sum vanish" guard would accept
    # this and silently return a garbage combined number. The real
    # precondition -- a 1-D cycle space -- must reject it.
    from crnl.networks.am_reversible import cycle_affinity
    net = ReactionNetwork(
        species=["A", "B", "C", "D"],
        reactions=[
            Reaction({"A": 1}, {"B": 1}, 1.0, name="f1:A->B"),
            Reaction({"B": 1}, {"A": 1}, 0.5, name="r1:B->A"),
            Reaction({"C": 1}, {"D": 1}, 2.0, name="f2:C->D"),
            Reaction({"D": 1}, {"C": 1}, 0.7, name="r2:D->C"),
        ],
        name="disjoint",
    )
    with pytest.raises(ValueError):
        cycle_affinity(net, reverse_pairing(net))


def test_cycle_affinity_is_order_independent():
    # false-rejection regression: picking "forward" by list-index convention
    # (rather than the actual cycle-space computation) can reject a
    # perfectly well defined network just because reactions were listed in a
    # different order. Reorder AM's reactions to [r1, f2, f3, f1, r2, r3]
    # (r1 listed before its forward f1) and require the same affinity as the
    # canonical order.
    from crnl.networks.am_reversible import cycle_affinity
    gamma = 0.3
    canonical = am_reversible(gamma)
    A_canonical = cycle_affinity(canonical, reverse_pairing(canonical))

    r1, f2, f3, f1, r2, r3 = (
        canonical.reactions[3],
        canonical.reactions[1],
        canonical.reactions[2],
        canonical.reactions[0],
        canonical.reactions[4],
        canonical.reactions[5],
    )
    reordered = ReactionNetwork(
        species=canonical.species,
        reactions=[r1, f2, f3, f1, r2, r3],
        name="am-reordered",
    )
    A_reordered = cycle_affinity(reordered, reverse_pairing(reordered))
    assert A_reordered == pytest.approx(A_canonical)
    assert A_reordered == pytest.approx(-3.0 * np.log(gamma))


def test_symmetric_point_is_a_fixed_point_at_every_gamma():
    # surprising but exact: at (1/3,1/3,1/3) the three net fluxes are equal and
    # the three stoichiometric vectors sum to zero, so S.v = 0 for every gamma.
    for gamma in (0.0, 0.2, 0.5, 0.8, 1.0):
        net = am_reversible(gamma)
        assert np.allclose(net.rhs([1 / 3, 1 / 3, 1 / 3]), 0.0, atol=1e-12)


def test_decision_mode_eigenvalue_matches_closed_form():
    """Both closed-form eigenvalues, against the engine's numeric Jacobian.

    The engine returns the FULL 3x3 Jacobian in (x, y, b). The closed forms are
    properties of the REDUCED system with b = 1 - x - y, so the reduction must be
    applied first: dF/dx|reduced = dF/dx - dF/db, since db/dx = db/dy = -1.
    Skipping that step happens to give the right antisymmetric eigenvalue by
    coincidence, but the symmetric one comes out a factor of 3 wrong
    (-1/3 instead of -1 at gamma=0).
    """
    from crnl.deterministic import jacobian
    from crnl.networks.am_reversible import lambda_antisym, lambda_sym
    for gamma in (0.0, 0.1, 0.25, 0.4, 0.5, 0.6):
        J = jacobian(am_reversible(gamma), [1 / 3, 1 / 3, 1 / 3])
        Jr = J[:2, :2] - J[:2, 2:3]        # eliminate b by the conservation law
        assert Jr[0, 0] - Jr[0, 1] == pytest.approx(lambda_antisym(gamma), abs=1e-9)
        assert Jr[0, 0] + Jr[0, 1] == pytest.approx(lambda_sym(gamma), abs=1e-9)


def test_restoring_gain_recovers_irreversible_am_at_gamma_zero():
    from crnl.networks.am_reversible import lambda_antisym
    # design.md 2.3's AM saddle eigenvalue is +1/3; it is the gamma->0 limit here
    assert lambda_antisym(0.0) == pytest.approx(1 / 3)
    assert lambda_antisym(GAMMA_C) == pytest.approx(0.0)
    assert lambda_antisym(0.6) < 0.0


def test_pure_corner_is_not_a_fixed_point():
    for gamma in (0.25, 0.7):
        net = am_reversible(gamma)
        r = net.rhs([1.0, 0.0, 0.0])
        assert r[0] == pytest.approx(-gamma)      # dx/dt = -gamma
        assert r[2] == pytest.approx(+gamma)


def test_attractors_and_blank_population():
    from crnl.networks.am_reversible import fixed_points
    for gamma in (0.05, 0.2, 0.3, 0.45):
        fps = fixed_points(gamma)
        assert len(fps) == 3                       # 3, not 4: all-blank left the simplex
        attractors = [f for f in fps if f["kind"] == "attractor"]
        assert len(attractors) == 2
        for f in attractors:
            assert f["b"] == pytest.approx(gamma / (1 + gamma))   # b* = g/(1+g)
            assert np.allclose(am_reversible(gamma).rhs([f["x"], f["y"], f["b"]]),
                               0.0, atol=1e-9)


def test_only_symmetric_point_above_gamma_c():
    from crnl.networks.am_reversible import fixed_points, delta_star
    for gamma in (0.5, 0.7):
        fps = fixed_points(gamma)
        assert len(fps) == 1
        assert fps[0]["kind"] == "symmetric"
        assert delta_star(gamma) == 0.0


def test_pitchfork_scaling_closed_form():
    from crnl.networks.am_reversible import delta_star
    # delta* / sqrt(gamma_c - gamma) -> 4*sqrt(2)/3 as gamma -> gamma_c
    target = 4 * np.sqrt(2) / 3
    ratio = delta_star(0.4999) / np.sqrt(GAMMA_C - 0.4999)
    assert ratio == pytest.approx(target, rel=1e-3)


def test_thresholds_are_reachable_for_every_gamma():
    # regression for a real design bug: a FIXED theta=0.5 is unreachable above
    # gamma~0.417 (delta*(0.49) = 0.187), which turns "decide" into "fluctuate
    # past the attractor" and would make a flip counter read zero forever.
    from crnl.networks.am_reversible import theta_decide, delta_star
    for gamma in (0.05, 0.1, 0.2, 0.3, 0.4, 0.45, 0.49):
        assert theta_decide(gamma) < delta_star(gamma)
    # and the fixed value really is unreachable up there, which is why we scale
    assert delta_star(0.49) < 0.5
    assert delta_star(0.4) > 0.5


def test_initial_counts_follow_project_convention():
    # design.md 4 / restoration_wall.py: B(0) = 0, committed split to realise
    # the bias. Total must be exactly omega.
    from crnl.networks.am_reversible import initial_counts, delta_star
    omega, gamma = 120, 0.3
    n0 = initial_counts(omega, gamma)
    assert n0.sum() == omega
    assert n0[2] == 0                                  # B starts empty
    realised = (n0[0] - n0[1]) / omega
    assert realised == pytest.approx(0.2 * delta_star(gamma), abs=2 / omega)
    assert n0[0] > n0[1]                               # X is the favoured symbol


def test_initial_counts_by_explicit_count_difference():
    """The gamma sweep must fix the integer count difference, not a fraction.

    delta_0 = 0.2*delta_star(gamma) is NOT representable on an integer lattice at
    any Omega, so a fraction-driven sweep makes the realised bias jitter
    non-monotonically with gamma -- and one molecule of bias is worth ~20 k_B of
    dissipation, which is enough to fake a fold-back in the headline curve. So
    `initial_counts` accepts an explicit count difference, and callers sweeping
    gamma use it.
    """
    from crnl.networks.am_reversible import initial_counts
    for omega in (60, 120, 180):
        for d in (4, 8, 12):
            n0 = initial_counts(omega, 0.3, count_diff=d)
            assert n0.sum() == omega
            assert n0[2] == 0
            assert int(n0[0]) - int(n0[1]) == d          # exact, by construction


def test_initial_counts_rejects_an_unreachable_bias():
    from crnl.networks.am_reversible import initial_counts
    with pytest.raises(ValueError):
        initial_counts(60, 0.3, count_diff=200)          # exceeds omega


def test_degenerate_counts_absorption_structure():
    """No absorbing state for N >= 2; EVERY state absorbing for N <= 1.

    All six reactions are bimolecular, so with 0 or 1 molecules nothing can fire.
    (At N = 2 the chain is also reducible -- three closed 2-state classes -- so
    experiments and other tests must avoid N <= 2.)
    """
    net = am_reversible(0.3)
    S = net.stoichiometry_matrix().astype(int)

    def total_propensity(n, N):
        cs = net.stochastic_constants(float(N))
        a = net.propensities(np.array(n), cs)
        # only count jumps that stay in the simplex
        return sum(a[j] for j in range(6)
                   if a[j] > 0 and (np.array(n) + S[:, j]).min() >= 0)

    # N <= 1: everything is absorbing
    for n in ([1, 0, 0], [0, 1, 0], [0, 0, 1]):
        assert total_propensity(n, 1) == 0.0

    # N >= 3: nothing is absorbing, including the corners
    for N in (3, 5, 7):
        for nX in range(N + 1):
            for nY in range(N + 1 - nX):
                n = [nX, nY, N - nX - nY]
                assert total_propensity(n, N) > 0.0, f"absorbing at {n}, N={N}"


def test_symmetric_point_stability_matches_the_bifurcation():
    # kind="symmetric" alone can't distinguish the saddle from the sole
    # attractor -- that's the whole content of the gamma_c bifurcation, so
    # `stable` must carry it.
    from crnl.networks.am_reversible import fixed_points
    fps_below = fixed_points(0.2)
    sym_below = next(f for f in fps_below if f["kind"] == "symmetric")
    assert sym_below["stable"] is False           # saddle below GAMMA_C

    fps_above = fixed_points(0.7)
    sym_above = next(f for f in fps_above if f["kind"] == "symmetric")
    assert sym_above["stable"] is True             # sole attractor above GAMMA_C


def test_landscape_functions_reject_invalid_gamma():
    # am_reversible already guards its input; the four landscape functions
    # must match, not silently return simplex-violating or inf/nan garbage.
    from crnl.networks.am_reversible import (
        lambda_antisym, lambda_sym, delta_star, fixed_points,
    )
    for fn in (lambda_antisym, lambda_sym, delta_star, fixed_points):
        with pytest.raises(ValueError):
            fn(-0.1)
        with pytest.raises(ValueError):
            fn(float("nan"))
