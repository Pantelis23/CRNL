"""Freeze-out as a time change (crnl/freezeout.py, FINDINGS.md 5.1).

Three groups:
  * the exact reduction -- the expanding SSA IS ordinary SSA truncated at
    internal time 1/H, checked bit-for-bit against an independent reference loop;
  * the fast cross-trial instrument -- agrees with the reference engine and
    keeps its bookkeeping straight;
  * the ingredients of the (3/2) ln Omega law -- d' = b*d exactly, b -> 1/3,
    and the transient deficit integral (2/3) ln(3/2).
"""

from __future__ import annotations

import numpy as np
import pytest

from crnl.deterministic import integrate
from crnl.expanding import gillespie_expanding
from crnl.cme import enumerate_states, generator
from crnl.freezeout import (
    am_generator,
    am_observables,
    am_order_exact,
    am_state_index,
    bimolecular_pairs,
    crossing_tau,
    deterministic_times,
    internal_clock_sweep,
    n_winner_observables,
)
from crnl.networks import approximate_majority, n_winner
from crnl.reactions import Reaction, ReactionNetwork
from crnl.stochastic import seed_for
from crnl.vectorized import Compiled, compile_network, propensities_fast


def _ssa_state_at(comp: Compiled, n0, rng, tau_max):
    """Independent reference: state of the ORDINARY chain at internal time tau_max.

    Deliberately not `gillespie_fast`, whose `while t < t_max` loop executes the
    step that *crosses* t_max and so returns the state one event too late.
    """
    S = comp.S
    n = np.array(n0, dtype=np.int64)
    t = 0.0
    while True:
        a = propensities_fast(comp, n)
        a0 = a.sum()
        if a0 <= 0.0:
            return n, "absorbed"
        tau = -np.log(rng.random()) / a0
        j = min(int(np.searchsorted(np.cumsum(a), rng.random() * a0)),
                comp.n_reactions - 1)
        if t + tau > tau_max:
            return n, "frozen"
        n = n + S[:, j]
        t += tau


@pytest.mark.parametrize("hubble", [0.05, 0.12, 0.3, 1.0, 3.0])
def test_expanding_ssa_is_ordinary_ssa_truncated_at_one_over_H(hubble):
    """The whole of FINDINGS 4-6 rests on this: same chain, finite time budget."""
    net = approximate_majority()
    comp = compile_network(net, 120)
    n0 = np.array([60, 60, 0])
    for trial in range(60):
        exp = gillespie_expanding(comp, n0, seed_for(120, trial, base=7),
                                  hubble=hubble, species=net.species)
        ref, status = _ssa_state_at(comp, n0, seed_for(120, trial, base=7),
                                    1.0 / hubble)
        assert np.array_equal(exp.n_final, ref), (
            f"H={hubble} trial={trial}: {exp.n_final} != {ref}")
        assert (exp.status == "frozen") == (status == "frozen")


def test_equivalence_holds_for_n_winner_too():
    """The reduction needs only uniform reaction order, not AM specifically."""
    net = n_winner(4)
    comp = compile_network(net, 80)
    n0 = np.zeros(len(net.species), dtype=np.int64)
    n0[:4] = 20
    for trial in range(30):
        exp = gillespie_expanding(comp, n0, seed_for(80, trial, base=3),
                                  hubble=0.2, species=net.species)
        ref, _ = _ssa_state_at(comp, n0, seed_for(80, trial, base=3), 5.0)
        assert np.array_equal(exp.n_final, ref)


def test_bimolecular_pairs_rejects_a_homodimer():
    net = ReactionNetwork(
        species=["A", "B"],
        reactions=[Reaction({"A": 2}, {"B": 2}, 1.0, name="dimer")],
        name="homodimer",
    )
    with pytest.raises(ValueError, match="coefficient 2"):
        bimolecular_pairs(compile_network(net, 10))


def test_bimolecular_pairs_reads_am():
    comp = compile_network(approximate_majority(), 40)
    sp_a, sp_b = bimolecular_pairs(comp)
    # X+Y->2B, B+X->2X, B+Y->2Y  with species order X, Y, B
    assert sorted(zip(sp_a.tolist(), sp_b.tolist())) == [(0, 1), (0, 2), (1, 2)]


def test_every_grid_point_receives_every_trial():
    """Bookkeeping guard: a fill bug would silently bias the mean, not crash."""
    comp = compile_network(approximate_majority(), 30)
    taus = np.geomspace(0.5, 40.0, 37)
    r = internal_clock_sweep(comp, [15, 15, 0], taus, 200,
                             np.random.default_rng(0),
                             lambda n: [np.ones(n.shape[0])])
    assert np.allclose(r["means"][0], 1.0), r["means"][0]


def test_instrument_matches_the_reference_expanding_engine():
    """Independent propensity and selection code; agreement is a real check."""
    from experiments.expansion import run_point

    omega = 60
    comp = compile_network(approximate_majority(), omega)
    hs = [0.08, 0.15, 0.3]
    taus = np.array(sorted(1.0 / h for h in hs))
    fast = internal_clock_sweep(comp, [30, 30, 0], taus, 8000,
                                np.random.default_rng(11), am_observables)
    for h in hs:
        i = int(np.argmin(np.abs(taus - 1.0 / h)))
        ref = run_point(omega, h, 3000, 0)
        sem = np.hypot(fast["sem0"][i], ref["order_sd"] / np.sqrt(3000))
        assert abs(fast["means"][0][i] - ref["order"]) < 4 * sem
        assert abs(fast["means"][1][i] - ref["relic"]) < 0.01


def test_n_winner_gap_convention_reduces_to_am_at_two_symbols():
    obs = n_winner_observables(2, convention="gap")
    n = np.array([[7, 3, 5], [4, 4, 2], [0, 0, 9]])
    assert np.allclose(obs(n)[0], am_observables(n)[0])


def test_n_winner_dominance_matches_expansion_radix():
    """Same formula as experiments/expansion_radix.py, so Sec.6 stays comparable."""
    obs = n_winner_observables(4, convention="dominance")
    n = np.array([[10, 2, 3, 1, 4], [5, 5, 5, 5, 0], [0, 0, 0, 0, 9]])
    got = obs(n)[0]
    for row, g in zip(n, got):
        committed = row[:4].astype(float)
        tot = committed.sum()
        want = 0.0 if tot == 0 else (committed.max()/tot - 0.25) / 0.75
        assert g == pytest.approx(want)
    assert np.allclose(obs(n)[1], [4.0, 4.0, 0.0])


def test_n_winner_observables_rejects_unknown_convention():
    with pytest.raises(ValueError, match="unknown convention"):
        n_winner_observables(3, convention="nope")


def test_crossing_tau_interpolates_and_reports_no_crossing():
    taus = np.array([1.0, 2.0, 4.0, 8.0])
    curve = np.array([0.1, 0.3, 0.7, 0.95])
    t = crossing_tau(taus, curve, 0.5)
    assert 2.0 < t < 4.0
    assert np.isnan(crossing_tau(taus, curve, 0.99))


def test_antisymmetric_mode_grows_at_exactly_b():
    """d = x - y obeys d' = b*d EXACTLY (not just to linear order).

    This is what makes the growth rate lambda = b* = 1/3 rather than a
    linearisation constant, and it is the whole basis of the (3/2) ln Omega law.
    """
    net = approximate_majority()
    S = net.stoichiometry_matrix()
    rng = np.random.default_rng(4)
    for _ in range(200):
        x = rng.random(3)
        x /= x.sum()
        dxdt = S @ net.fluxes(x)
        assert dxdt[0] - dxdt[1] == pytest.approx(x[2] * (x[0] - x[1]), abs=1e-14)


def test_symmetric_transient_deficit_is_two_thirds_ln_three_halves():
    """int_0^inf (1/3 - b) dt = (2/3) ln(3/2) from b(0) = 0.

    The exponential growth of d is exp(int b dt) = exp(t/3 - 0.2703), so this
    constant is the O(1) offset the (3/2) ln Omega law absorbs into its intercept.
    """
    net = approximate_majority()
    traj = integrate(net, [0.5, 0.5, 0.0], t_span=(0.0, 60.0), n_eval=60001)
    b = traj.x[2]
    deficit = np.trapezoid(1.0 / 3.0 - b, traj.t)
    assert deficit == pytest.approx((2.0 / 3.0) * np.log(1.5), rel=2e-4)
    assert b[-1] == pytest.approx(1.0 / 3.0, abs=1e-9)


@pytest.mark.parametrize("total", [8, 30])
def test_fast_am_generator_matches_the_network_agnostic_cme(total):
    """The Omega^3 exact route is only trustworthy if its generator is the same one."""
    net = approximate_majority()
    fast, xs, ys = am_generator(total)
    ref = generator(net, total, total)
    states, _ = enumerate_states(3, total)
    assert np.array_equal(states[:, 0], xs)
    assert np.array_equal(states[:, 1], ys)
    assert abs(fast - ref).max() < 1e-12
    assert np.allclose(np.asarray(fast.sum(axis=1)).ravel(), 0.0, atol=1e-12)


def test_am_state_index_is_a_bijection_on_the_simplex():
    xs, ys, index_of = am_state_index(12)
    assert len(xs) == 13 * 14 // 2
    assert np.array_equal(index_of(xs, ys), np.arange(len(xs)))
    assert (xs + ys <= 12).all()


def test_exact_cme_curve_agrees_with_the_sampled_one():
    """Two instruments, no shared code path: sparse linear algebra vs SSA."""
    omega = 90
    taus, curve = am_order_exact(omega, 26.0, 261)
    comp = compile_network(approximate_majority(), omega)
    fast = internal_clock_sweep(comp, [45, 45, 0], np.geomspace(2.0, 25.0, 40),
                                12000, np.random.default_rng(5), am_observables)
    sampled = fast["means"][0]
    exact_at = np.interp(fast["taus"], taus, curve)
    z = (sampled - exact_at) / np.maximum(fast["sem0"], 1e-9)
    assert np.abs(z).max() < 5.0, np.abs(z).max()
    t_exact = crossing_tau(taus[1:], curve[1:], 0.5)
    t_ssa = crossing_tau(fast["taus"], sampled, 0.5)
    assert abs(t_exact - t_ssa) < 0.15, (t_exact, t_ssa)


def test_deterministic_slopes_converge_to_three_halves_and_five_halves():
    """The ODE route reproduces both logarithms with no stochastics involved."""
    omegas = [40960, 163840, 655360]
    ts = [deterministic_times(w) for w in omegas]
    for lo, hi in zip(ts, ts[1:]):
        d = np.log(4.0)
        assert (hi["t_level"] - lo["t_level"]) / d == pytest.approx(1.5, abs=0.01)
        assert (hi["t_clear"] - lo["t_clear"]) / d == pytest.approx(2.5, abs=0.01)
    # Both approach from ABOVE. How FAR above at small Omega depends on the
    # quenched start convention, so this route pins the LIMITS, not the size of
    # the finite-Omega excess -- see FINDINGS.md Sec.5.1.
    small = [deterministic_times(w) for w in (40, 160)]
    assert (small[1]["t_level"] - small[0]["t_level"]) / np.log(4.0) > 1.5
    assert (small[1]["t_clear"] - small[0]["t_clear"]) / np.log(4.0) > 2.6


def test_deterministic_seed_scale_only_shifts_times():
    """seed_scale multiplies delta0, so it adds 3 ln(scale) and changes no slope."""
    a = deterministic_times(10240)
    b = deterministic_times(10240, seed_scale=2.0)
    shift = 3.0 * np.log(2.0)
    assert b["t_level"] - a["t_level"] == pytest.approx(-shift, abs=0.02)
    assert b["t_clear"] - a["t_clear"] == pytest.approx(-shift, abs=0.02)
    with pytest.raises(ValueError, match="seed_scale must be positive"):
        deterministic_times(100, seed_scale=0.0)
