"""SSA instrumentation: the loop, its counters, and the reversible network's
first live homodimer path.

The instrumented loop exists because gillespie_fast returns only a final state --
no trajectory, no stopping predicate -- so Parts A-C need new loop machinery
anyway. The price of a second loop is DRIFT, so the first test pins them
bit-for-bit on the same seed. That is the whole reason the fast loop was left
unmodified.
"""

import numpy as np
import pytest

from crnl.cme import first_passage
from crnl.networks.am_reversible import (
    am_reversible, cycle_affinity, delta_star, initial_counts, reverse_pairing,
)
from crnl.stochastic import seed_for
from crnl.thermo import (
    FlipCounter, decompose, gillespie_instrumented, ln_multinomial,
)
from crnl.vectorized import compile_network, gillespie_fast


def _setup(gamma=0.3, omega=60):
    net = am_reversible(gamma)
    return net, reverse_pairing(net), compile_network(net, float(omega)), omega


# -- the instrumented loop IS the fast loop (spec test 13a) -----------------

def test_instrumented_matches_fast_bit_for_bit():
    """Same seed, same draws, same chain -- to the last digit of t_final.

    If this fails, do NOT relax it to approx: it means the draw order diverged,
    which invalidates every SSA number in this plan.
    """
    net, pairing, comp, omega = _setup()
    n0 = initial_counts(omega, 0.3, count_diff=4)
    a = gillespie_fast(comp, n0, seed_for(omega, 7, base=11),
                       max_steps=3000, species=list(net.species))
    b = gillespie_instrumented(comp, n0, seed_for(omega, 7, base=11), pairing,
                               max_steps=3000, species=list(net.species))
    assert b.steps == a.steps
    assert b.t_final == a.t_final
    assert np.array_equal(b.n_final, a.n_final)
    assert b.absorbed == a.absorbed


# -- the counter against an INDEPENDENT oracle (spec test 13b) --------------

def test_net_firings_reproduces_the_exact_cme_expectation():
    """The sampled mean <M> must match cme.first_passage's exact value.

    This is the test that has content. Checking `decompose(...)` against
    `lnW + M*ln(1/gamma)` instead is TAUTOLOGICAL -- both sides consume the same
    net_firings, so it passes even with the forward/reverse sign inverted
    (verified: flipping the sign gives -313.5 == -313.5, test still green). The
    CME value is computed by a completely different route, so it discriminates.
    """
    gamma, omega, trials = 0.3, 60, 500
    net, pairing, comp, _ = _setup(gamma, omega)
    theta = max(2, int(round(0.7 * delta_star(gamma) * omega)))
    n0 = initial_counts(omega, gamma,
                        count_diff=max(1, int(round(0.2 * delta_star(gamma) * omega))))

    def stop(n):
        return abs(int(n[0]) - int(n[1])) >= theta

    exact = first_passage(net, omega, float(omega), n0, stop, pairing)
    assert exact["valid"]

    firings = []
    for t in range(trials):
        r = gillespie_instrumented(comp, n0, seed_for(omega, t, base=5), pairing,
                                   stop=stop, max_steps=5_000_000,
                                   species=list(net.species))
        assert r.stopped
        firings.append(r.net_firings)
    mean = float(np.mean(firings))
    sem = float(np.std(firings, ddof=1) / np.sqrt(trials))
    assert abs(mean - exact["net_reaction_firings"]) < 4.0 * sem, (
        f"sampled <M>={mean:.1f}+-{sem:.1f} vs exact "
        f"{exact['net_reaction_firings']:.1f}")


def test_stop_predicate_halts_and_reports_it():
    gamma = 0.3
    net, pairing, comp, omega = _setup(gamma)
    n0 = initial_counts(omega, gamma, count_diff=2)
    theta = int(round(0.7 * delta_star(gamma) * omega))
    r = gillespie_instrumented(
        comp, n0, seed_for(omega, 1, base=2), pairing,
        stop=lambda n: abs(int(n[0]) - int(n[1])) >= theta,
        max_steps=200_000, species=list(net.species))
    assert r.stopped
    assert abs(int(r.n_final[0]) - int(r.n_final[1])) >= theta


def test_incomplete_pairing_is_rejected():
    """An unpaired reaction must RAISE, not be silently counted as a reverse.

    `pairing[j] > j` maps pairing[j] == -1 to False, i.e. to "reverse". With r3
    removed the loop otherwise runs happily and returns a wrong integer, while
    thermo.entropy_step raises on the same input. The guard must enforce the
    precondition it claims.
    """
    from crnl.reactions import ReactionNetwork

    net = am_reversible(0.3)
    trimmed = ReactionNetwork(species=list(net.species),
                              reactions=[r for i, r in enumerate(net.reactions)
                                         if i != 2],
                              name="am_reversible_broken")
    comp = compile_network(trimmed, 60.0)
    bad = np.array([3, -1, 0, 1, -1], dtype=np.int64)[:trimmed.n_reactions]
    with pytest.raises(ValueError, match="pairing"):
        gillespie_instrumented(comp, np.array([31, 29, 0]),
                               seed_for(60, 0, base=1), bad, max_steps=10)


# -- flip detection (spec test 15) ------------------------------------------

def test_flip_counter_counts_a_clean_flip_once():
    fc = FlipCounter(arm=0.3)
    for d in [0.4, 0.5, 0.45, -0.4, -0.5]:
        fc.update(d)
    assert fc.flips == 1
    assert fc.side == -1


def test_flip_counter_ignores_noise_recrossings():
    """Wobbling across zero without reaching the far arm is not a flip."""
    fc = FlipCounter(arm=0.3)
    for d in [0.4, 0.1, -0.05, 0.2, -0.2, 0.35, 0.05, -0.29]:
        fc.update(d)
    assert fc.flips == 0
    assert fc.side == +1


def test_flip_counter_counts_each_one_way_crossing():
    """A round trip is TWO flips. This is the convention that fixes the
    flips/T -> 1/tau identity used by Part B; 1/(2 tau) is the round-trip rate
    and using it makes tau_SSA exactly half of tau_CME."""
    fc = FlipCounter(arm=0.3)
    for d in [0.4, -0.4, 0.4]:
        fc.update(d)
    assert fc.flips == 2


def test_flip_counter_starts_unarmed_and_arms_without_counting():
    fc = FlipCounter(arm=0.3)
    fc.update(0.0)
    assert (fc.flips, fc.side) == (0, 0)
    fc.update(0.4)
    assert (fc.flips, fc.side) == (0, +1)


def test_instrumented_loop_counts_flips():
    gamma, omega = 0.45, 20
    net, pairing, comp, _ = _setup(gamma, omega)
    n0 = initial_counts(omega, gamma, count_diff=2)
    r = gillespie_instrumented(comp, n0, seed_for(omega, 0, base=1), pairing,
                               flip_arm=0.7 * delta_star(gamma), omega=omega,
                               max_steps=200_000, species=list(net.species))
    assert r.flips > 0


# -- SSA -> ODE for the FULL reversible network (spec test 3) ---------------

def test_reversible_ssa_converges_to_the_ode():
    """design.md 6's "single best test that the units convention is right", on
    the engine's first live homodimer path (all three reverses 2B->X+Y, 2X->B+X,
    2Y->B+Y are homodimers; irreversible AM had none).

    Asserted per-Omega against the Monte-Carlo SEM, NOT as errs[1] < errs[0]:
    both errors are dominated by sampling noise at reachable trial counts, so
    the ordering assertion is a coin flip (measured: it FAILS on seed base 18
    and passes on 17 and 19 with no code change).
    """
    from crnl.deterministic import integrate

    gamma = 0.7                                  # > GAMMA_C: single stable state
    net = am_reversible(gamma)
    pairing = reverse_pairing(net)
    x0 = np.array([0.6, 0.2, 0.2])
    t_end = 12.0
    target = integrate(net, x0, t_span=(0.0, t_end)).x[:, -1]

    for omega, trials in ((60, 400), (240, 400)):
        comp = compile_network(net, float(omega))
        n0 = np.array([int(round(f * omega)) for f in x0], dtype=np.int64)
        n0[2] = omega - n0[0] - n0[1]
        finals = np.array([
            gillespie_instrumented(comp, n0, seed_for(omega, t, base=17),
                                   pairing, t_max=t_end,
                                   species=list(net.species)).n_final / omega
            for t in range(trials)])
        mean = finals.mean(axis=0)
        sem = finals.std(axis=0, ddof=1) / np.sqrt(trials)
        dev = np.abs(mean - target)
        assert (dev < 4.0 * sem + 1.5 / omega).all(), (
            f"Omega={omega}: |mean-ode|={dev}, sem={sem}")


# -- gamma -> 0 recovers irreversible AM statistically (spec test 5) --------

def test_small_gamma_matches_irreversible_am_error_rate():
    """Conditional-on-decision error fraction, because the reversible model has
    no all-blank bin for gamma > 0.

    Omega and the bias are chosen so the all-blank bin is actually POPULATED in
    the irreversible arm -- at Omega=60 with count_diff=6 it is empty, which
    makes the 'conditional' denominator identical to the unconditional one and
    the test silently stops testing its own headline.
    """
    from crnl.networks import approximate_majority

    omega, trials, diff = 12, 3000, 2
    n_x = (omega + diff) // 2
    n0 = np.array([n_x, omega - n_x, 0], dtype=np.int64)

    comp_irr = compile_network(approximate_majority(), float(omega))
    wrong = decided = blank = 0
    for t in range(trials):
        r = gillespie_fast(comp_irr, n0, seed_for(omega, t, base=23))
        x, y = int(r.n_final[0]), int(r.n_final[1])
        if x == 0 and y == 0:
            blank += 1
            continue
        decided += 1
        wrong += y > x
    assert blank > 0, "all-blank bin empty: the conditional denominator is untested"
    p_irr = wrong / decided

    gamma = 1e-4
    net = am_reversible(gamma)
    pairing = reverse_pairing(net)
    comp_rev = compile_network(net, float(omega))
    theta = int(round(0.9 * delta_star(gamma) * omega))
    wrong_rev = 0
    for t in range(trials):
        r = gillespie_instrumented(
            comp_rev, n0, seed_for(omega, t, base=23), pairing,
            stop=lambda n: abs(int(n[0]) - int(n[1])) >= theta,
            max_steps=2_000_000, species=list(net.species))
        assert r.stopped
        wrong_rev += int(r.n_final[1]) > int(r.n_final[0])
    p_rev = wrong_rev / trials

    sem = np.sqrt(p_irr * (1 - p_irr) / decided + p_rev * (1 - p_rev) / trials)
    assert abs(p_rev - p_irr) < 4.0 * sem, f"{p_rev=} vs {p_irr=} sem={sem:.4f}"
