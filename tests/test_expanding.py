"""Tests for the expanding-volume (freeze-out) SSA (crnl/expanding.py).

The exact waiting-time math is the fragile core, so it is tested directly
against the analytic decaying-rate Poisson distribution, plus the H->0 reduction
to ordinary Gillespie and the qualitative freeze-out behavior.
"""

from __future__ import annotations

import numpy as np
import pytest

from crnl.networks import approximate_majority
from crnl.reactions import Reaction, ReactionNetwork
from crnl.vectorized import compile_network, gillespie_fast
from crnl.expanding import (
    next_event_time,
    gillespie_expanding,
    classify_freeze,
    common_order,
)
from crnl import seed_for


# --------------------------------------------------------------------------- #
# The exact waiting-time sampler                                              #
# --------------------------------------------------------------------------- #

def test_next_event_time_h0_is_ordinary_exponential():
    # lam <= 0 must reduce to the standard Gillespie wait -ln(u)/a0
    for u in (0.1, 0.5, 0.9):
        assert next_event_time(3.0, 0.0, u) == pytest.approx(-np.log(u) / 3.0)


def test_next_event_time_freeze_probability_matches_analytics():
    # For a rate a0*exp(-lam s): P(no event ever) = exp(-a0/lam).
    rng = np.random.default_rng(0)
    for a0, lam in [(5.0, 1.0), (2.0, 3.0), (10.0, 2.0)]:
        N = 200_000
        frozen = sum(next_event_time(a0, lam, rng.random()) is None for _ in range(N))
        p_emp = frozen / N
        p_ana = np.exp(-a0 / lam)
        assert p_emp == pytest.approx(p_ana, abs=0.005)


def test_next_event_time_cdf_matches_analytics():
    # Unconditional CDF: P(T <= tau) = 1 - exp(-(a0/lam)(1 - exp(-lam tau))).
    rng = np.random.default_rng(1)
    a0, lam = 5.0, 1.0
    N = 200_000
    taus = np.array([next_event_time(a0, lam, rng.random()) or np.inf for _ in range(N)])
    for tau in (0.1, 0.3, 0.8):
        emp = (taus <= tau).mean()
        ana = 1.0 - np.exp(-(a0 / lam) * (1.0 - np.exp(-lam * tau)))
        assert emp == pytest.approx(ana, abs=0.005)


# --------------------------------------------------------------------------- #
# The expanding SSA on AM                                                     #
# --------------------------------------------------------------------------- #

def test_common_order_rejects_mixed_order_network():
    # a0(t) is a single exponential only for uniform reaction order
    net = ReactionNetwork(
        species=["A", "B"],
        reactions=[Reaction({"A": 1}, {"B": 1}, 1.0),          # unimolecular
                   Reaction({"A": 1, "B": 1}, {"A": 2}, 1.0)],  # bimolecular
    )
    with pytest.raises(ValueError):
        common_order(compile_network(net, 10.0))


def test_am_common_order_is_two():
    assert common_order(compile_network(approximate_majority(), 50.0)) == 2


def test_h_to_zero_matches_ordinary_gillespie():
    # design.md-style check: tiny H must reproduce ordinary consensus (always
    # absorbs, ~50/50, no undecided).
    net = approximate_majority()
    omega = 80
    n0 = np.array([40, 40, 0])
    comp = compile_network(net, omega)
    trials = 2000
    names = list(net.species)

    ord_x = 0
    for t in range(trials):
        r = gillespie_fast(comp, n0, seed_for(omega, t), species=names)
        ord_x += r.n_final[0] > 0 and r.n_final[1] == 0

    exp_x = undecided = absorbed = 0
    for t in range(trials):
        r = gillespie_expanding(comp, n0, seed_for(omega, t), hubble=1e-9, species=names)
        absorbed += r.status == "absorbed"
        c = classify_freeze(r)
        undecided += c == "undecided"
        exp_x += c == "X"

    assert undecided == 0                       # nothing freezes at H->0
    assert absorbed == trials                   # everything resolves
    assert exp_x / trials == pytest.approx(ord_x / trials, abs=0.04)


def test_fast_expansion_freezes_undecided():
    net = approximate_majority()
    omega = 80
    n0 = np.array([40, 40, 0])
    comp = compile_network(net, omega)
    frozen = undecided = 0
    for t in range(500):
        r = gillespie_expanding(comp, n0, seed_for(omega, t), hubble=5.0,
                                species=list(net.species))
        frozen += r.frozen
        undecided += classify_freeze(r) == "undecided"
    assert frozen == 500                        # fast expansion always freezes
    assert undecided / 500 > 0.9                # ...mid-decision


def test_order_parameter_decreases_with_expansion():
    net = approximate_majority()
    omega = 120
    n0 = np.array([60, 60, 0])
    comp = compile_network(net, omega)
    names = list(net.species)

    def order(H):
        vals = []
        for t in range(400):
            r = gillespie_expanding(comp, n0, seed_for(omega, t, base=int(H * 100) + 1),
                                    hubble=H, species=names)
            x, y = int(r.n_final[0]), int(r.n_final[1])
            vals.append(abs(x - y) / (x + y) if (x + y) else 1.0)
        return np.mean(vals)

    o = [order(H) for H in (0.01, 0.1, 0.5, 2.0)]
    # monotonically non-increasing, spanning decided -> undecided
    assert all(o[i] >= o[i + 1] - 0.02 for i in range(len(o) - 1))
    assert o[0] > 0.95 and o[-1] < 0.2


def test_negative_hubble_is_rejected():
    # contraction (H<0) makes the rate grow; the closed-form inversion here is
    # only valid for expansion, so it must be refused rather than silently wrong.
    net = approximate_majority()
    comp = compile_network(net, 60.0)
    with pytest.raises(ValueError):
        gillespie_expanding(comp, np.array([30, 30, 0]), seed_for(60, 0), hubble=-0.5,
                            species=list(net.species))


def test_critical_h_interpolates_crossing():
    from experiments.expansion_radix import critical_h
    rows = [
        {"hubble": 0.01, "D": 1.0},
        {"hubble": 0.1, "D": 0.6},
        {"hubble": 1.0, "D": 0.2},   # crossing 0.5 between H=0.1 and H=1.0
    ]
    hs = critical_h(rows, level=0.5)
    assert 0.1 < hs < 1.0
    # returns nan if the level is never crossed
    assert np.isnan(critical_h([{"hubble": 0.1, "D": 0.9},
                                {"hubble": 1.0, "D": 0.8}], level=0.5))


def test_n_winner_bigger_alphabet_freezes_easier():
    # ties radix + expansion: at a fixed intermediate H, a larger alphabet is
    # more frozen (lower winner-dominance D) than a smaller one.
    from experiments.expansion_radix import run_point
    d2 = run_point(2, 160, 0.08, 800, 0)["D"]
    d8 = run_point(8, 160, 0.08, 800, 0)["D"]
    assert d2 - d8 > 0.1       # n=8 is clearly more frozen (less decided) than n=2


def test_expanding_conserves_count():
    net = approximate_majority()
    omega = 100
    comp = compile_network(net, omega)
    for H in (0.0, 0.2, 1.0):
        r = gillespie_expanding(comp, np.array([50, 50, 0]), seed_for(omega, 1),
                                hubble=H, species=list(net.species))
        assert int(r.n_final.sum()) == omega     # every reaction is 2->2
