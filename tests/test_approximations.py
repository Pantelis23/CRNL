"""Tests for the coarse-graining hierarchy (FINDINGS 21).

The two that matter are conservation and the small-step limit. Both integrators
must preserve the conserved totals EXACTLY -- a step is always a whole combination
of stoichiometry columns -- and tau-leaping at a small enough window must agree
with the exact SSA, since that is the claim the whole hierarchy is calibrated on.
"""
from __future__ import annotations

import numpy as np
import pytest

from crnl.approximations import cle_run, tau_leap_run
from crnl.networks.am_reversible import am_reversible, reverse_pairing
from crnl.thermo import gillespie_instrumented
from crnl.vectorized import compile_network

N0 = np.array([40, 20, 30])
TOTAL = 90


def _comp(gamma=0.3):
    return compile_network(am_reversible(gamma), 90.0)


@pytest.mark.parametrize("runner,kw", [(cle_run, {"dt": 0.02}),
                                       (tau_leap_run, {"tau": 0.05})])
def test_conservation_is_exact(runner, kw):
    """Not approximate: every step is a combination of conserving columns."""
    for seed in range(12):
        r = runner(_comp(), N0, np.random.default_rng(seed), t_max=20.0, **kw)
        assert float(r.n_final.sum()) == pytest.approx(TOTAL, abs=1e-9), r.n_final


@pytest.mark.parametrize("runner,kw", [(cle_run, {"dt": 0.02}),
                                       (tau_leap_run, {"tau": 0.05})])
def test_counts_never_go_negative(runner, kw):
    """Negativity is handled by halving the step, never by clipping a species.

    Clipping would break conservation and would hand the minority species a floor
    it did not earn -- the failure mode behind three withdrawn results here.
    """
    for seed in range(12):
        r = runner(_comp(), N0, np.random.default_rng(seed), t_max=20.0, **kw)
        assert (r.n_final >= 0).all(), r.n_final


def test_tau_leap_keeps_integer_counts():
    r = tau_leap_run(_comp(), N0, np.random.default_rng(0), tau=0.05, t_max=10.0)
    assert np.allclose(r.n_final, np.round(r.n_final))


def test_cle_does_not_quantise():
    """Rounding CLE counts to integers would smuggle in a discreteness the CLE
    does not have, which is the very thing being measured."""
    r = cle_run(_comp(), N0, np.random.default_rng(0), dt=0.02, t_max=10.0)
    assert not np.allclose(r.n_final, np.round(r.n_final))


def test_small_tau_agrees_with_exact_ssa():
    """The calibration the hierarchy rests on: tau -> 0 recovers the SSA.

    Compares the splitting probability from a biased start. 3000 trials at
    p ~ 0.2 gives a standard error near 0.007 per arm, so 0.03 absolute is a
    ~3 sigma band on the difference.
    """
    net = am_reversible(0.3)
    comp = compile_network(net, 90.0)
    pair = reverse_pairing(net)
    thr = 40
    stop_i = lambda n: abs(int(n[0]) - int(n[1])) >= thr
    stop_f = lambda n: abs(n[0] - n[1]) >= thr

    rng = np.random.default_rng(1)
    wrong = 0
    for _ in range(3000):
        r = gillespie_instrumented(comp, N0, rng, pair, stop=stop_i,
                                   max_steps=2_000_000)
        wrong += int(r.n_final[0] <= r.n_final[1])
    p_ssa = wrong / 3000

    rng = np.random.default_rng(2)
    wrong = ok = 0
    for _ in range(3000):
        r = tau_leap_run(comp, N0, rng, tau=0.01, stop=stop_f, t_max=400.0)
        if r.hit_budget:
            continue
        ok += 1
        wrong += int(r.n_final[0] <= r.n_final[1])
    p_tau = wrong / ok

    assert ok > 2900, ok
    assert abs(p_tau - p_ssa) < 0.03, (p_ssa, p_tau)


def test_retries_are_reported_not_hidden():
    """A run that constantly retries is reporting a bad step, not a result."""
    r = tau_leap_run(_comp(), np.array([1, 1, 88]), np.random.default_rng(0),
                     tau=5.0, t_max=50.0)
    assert r.retries > 0
