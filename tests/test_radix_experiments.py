from __future__ import annotations

import numpy as np
import pytest


def test_champion_counts_sum_and_anchor():
    from experiments.radix_wall import champion_counts
    # sums to omega, has a blank slot of 0, champion is index 0
    c = champion_counts(4, 100, 0.10)
    assert c.sum() == 100
    assert c[-1] == 0                      # blank starts empty
    assert c[0] > c[1:-1].max()            # champion STRICTLY leads over rivals


def test_champion_strictly_leads_when_delta_omega_small():
    # regression (fable-5): delta*omega small used to round to a tie [3,3,3,2].
    from experiments.radix_wall import champion_counts
    for n, omega, delta in [(4, 11, 0.05), (5, 10, 0.05), (3, 11, 0.05)]:
        c = champion_counts(n, omega, delta)
        assert c.sum() == omega
        assert c[0] > c[1:-1].max(), (n, omega, delta, c)  # strict lead over rivals


def test_champion_counts_reduces_to_55_45_at_n2():
    from experiments.radix_wall import champion_counts
    c = champion_counts(2, 100, 0.10)
    # fixed pairwise margin delta=0.10 -> champion 55, rival 45, blank 0
    assert list(c) == [55, 45, 0]


def test_champion_pairwise_margin_is_delta():
    from experiments.radix_wall import champion_counts
    # champion fraction minus each rival fraction == delta (to rounding)
    n, omega, delta = 5, 10_000, 0.10
    c = champion_counts(n, omega, delta)
    f = c[:-1] / omega
    lead = f[0] - f[1]
    assert lead == pytest.approx(delta, abs=2e-3)


def test_fit_c_recovers_known_slope():
    from experiments.radix_wall import fit_c
    # synthetic: champ_loss = exp(-c*Omega), c=0.02, plenty of events
    c_true = 0.02
    omegas = np.array([20, 40, 60, 80, 100, 120])
    trials = 100_000
    loss = np.exp(-c_true * omegas)
    champ_loss_counts = (loss * trials).astype(int)
    fit = fit_c(omegas, champ_loss_counts, trials)
    assert fit is not None
    assert fit["c"] == pytest.approx(c_true, rel=0.05)
