"""§105 — the corrected head-start algebra, and the derived head start."""

import numpy as np
import pytest

from experiments.predicting_transmission import MEASURED, descent_rate
from experiments.the_head_start import (
    bistability_edge,
    conditional_traversal,
    p_headstart,
    p_longer_window,
    rate_weighted_head_start,
)

K, _ = descent_rate(30)


def test_the_two_forms_agree_only_at_zero_delay():
    """§105's whole point: a head start is not a longer window."""
    assert p_headstart(K, 2.0, 0.0) == pytest.approx(p_longer_window(K, 2.0, 0.0), abs=1e-12)
    for d in (0.05, 0.2, 1.0):
        assert abs(p_headstart(K, 2.0, d) - p_longer_window(K, 2.0, d)) > 1e-3


def test_required_delay_is_near_constant_under_the_correct_form():
    """Under §104's form it spans 18.9x; under the correct one, 2.2x."""
    from scipy.optimize import brentq
    dc = [brentq(lambda D: p_headstart(K, t, D) - m, 0.0, 50.0)
          for t, m in sorted(MEASURED.items())]
    dw = [brentq(lambda D: p_longer_window(K, t, D) - m, 0.0, 500.0)
          for t, m in sorted(MEASURED.items())]
    assert max(dc) / min(dc) < 3.0
    assert max(dw) / min(dw) > 10.0
    # and §104's version grows monotonically with the window -- which is what misled it
    assert dw == sorted(dw)


def test_downstream_loses_its_rail_above_the_upstream_saddle():
    """The premise, and it is a root count rather than a fit."""
    x = bistability_edge()
    assert x == pytest.approx(1.5795, abs=1e-3)
    assert x > 1.0, "stage 2 must start sliding before stage 1 formally crosses"


def test_naive_head_start_is_refuted_and_inconsistent():
    """§105.1: raw traversal is 9x too big AND longer than stage 2's whole descent."""
    tau, hprob = conditional_traversal(30, bistability_edge())
    assert tau == pytest.approx(0.4464, abs=2e-3)
    assert 0.0 < hprob < 1.0
    assert tau / 0.0486 > 5.0, "P2 fails"
    assert tau > 1.0 / K, "P3 fails: longer than the descent it is supposed to precede"
    # and it would put p_transmit near 1 where the measurement says 0.7254
    assert p_headstart(K, 0.5, tau) > 0.95
    assert MEASURED[0.5] < 0.75


def test_critical_slowing_at_the_saddle_node():
    """§105.2's mechanism: the descent is ~9x slower at the edge than at the low rail."""
    from experiments.escape_accounts_for_it import escape_rate
    edge = escape_rate(30, bistability_edge())
    low = escape_rate(30, 0.15)
    assert low / edge > 5.0


def test_rate_weighted_head_start_is_derived_not_fitted():
    """§105.2: P4, inside its pre-registered factor-of-two gate."""
    d_eff = rate_weighted_head_start(30, bistability_edge())
    assert d_eff == pytest.approx(0.0686, abs=3e-3)
    ratio = d_eff / 0.0486
    assert 0.5 < ratio < 2.0
    # and it is far below the raw traversal it corrects
    tau, _ = conditional_traversal(30, bistability_edge())
    assert d_eff < tau / 5.0


def test_derived_delay_beats_zero_delay_with_nothing_fitted():
    """The parameter-free curve must improve, and the one-signed bias must go."""
    d = rate_weighted_head_start(30, bistability_edge())
    res0 = {t: (p_headstart(K, t, 0.0) - m) / m for t, m in MEASURED.items()}
    resd = {t: (p_headstart(K, t, d) - m) / m for t, m in MEASURED.items()}
    assert all(r < 0 for r in res0.values()), "Delta=0 is one-signed low"
    assert max(abs(r) for r in resd.values()) < max(abs(r) for r in res0.values())
    assert not all(r < 0 for r in resd.values()), "the one-signed bias must be gone"
    assert abs(resd[2.0]) < 0.01


def test_chain_closure_survives_every_variant():
    """§103's gate holds under measured, Delta=0 and derived-Delta p_transmit."""
    from experiments.chain_without_a_joint_solve import (
        MEASURED as CH, chain_operating_points, split_from,
    )
    d = rate_weighted_head_start(30, bistability_edge())
    for pt in (0.9376, p_headstart(K, 2.0, 0.0), p_headstart(K, 2.0, d)):
        for (om, D), ref in CH.items():
            mus, _ = chain_operating_points(om, D)
            _, c, p = split_from(om, mus, 2.0, p_transmit=pt, legacy=True)
            ratio = (c / p) / (ref["contam"] / ref["pure"])
            assert 0.5 < ratio < 2.0, f"p_t={pt} Om={om} D={D}: {ratio}"
