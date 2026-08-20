"""§102 — the contaminated channel against escape rates at the operating points."""

import numpy as np
import pytest

from experiments.depth_compounding import R3
from experiments.escape_accounts_for_it import (
    P_TRANSMIT,
    escape_rate,
    operating_points,
    predict,
    rate_limits,
)
from experiments.what_reflection_costs import spectral_gap


@pytest.mark.parametrize("om", (14, 30))
def test_pinned_at_the_rail_is_stage_one(om):
    """h(r3) = 1, so a pinned downstream stage IS stage 1 -- exactly, not approximately."""
    k = escape_rate(om, R3)
    g, _ = spectral_gap(om, False)
    assert k == pytest.approx(g, rel=1e-12)


@pytest.mark.parametrize("om", (14, 30))
def test_escape_rate_rises_steeply_as_the_operating_point_degrades(om):
    xs = [R3, 3.0, 2.8, 2.6, 2.4, 2.2, 2.0]
    ks = [escape_rate(om, x) for x in xs]
    assert ks == sorted(ks), "must be monotone in the operating point"
    assert ks[-1] / ks[0] > 2.0
    # steeper at the deeper barrier -- §102.1 needs this asymmetry
    if om == 30:
        assert ks[-1] / ks[0] > 5.0


@pytest.mark.parametrize("om,D,want", ((14, 2, 1.159), (14, 3, 1.040), (30, 2, 1.359)))
def test_suspect_predicts_the_split_within_the_registered_gate(om, D, want):
    """P3's gate was a factor of two, pre-registered and deliberately loose."""
    mus, tot, pure, contam = operating_points(om, D, 2.0)
    _, cp, pp = predict(om, mus, 2.0, legacy=True)   # §102's published values
    ratio = (cp / pp) / (contam / pure)
    assert ratio == pytest.approx(want, abs=0.02)
    assert 0.5 < ratio < 2.0, "the pre-registered gate"


def test_operating_points_degrade_down_the_chain():
    mus, *_ = operating_points(14, 3, 2.0)
    assert mus == sorted(mus, reverse=True), f"each stage should sit lower: {mus}"


@pytest.mark.parametrize("om,D", ((14, 2), (14, 3), (30, 2)))
def test_effective_rate_is_bracketed_by_the_two_averaging_limits(om, D):
    """§102.1: k(<x>) <= k_eff <= <k(x)>, and the position sits near the fast end."""
    k_mean, k_avg, k_eff, pos = rate_limits(om, D, 2.0)
    assert k_mean < k_eff < k_avg, "the measurement must lie between the limits"
    assert 0.0 < pos < 0.5, f"expected the fast end, got {pos}"


def test_the_two_limits_separate_more_at_the_deeper_barrier():
    """Why P3's residual is worst at Omega=30: the rate curve is steeper there."""
    sep = {}
    for om in (14, 30):
        k_mean, k_avg, _, _ = rate_limits(om, 2, 2.0)
        sep[om] = k_avg / k_mean
    assert sep[30] > sep[14]
    assert sep[14] == pytest.approx(2.435, abs=0.02)
    assert sep[30] == pytest.approx(4.759, abs=0.02)


def test_p4_as_written_is_ill_posed_across_depths():
    """The model's own formula sums over upstream stages, so contam/pure depends on D.

    Guard the reasoning, not just the number: k1/k_last is near-constant across the three
    cells while the measured ratio varies by 3x, so P4 read literally would kill a suspect
    that P3 confirms.
    """
    k_ratios, measured = [], []
    for om, D in ((14, 2), (14, 3), (30, 2)):
        mus, tot, pure, contam = operating_points(om, D, 2.0)
        ks, _, _ = predict(om, mus, 2.0, legacy=True)   # §102's P4 as written
        k_ratios.append(ks[0] / ks[-1])
        measured.append(contam / pure)
    assert max(k_ratios) / min(k_ratios) < 1.10, "k1/k_last barely moves"
    assert max(measured) / min(measured) > 3.0, "the measured ratio moves a lot"
