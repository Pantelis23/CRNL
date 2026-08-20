"""§107 — the closure with §106's corrections folded in."""

import json
import pathlib

import numpy as np
import pytest

from experiments.depth_compounding import R3
from experiments.the_corrected_closure import closure, geometric_rate, pi_low
from experiments.what_reflection_costs import spectral_gap

RESULTS = pathlib.Path("results/the_corrected_closure.json")
ALL_ON = dict(indexing=True, two_state=True, geometric=True, derived_pt=True)


@pytest.fixture(scope="module")
def rows():
    return json.loads(RESULTS.read_text())


@pytest.mark.parametrize("om", (14, 30, 70))
def test_stage_one_rate_is_the_gap_and_averaging_leaves_it_alone(om):
    """P1: stage 1 has no upstream, so correction (3) must not touch it."""
    g, _ = spectral_gap(om, False)
    for kw in ({"indexing": True}, {"indexing": True, "geometric": True}):
        _, ks, _ = closure(om, **kw)
        assert ks[0] == pytest.approx(g, rel=1e-9)


@pytest.mark.parametrize("om", (14, 30, 55))
def test_geometric_mean_exceeds_the_rate_at_the_mean(om):
    """ln k is convex in the input, so the action average beats the rate at the mean."""
    from experiments.escape_accounts_for_it import escape_rate
    from experiments.chain_without_a_joint_solve import chain_operating_points
    mus, _ = chain_operating_points(om, 2)
    assert geometric_rate(om, mus[0]) > escape_rate(om, mus[0])


@pytest.mark.parametrize("om", (14, 30, 70))
def test_pi_low_is_below_one_and_falls_with_omega(om):
    p = pi_low(om, R3)
    assert 0.0 < p < 1.0


def test_pi_low_falls_with_omega():
    ps = [pi_low(om, R3) for om in (14, 30, 70)]
    assert ps == sorted(ps, reverse=True)
    assert ps[0] / ps[-1] > 1.5


def test_ablation_credits_the_geometric_correction(rows):
    """P3: geometric carries the improvement; two-state is nearly invisible in the ratio."""
    r30 = next(r for r in rows if r["omega"] == 30)["ratios"]
    unc, idx, two, geo, pt = r30
    assert abs(two - idx) < 0.05, "two-state must nearly cancel in the ratio"
    assert (idx - geo) > 4 * (idx - two), "geometric must carry most of the shift"
    assert geo == pytest.approx(1.0138, abs=5e-3)


def test_corrections_recentre_the_model(rows):
    """The corrected model crosses 1 instead of sitting above it everywhere."""
    unc = [r["ratios"][0] for r in rows]
    cor = [r["ratios"][-1] for r in rows]
    assert all(v > 1.0 for v in unc), "uncorrected is above 1 everywhere"
    assert min(cor) < 1.0 < max(cor), "corrected must cross 1"


def test_the_drift_is_not_removed(rows):
    """P2 refuted on the half that matters: still monotone, span no better."""
    cor = [r["ratios"][-1] for r in rows]
    unc = [r["ratios"][0] for r in rows]
    assert cor == sorted(cor), "still monotone in Omega"
    assert max(cor) / min(cor) >= max(unc) / min(unc), "the span did not improve"


def test_every_cell_still_inside_the_registered_gate(rows):
    """What survives: §102's factor-of-two gate holds across the whole swept range."""
    for r in rows:
        assert 0.5 < r["ratios"][-1] < 2.0, f"Omega={r['omega']}: {r['ratios'][-1]}"
