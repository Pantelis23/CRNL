"""§108 — the drift is two mechanisms, and one is provably not an averaging effect."""

import json
import pathlib

import numpy as np
import pytest

from experiments.two_mechanisms import (
    OMEGAS,
    candidate_averages,
    effective_rate,
    joint_absorbing,
    ln_k_curvature,
)

RESULTS = pathlib.Path("results/two_mechanisms.json")


@pytest.fixture(scope="module")
def data():
    return json.loads(RESULTS.read_text())


@pytest.mark.parametrize("om", (14, 30, 70))
def test_ln_k_is_convex_in_the_input(om):
    """The linchpin: convexity puts every average at or above k(<x>) by Jensen."""
    d2 = ln_k_curvature(om)
    assert (d2 > 0).all(), f"Omega={om}: {d2.min()}"


def test_curvature_grows_with_omega(data):
    """Why the averaging enhancement strengthens at large Omega."""
    c = data["ln_k_curvature"]
    lows = [c[str(om)][0] for om in OMEGAS]
    assert lows == sorted(lows), "minimum curvature must increase with Omega"


def test_candidate_averages_are_ordered(data):
    """k(<x>) <= exp(<ln k>) <= <k>, which is what makes §108's proof bite."""
    for r in data["cells"]:
        assert r["k_at_mean"] <= r["k_geometric"] <= r["k_arithmetic"]


def test_measured_rate_falls_below_every_average_at_small_omega(data):
    """§108.3: outside the reach of ANY averaging prescription, by an inequality."""
    cells = {r["omega"]: r for r in data["cells"]}
    assert cells[14]["true_over_mean"] == pytest.approx(0.7880, abs=5e-3)
    assert cells[14]["true_over_mean"] < 1.0
    assert cells[20]["true_over_mean"] < 1.0
    # and the same cells are below the geometric mean too, a fortiori
    assert cells[14]["lambda_true"] < cells[14]["k_geometric"]


def test_measured_rate_crosses_upward_with_omega(data):
    """Two mechanisms: suppressing at small Omega, enhancing at large."""
    ratios = [r["true_over_mean"] for r in data["cells"]]
    assert ratios == sorted(ratios), "must rise monotonically"
    assert min(ratios) < 1.0 < max(ratios), "and must cross 1"


def test_removing_return_trips_does_not_flatten_the_drift(data):
    """The pre-registered candidate (f), refuted."""
    cells = data["cells"]
    mus = [r["pure_absorbing"] for r in cells]
    assert all(m > 0 for m in mus)
    ratios = [r["true_over_mean"] for r in cells]
    assert max(ratios) / min(ratios) > 1.7, "the drift survives absorption"


def test_absorbing_chain_is_cheaper_than_the_free_one_at_stage_one():
    """Wiring: absorbing stage 1 can only reduce P(s2 low, s1 high) vs a free upstream."""
    from experiments.free_upstream_depth import channel_split, solve
    om, t = 30, 2.0
    pa = joint_absorbing(om, t)
    p, ref, dims, strides, walled = solve(om, 2, 0, t)
    _, pure_free, _ = channel_split(p, om, dims, strides, ref, walled)
    assert pa < pure_free, "removing return trips must lower pure"


def test_effective_rate_inversion_is_consistent():
    """A round trip through the two-state form."""
    om = 30
    pa = joint_absorbing(om)
    lam = effective_rate(om, pa)
    kmean, kgeo, kari = candidate_averages(om)
    assert kmean < lam < kari, "the Omega=30 cell sits inside the bracket"
    assert lam == pytest.approx(4.4897e-3, rel=2e-2)
