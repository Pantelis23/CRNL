"""§109 — the seed carried §108's effect; both mismatches fixed."""

import json
import pathlib

import numpy as np
import pytest

from experiments.chain_without_a_joint_solve import chain_operating_points
from experiments.depth_compounding import R3
from experiments.escape_accounts_for_it import escape_rate
from experiments.it_was_the_seed import (
    OMEGAS,
    joint_absorbing_seeded,
    one_stage_low,
)
from experiments.the_corrected_closure import pi_low
from experiments.what_reflection_costs import spectral_gap

RESULTS = pathlib.Path("results/it_was_the_seed.json")


@pytest.fixture(scope="module")
def data():
    return json.loads(RESULTS.read_text())


def test_delta_seed_carries_a_huge_transient_with_no_averaging_available(data):
    """§109.1: a SINGLE stage, upstream pinned -- so this cannot be an averaging effect."""
    rows = data["p1"]["rows"]
    d = [r["delta"] for r in rows]
    q = [r["qsd"] for r in rows]
    assert max(d) / min(d) > 20.0, "the delta seed's transient must be enormous"
    assert max(q) / min(q) < 1.2, "the QSD seed must be nearly flat"
    assert d == sorted(d), "and the delta transient must rise monotonically"


def test_qsd_seed_converges_to_the_two_state_forward_rate(data):
    """P1: P(low)/t -> pi_low * lambda."""
    target = data["p1"]["target"]
    last = data["p1"]["rows"][-1]["qsd"]
    assert last / target > 0.95
    om = 30
    mus, _ = chain_operating_points(om, 2)
    assert target == pytest.approx(pi_low(om, mus[0]) * escape_rate(om, mus[0]), rel=1e-9)


def test_matched_seeding_closes_t_casc_z(data):
    """§109.2: the flat 21.5% was the reflected-vs-QSD seed mismatch."""
    p2 = data["p2"]
    refl = [r["refl_over_qsd"] for r in p2]
    ratios = [r["qsd_over_model"] for r in p2]
    assert all(1.2 < v < 1.3 for v in refl), "the reflected law escapes ~25% faster"
    assert max(refl) / min(refl) < 1.05, "and does so flatly in Omega"
    assert max(ratios) / min(ratios) < 1.03, "matched seeding leaves a FLAT residual"
    assert all(0.93 < v < 1.0 for v in ratios), "of only 3-5%"


def test_section_108_premise_does_not_survive_matched_seeding(data):
    """The withdrawal: no cell is materially below k(<x>) once the seed is matched."""
    frac = ([r for r in data["p1"]["rows"] if r["t"] == 2.0][0]["qsd"]
            / data["p1"]["target"])
    assert frac < 0.96, "the QSD instrument itself under-reads at t=2"
    corrected = [r["qsd"] / frac for r in data["p3"]]
    assert all(v > 1.0 for v in corrected), f"none may be below k(<x>): {corrected}"
    # and the delta-seeded version -- §108's -- genuinely was below
    assert min(r["delta"] for r in data["p3"]) < 0.8


def test_the_drift_survives_both_seed_fixes(data):
    """P4: fixing a seed cannot remove a trend across five barrier depths."""
    q = [r["qsd"] for r in data["p3"]]
    assert q == sorted(q), "still monotone"
    assert max(q) / min(q) > 1.9, "and the span is undiminished"


def test_system_traverses_and_exits_the_bracket(data):
    """§109.3: from the fast limit at Omega=14 to past the frozen limit at Omega=70."""
    cells = {r["omega"]: r for r in data["p3"]}
    assert cells[14]["qsd"] < 1.1, "at the fast limit"
    assert cells[70]["qsd"] > cells[70]["bracket_top"], "past the frozen limit"


@pytest.mark.parametrize("om", (14, 70))
def test_reflected_seed_escapes_faster_than_qsd(om):
    """Computed live: the reflected law is not depleted near the saddle."""
    r = one_stage_low(om, 0.0, 2.0, "reflected", first=True)
    q = one_stage_low(om, 0.0, 2.0, "qsd", first=True)
    assert r > q


def test_transient_fraction_is_omega_dependent():
    """§109's first correction extrapolated one fraction from Omega=30. It is not constant."""
    from experiments.escape_accounts_for_it import escape_rate
    fracs = {}
    for om in (14, 30):
        mus, _ = chain_operating_points(om, 2)
        x = mus[0]
        target = pi_low(om, x) * escape_rate(om, x)
        fracs[om] = (one_stage_low(om, x, 2.0, "qsd") / 2.0) / target
    assert fracs[14] == pytest.approx(0.8787, abs=5e-3)
    assert fracs[30] == pytest.approx(0.9471, abs=5e-3)
    assert fracs[30] / fracs[14] > 1.06, "the shortfall at Omega=14 is more than twice Omega=30's"


def test_correcting_per_cell_makes_the_withdrawal_unmarginal(data):
    """With the per-Omega fraction the Omega=14 cell is 12.6% above k(<x>), not 4.5%."""
    from experiments.escape_accounts_for_it import escape_rate
    cells = {r["omega"]: r for r in data["p3"]}
    om = 14
    mus, _ = chain_operating_points(om, 2)
    x = mus[0]
    frac = (one_stage_low(om, x, 2.0, "qsd") / 2.0) / (pi_low(om, x) * escape_rate(om, x))
    corrected = cells[om]["qsd"] / frac
    assert corrected > 1.10, f"comfortably above k(<x>), not marginal: {corrected}"


def test_the_head_start_algebra_against_quadrature():
    """§109.4 Class C: verify §104's error is arithmetic, independently of my derivation."""
    from scipy.integrate import quad
    from experiments.the_head_start import p_headstart, p_longer_window
    k, D, t = 5.5169, 0.0686, 2.0
    num = quad(lambda s: 1 - np.exp(-k * (t - s + D)), 0, t)[0] / t
    assert p_headstart(k, t, D) == pytest.approx(num, abs=1e-12)
    assert abs(p_longer_window(k, t, D) - num) > 1e-2
