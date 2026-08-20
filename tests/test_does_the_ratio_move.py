"""§110 — the position is anti-correlated with the timescale ratio."""

import json
import pathlib

import pytest

from experiments.does_the_ratio_move import (
    OMEGAS,
    downstream_crossing,
    upstream_relaxation,
)

RESULTS = pathlib.Path("results/does_the_ratio_move.json")


@pytest.fixture(scope="module")
def rows():
    return json.loads(RESULTS.read_text())


def test_crossing_time_rises_with_omega(rows):
    """P2: a sharper instanton lingers longer near the saddle."""
    xs = [r["tau_cross"] for r in rows]
    assert xs == sorted(xs)
    assert max(xs) / min(xs) > 1.8


def test_upstream_clock_is_not_quite_a_constant(rows):
    """P1 failed, mildly -- recorded rather than smoothed over."""
    ups = [r["tau_up"] for r in rows]
    assert max(ups) / min(ups) > 1.2, "the prediction of <20% variation was wrong"
    assert ups != sorted(ups), "and it is non-monotone"


def test_position_and_ratio_move_in_opposite_directions(rows):
    """P3: the refutation, and P4 -- every cell, not the endpoints."""
    rt = [r["ratio"] for r in rows]
    ps = [r["position"] for r in rows]
    assert rt == sorted(rt, reverse=True), "the ratio falls monotonically"
    assert ps == sorted(ps), "the position rises monotonically"
    # anti-correlated in every adjacent pair, not just end to end
    for i in range(len(rows) - 1):
        assert (rt[i + 1] - rt[i]) * (ps[i + 1] - ps[i]) < 0, f"cell {i}"


def test_position_traverses_and_exits_the_bracket(rows):
    ps = [r["position"] for r in rows]
    assert ps[0] < 0.05, "starts at the fast limit"
    assert ps[-1] > 1.0, "and ends past the frozen limit"


def test_the_falling_ratio_is_driven_by_the_crossing_time(rows):
    """It is tau_cross rising, not tau_up moving, that turns the ratio over."""
    first, last = rows[0], rows[-1]
    assert last["tau_cross"] / first["tau_cross"] > 1.8
    assert 0.8 < last["tau_up"] / first["tau_up"] < 1.4


@pytest.mark.parametrize("om", (14, 70))
def test_clocks_recompute_live(om):
    """Guard the stored numbers against the functions that produced them."""
    stored = {r["omega"]: r for r in json.loads(RESULTS.read_text())}[om]
    tau_up, l1, l2 = upstream_relaxation(om)
    tau_x, _, _ = downstream_crossing(om)
    assert tau_up == pytest.approx(stored["tau_up"], rel=1e-6)
    assert tau_x == pytest.approx(stored["tau_cross"], rel=1e-6)
    assert l1 < l2, "the escape must be the slowest mode"
