"""§113 — the box audit: what the threading must not change, and what the box does and does not move."""

import json
import math
import pathlib

import numpy as np
import pytest

from experiments.what_the_box_was_doing import _converged


# ------------------------------------------------------------------ the threading must be inert

def test_spectral_gap_default_reproduces_the_published_reference():
    """§102's P1 pins this to solver tolerance; threading cap_mult must not touch it."""
    from experiments.what_reflection_costs import spectral_gap

    assert spectral_gap(14, False)[0] == pytest.approx(7.29805e-02, rel=1e-6)
    assert spectral_gap(30, False)[0] == pytest.approx(3.80268e-03, rel=1e-6)
    # and the explicit default must equal the implicit one
    assert spectral_gap(14, False)[0] == spectral_gap(14, False, cap_mult=1.25)[0]


def test_raising_the_dense_threshold_changed_no_published_cell():
    """§113 raised spectral_gap's dense/ARPACK switch from 400 to 1500 states. Every published
    cell runs Omega <= 55 at cap_mult 1.25, i.e. <= 220 states -- already dense either way."""
    from experiments.what_reflection_costs import stage1_generator

    for om in (14, 30, 55):
        for refl in (False, True):
            Q, _ = stage1_generator(om, refl, cap_mult=1.25)
            assert Q.shape[0] < 400, \
                f"Omega={om} refl={refl}: {Q.shape[0]} states was NOT already dense"


def test_rate_limits_default_reproduces_section_102():
    from experiments.escape_accounts_for_it import rate_limits

    k_mean, k_avg, k_eff, pos = rate_limits(14, 2, 2.0)
    assert k_avg / k_mean == pytest.approx(2.435, abs=0.02)
    assert 0.0 < pos < 0.5
    assert (k_mean, k_avg, k_eff, pos) == rate_limits(14, 2, 2.0, cap_mult=1.25)


def test_the_two_clocks_default_reproduces_section_110s_stored_table():
    from experiments.does_the_ratio_move import downstream_crossing, upstream_relaxation

    src = json.loads(pathlib.Path("results/does_the_ratio_move.json").read_text())
    for r in src:
        om = r["omega"]
        assert upstream_relaxation(om)[0] == pytest.approx(r["tau_up"], rel=1e-9)
        assert downstream_crossing(om)[0] == pytest.approx(r["tau_cross"], rel=1e-9)


def test_solve_threading_is_inert():
    """free_upstream_depth.solve gained cap_mult; the default must be bit-identical."""
    from experiments.free_upstream_depth import solve

    a = solve(14, 2, 0, 2.0)[0]
    b = solve(14, 2, 0, 2.0, cap_mult=1.25)[0]
    assert np.array_equal(a, b)


# ------------------------------------------------- what the box provably cannot touch

def test_tau_cross_ignores_the_box_by_construction_not_by_measurement():
    """Its h-transform runs on [saddle, rail], both ABSORBING, so the wall is outside the
    domain and `downstream_crossing` never reads cap_mult at all.

    §113's first draft SWEPT this argument and reported the five identical results as a 0.00%
    convergence measurement. That is a number the harness produced and the chemistry did not
    (rule 10). The structural fact is asserted here instead; a sweep would prove nothing, and
    re-running the same linear solve five times does not even give bit-identical answers,
    because BLAS threading moves the last ulp -- which is how the vacuous sweep was caught.
    """
    import inspect

    from experiments.does_the_ratio_move import downstream_crossing

    body = inspect.getsource(downstream_crossing).split('"""')[2]
    assert body.count("cap_mult") == 0, "it now uses the box; §113's P5 must actually sweep it"
    for om in (14, 30, 55):
        vals = [downstream_crossing(om, cap_mult=c)[0] for c in (1.25, 4.0)]
        assert vals[0] == pytest.approx(vals[1], rel=1e-12)


@pytest.mark.parametrize("om", (14, 30, 55))
def test_descent_rate_is_box_independent_but_only_measured(om):
    """Not provable the same way: its domain DOES reach the wall. It is flat because the
    upstream is pinned at the LOW rail, which removes the downstream's high state entirely,
    so upward excursions are exponentially unlikely. Measured, not derived."""
    from experiments.predicting_transmission import descent_rate

    vals = [descent_rate(om, cap_mult=c)[0] for c in (1.25, 2.0, 4.0)]
    assert max(vals) / min(vals) - 1 < 1e-9, vals


# ------------------------------------------------- what the box does move, and in which direction

@pytest.mark.parametrize("om,worst", ((14, 0.02), (30, 0.005), (55, 0.001)))
def test_spectral_gap_carries_a_one_signed_box_error_that_falls_with_volume(om, worst):
    """A wider box can only LENGTHEN an escape, so it can only LOWER the rate."""
    from experiments.what_reflection_costs import spectral_gap

    vals = [spectral_gap(om, False, cap_mult=c)[0] for c in (1.25, 1.6, 2.0, 3.0, 4.0)]
    assert vals[0] > vals[-1], "the converged gap must be the smaller one"
    # One-signed only while the differences are above the solver's noise floor. Demanding a
    # strictly decreasing sequence all the way to cap_mult 4.0 fails on eigenvalue noise at the
    # 1e-11 level -- a fixed gate on a converged quantity, which is what rule 20 forbids.
    head = abs(vals[1] - vals[0])
    for i in range(len(vals) - 1):
        step = vals[i] - vals[i + 1]
        assert step > 0 or abs(step) < 1e-3 * head, f"not one-signed above the floor: {vals}"
    assert abs(vals[0] / vals[-1] - 1) > worst, "this Omega should still show a visible error"


def test_the_box_error_shrinks_as_volume_grows():
    from experiments.what_reflection_costs import spectral_gap

    errs = []
    for om in (14, 30, 55):
        v = [spectral_gap(om, False, cap_mult=c)[0] for c in (1.25, 4.0)]
        errs.append(abs(v[0] / v[1] - 1))
    assert errs == sorted(errs, reverse=True), errs


# ------------------------------------------------- the published conclusions

def test_section_110s_refutation_survives_every_box():
    """Rule 14: §110 published a WITHDRAWAL, so it gets verified as carefully as an assertion.

    Its claim is that tau_up/tau_cross FALLS across Omega while the position RISES. If a box
    flipped that, the refutation would be an artifact.
    """
    from experiments.does_the_ratio_move import downstream_crossing, upstream_relaxation

    oms = (14, 20, 30, 55, 70)
    for c in (1.25, 2.0, 4.0):
        series = [upstream_relaxation(om, cap_mult=c)[0] / downstream_crossing(om, cap_mult=c)[0]
                  for om in oms]
        assert series == sorted(series, reverse=True), f"cap={c}: {series}"
        assert series[-1] < series[0]


def test_section_102s_position_survives_the_box():
    """§102 published 0 < pos < 0.5 ("near the fast end") and a test asserts it."""
    from experiments.escape_accounts_for_it import rate_limits

    pos = [rate_limits(14, 2, 2.0, cap_mult=c)[3] for c in (1.25, 2.0, 4.0)]
    assert all(0.0 < p < 0.5 for p in pos), pos
    assert abs(pos[0] - pos[-1]) < 0.01, f"moved by {abs(pos[0]-pos[-1])}"


def test_the_bracket_width_moves_more_than_predicted():
    """P4 predicted the width would cancel to under 2% at Omega=14. It does not -- 3.07%.

    Pinned so the miss stays visible: the one-signed bias cancels in the ratio by a factor of
    about 2, not to nothing, and the published conclusion survives for a different reason (the
    position's numerator moves with its denominator).
    """
    from experiments.escape_accounts_for_it import rate_limits

    w = [rate_limits(14, 2, 2.0, cap_mult=c)[1] / rate_limits(14, 2, 2.0, cap_mult=c)[0]
         for c in (1.25, 4.0)]
    moved = abs(w[0] / w[1] - 1) * 100
    assert moved > 2.0, f"the prediction would have held after all: {moved:.2f}%"
    assert moved < 5.0, f"but not by much: {moved:.2f}%"


# ------------------------------------------------- the convergence criterion itself

def test_convergence_criterion_accepts_a_settled_sequence_with_noise_at_the_floor():
    """The first version demanded strictly-shrinking differences and rejected this."""
    ok, _ = _converged([1.0, 0.9, 0.8999, 0.89990000001, 0.8999])
    assert ok


def test_convergence_criterion_rejects_a_sequence_still_moving():
    ok, _ = _converged([1.0, 0.9, 0.8, 0.7, 0.6])
    assert not ok


def test_convergence_criterion_needs_three_points():
    ok, _ = _converged([1.0, 0.9])
    assert not ok
    ok, _ = _converged([1.0, 0.9, float("nan"), float("nan"), float("nan")])
    assert not ok


def test_convergence_criterion_is_scale_free():
    """It compares the tail step to the head step, so multiplying everything changes nothing.

    That is what keeps it from being the fixed tolerance rule 20 forbids.
    """
    seq = [1.0, 0.5, 0.4, 0.399, 0.3990001]
    for scale in (1e-9, 1.0, 1e9):
        assert _converged([v * scale for v in seq])[0] == _converged(seq)[0]
