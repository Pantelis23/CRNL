"""§103 — the chain predicted from single-element quantities, no joint generator."""

import numpy as np
import pytest

from experiments.chain_without_a_joint_solve import (
    MEASURED,
    chain_operating_points,
    split_from,
)
from experiments.depth_compounding import R3
from experiments.jensen_shift import F

GATE = (0.5, 2.0)                    # §102's pre-registered factor of two


def test_the_rail_is_a_fixed_point_of_the_transfer_map():
    """§103.1's root cause: neutrality at the rail makes r3 a fixed point of F."""
    assert F(R3) == pytest.approx(R3, abs=1e-9)


def test_transfer_map_alone_cannot_degrade_a_chain():
    """Iterating <F> from near the rail converges TOWARD it -- the withdrawn behaviour."""
    x = 2.85
    seq = [x]
    for _ in range(4):
        x = F(x)
        seq.append(x)
    assert seq == sorted(seq), f"F alone must climb toward the rail: {seq}"
    assert seq[-1] < R3 + 1e-9


def test_operating_points_degrade_with_depth():
    """With the intrinsic term restored, the predicted chain falls away from the rail."""
    mus, _ = chain_operating_points(30, 3)
    assert mus == sorted(mus, reverse=True), f"must degrade: {mus}"
    assert mus[0] < R3, "even stage 1 sits below its deterministic rail"


@pytest.mark.parametrize("om,D", sorted(MEASURED))
def test_predicted_split_is_inside_the_registered_gate(om, D):
    mus, _ = chain_operating_points(om, D)
    _, cp, pp = split_from(om, mus, 2.0, legacy=True)   # §103's published gate
    ref = MEASURED[(om, D)]
    ratio = (cp / pp) / (ref["contam"] / ref["pure"])
    assert GATE[0] < ratio < GATE[1], f"Om={om} D={D}: {ratio}"


@pytest.mark.parametrize("om,D,want", ((30, 2, 1.3702), (30, 3, 1.1689),
                                       (14, 2, 1.1844), (14, 3, 1.0585)))
def test_predicted_split_reproduces_the_recorded_values(om, D, want):
    mus, _ = chain_operating_points(om, D)
    _, cp, pp = split_from(om, mus, 2.0, legacy=True)   # §103's published values
    ref = MEASURED[(om, D)]
    assert (cp / pp) / (ref["contam"] / ref["pure"]) == pytest.approx(want, abs=2e-3)


def test_operating_points_are_accurate_where_the_expansion_is_healthy():
    """P1/P2 held at Omega=30 and failed at Omega=14; both are recorded, not just the win."""
    mus30, _ = chain_operating_points(30, 2)
    for pred, meas in zip(mus30, MEASURED[(30, 2)]["mus"]):
        assert abs(pred / meas - 1) < 0.01
    mus14, _ = chain_operating_points(14, 3)
    errs = [pred / meas - 1 for pred, meas in zip(mus14, MEASURED[(14, 3)]["mus"])]
    assert max(errs) > 0.02, "Omega=14 must miss the 2% prediction -- P2 was half-refuted"
    assert all(e > 0 for e in errs), "and every miss is one-signed, predicting less damage"


def test_conclusion_survives_the_whole_measured_range_of_p_transmit():
    """The one empirical input cannot carry the result (§100.2 gives 0.7254-0.9830)."""
    for pt in (0.7254, 0.8860, 0.9376, 0.9830):
        for (om, D), ref in MEASURED.items():
            mus, _ = chain_operating_points(om, D)
            _, cp, pp = split_from(om, mus, 2.0, p_transmit=pt, legacy=True)
            ratio = (cp / pp) / (ref["contam"] / ref["pure"])
            assert GATE[0] < ratio < GATE[1], f"p_t={pt} Om={om} D={D}: {ratio}"


def test_prediction_side_builds_no_joint_generator():
    """The claim is structural: every predicted quantity is 1-D. Guard the state-space size.

    A joint generator at Omega=30, D=3 is 121**3 = 1,771,561 states; every object the
    prediction touches is at most a few hundred.
    """
    import experiments.chain_without_a_joint_solve as mod
    from experiments.margin_law import upstream_qsd
    xs, px = upstream_qsd(30)
    assert len(xs) < 200
    assert len(mod.lattice(30)) < 200
