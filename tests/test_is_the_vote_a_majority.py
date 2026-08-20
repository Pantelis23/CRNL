"""§111 — the pooled vote is a majority only if the rails are symmetric."""

import json
import math
import pathlib

import pytest

from experiments.depth_compounding import R1, R2, R3
from experiments.is_the_vote_a_majority import (
    MFPT_CEILING,
    hold_lifetime,
    m_of_k,
    merge_commit,
    merge_fraction,
    predicted_crossover,
    x_merged,
)

RESULTS = pathlib.Path("results/is_the_vote_a_majority.json")


@pytest.fixture(scope="module")
def data():
    return json.loads(RESULTS.read_text())


def test_merge_fraction_is_not_one_half():
    """AM's symmetric rails would give 0.5; Schlogl's give 0.72."""
    assert merge_fraction() == pytest.approx((R3 - R2) / (R3 - R1))
    assert merge_fraction() == pytest.approx(0.71972, abs=1e-4)
    assert merge_fraction() > 0.6, "asymmetry must be substantial, not marginal"


@pytest.mark.parametrize("k,m,maj", ((3, 3, 2), (5, 4, 3), (7, 6, 4)))
def test_threshold_is_not_the_majority(k, m, maj):
    assert m_of_k(k) == m
    assert m_of_k(k) != maj, "a majority merge would be the AM answer"


def test_k_three_requires_unanimity():
    """The sharp case: m = k, so all three tanks must fail."""
    assert m_of_k(3) == 3


def test_hold_lifetime_is_exponential_in_volume(data):
    """P1: ln T = cN + a, R2 ~ 1."""
    assert data["r2"] > 0.999
    assert data["c"] == pytest.approx(0.183476, abs=1e-4)


def test_mfpt_solve_refuses_to_return_a_negative_lifetime():
    """§111.2(a): the silent failure above N~150 is now an exception, not a number."""
    assert hold_lifetime(100) > 0
    with pytest.raises(FloatingPointError):
        hold_lifetime(220)


def test_every_crossover_volume_respects_the_ceiling(data):
    for xo in data["crossovers"]:
        for N in xo["volumes"]:
            assert xo["k"] * N <= MFPT_CEILING


def test_k_three_ratio_is_volume_independent_and_matches_the_closed_form(data):
    """§111.1: m = k makes both protocols grow at 3c, so the ratio is a constant."""
    xo = next(x for x in data["crossovers"] if x["k"] == 3)
    lo = xo["ln_ratio"]
    predicted = 2 * (math.log(0.15106815) - data["a"])
    assert predicted == pytest.approx(-3.4131, abs=5e-3)
    assert lo[0] == pytest.approx(predicted, rel=0.01), "0.31% at the smallest volume"
    assert all(v < 0 for v in lo), "re-merging wins at every volume tested"
    assert max(lo) - min(lo) < 0.5, "and the ratio is near-constant, not crossing"


def test_the_crossover_is_strongly_k_dependent(data):
    """P3: the opposite of §34's k-independence on AM."""
    preds = {x["k"]: x["predicted"] for x in data["crossovers"]}
    assert not math.isfinite(preds[3]), "k=3 has no crossover of this form"
    assert preds[7] / preds[5] > 1.5, "the finite ones differ by a large factor"


def test_absolute_crossover_test_at_k_five(data):
    """P4: predicted from the hold protocol alone, nothing fitted to a crossover."""
    xo = next(x for x in data["crossovers"] if x["k"] == 5)
    assert xo["measured"] == pytest.approx(17.904, abs=0.1)
    assert xo["predicted"] == pytest.approx(19.131, abs=0.1)
    assert abs(xo["ratio"] - 1) < 0.10


def test_merge_threshold_is_soft_near_the_saddle():
    """§111.2(b): a pool landing within a fluctuation of the saddle is a coin flip."""
    assert abs(x_merged(7, 5) - R2) < 0.05, "k=7 j=5 lands just above the saddle"
    p = merge_commit(7, 5, 14)
    assert 0.3 < p < 0.7, f"must be soft, not a step: {p}"
    # while a cell far from the saddle is sharp
    assert merge_commit(7, 3, 14) > 0.99
    assert merge_commit(7, 6, 14) < 0.01


def test_predicted_crossover_is_infinite_when_m_equals_k():
    assert not math.isfinite(predicted_crossover(0.18, -0.18, 3, 0.15))
    assert math.isfinite(predicted_crossover(0.18, -0.18, 5, 0.15))
