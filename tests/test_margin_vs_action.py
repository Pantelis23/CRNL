"""§99.1 — the post-hoc discrimination between margin/sigma and A*Omega."""

import numpy as np
import pytest

from experiments import margin_vs_action as mva


@pytest.fixture(scope="module")
def data():
    return mva.load()


def test_ratios_differ_across_elements():
    """If the two axes were one axis this test would be void; they are 22% apart."""
    r_cal, r_new = mva.ratios()
    assert r_cal == pytest.approx(1.7931, abs=1e-3)
    assert r_new == pytest.approx(2.1945, abs=1e-3)
    assert abs(r_new / r_cal - 1) > 0.02


def test_both_candidates_are_inside_the_traced_range(data):
    """§59/§98.1: an extrapolated np.interp silently returns an endpoint. Guard it."""
    for r in mva.discriminate(data):
        assert r["inside"], f"{r['variable']} left the traced range"


def test_action_transfers_better(data):
    res = {r["variable"]: r for r in mva.discriminate(data)}
    assert res["A*Omega"]["predicted"] == pytest.approx(0.6383, abs=2e-3)
    assert res["(margin/sigma)^2"]["predicted"] == pytest.approx(0.5808, abs=2e-3)
    assert abs(res["A*Omega"]["rel_err"]) == pytest.approx(0.066, abs=3e-3)
    assert abs(res["(margin/sigma)^2"]["rel_err"]) == pytest.approx(0.150, abs=3e-3)
    assert abs(res["A*Omega"]["rel_err"]) < abs(res["(margin/sigma)^2"]["rel_err"])


def test_both_predictions_undershoot(data):
    """Both miss low, so the sign carries no information about which is right."""
    assert all(r["rel_err"] < 0 for r in mva.discriminate(data))


def test_x_increasing_assertion_is_live(data):
    """§98.1's bug was decreasing x passed to np.interp. The guard must actually fire."""
    rows = sorted(data["sweep"], key=lambda r: r["AOm"])
    aoms = np.array([r["AOm"] for r in rows], float)[::-1]
    with pytest.raises(AssertionError):
        assert np.all(np.diff(aoms) > 0), "A*Omega must increase for np.interp (§98.1)"
