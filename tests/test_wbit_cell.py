"""§125 -- W-bit's analog cell as a CRNL element: faithful, exact, and Arrhenius in 1/sigma^2."""

import json
import pathlib

import pytest

from experiments.wbit_cell import WBIT_SRC, fixed_points, flip_rate

STORE = pathlib.Path("results/wbit_cell.json")


def _d():
    return json.loads(STORE.read_text())


def test_reduced_map_is_wbit_itself():
    """W1: one-step law matches W-bit's own AnalogWbitNetwork within 3 standard errors."""
    for row in _d()["W1"]:
        assert abs(row["dmean_over_se"]) < 3 and abs(row["dsd_rel"]) < 0.03


def test_exact_rate_converged_and_matches_monte_carlo():
    w = _d()["W2"]
    r = [w["rates_M"][m] for m in ("400", "800", "1600")]
    assert abs(r[2] / r[1] - 1) < abs(r[1] / r[0] - 1) < 1e-3
    assert abs(w["mc"]["rate"] - r[2]) < 2 * w["mc"]["se"]


def test_arrhenius_slope_approaches_the_minimum_action():
    """W3: the local slope rises toward S (rule 20: approach, not a fixed band)."""
    rows = _d()["W3"]["rows"]
    tail = [r["local_slope_over_S"] for r in rows[-5:]]
    assert tail == sorted(tail) and 0.99 < tail[-1] < 1.0


def test_rate_reproduces_from_code():
    assert abs(flip_rate(800, 1.5, 0.5, 0.5, 0.4) / _d()["W2"]["rates_M"]["800"] - 1) < 1e-12
    assert len(fixed_points(1.5, 0.5, 0.5)) == 3


@pytest.mark.skipif(not WBIT_SRC.exists(), reason="sibling W-bit repo not present")
def test_wbit_class_still_matches():
    """The seam: if W-bit's step() changes, this fails (one draw-set at the high well)."""
    import math

    import numpy as np
    from scipy.special import expit

    from experiments.wbit_cell import wbit_one_step
    p0 = fixed_points(1.5, 0.5, 0.5)[-1]
    w = wbit_one_step(1.5, 0.5, 0.5, 0.4, p0, 5000, seed=11)
    r = 0.5 * p0 + 0.5 * expit((1.5 * (2 * p0 - 1) + np.random.default_rng(2).normal(0, 0.4 * math.sqrt(2), 5000)) / 0.5)
    assert abs(w.mean() - r.mean()) < 4 * math.hypot(w.std(), r.std()) / math.sqrt(5000)


def test_min_action_is_finite_and_reproduces():
    """The stored S once came out inf (a convergence check that skipped the first update)."""
    from experiments.wbit_cell import min_action
    S = min_action(1.5, 0.5, 0.5)
    assert 0.3 < S < 0.36 and abs(S / _d()["W3"]["S"] - 1) < 1e-12
