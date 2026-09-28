"""§119 — a weaker upstream: the downstream's error floor is slaved to the upstream's dip rate."""

import json
import pathlib

from experiments.weak_upstream import ARMS, edge, element

STORE = pathlib.Path("results/weak_upstream.json")


def _d():
    return json.loads(STORE.read_text())


def _last(key):
    return _d()[key]["series"][-1]


def _floor(key):
    v = _d()[key]
    s = v["series"][-1]
    return s["p_free"] / (v["k1"] * s["t"])


def test_zone_widths_from_the_cubic():
    w = {k: edge(element(b)) - element(a).r2 for k, (a, b) in ARMS.items()}
    assert w["N"] < 0 < w["c"] < w["Z"] < w["b"] < w["X"] < w["S"]
    assert abs(w["a"] - w["Z"]) < 0.01, "a and Z: same width, no shared element"


def test_identical_elements_stay_flat():
    """P3: §118's shape on a second element."""
    E = [_last(f"S:{om}")["E"] for om in (30, 40, 55)]
    assert max(E) / min(E) < 1.5


def test_the_k1_over_k2_law_failed():
    """P2 refuted: E * k2/k1 falls ~80x over Omega = 30 -> 55 in arm X; E saturates."""
    d = _d()
    prod = [d[f"X:{om}"]["series"][-1]["E"] * d[f"X:{om}"]["k2"] / d[f"X:{om}"]["k1"]
            for om in (30, 55)]
    assert prod[0] / prod[1] > 30
    assert _last("X:55")["E"] > 10


def test_zone_is_necessary():
    """P4: same stage 1, no zone -> p_d an order of magnitude down, and the gap widens."""
    r = [_last(f"N:{om}")["p_d"] / _last(f"Z:{om}")["p_d"] for om in (40, 55)]
    assert r[0] < 0.1 and r[1] < r[0]


def test_zone_width_alone_does_not_set_p_d():
    """Q1 FAILED: a and Z share their width and differ ~5x."""
    for om in (40, 55):
        assert _last(f"Z:{om}")["p_d"] / _last(f"a:{om}")["p_d"] > 1.5


def test_p_d_monotone_in_width():
    """Q2 held at both Omega."""
    for om in (40, 55):
        p = {k: _last(f"{k}:{om}")["p_d"] for k in "NdcZbX"}
        assert p["N"] < p["d"] < p["c"] < p["Z"] < p["b"] < p["X"]


def test_the_floor_out_of_sample():
    """Q3, the result: P(stage 2 low | stage 1 high)/(k1 t) moves < 1.5x from Omega = 40 to 55
    on the three pairs never run before the prediction, while the downstream's own rate does not
    enter -- it moves by orders of magnitude where it is resolvable."""
    d = _d()
    for arm in "abc":
        f = [_floor(f"{arm}:{om}") for om in (40, 55)]
        assert max(f) / min(f) < 1.5
    k2 = [d[f"c:{om}"]["k2"] for om in (40, 55)]
    assert k2[0] / k2[1] > 100


def test_conditioned_models_are_off_by_orders_of_magnitude():
    """E where p0 is above the solver floor (C2): 15x with the weakest gap, 1e4x on arm c."""
    assert _last("X:55")["E"] > 10
    assert _last("c:55")["E"] > 5e3
    assert _last("c:55")["p0"] > 1e-9


def test_instrument_converged():
    """C1: the box moves nothing by 5%; C2: halving the step moves nothing above 1e-5."""
    for key, r in _d()["convergence"].items():
        lim = 0.05 if key.endswith("box2.0") else 1e-5
        assert all(abs(v) < lim for v in r["rel"].values()), key
