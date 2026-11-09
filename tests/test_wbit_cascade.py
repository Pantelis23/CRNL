"""§126 -- two W-bit cells in series: the mechanism's sign transfers; per-cell reliability does not compose."""

import json
import pathlib

STORE = pathlib.Path("results/wbit_cascade.json")
S = (0.34, 0.3, 0.27, 0.25)


def _d():
    return json.loads(STORE.read_text())


def test_grid_converged():
    d = _d()
    for k in "ZN":
        for s in S:
            assert abs(d[f"{k}:{s}:160"]["R"] / d[f"{k}:{s}:120"]["R"] - 1) < 0.05
            assert abs(d[f"{k}:{s}:160"]["mass_check"] - 1) < 1e-10


def test_zone_rises_no_zone_falls():
    d = _d()
    z = [d[f"Z:{s}:160"]["R"] for s in S]
    n = [d[f"N:{s}:160"]["R"] for s in S]
    assert z == sorted(z) and n == sorted(n, reverse=True)


def test_per_cell_reliability_does_not_compose():
    d = _d()
    r = d["Z:0.25:160"]
    assert r["R"] / (r["cell2_own_rate"] / r["dip_rate"]) > 10


def test_identical_healthy_cell_two():
    from experiments.wbit_cascade import HI1, PAIRS
    drives = [b2 + c * (2 * HI1 - 1) for c, b2 in PAIRS.values()]
    assert abs(drives[0] - drives[1]) < 1e-12
