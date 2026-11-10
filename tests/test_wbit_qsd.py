"""§129: the positive-arithmetic exit solver agrees with GMRES, and N's stored slopes rise monotonically."""

import json
import math
from pathlib import Path

import pytest

from experiments.wbit_cascade import PAIRS
from experiments.wbit_qsd import exits_pos
from experiments.wbit_turn import exits_mf

RES = Path(__file__).resolve().parents[1] / "results" / "wbit_qsd.json"
SIG = [0.18, 0.16, 0.15, 0.14, 0.13, 0.12, 0.11, 0.1]


@pytest.mark.parametrize("name", ["Z", "N"])
def test_positive_solver_matches_gmres(name):
    c, b2 = PAIRS[name]
    p, g = exits_pos(60, 0.30, c, b2), exits_mf(60, 0.30, c, b2)
    assert p["R"] == pytest.approx(g["R"], rel=1e-9)
    assert abs(p["mass_check"] - 1) < 1e-12


def _slopes(rep, name, M):
    R = [rep[f"{name}:{s}:{M}"]["R"] for s in SIG]
    return [(math.log(R[i]) - math.log(R[i - 1])) / (1 / SIG[i] ** 2 - 1 / SIG[i - 1] ** 2) for i in range(1, len(R))]


def test_n_slope_rises_every_step_and_stays_below_gamma():
    rep = json.loads(RES.read_text())
    for M in (640, 960):
        n = _slopes(rep, "N", M)
        assert all(b > a for a, b in zip(n, n[1:]))
        assert max(n) < -0.1066
        z = _slopes(rep, "Z", M)
        k = z.index(max(z))
        assert 0 < k < len(z) - 1 and all(b < a for a, b in zip(z[k:], z[k + 1:]))
        assert min(z[k:]) > 0.0125


def test_grid_and_mass():
    rep = json.loads(RES.read_text())
    for v in rep.values():
        assert abs(v["mass_check"] - 1) < 1e-12
    for n in "ZN":
        for s in (0.12, 0.11, 0.1):
            assert rep[f"{n}:{s}:1280"]["R"] == pytest.approx(rep[f"{n}:{s}:960"]["R"], rel=2e-3)
