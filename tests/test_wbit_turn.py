"""§128: the matrix-free exit solver agrees with §126's direct solve, and the stored slopes turn."""

import json
import math
from pathlib import Path

import pytest

from experiments.wbit_cascade import PAIRS, exits
from experiments.wbit_turn import exits_mf

RES = Path(__file__).resolve().parents[1] / "results" / "wbit_turn.json"


@pytest.mark.parametrize("name", ["Z", "N"])
def test_matrix_free_matches_direct(name):
    c, b2 = PAIRS[name]
    d, m = exits(60, 0.30, c, b2), exits_mf(60, 0.30, c, b2)
    assert m["R"] == pytest.approx(d["R"], rel=1e-9)
    assert abs(m["mass_check"] - 1) < 1e-10


def _slopes(rep, name, sig, M):
    R = [rep[f"{name}:{s}:{M}"]["R"] for s in sig]
    return [(math.log(R[i]) - math.log(R[i - 1])) / (1 / sig[i] ** 2 - 1 / sig[i - 1] ** 2) for i in range(1, len(R))]


def test_both_slope_sequences_reverse_in_readable_cells():
    rep = json.loads(RES.read_text())
    for M in (480, 640, 960):
        sig = [0.18, 0.16, 0.15] if M == 960 else [0.25, 0.22, 0.20, 0.18, 0.16, 0.15]
        n = _slopes(rep, "N", sig, M)
        for s in sig:                                   # readability rule: mass deficit < 1% of q_fall
            r = rep[f"N:{s}:{M}"]
            assert abs(r["mass_check"] - 1) < 0.01 * r["q_fall"]
        assert n[-1] > n[-2] and min(n) > -0.267        # N turns back, never past cell 2's own route
    z = _slopes(rep, "Z", [0.34, 0.30, 0.27, 0.25, 0.22, 0.20, 0.18], 640)
    assert z.index(min(z)) not in (0, len(z) - 1)        # Z's slope has an interior minimum
