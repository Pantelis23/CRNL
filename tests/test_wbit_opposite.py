"""§132: the opposite sweep's stored verdict, and its construction."""

import json
from pathlib import Path

import pytest

from experiments.wbit_opposite import C, HS, b2_of
from experiments.wbit_cascade import HI1

RES = Path(__file__).resolve().parents[1] / "results" / "wbit_opposite.json"


def test_c_is_held_fixed():
    for h in HS:
        assert h - b2_of(h) == pytest.approx(C * (2 * HI1 - 1), abs=1e-12)


def test_stored_verdict_kills_window_as_stated():
    rep = json.loads(RES.read_text())
    for s in (0.15, 0.12):
        v = rep[f"verdict:{s}"]
        assert v["clean"] == list(HS)
        assert all(x < 0 for x in v["steps"])                        # O2: monotone
        frac = v["drop"] / v["predicted"]
        assert 0 < frac < 0.5                                         # right sign, below the factor-2 band
        assert abs(v["drop"]) > 0.5                                   # not flat either
