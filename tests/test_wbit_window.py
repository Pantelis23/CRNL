"""§131: the positive 1-D rate, the healthy-drive parameter, and the stored c-dependence at fixed Gamma."""

import json
from pathlib import Path

import pytest

from experiments.wbit_cascade import PAIRS
from experiments.wbit_cell import flip_rate
from experiments.wbit_qsd import exits_pos
from experiments.wbit_window import B2, h_of, rate_pos

RES = Path(__file__).resolve().parents[1] / "results" / "wbit_window.json"


def test_positive_rate_matches_direct_where_readable():
    assert rate_pos(800, 0.2, 0.0) == pytest.approx(flip_rate(800, 1.5, 0.5, 0.5, 0.2, 0.0), rel=2e-3)


def test_h_default_is_unchanged():
    c, b2 = PAIRS["N"]
    a = exits_pos(60, 0.30, c, b2)
    b = exits_pos(60, 0.30, c, b2, h=0.1)
    assert a["R"] == b["R"]


def test_b2_is_held_fixed():
    for c in (0.30, 0.34, 0.365):
        assert h_of(c) - c * (2 * 0.929284 - 1) == pytest.approx(B2, abs=1e-12)


def test_stored_clean_fits():
    rep = json.loads(RES.read_text())
    for s in (0.15, 0.12):
        fit = rep[f"ext_fit:{s}"]
        assert -20 < fit["slope"] < -5
        mags = [abs(x) for x in fit["local"]]
        assert max(mags) / min(mags) < 1.5
        for c in fit["clean"]:
            assert rep[f"{c}:{s}:320"]["own_over_R"] < 0.1
    assert abs(rep["ext_fit:0.12"]["slope"]) >= abs(rep["ext_fit:0.15"]["slope"])
