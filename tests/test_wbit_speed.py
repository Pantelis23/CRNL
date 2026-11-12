"""§130: the design window's closed-form exponent, the saddle-node latency law, and the trade-off table."""

import json
from pathlib import Path

import numpy as np
import pytest

import experiments.wbit_speed as ws
from experiments.wbit_action import S_DIP, escape_action

RES = Path(__file__).resolve().parents[1] / "results" / "wbit_speed.json"


def test_window_bounds():
    assert ws.C_MIN == pytest.approx(0.17910, abs=1e-5)
    assert ws.C_HALF == pytest.approx(0.35821, abs=1e-5)
    assert ws.edge_of(ws.C_HALF) == pytest.approx(0.5, abs=1e-12)
    assert ws.H - ws.C_MIN_TRUE * ws.K == pytest.approx(-ws.D_C_TRUE, abs=1e-15)


def test_closed_form_gamma_inside_window():
    """Gamma(c) = -S2(b2(c)) for an edge below 1/2, by the fast collocation route at long L."""
    c = 0.26
    rows = ws.gamma_colloc_long(c, (100,))
    assert rows[-1]["feasible"]
    assert S_DIP - rows[-1]["W_star"] == pytest.approx(-escape_action(ws.b2_of(c)), abs=2e-5)


def test_fast_collocation_reproduces_section_127():
    rows = ws.gamma_colloc_long(0.6707, (40,))
    assert rows[-1]["W_star"] == pytest.approx(0.31646, abs=2e-5)


def test_latency_follows_the_normal_form():
    prod = [ws.det_latency(ws.C_MIN_TRUE + e) * e ** 0.5 for e in (1e-4, 1e-5, 1e-6)]
    gaps = [abs(p - ws.NF_LIMIT) for p in prod]
    assert gaps[0] > gaps[1] > gaps[2] and gaps[2] < 0.01


def test_stored_tradeoff_is_monotone():
    rep = json.loads(RES.read_text())
    rows = sorted(rep["P4"], key=lambda r: r["c"])
    R = [r["R"] for r in rows]
    tau = [r["tau_noisy"] for r in rows]
    assert all(b > a for a, b in zip(R, R[1:]))           # R rises with c (edge moves up)
    assert all(b < a for a, b in zip(tau, tau[1:]))       # latency falls with c
    for r in rep["P4"]:
        assert abs(r["mass_check"] - 1) < 1e-12
    for k, v in rep.items():
        if k.startswith("P1long:") and "diff" in v:
            assert abs(v["diff"]) < 1e-4
