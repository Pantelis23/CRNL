"""§134: the one-calibration spread-exit model against the exact totals (stored)."""

import json
import math
from pathlib import Path

import numpy as np

from experiments.sigma_total import GAMMA, eps_fn

RES = Path(__file__).resolve().parents[1] / "results" / "sigma_flux.json"


def test_model_reproduces_stored_and_tracks_the_exact_total():
    s = json.loads(RES.read_text())
    eps = eps_fn(s)
    bins = np.array(s["XK:800:sigma"]["flux_bins_x1"])
    F = np.array(s["XK:800:sigma"]["flux_density"])
    C = F / (800 ** 0.5 * np.exp((GAMMA - eps(bins)) * 800))
    model = lambda o: float((C * o ** 0.5 * np.exp((GAMMA - eps(bins)) * o)).sum())
    offs = []
    for w in ("800-1600", "1600-2000", "2000-2800", "2800-3400"):
        o1, o2 = map(int, w.split("-"))
        a_mod = (math.log(model(o2) / model(o1)) - GAMMA * (o2 - o1)) / math.log(o2 / o1)
        rec = s["S134"][w]
        assert abs(a_mod - rec["a_model"]) < 1e-9
        offs.append(rec["a_model"] - rec["a_measured"])
    assert all(0 < b < a for a, b in zip(offs, offs[1:]))        # converging, window by window
    assert offs[-1] < 0.05


def test_cap_and_conservation():
    s = json.loads(RES.read_text())
    assert abs(s["XK:1600:sigma:cap1.1"]["total"] / s["XK:1600:sigma"]["total"] - 1) < 1e-6
    for om in (2000, 2800):
        a, b = s[f"XK:{om}:sad2:cap1.1"]["total"], s[f"XK:{om}:sigma:cap1.1"]["total"]
        assert abs(b / a - 1) < 1e-9
