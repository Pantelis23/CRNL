"""§133: the general-cut flux reproduces fall_profile, conserves the total, and the stored geometry reading."""

import json
from pathlib import Path

import numpy as np
import pytest

import experiments.episode_committor as ec
from experiments.sigma_flux import cut_flux, sigma_table
from experiments.weak_upstream import ARMS, element

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def xk():
    return element(ARMS["XK"][0]), element(ARMS["XK"][1])


def test_saddle_cut_reproduces_fall_profile(xk):
    A, B = xk
    a = np.array(cut_flux(A, B, 120, "sad2")["flux_density"])
    b = np.array(ec.fall_profile(A, B, 120)["flux_density"])
    m = b > 1e-6 * b.max()
    assert np.max(np.abs(a[m] / b[m] - 1)) < 1e-9


def test_totals_conserved_across_cuts(xk):
    A, B = xk
    sig = sigma_table(A, B)
    t1 = cut_flux(A, B, 120, "sad2")["total"]
    t2 = cut_flux(A, B, 120, "sigma", sig)["total"]
    assert t2 == pytest.approx(t1, rel=1e-9)


def test_stored_geometry_reading():
    s = json.loads((ROOT / "results" / "sigma_flux.json").read_text())
    gaps = [s[f"XK:{om}:sad2"]["mean_x1"] - s[f"XK:{om}:sigma"]["mean_x1"] for om in (200, 400, 800, 1600)]
    assert all(b > a for a, b in zip(gaps, gaps[1:])) and 0.05 < gaps[-1] < 0.25
    sd = [s[f"XK:{om}:sigma"]["sd_x1"] for om in (400, 800, 1600)]
    assert sd[-1] > sd[0]                                   # on Sigma the exit widens


def test_action_along_sigma_falls_toward_the_saddle():
    r = np.array(json.loads((ROOT / "results" / "sigma_flux.json").read_text())["ray_profile_XK"]["rows"])
    d = json.loads((ROOT / "results" / "joint_action.json").read_text())
    ws = d["V1_s1"]["E-070"] - d["XK"]["Gamma"]
    x, S = r[:, 1], r[:, 3] - ws
    env = [S[(x >= a) & (x < a + 0.1)].min() for a in (1.3, 1.6, 2.0, 2.4, 2.8)]
    assert all(b < a for a, b in zip(env, env[1:])) and min(env) > 0


def test_pin_scale_cancels_exactly_and_underflow_fails_loudly(xk, monkeypatch):
    """§134.6: the pin is a power of two, so totals are bit-identical to pin = 1 where nothing underflows; a pin
    small enough to underflow dip must raise rather than divide by zero (the Omega = 5000 failure)."""
    A, B = xk
    tot = {}
    for s in (1.0, 2.0 ** 900):
        monkeypatch.setattr(ec, "PIN_SCALE", s)
        tot[s] = cut_flux(A, B, 120, "sad2")
    assert tot[1.0]["total"] == tot[2.0 ** 900]["total"] and tot[1.0]["sd_x1"] == tot[2.0 ** 900]["sd_x1"]
    monkeypatch.setattr(ec, "PIN_SCALE", 2.0 ** -1060)
    with pytest.raises(FloatingPointError):
        cut_flux(A, B, 120, "sad2")
