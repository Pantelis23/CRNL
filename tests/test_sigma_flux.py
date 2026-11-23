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


def test_raw_flux_rebins_to_the_binned_cell(xk):
    """§134.9: cut_flux(raw=True) is the same flux, unbinned: re-binned it reproduces the binned density to rounding."""
    A, B = xk
    sig = sigma_table(A, B)
    x1, f = cut_flux(A, B, 120, "sigma", sig, raw=True)
    d = np.array(cut_flux(A, B, 120, "sigma", sig)["flux_density"])
    rb = np.array([f[(x1 >= 1.0 + 0.05 * k) & (x1 < 1.05 + 0.05 * k)].sum() / 0.05 for k in range(44)])
    m = d > 1e-12 * d.max()
    assert np.max(np.abs(rb[m] / d[m] - 1)) < 1e-12


def test_core_width_scales_as_the_exact_gaussian():
    """§134.10: the stored core widths shrink with nu_eff in [0.43, 0.5] on both kernels over 1600 -> 4000 (never in
    the cubic band), and track the exact Gaussian e^{-h Omega d^2 / 2}, h = 0.1502, at a steady ~1.18x."""
    import math
    c = json.loads((ROOT / "results" / "sigma_core.json").read_text())
    for h in ("h0.01", "h0.02"):
        nu = math.log(c["1600"][h]["fwhm"] / c["4000"][h]["fwhm"]) / math.log(2.5)
        assert 0.43 <= nu <= 0.5
    r = []
    for om in ("1600", "2800", "4000"):
        f = c[om]["h0.01"]["fwhm"]
        r.append(math.sqrt(f * f - (2.3548 * 0.01) ** 2) / (2 * math.sqrt(math.log(2) / (0.5 * 0.1502 * int(om)))))
    assert max(r) - min(r) < 0.05 and 1.1 < min(r) < 1.25


def test_maier_stein_class_and_skew_do_not_explain_the_core():
    """§134.11: XK's saddle has mu > 1, its wedge-side rate is h/2, and the skew (c = 0.20) is far too small to
    produce the measured 18% excess width."""
    m = json.loads((ROOT / "results" / "ms_core_xk.json").read_text())
    assert m["mu"] > 1 and abs(m["k_wedge_x1"] - 0.1502 / 2) < 2e-4
    assert abs(m["k_wedge_x1"] / m["k_other_x1"] - 1) < 0.01
    assert all(1.15 < c["ratio_meas_over_pred"] < 1.22 for c in m["cells"].values())


def test_density_is_wkb_regular_along_the_path_and_anomalous_at_the_saddle():
    """§134.12: ln P + Omega W converges along the optimal path (last-window slope <= 0.01 in magnitude away from H),
    while at H itself it keeps growing (slope > 0.3 in every window)."""
    import math
    d = json.loads((ROOT / "results" / "density_saddle.json").read_text())
    ks = sorted(d, key=int)
    assert all((d[b]["L"] - d[a]["L"]) / math.log(int(b) / int(a)) > 0.3 for a, b in zip(ks, ks[1:]))
    p = json.loads((ROOT / "results" / "density_along_path.json").read_text())
    for x2t in ("2.5", "2.0", "1.6", "1.3", "1.0", "0.8", "0.7"):
        s = (p["1600"][x2t]["L"] - p["1200"][x2t]["L"]) / math.log(1600 / 1200)
        assert abs(s) <= 0.01


def test_ridge_test_premise_and_convergence():
    """§134.13: X's ridge differs from XK's by > 2x at Omega = 1600 (the test's premise), and the two pairs' saddle
    excesses converge to within 0.01 in the last window while differing by > 0.2 in the first."""
    import math
    x = json.loads((ROOT / "results" / "ridge_x.json").read_text())
    k = json.loads((ROOT / "results" / "density_saddle.json").read_text())
    c = json.loads((ROOT / "results" / "sigma_core.json").read_text())
    assert c["1600"]["share_within_0.05"] / x["1600"]["core_share"] > 2
    E = lambda d, a, b: (d[b]["L"] - d[a]["L"]) / math.log(int(b) / int(a))
    assert abs(E(x, "1200", "1600") - E(k, "1200", "1600")) < 0.01
    assert abs(E(x, "200", "400") - E(k, "200", "400")) > 0.2


def test_saddle_density_excess_levels_off():
    """§134.14: the excess at H rises in every window but its increments shrink, and the two largest windows sit in
    [0.35, 0.42] within 0.005 of each other."""
    import math
    d = json.loads((ROOT / "results" / "density_saddle.json").read_text())
    ks = sorted(d, key=int)
    E = [(d[b]["L"] - d[a]["L"]) / math.log(int(b) / int(a)) for a, b in zip(ks, ks[1:])]
    inc = [y - x for x, y in zip(E, E[1:])]
    assert all(v > 0 for v in inc) and all(b < a for a, b in zip(inc, inc[1:]))
    assert all(0.35 <= e <= 0.42 for e in E[-2:]) and abs(E[-1] - E[-2]) < 0.005


def test_pair_s_bvp_and_density_excess():
    """§134.15: S's W* is steady along the BVP continuation (Schur-projected A), and its saddle excess is still
    climbing at Omega = 1600, below XK's."""
    import math
    c = json.loads((ROOT / "results" / "bvp_s_tol1e-06_drop0.json").read_text())["cells"]
    W = [x["W"] for x in c if x.get("axis") != "tol" and x["status"] == 0 and x["offset_A"] <= 1e-5
         and x["offset_B"] <= 1e-4]
    assert len(W) >= 3 and max(W) - min(W) < 1e-9 and abs(W[-1] - 0.14457152566) < 2e-9
    d = json.loads((ROOT / "results" / "density_mu_S.json").read_text())
    ks = sorted(d, key=int)
    E = [(d[b]["L"] - d[a]["L"]) / math.log(int(b) / int(a)) for a, b in zip(ks, ks[1:])]
    assert all(y > x for x, y in zip(E, E[1:])) and E[-1] < 0.368 - 0.05


def test_saddle_exponent_moves_with_mu():
    """§134.16: at 2800 -> 4000 pair S (mu = 3.33) sits >= 0.05 below XK (mu = 5.21) with a last increment < 0.02."""
    import math
    E = lambda d, a, b: (d[b]["L"] - d[a]["L"]) / math.log(int(b) / int(a))
    s = json.loads((ROOT / "results" / "density_mu_S.json").read_text())
    k = json.loads((ROOT / "results" / "density_saddle.json").read_text())
    es, ek = E(s, "2800", "4000"), E(k, "2800", "4000")
    assert ek - es >= 0.05 and abs(es - E(s, "1600", "2800")) < 0.02


def test_linear_saddle_is_wkb_regular_at_both_mu():
    """§134.17: the exact linear inner saddle gives beta = 0 (|beta| <= 0.002 over r = 5-200) at XK's and S's (mu, c)."""
    import numpy as np
    for tag in ("xk", "s"):
        run = json.loads((ROOT / "results" / f"saddle_inner_{tag}.json").read_text())["runs"]["1e+06"]
        r, sl = np.array(run["r"]), np.array(run["slope"]) / 2
        m = (r >= 5) & (r <= 200)
        assert np.all(np.abs(sl[m]) <= 0.002)


def test_gaussian_inner_saddle_verified():
    """§134.19: the closed-form Green's function integrates to 1, and Monte Carlo agrees with the time-integrated
    density to within 1% at the finer dt (|z| < 2)."""
    c = json.loads((ROOT / "results" / "saddle_inner_check_xk.json").read_text())
    assert all(abs(v - 1) < 1e-6 for v in c["W1"].values())
    fine = c["W2"]["0.0005"]
    assert all(abs(r["ratio"] - 1) < 0.01 and abs(r["z"]) < 2 for r in fine.values())


def test_excess_profile_rises_up_the_incoming_side():
    """§134.21: at fixed stretched distance along x2 through H, the excess's last-window slope rises monotonically from
    below H to above it (r = -3 .. +10)."""
    import math
    d = json.loads((ROOT / "results" / "excess_profile.json").read_text())
    s = [(d["1600"][r]["L_flat"] - d["1200"][r]["L_flat"]) / math.log(1600 / 1200) for r in ("-3", "-1", "0", "1", "3", "10")]
    assert all(b > a for a, b in zip(s, s[1:])) and s[0] < 0.15 and s[-1] > 1.0


def test_depletion_zone_along_the_incoming_path():
    """§134.25: on the orbit at fixed stretched distance, r = 5 and 10 have negative last-window slopes while r = 1 is
    positive (the excess sits only near H)."""
    import math
    d = json.loads((ROOT / "results" / "excess_along_path.json").read_text())
    s = {r: (d["1600"][r]["L"] - d["1200"][r]["L"]) / math.log(1600 / 1200) for r in ("1", "5", "10")}
    assert s["1"] > 0.2 and s["5"] < -0.3 and s["10"] < -0.3


def test_depletion_deeper_for_the_closer_path():
    """§134.26: pair X (path 5x closer to Sigma) has a deeper depletion than XK at r = 5, 10 and 20."""
    import math
    x = json.loads((ROOT / "results" / "excess_along_path_X.json").read_text())
    k = json.loads((ROOT / "results" / "excess_along_path.json").read_text())
    sl = lambda d, r: (d["1600"][r]["L"] - d["1200"][r]["L"]) / math.log(1600 / 1200)
    assert all(sl(x, r) < sl(k, r) for r in ("5", "10", "20")) and min(sl(x, r) for r in ("5", "10", "20")) < -0.782
