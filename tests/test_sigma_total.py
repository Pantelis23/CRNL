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


def test_saddle_climb_bound_and_the_envelope_violating_it():
    """§134.2: h along Sigma at the joint saddle, and the ray envelope above 0.5 h delta^2 near the saddle."""
    import os
    import subprocess
    import sys
    root = Path(__file__).resolve().parents[1]
    code = (
        "import numpy as np\n"
        "from scipy.linalg import solve_continuous_lyapunov\n"
        "import model\nfrom model import *\n"
        "xs1, xs2 = 3.0, R2_2\n"
        "f1p = 2*A1*xs1 - (3*xs1*xs1 + C1); hh, dh = model.hill(xs1), model.dhill(xs1)\n"
        "f2p = 2*hh*A2*xs2 - (3*xs2*xs2 + C2)\n"
        "A = np.array([[f1p,0.0],[dh*A2*xs2*xs2, f2p]]); D = np.diag([2*(A1*xs1*xs1+B1), 2*(hh*A2*xs2*xs2+B2)])\n"
        "Z = np.linalg.inv(solve_continuous_lyapunov(A, -D))\n"
        "w, v = np.linalg.eig(A); es = np.real(v[:, np.argmin(np.real(w))]); es = es/es[0]\n"
        "print(float(es@Z@es))\n")
    out = subprocess.run([sys.executable, "-c", code], cwd=root / "experiments" / "joint_action",
                         env={**os.environ, "HILL_K": "1.2"}, capture_output=True, text=True, check=True)
    h = float(out.stdout.strip())
    assert abs(h - 0.150044) < 1e-5
    s = json.loads(RES.read_text())
    eps = eps_fn(s)
    assert eps(2.991) > 0.5 * h * 0.009 ** 2                   # the envelope is not W there


def test_omega_4000_registered_windows():
    """§134.5: M6 and T-CASC-au's window, both written before Omega = 4000 existed, and its conservation."""
    s = json.loads(RES.read_text())
    a, b = s["XK:4000:sad2:cap1.1"]["total"], s["XK:4000:sigma:cap1.1"]["total"]
    assert abs(b / a - 1) < 1e-9
    m6 = s["S134"]["3400-4000"]
    assert abs(m6["a_measured"] - m6["a_model"]) <= 0.01
    tc = s["S134"]["2800-4000"]
    assert abs(tc["a_measured"] - tc["a_model"]) < 0.03 and tc["a_measured"] <= 0.01
    seq = [s["S134"][w]["a_measured"] for w in ("800-1600", "1600-2000", "2000-2800", "2800-3400", "3400-4000")]
    assert all(y > x for x, y in zip(seq, seq[1:]))       # measured a rises at every step


def test_omega_5000_registered_and_its_gamma_dependence():
    """§134.7: M7 (refuted) and M8 (held), written before Omega = 5000 existed; and M7's verdict read across every
    Gamma the two independent methods allow, since the measured a shifts by dGamma * Omega."""
    s = json.loads(RES.read_text())
    t4, t5 = s["XK:4000:sad2:cap1.1"]["total"], s["XK:5000:sad2:cap1.1"]["total"]
    assert abs(s["XK:5000:sigma:cap1.1"]["total"] / t5 - 1) < 1e-9          # conservation across cuts
    assert s["XK:5000:sigma:cap1.1"]["sd_x1"] < 0.665                          # M8
    a = lambda G: (math.log(t5 / t4) - G * 1000) / math.log(5000 / 4000)
    assert a(GAMMA) >= 0                                                       # M7's band: refuted
    ja = json.loads((RES.parent / "joint_action.json").read_text())["XK"]
    g_mam = max(v["Gamma"] for v in ja["mam"].values())
    g_ray = max(ja["convergence_rays"].values())
    assert a(g_ray) >= 0 and a(g_mam) >= 0                                     # holds for both methods' Gamma
    assert a(7.345e-3) < 0                                                     # ...and would flip at the free fit's


def test_xk_minimum_action_rises_and_stays_below_the_flip():
    """§134.7: Gamma_MAM rises with N (G1) and never reaches 7.3410e-3, the rays' value plus their stated error."""
    ja = json.loads((RES.parent / "joint_action.json").read_text())["XK"]["mam"]
    g = [ja[k]["Gamma"] for k in sorted(ja, key=int)]
    assert all(y > x for x, y in zip(g, g[1:])) and max(g) < 7.3410e-3
