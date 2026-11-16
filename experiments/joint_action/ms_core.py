"""§134.11: the saddle core's width against Maier & Stein's skewed exit distribution (mu > 1), in absolute terms.

Literature (rule 21). Maier & Stein, "Limiting exit location distributions in the stochastic exit problem", SIAM J.
Appl. Math. 57 (1997), arXiv adap-org/9407003. With mu = |lambda_s(H)| / lambda_u(H) > 1 the exit location density on
the separatrix lives on the eps^{1/2} lengthscale and is generically two-sided Gaussian (their eq. 78), in coordinates
with D(H) = I, lambda_u = 1 and B(H) = [[1, 0], [c, -mu]] (boundary along e_s):
    wedge side:      p ~ exp(-mu s^2 / eps)
    other side:      p ~ exp(-mu (mu-1)^2 / (c^2 + (mu-1)^2) * s^2 / eps)
(read from extracted PDF text; consistent with both limits stated there: equal rates at c = 0, and the non-wedge rate
-> 0 as mu -> 1). Here eps = 1/Omega and D = diag(lambda_i + mu_i) at H for the CME.

§134.10 measured the core's FWHM at 1.17 / 1.18 / 1.19 x the one-sided Gaussian e^{-h Omega d^2 / 2}, h = 0.1502
(which is exactly the wedge-side rate: lambda_s / D_11 = 0.0751 = h/2).

Predictions, written before computing c:
  MS1 (absolute) The two-sided Gaussian with the computed c predicts FWHM = sqrt(ln2/(k_w Omega)) + sqrt(ln2/(k_o Omega))
      (k_w, k_o the rates in x1 units). Its ratio to the measured FWHM (h = 0.01, kernel removed) lies in [0.9, 1.1]
      at each of Omega = 1600, 2800, 4000 -- i.e. Maier-Stein accounts for the 18% excess. Outside that band at
      any Omega: the excess is NOT the wedge asymmetry.
  MS2 The broader side is the non-wedge side. (Which side of x1 = 3 is broader was already seen in §134.10's
      w_left / w_right -- left -- so only the agreement of the geometry with it is tested, not the side itself.)
  My expectation: MS1 holds.

Usage (from this directory): HILL_K=1.2 python ms_core.py
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from model import lam1, mu1, lam2, mu2, dlam1, dmu1, hill, dhill, A2, C2, R2_2   # noqa: E402

H = np.array([3.0, R2_2])
x1, x2 = H
Bm = np.array([[dlam1(x1) - dmu1(x1), 0.0],
               [dhill(x1) * A2 * x2 * x2, hill(x1) * A2 * 2 * x2 - (3 * x2 * x2 + C2)]])
D = np.diag([lam1(x1) + mu1(x1), lam2(x1, x2) + mu2(x2)])
Dh = np.sqrt(D)
By = np.linalg.inv(Dh) @ Bm @ Dh                                  # y = D^{-1/2} x: unit diffusion
w, V = np.linalg.eig(By)
iu, is_ = int(np.argmax(w.real)), int(np.argmin(w.real))
lu, ls = w[iu].real, -w[is_].real
mu = ls / lu
es = V[:, is_].real; es /= np.linalg.norm(es)                       # boundary direction (Sigma's tangent)
n = np.array([es[1], -es[0]])                                       # orthonormal normal
Q = np.column_stack([n, es])                                        # new axes: x1' normal, x2' along the boundary
Bq = Q.T @ By @ Q / lu                                              # time rescaled: lambda_u = 1
c = Bq[1, 0]
kw = ls                                                             # exp(-Omega * ls * s_y^2): wedge side
ko = ls * (mu - 1) ** 2 / (c * c + (mu - 1) ** 2)
scale = (Dh @ es)[0]                                                # x1 displacement per unit s_y
kw_x, ko_x = kw / scale ** 2, ko / scale ** 2
print(f"lambda_u={lu:.5f} lambda_s={ls:.5f} mu={mu:.4f} c={c:+.4f} (check Bq[0,1]={Bq[0, 1]:+.2e}, "
      f"Bq diag {Bq[0, 0]:+.4f} {Bq[1, 1]:+.4f})")
print(f"rates in x1 units: wedge k_w={kw_x:.5f} (h/2 = 0.0751), other k_o={ko_x:.5f}, width ratio "
      f"{math.sqrt(kw_x / ko_x):.3f}")
core = json.loads((HERE.parents[1] / "results" / "sigma_core.json").read_text())
out = {"mu": mu, "c": c, "k_wedge_x1": kw_x, "k_other_x1": ko_x, "cells": {}}
for om in ("1600", "2800", "4000"):
    O = int(om)
    pred = math.sqrt(math.log(2) / (kw_x * O)) + math.sqrt(math.log(2) / (ko_x * O))
    f = core[om]["h0.01"]["fwhm"]
    meas = math.sqrt(f * f - (2.3548 * 0.01) ** 2)
    wl, wr = core[om]["h0.01"]["w_left"], core[om]["h0.01"]["w_right"]
    out["cells"][om] = {"fwhm_pred": pred, "fwhm_meas": meas, "ratio_meas_over_pred": meas / pred,
                        "w_left": wl, "w_right": wr}
    print(f"Omega={O}: FWHM predicted {pred:.4f}, measured {meas:.4f}, measured/predicted {meas / pred:.3f}; "
          f"w_left {wl:.4f} w_right {wr:.4f}")
(HERE.parents[1] / "results" / "ms_core_xk.json").write_text(json.dumps(out, indent=1))

# OUTCOME (§134.11): mu = 5.2083 (the mu > 1 class: eps^{1/2} lengthscale -- matches §134.10's nu ~ 1/2), c = +0.2013,
# so the wedge asymmetry is negligible: k_w = 0.07519 (= h/2, as it must), k_o = 0.07502, width ratio 1.001.
# MS1 FAILED: measured/predicted FWHM = 1.173 / 1.179 / 1.191 at Omega = 1600 / 2800 / 4000 -- the 18% excess is NOT the
# Maier-Stein skew. MS2 moot: predicted asymmetry 0.1%; measured w_left ~ w_right at 1600 and 4000 (2800's was noise).
# My expectation (MS1 holds) was WRONG.
