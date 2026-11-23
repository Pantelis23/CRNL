"""§134.29 (T-CASC-bb): the grazing limit's closed form, E0 = (mu - 1) / (2 mu), checked on the exact linear saddle.

Derivation (written before running). Source on Sigma (g = 0) at stretched distance lam L upstream; the deterministic flow
carries it into H along the stable axis, X2 = lam L e^{-mu t}, reaching the inner scale at t* ~ ln(lam L) / mu. After
that the Green's function at H decays as det Sigma_t^{-1/2} ~ e^{-t} (the unstable direction spreads as e^{2t}), so
v(H) ~ int_{t*}^inf e^{-t} dt = (lam L)^{-1/mu}. The reference point on the incoming path is crossed at a speed ~ lam,
so v(y_ref) ~ lam^{-1}. Both actions vanish (the path is deterministic). Hence Q ~ (1 - 1/mu) ln lam and, with
lam = sqrt(Omega lambda_u), E = dQ/dlnOmega -> (mu - 1)/(2 mu). The coupling c does not enter the exponent
(X1's dynamics are independent of X2, and X1 stays 0 on Sigma).

Predictions:
  Z1 with g0 = 0 the linear model's E approaches (mu - 1)/(2 mu) as lam grows -- read as the sequence over windows
     (rule 20), at mu = 2, 3.333, 5.208, 10, with c = 0 and with c = 0.2: the last window within 0.005 of the
     formula, and the error shrinking window to window. Any mu where the sequence settles elsewhere kills the formula.
  Z2 E0 does not depend on c (the c = 0 and c = 0.2 sequences agree to within 0.005 in the last window).
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import saddle_inner as si                                                      # noqa: E402


def lnv_S(B, y, x0):
    rho = max(np.linalg.norm(x0), 2.0)
    ts = np.concatenate([np.geomspace(1e-6, 1.0, 800), np.linspace(1.0, 2.0 * math.log(rho) + 40, 40000)[1:]])
    return si.ln_v_and_S(B, np.asarray(y, float), np.asarray(x0, float), ts)


if __name__ == "__main__":
    out = {}
    L0 = 0.1264
    lams = [math.sqrt(o) for o in (400, 1600, 6400, 25600, 102400, 409600)]
    for mu in (2.0, 3.333, 5.208, 10.0):
        for c in (0.0, 0.2):
            B = np.array([[1.0, 0.0], [c, -mu]])
            x0p = np.array([0.0, -L0])
            tt = np.linspace(0, 5.0 / mu, 20001)
            traj = np.array([expm(B * t) @ x0p for t in tt])
            yref = traj[int(np.argmin(np.abs(np.abs(traj[:, 1]) - L0 / 2)))]
            Q = []
            for lam in lams:
                l0, s0 = lnv_S(B, np.zeros(2), lam * x0p)
                lr, sr = lnv_S(B, lam * yref, lam * x0p)
                Q.append((l0 + s0) - (lr + sr))
            E = [(Q[i + 1] - Q[i]) / (2 * math.log(lams[i + 1] / lams[i])) for i in range(len(lams) - 1)]
            f = (mu - 1) / (2 * mu)
            key = f"mu{mu:g}_c{c:g}"
            out[key] = {"mu": mu, "c": c, "E": E, "formula": f}
            print(f"{key}: formula {f:.4f}  E per window " + " ".join(f"{e:+.4f}" for e in E), flush=True)
    (HERE.parents[1] / "results" / "saddle_graze_limit.json").write_text(json.dumps(out, indent=1))

# OUTCOME (§134.29). First grid (x1): sequences fell away at large lam (mu = 5.208: ..., 0.366, 0.277) -- a TIME-GRID
# ARTEFACT (v(y_ref)'s integrand in t is ~0.005 wide at lam = 640, against a grid step ~0.001): x4 and x16 agree,
# 0.3947 / 0.3970 / 0.3986. Converged grid (x4), windows to Omega = 409600:
#   mu=2: 0.2461 0.2466 0.2475 0.2467 (formula 0.2500) | mu=3.333: 0.3417 0.3430 0.3439 0.3410 (0.3500)
#   mu=5.208: 0.3928 0.3947 0.3969 0.3985 (0.4040)     | mu=10: 0.4370 0.4388 0.4407 0.4423 (0.4500)
# Z1 AS WRITTEN FAILS: a fixed 0.005 gate on a slowly converging quantity (rule 20, against my own criterion). All four
# approach the formula from below to within 0.3-2%; mu = 5.2 and 10 still rising monotonically; mu = 2 and 3.33 dip in
# the last window (lam = 640; 16x grid not checked there). Z2 HELD: c = 0 and c = 0.2 agree to 4 decimals.
# §134.28 re-verified on the x4 grid: changes <= 0.002 (results/saddle_graze_grid4.json) -- it stands.
