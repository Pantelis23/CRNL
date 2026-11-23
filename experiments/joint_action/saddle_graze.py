"""§134.28: the exact linear saddle fed along a GRAZING path -- suspect 6's isolating test, with the gap as the only lever.

§134.17 fed the linear saddle (Maier-Stein normal form B = [[1, 0], [c, -mu]], D = I) from far along the asymptotic
incoming direction e_in and found no excess (beta = 0). But the chain's path runs PARALLEL to Sigma at a small gap before
entering H (§134.18). Here the source sits at physical normal-form position (g0, L0): normal offset g0 = gap / sqrt(D22)
from Sigma (the X2 axis) and L0 = (3 - 1.9) / sqrt(D11) upstream along it. The chain's Omega-scaling is a uniform
stretch of every coordinate by lam = sqrt(Omega lambda_u). With y_ref the deterministic trajectory point from the source
where X2 has fallen to half (a fixed physical point on the incoming path),
    Q(Omega) = [ln v(0) + S(0)] - [ln v(lam y_ref) + S(lam y_ref)],   E_lin = dQ / dln Omega,
which is the linear model's analogue of the chain's excess at H (L(H) relative to the regular outer path).

Geometries: XK-like (its mu, c; gap 0.0107 at x1 = 2.7), X-like (gap 0.0024), and a non-grazing control (gap 0.10, the
same as XK otherwise). Omega = 400, 800, 1600, 2800, 4000.

Predictions, written before running:
  G1 (the test) the XK-like geometry gives E_lin in [0.30, 0.45] over Omega = 1600 -> 4000, matching the chain's
     0.37 within its spread -> suspect 6 (grazing geometry makes the excess) SURVIVES, now isolated in an exact
     model. E_lin within +-0.05 of 0 there -> suspect 6 as a local-geometry mechanism DIES. Otherwise UNRESOLVED.
  G2 (control) the non-grazing control gives |E_lin| < 0.05 over the same Omega.
  My expectation: G1 holds and G2 holds.

Usage (from this directory): python saddle_graze.py
"""
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def params(hk, r2="0.606"):
    """(mu, c, lambda_u, D11, D22) for the pair, computed in a fresh process (model.py reads HILL_K at import)."""
    code = ("import sys, json, numpy as np; sys.path.insert(0, '.')\n"
            "import saddle_inner as si\n"
            "from model import lam1, mu1, lam2, mu2, dlam1, dmu1, hill, dhill, A2, C2, R2_2\n"
            "mu, c, _ = si.normal_form()\n"
            "x2 = R2_2\n"
            "lu = hill(3.0) * A2 * 2 * x2 - (3 * x2 * x2 + C2)\n"
            "print(json.dumps([mu, c, lu, lam1(3.0) + mu1(3.0), lam2(3.0, x2) + mu2(x2)]))\n")
    out = subprocess.run([sys.executable, "-c", code], cwd=HERE, env=dict(os.environ, HILL_K=hk, R2_STAGE2=r2),
                         capture_output=True, text=True, check=True).stdout.strip().splitlines()[-1]
    return json.loads(out)


def lnv_S(B, y, x0):
    import saddle_inner as si
    rho = max(np.linalg.norm(x0), 2.0)
    ts = np.concatenate([np.geomspace(1e-6, 1.0, 800), np.linspace(1.0, 2.0 * math.log(rho) + 30, 30000)[1:]])
    return si.ln_v_and_S(B, np.asarray(y, float), np.asarray(x0, float), ts)


if __name__ == "__main__":
    sys.path.insert(0, str(HERE))
    out = {}
    cases = {"XK_like": ("1.2", 0.0107, "0.606"), "X_like": ("1.0", 0.0024, "0.606"),
             "control_gap0.10": ("1.2", 0.10, "0.606"), "S_like": ("1.0", 0.0122, "1.005")}
    only = os.environ.get("ONLY")
    for name, (hk, gap, r2) in cases.items():
        if only and name != only:
            continue
        mu, c, lu, D11, D22 = params(hk, r2)
        B = np.array([[1.0, 0.0], [c, -mu]])
        g0, L0 = gap / math.sqrt(D22), (3.0 - 1.9) / math.sqrt(D11)
        x0p = np.array([g0, -L0])                                      # upstream along Sigma (the X2 axis)
        # the deterministic trajectory from the source: X(t) = e^{Bt} x0; y_ref where X2 has fallen to half
        from scipy.linalg import expm
        tt = np.linspace(0, 5.0 / mu, 20001)
        traj = np.array([expm(B * t) @ x0p for t in tt])
        k = int(np.argmin(np.abs(np.abs(traj[:, 1]) - L0 / 2)))
        yref = traj[k]
        rec = {"mu": mu, "c": c, "g0": g0, "L0": L0, "yref": yref.tolist(), "cells": {}}
        for om in (400, 800, 1600, 2800, 4000):
            lam = math.sqrt(om * lu)
            l0, s0 = lnv_S(B, np.zeros(2), lam * x0p)
            lr, sr = lnv_S(B, lam * yref, lam * x0p)
            rec["cells"][str(om)] = {"Q": (l0 + s0) - (lr + sr), "S0": s0}
        oms = sorted(rec["cells"], key=int)
        rec["E"] = {f"{a}-{b}": (rec["cells"][b]["Q"] - rec["cells"][a]["Q"]) / math.log(int(b) / int(a))
                    for a, b in zip(oms, oms[1:])}
        out[name] = rec
        print(f"{name}: mu={mu:.3f} c={c:+.3f} g0={g0:.5f} L0={L0:.4f}  E_lin per window " +
              " ".join(f"{w}:{v:+.3f}" for w, v in rec["E"].items()), flush=True)
    res = HERE.parents[1] / "results" / "saddle_graze.json"
    prev = json.loads(res.read_text()) if res.exists() else {}
    prev.update(out)
    res.write_text(json.dumps(prev, indent=1))

# OUTCOME (§134.28): E_lin per window 400-800 / 800-1600 / 1600-2800 / 2800-4000:
#   XK_like +0.397 / +0.388 / +0.384 / +0.380 (chain E at H: +0.373 / +0.374) -> G1 HOLDS, an absolute match within 0.01;
#   X_like  +0.400 / +0.393 / +0.391 / +0.390 (chain X +0.393 -> +0.359, still converging);
#   control +0.204 / +0.132 / +0.077 / +0.047 -- decaying toward 0; G2 as written FAILS at 1600-2800 (0.077 vs |E| < 0.05)
#   and holds at 2800-4000: a fixed tolerance on a converging quantity (rule 20), recorded against the criterion.
# REGISTERED BEFORE THE S RUN: S_like (mu = 3.33, S's gap 0.0122 at x1 = 2.7) must give E_lin(2800-4000) within 0.03 of
# the chain's S value +0.318 (§134.16); outside it, the linear grazing model does not carry the mu-dependence.
# -> S_like: +0.346 / +0.340 / +0.339 / +0.339; chain S +0.318 at 2800-4000 (still rising; geometric-tail estimate ~0.337).
#    |0.339 - 0.318| = 0.021 <= 0.03 -> HOLDS. The linear grazing model carries the mu-dependence.
