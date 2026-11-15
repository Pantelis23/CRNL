"""§134.8 suspect test: is XK's W* (BVP 0.1484037485, §134.8) the GLOBAL minimum, or does a second route exist?

Why. With Gamma fixed at the converged BVP value 7.338162e-3, the exact totals' local slope d lnR/dOmega exceeds Gamma
by +1.0 / +2.1 / +3.0 e-6 at Omega_eff ~ 3100 / 3700 / 4500 -- growing. A prefactor C Omega^a (with or without 1/Omega
corrections) does not fit, and the spread-exit picture can only approach Gamma from BELOW. An effective exponent above
Gamma and climbing is what a second, lower-action route would produce (subdominant at small Omega). MAM and the BVP
both started from X's ray path, so both could be sitting in one basin.

How. mam.py's action at N = 400 (where the X-path start gives W = 0.148407329783, §134.7), from starts through
waypoints (x1w, x2w) spread over the region x1 in [s1, 3], x2 in [0.6, 3], in both orders (stage 1 down first / stage 2
down first), with x1 >= s1 enforced as in mam.py.

Predictions, written before running:
  S1  Every start ends at W >= 0.148407329783 - 1e-8 (the same route), i.e. no second route at this resolution.
      This is the kill test for the "second route" suspect: ONE start below by more than 1e-7 kills W*'s globality
      (rule 19: a single counterexample decides it), and then the BVP is re-run from that start.
  My expectation: S1 holds, and the suspect dies -- which would leave the totals' growing excess slope unexplained.

Usage (from the repo root): HILL_K=1.2 PYTHONPATH=. python experiments/joint_action/mam_multistart.py
"""
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.argv = ["mam.py", "400"]
sys.path.insert(0, str(HERE))
os.chdir(HERE)
import numpy as np                                  # noqa: E402
from scipy.optimize import minimize                 # noqa: E402
import mam                                          # noqa: E402
from model import V1, S1                            # noqa: E402

N = 400
BASE = 0.148407329783


def path_through(w, first):
    """A -> corner -> waypoint -> B, smoothed; `first` = which stage drops first."""
    A, B = mam.XA, mam.XB
    c = np.array([w[0], A[1]]) if first == 1 else np.array([A[0], w[1]])
    pts = np.array([A, c, w, np.array([B[0], w[1]]) if first == 1 else np.array([w[0], B[1]]), B])
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(seg)]); s /= s[-1]
    q = 0.5 - 0.5 * np.cos(np.pi * np.linspace(0, 1, N + 1))
    X = np.column_stack([np.interp(q, s, pts[:, 0]), np.interp(q, s, pts[:, 1])])
    X[:, 0] = np.maximum(X[:, 0], S1 + 1e-6)
    return X


def run(X):
    mam.N = N
    dX = np.diff(X, axis=0)
    u0 = np.log(np.maximum(np.linalg.norm(dX, axis=1), 1e-12) / 0.5)
    z = np.concatenate([X[1:-1].ravel(), u0])
    bnds = [(S1 + 1e-9, 3.6), (0.05, 3.6)] * (N - 1) + [(-40, 12)] * N
    best = None
    for _ in range(8):
        r = minimize(mam.action, z, jac=True, method="L-BFGS-B", bounds=bnds,
                     options={"maxiter": 400000, "maxfun": 800000, "ftol": 1e-16, "gtol": 1e-13})
        z = r.x
        if best is not None and abs(best - r.fun) < 1e-13:
            break
        best = r.fun
    Xo = np.vstack([mam.XA, z[:2 * (N - 1)].reshape(N - 1, 2), mam.XB])
    return float(r.fun), Xo


if __name__ == "__main__":
    out = []
    for first in (1, 2):
        for x1w in (1.02, 1.15, 1.4, 1.8, 2.4):
            for x2w in (0.7, 1.2, 2.0):
                W, Xo = run(path_through(np.array([x1w, x2w]), first))
                k = np.argmin(Xo[:, 0])
                rec = {"first": first, "waypoint": [x1w, x2w], "W": W, "dW": W - BASE, "Gamma": float(V1(S1) - W),
                       "min_x1": float(Xo[k, 0]), "x2_at_min_x1": float(Xo[k, 1])}
                out.append(rec)
                print(f"first={first} w=({x1w},{x2w}) W={W:.12f} dW={W - BASE:+.2e} min_x1={Xo[k, 0]:.4f} "
                      f"at x2={Xo[k, 1]:.3f}", flush=True)
                (HERE.parents[1] / "results" / "mam_multistart_xk.json").write_text(json.dumps(out, indent=1))

# OUTCOME (§134.8): all 30 starts end on ONE route (min x1 = 1.2704 at x2 ~ 1.18-1.19), all above the BVP's converged
# W* = 0.1484037485 -- the "second route" suspect is dead within the searched region. But S1 AS WRITTEN would have
# printed "killed": 8 starts ended 1e-7..5.9e-7 below the X-path start -- the same route, MAM's N = 400 node placement
# scattering W by ~+-6e-7. The 1e-7 threshold sat inside the method's own noise (rule 19). The deciding test is
# "a different path, or a W below the converged W*"; that is the one recorded.
