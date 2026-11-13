"""§134.7 kill test: XK's Gamma by direct minimum action (mam.py's action), converged in N, independent of the rays.

Why. XK's Gamma = 7.340e-3 (ray shooting; EPS 1e-6/1e-7/1e-8 spread 8e-7) was never checked by MAM -- X's was, to
1.3e-6 (results/joint_action.json). With Omega = 5000 in, the exact totals 800..5000 are fitted 20-40x better by a
CONSTANT prefactor exponent a ~ -0.02 and Gamma ~ 7.345e-3 than by Gamma = 7.340e-3 and a rising a. A Gamma error of
5e-6 alone produces the "a rising through zero" of §134.2-§134.5, because the measured a shifts by dGamma * Omega_eff.

How. mam.py's discretised action (free segment times, L-BFGS-B) from (3,3) to the joint saddle (3, 0.606), HILL_K = 1.2,
N = 100, 200, 400, 800, 1600, each N started from the previous path resampled (the first from X's ray path). W_MAM
approaches W* from above in N (X: 0.15407681 / 0.15406076 / 0.15405679 / 0.15405548 -> 0.15405354), so
Gamma_MAM(N) = V1(s1) - W_MAM(N) approaches Gamma from below.

Predictions, written before running:
  G1  Gamma_MAM(N) rises monotonically with N (as X's did).
  G2  The verdict is read from the sequence, not from a tolerance (rule 20), and the bound direction is used as a bound:
        - if any Gamma_MAM(N) exceeds 7.3410e-3 (the ray value + its stated 1e-6), the ray Gamma is too LOW and the
          free-fit reading (Gamma ~ 7.345e-3, a constant) becomes the leading one. One such N settles it (rule 19);
        - if the sequence's increments shrink geometrically and every candidate N -> inf extrapolation (Richardson in
          1/N and 1/N^2; rule 15: all reported) lands within 1e-6 of 7.340e-3, the ray value stands and the free-fit
          preference must come from curvature in the prefactor, not from Gamma;
        - otherwise UNRESOLVED, and written as such.
  My expectation, stated so it can be wrong: Gamma_MAM ends above 7.341e-3 (the data's pull is 11-17x the ray spread
  and in the direction a ray that stops short of the saddle would err).

Usage (from the repo root): HILL_K=1.2 PYTHONPATH=. python experiments/joint_action/mam_converge.py 100,200,400,800,1600
"""
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
NS = [int(v) for v in sys.argv[1].split(",")]
sys.argv = ["mam.py", str(NS[0])]
sys.path.insert(0, str(HERE))
os.chdir(HERE)
import numpy as np                                  # noqa: E402
from scipy.optimize import minimize                 # noqa: E402
import mam                                          # noqa: E402
from model import V1, S1                            # noqa: E402


def resample(X, N):
    seg = np.linalg.norm(np.diff(X, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(seg)]); s /= s[-1]
    q = 0.5 - 0.5 * np.cos(np.pi * np.linspace(0, 1, N + 1))
    return np.column_stack([np.interp(q, s, X[:, 0]), np.interp(q, s, X[:, 1])])


def run(N, X):
    mam.N = N
    dX = np.diff(X, axis=0)
    u0 = np.log(np.maximum(np.linalg.norm(dX, axis=1), 1e-12) / 0.5)
    z0 = np.concatenate([X[1:-1].ravel(), u0])
    bnds = [(S1 + 1e-9, 3.6), (0.05, 3.6)] * (N - 1) + [(-40, 12)] * N
    best = None
    for rnd in range(8):
        r = minimize(mam.action, z0, jac=True, method="L-BFGS-B", bounds=bnds,
                     options={"maxiter": 400000, "maxfun": 800000, "ftol": 1e-16, "gtol": 1e-13})
        z0 = r.x
        if best is not None and abs(best - r.fun) < 1e-13:
            break
        best = r.fun
    Xo = np.vstack([mam.XA, r.x[:2 * (N - 1)].reshape(N - 1, 2), mam.XB])
    return float(r.fun), Xo, rnd


if __name__ == "__main__":
    X = mam.init_from("optpath_eps1e-7_rtol1e-13.npy")
    out = {"HILL_K": os.environ.get("HILL_K"), "V1_s1": float(V1(S1)), "cells": {}}
    for N in NS:
        t = time.time()
        W, X, rounds = run(N, resample(X, N))
        G = float(V1(S1) - W)
        out["cells"][str(N)] = {"W": W, "Gamma": G, "min_x1": float(X[:, 0].min()), "rounds": rounds + 1,
                                "seconds": time.time() - t}
        print(f"N={N} W={W:.12f} Gamma={G:.9e} min_x1={X[:, 0].min():.5f} rounds={rounds + 1} "
              f"({time.time() - t:.0f} s)", flush=True)
        Path(HERE.parents[1] / "results" / "mam_xk.json").write_text(json.dumps(out, indent=1))

# OUTCOME (§134.7): G1 HELD -- Gamma_MAM rose at every N: 7.311018 / 7.329363 / 7.334580 / 7.335956 / 7.336251 /
# 7.336740 / 7.336913 e-3 for N = 100 .. 6400. G2 UNRESOLVED: no N exceeded 7.3410e-3, so the rays are not refuted, but
# the increments are not geometric (1600 -> 3200 exceeded 800 -> 1600), so the optimiser was not fully converged at
# some N; the extrapolations are 7.33697-7.33708e-3, not within 1e-6 of 7.340e-3. My stated expectation (above
# 7.341e-3) was WRONG: neither method supports the free fit's 7.345e-3. N = 6400 took 5.4 h.
