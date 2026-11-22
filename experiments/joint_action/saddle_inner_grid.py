"""§134.19: re-verify §134.17 (beta = 0 from the linear saddle) by an INDEPENDENT solve (rule 14: a withdrawal is a claim).

§134.17 integrated closed-form Gaussian Green's functions over time. Here the same problem -- dX = B X dt + dW,
B = [[1, 0], [c, -mu]], D = I/2 per unit (generator (1/2) Laplacian), point source at x0 = rho e_in -- is solved as the
stationary Fokker-Planck equation on a grid: finite volumes with Scharfetter-Gummel fluxes (stable at any Peclet
number), absorbing outer boundary far away, sparse direct solve. Nothing is shared with saddle_inner.py except (mu, c,
e_in) and the source position.

Predictions, written before running (XK's mu, c; rho = 30, grid h = 0.05 and 0.025):
  V1 (absolute, rule 16) v_grid / v_gauss at y = 0 and y = r e_in (r = 3, 5, 10) lies within 2% of 1 on the finer grid,
     and moves toward 1 from h = 0.05 to 0.025 (rule 20: convergence, not a fixed gate).
  V2 R_grid(r) = [ln v(0) + S(0)] - [ln v(r e_in) + S(r e_in)] has slope 2 beta <= 0.02 over r = 3-10 (S from the Gaussian
     action, the same for both). A slope near 0.74 (beta = 0.37) would withdraw §134.17 and return suspect 5.
  My expectation: both hold -- §134.17 stands.
"""
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import saddle_inner as si                                       # noqa: E402  (only mu, c, e_in, and v_gauss for comparison)


def bern(z):
    z = np.asarray(z, float)
    out = np.empty_like(z)
    small = np.abs(z) < 1e-6
    out[small] = 1 - z[small] / 2
    out[~small] = z[~small] / np.expm1(z[~small])
    return out


def solve_grid(mu, c, x0, h, L1=45.0, L2=10.0):
    x1 = np.arange(-L1, L1 + h / 2, h); x2 = np.arange(-L2, L2 + h / 2, h)
    n1, n2 = len(x1), len(x2)
    idx = lambda i, j: i * n2 + j
    D = 0.5
    rows, cols, vals = [], [], []
    diag = np.zeros(n1 * n2)
    I, J = np.meshgrid(np.arange(n1), np.arange(n2), indexing="ij")
    I, J = I.ravel(), J.ravel()
    # x1-faces between (i, j) and (i+1, j): drift a1 = x1 at the face
    for di, dj, axis in ((1, 0, 1), (0, 1, 2)):
        ok = (I + di < n1) & (J + dj < n2)
        i, j = I[ok], J[ok]
        if axis == 1:
            a = 0.5 * (x1[i] + x1[i + 1])
        else:
            xm2 = 0.5 * (x2[j] + x2[j + 1])
            a = c * x1[i] - mu * xm2
        pe = a * h / D
        bm, bp = bern(-pe), bern(pe)
        # flux from p=(i,j) to q=(i+di,j+dj): F = (D/h) [B(-pe) v_p - B(pe) v_q]; divergence /h
        p, q = idx(i, j), idx(i + di, j + dj)
        k = D / h / h
        # node p loses F, node q gains F
        rows += [p, p, q, q]; cols += [p, q, p, q]
        vals += [-k * bm, k * bp, k * bm, -k * bp]
    # absorbing boundary: flux out of the outermost faces to v = 0 ghost nodes
    for i_b, j_b, a_fn in ():
        pass
    A = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n1 * n2, n1 * n2)).tocsr()
    # boundary nodes: impose v = 0 (absorbing) by overwriting their rows with the identity
    bd = (I == 0) | (I == n1 - 1) | (J == 0) | (J == n2 - 1)
    A = A.tolil()
    for p in np.where(bd)[0]:
        A.rows[p] = [p]; A.data[p] = [1.0]
    A = A.tocsc()
    b = np.zeros(n1 * n2)
    i0, j0 = int(round((x0[0] + L1) / h)), int(round((x0[1] + L2) / h))
    b[idx(i0, j0)] = -1.0 / (h * h)                              # unit point source (A v + delta = 0)
    v = spla.spsolve(A, b)
    return x1, x2, v.reshape(n1, n2), (x1[i0], x2[j0])


def interp(x1, x2, V, y):
    from scipy.interpolate import RegularGridInterpolator
    return float(RegularGridInterpolator((x1, x2), np.log(np.maximum(V, 1e-300)))(np.atleast_2d(y))[0])


if __name__ == "__main__":
    mu, c, _ = si.normal_form()
    B = np.array([[1.0, 0.0], [c, -mu]])
    e_in = si.incoming_direction(B)
    rho = 30.0
    out = {"mu": mu, "c": c, "rho": rho, "cells": {}}
    for h in (0.05, 0.025):
        x1, x2, V, xs = solve_grid(mu, c, rho * e_in, h)
        x0g = np.array(xs)                                       # the source node actually used
        ts = np.concatenate([np.geomspace(1e-4, 1.0, 400), np.linspace(1.0, 2.0 * math.log(rho) + 25, 20000)[1:]])
        pts = {"0": np.zeros(2)} | {f"r{r}": r * e_in for r in (3, 5, 10)}
        rec = {}
        for k, y in pts.items():
            lg = interp(x1, x2, V, y)
            lgs, S = si.ln_v_and_S(B, y, x0g, ts)
            rec[k] = {"ln_v_grid": lg, "ln_v_gauss": lgs, "S": S, "ratio": math.exp(lg - lgs)}
        R = {k: (rec["0"]["ln_v_grid"] + rec["0"]["S"]) - (rec[k]["ln_v_grid"] + rec[k]["S"]) for k in ("r3", "r5", "r10")}
        slope = (R["r10"] - R["r3"]) / math.log(10 / 3)
        out["cells"][f"{h:g}"] = {"points": rec, "R": R, "two_beta_3_10": slope}
        print(f"h={h}: ratios grid/gauss " + " ".join(f"{k}:{v['ratio']:.4f}" for k, v in rec.items())
              + f"  2beta(3-10) grid = {slope:+.4f}", flush=True)
    (HERE.parents[1] / "results" / "saddle_inner_grid_xk.json").write_text(json.dumps(out, indent=1))

# OUTCOME (§134.19, first design): V1 and V2 were INFEASIBLE AS REGISTERED. With rho = 30 the density at H is e^{-900}
# (Gaussian: ln v(0) = -900.8), far below double range; a direct sparse solve resolves only ~1e-16 of its peak, so the
# grid's v(0) ~ 1e-76 was round-off (212,625 negative values). The design could only print noise (rule 19).
#
# REVISED BEFORE RUNNING (§134.19 revision): validate the Gaussian computation itself, which is exact for any rho, so a
# check at small rho carries to rho = 1e6.
#   W1 normalisation: int G_t(y) dy = 1 to <= 1e-6 at t = 0.3, 1, 3, 10 (quadrature on a box around the mean).
#   W2 Monte Carlo: dX = B X dt + dW (Euler-Maruyama, dt = 1e-3 and 5e-4), source rho = 2 along e_in, paths stopped at
#      |X| > 12; occupation density in boxes of side 0.2 at y = 0 and y = 1.0 e_in. MC / Gaussian within 3 MC sigma of 1
#      at both points and both dt -- or its trend in dt points to 1.
#   Either failing withdraws §134.17's beta = 0 as unverified.
