"""§134.17 (T-CASC-ay, exact part): beta(mu, c) from the LINEAR saddle alone -- an exact Gaussian computation, no fit.

Why local. Along the MPEP the WKB transport integrand -[tr B + (1/2) tr(D Z)] vanishes exactly at H (on the Lagrangian plane
spanned by the two noise modes, tr(B + D Z) = lambda_s - lambda_u), so the outer prefactor K stays finite; any power of
Omega at H must come from the O(Omega^{-1/2}) inner region, which is a linear problem.

Model. Maier-Stein normal form near H: dX = B X dt + dW, B = [[1, 0], [c, -mu]] (time in units of 1/lambda_u, D = I,
stretched coordinates X = x / eps^{1/2}, so eps has dropped out). Feed it from a point source far up the incoming MPEP,
x0 = rho * e_in, where e_in is the x-part of the Hamiltonian flow's noise-stable eigenvector (eigenvalue -1), the
direction the MPEP arrives along (§134.8: offset ~ e^{-lambda_u t}). The stationary density from the source is exact:
    v(y) = int_0^inf G_t(y | x0) dt,   G_t Gaussian with mean e^{Bt} x0 and covariance Sigma_t = int_0^t e^{Bs} e^{B^T s} ds.
Its WKB action from the source is S(y) = inf_t (1/2) (y - m_t)^T Sigma_t^{-1} (y - m_t). Define, along the incoming ray at
stretched distance r from H (1 << r << rho),
    R(r) = [ln v(0) + S(0)] - [ln v(r e_in) + S(r e_in)].
If the inner density at H carries Omega^beta relative to the outer K (stretched r ~ Omega^{1/2} x), then dR/d ln r -> 2 beta.

Predictions, written before running:
  I1  With XK's (mu, c) = (5.2083, +0.2013), beta from the plateau of dR/dlnr lies within +-0.02 of the measured 0.374
      (§134.14). With S's (mu = 3.333, c from ms_core's transform) it lies in S's measured band [0.32, 0.37).
      Both must hold; either failing means the linear inner problem does not by itself set beta (the outer data, or the
      jump-process corrections, matter).
  I2  (instrument) the plateau is flat to +-0.01 over a decade of r, and does not move when rho is raised 10x.
  My expectation, stated so it can be wrong: I1 holds for XK; S within its band.

Usage (from this directory): HILL_K=1.2 python saddle_inner.py  (and R2_STAGE2=1.005 HILL_K=1.0 for S)
"""
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from model import lam1, mu1, lam2, mu2, dlam1, dmu1, hill, dhill, A2, C2, R2_2   # noqa: E402


def normal_form():
    """(mu, c) exactly as ms_core.py computes them (D(H) = I, lambda_u = 1, boundary along e_s)."""
    x1, x2 = 3.0, R2_2
    Bm = np.array([[dlam1(x1) - dmu1(x1), 0.0],
                   [dhill(x1) * A2 * x2 * x2, hill(x1) * A2 * 2 * x2 - (3 * x2 * x2 + C2)]])
    Dh = np.sqrt(np.diag([lam1(x1) + mu1(x1), lam2(x1, x2) + mu2(x2)]))
    By = np.linalg.inv(Dh) @ Bm @ Dh
    w, V = np.linalg.eig(By)
    iu, is_ = int(np.argmax(w.real)), int(np.argmin(w.real))
    lu, ls = w[iu].real, -w[is_].real
    es = V[:, is_].real; es /= np.linalg.norm(es)
    n = np.array([es[1], -es[0]])
    Q = np.column_stack([n, es])
    Bq = Q.T @ By @ Q / lu
    return ls / lu, Bq[1, 0], Bq


def incoming_direction(B):
    """x-part of the Hamiltonian flow's eigenvector with eigenvalue -1 (H = p.p/2 + p.Bx): J = [[B, I], [0, -B^T]]."""
    J = np.block([[B, np.eye(2)], [np.zeros((2, 2)), -B.T]])
    w, V = np.linalg.eig(J)
    k = int(np.argmin(np.abs(w + 1.0)))
    e = V[:2, k].real
    return e / np.linalg.norm(e)


def green_terms(B, y, x0, ts):
    """ln G_t(y | x0) and (1/2) d^T Sigma^{-1} d, in closed form. B = [[1, 0], [c, -mu]]. In B's eigencoordinates
    Y1 = X1, Y2 = X2 - k X1 (k = c / (1 + mu)) the drift is diag(1, -mu) and the noise covariance is [[1, -k], [-k, 1 + k^2]],
    so Sigma_Y has closed-form entries and only Sigma_11 grows; the quadratic form goes through a Schur complement with no
    cancellation (a generic solve loses everything once Sigma's eigenvalues differ by ~1e19)."""
    mu, c = -B[1, 1], B[1, 0]
    k = c / (1 + mu)
    t = np.asarray(ts)
    s11 = np.expm1(2 * t) / 2
    s12 = -k * (-np.expm1((1 - mu) * t)) / (mu - 1)
    s22 = (1 + k * k) * (-np.expm1(-2 * mu * t)) / (2 * mu)
    S22 = s22 - s12 * s12 / s11
    toY = lambda v: np.array([v[0], v[1] - k * v[0]])
    yY, x0Y = toY(y), toY(x0)
    m1, m2 = np.exp(t) * x0Y[0], np.exp(-mu * t) * x0Y[1]
    d1, d2 = yY[0] - m1, yY[1] - m2
    act = 0.5 * (d1 * d1 / s11 + (d2 - s12 / s11 * d1) ** 2 / S22)
    lndet = np.log(s11) + np.log(S22)                # the Jacobian of the X -> Y map is 1
    return -math.log(2 * math.pi) - 0.5 * lndet - act, act


def ln_v_and_S(B, y, x0, ts):
    lnG, act = green_terms(B, y, x0, ts)
    dt = np.gradient(ts)
    m = lnG.max()
    return m + math.log(np.sum(np.exp(lnG - m) * dt)), act.min()


if __name__ == "__main__":
    mu, c, Bq = normal_form()
    B = np.array([[1.0, 0.0], [c, -mu]])
    e_in = incoming_direction(B)
    tag = os.environ.get("INNER_TAG", f"mu{mu:.3f}")
    out = {"mu": mu, "c": c, "e_in": e_in.tolist(), "runs": {}}
    print(f"mu={mu:.4f} c={c:+.4f} e_in={e_in.round(4)}", flush=True)
    for rho in (1e5, 1e6):
        x0 = rho * e_in
        ts = np.concatenate([np.geomspace(1e-4, 1.0, 400), np.linspace(1.0, 2.0 * math.log(rho) + 25, 20000)[1:]])
        l0, s0 = ln_v_and_S(B, np.zeros(2), x0, ts)
        rs = np.geomspace(2, 3e3, 25)
        R = []
        for r in rs:
            lr, sr = ln_v_and_S(B, r * e_in, x0, ts)
            R.append((l0 + s0) - (lr + sr))
        slope = np.gradient(np.array(R), np.log(rs))
        out["runs"][f"{rho:g}"] = {"r": rs.tolist(), "R": R, "slope": slope.tolist()}
        print(f"rho={rho:g}: beta(r) = " + " ".join(f"{r:.0f}:{s / 2:+.4f}" for r, s in zip(rs[::3], slope[::3])),
              flush=True)
    (HERE.parents[1] / "results" / f"saddle_inner_{tag}.json").write_text(json.dumps(out, indent=1))

# OUTCOME (§134.17): beta = 0 for BOTH pairs. Plateau over r = 5-200 at rho = 1e6: |beta| <= 0.0004 for XK
# (mu = 5.208, c = +0.201) and S (mu = 3.333, c = +0.114); the large-r tail is the source's own influence. I1 FAILED:
# the linear inner saddle is WKB-regular for mu > 1 (as Maier-Stein's O(1) prefactor says) and makes no Omega^0.37.
# My expectation was WRONG. Development: a generic solve lost all precision once Sigma's eigenvalues differed by ~1e19;
# the closed form in B's eigencoordinates (Schur complement) fixed it.
