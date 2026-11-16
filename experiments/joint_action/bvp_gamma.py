"""§134.8 (T-CASC-av's kill test): W* as a heteroclinic boundary-value problem, a third method independent of the rays
and of MAM, whose numerical parameters can each be swept to convergence.

The minimum-action path from the joint stable point A = (3, 3) to the joint saddle B = (3, R2_2) is the heteroclinic
orbit of the WKB Hamiltonian flow from (A, p = 0) to (B, p = 0). Solved on [0, T] with projection boundary conditions
(Beyn): y(0) - A in the unstable subspace of the linearised flow at A, y(T) - B in the stable subspace at B. The
action carries its own quadratic end corrections, so it is exact to O(|delta|^3) at both ends:
    W* = S(T) - (1/2) (x(T) - B) . p(T),   with S(0) = (1/2) (x(0) - A) . p(0),   S' = p . x'.
(Near a hyperbolic point W = W(pt) + (1/2) dx^T Z dx along its Lagrangian manifold, and p = Z dx there.)
Numerical axes: the time span T (sets the end offsets) and solve_bvp's tol. The initial guess comes from mam.py's
path at N = 400.

Predictions, written before running:
  B1 (validation, X: HILL_K = 1): W* converges in T and tol, and agrees with the rays' 0.15405354 to within the rays'
     own stated Gamma_err (2e-6). If it does not, the BVP is not trusted for XK and nothing is read from it.
  B2 (XK: HILL_K = 1.2): the converged Gamma, read from the sequence over T and tol (rule 20), decides T-CASC-av:
        Gamma >= 7.344e-3              -> a constant a ~ -0.02, and §134's "a rising" was Gamma (rule 21);
        7.3369e-3 <= Gamma <= 7.3411e-3 -> a really rises, and M7's "refuted" stands;
        and, separately, M6 (§134.5) held iff Gamma >= 7.3388e-3.
     If the sequence over T or tol has not settled to well below these gaps, the verdict is UNRESOLVED.
  My expectation, stated so it can be wrong: XK's Gamma lands within 1e-6 of the rays' 7.3402e-3, MAM being slow
  from below (as §134.7's increments suggest), not wrong.

Usage (from this directory): HILL_K=1.2 python bvp_gamma.py
"""
import os
import sys
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp

HERE = Path(__file__).resolve().parent
os.chdir(HERE)
sys.argv = ["mam.py", "400"]
sys.path.insert(0, str(HERE))
import mam                                                # noqa: E402
from model import V1, S1, R2_2, ham_rhs, KH               # noqa: E402
from scipy.optimize import minimize                       # noqa: E402

A = np.array([3.0, 3.0, 0.0, 0.0])
B = np.array([3.0, R2_2, 0.0, 0.0])


def f4(y):
    return np.array(ham_rhs(0.0, np.concatenate([y, np.zeros((1,) + y.shape[1:])]))[:4])


def jac(pt, h=1e-6):
    J = np.zeros((4, 4))
    for k in range(4):
        e = np.zeros(4); e[k] = h
        J[:, k] = (f4((pt + e)[:, None])[:, 0] - f4((pt - e)[:, None])[:, 0]) / (2 * h)
    return J


def projector_rows(J, keep):
    """Rows of the left eigenvectors whose eigenvalues are NOT in `keep` ('unstable' / 'stable')."""
    w, V = np.linalg.eig(J)
    W = np.linalg.inv(V)                                    # rows: left eigenvectors
    drop = w.real < 0 if keep == "unstable" else w.real > 0
    return np.real_if_close(W[drop]).real


def schur_complement_rows(J, keep):
    """Rows spanning the orthogonal complement of J's unstable ("rhp") or stable ("lhp") invariant subspace, from a
    sorted real Schur form -- stable when J is defective (pair S: identical stages make A's Jacobian a Jordan block,
    where eig's eigenvectors are near-parallel and their inverse is garbage; §134.15)."""
    from scipy.linalg import schur
    T, Z, k = schur(J, output="real", sort=keep)
    return Z[:, k:].T


LA = schur_complement_rows(jac(A), "rhp")                   # (y0 - A) lies in A's unstable subspace
LB = projector_rows(jac(B), "stable")                       # (yT - B) has no unstable-mode component


def mam_guess(N=400):
    mam.N = N
    X = mam.init_from("optpath_eps1e-7_rtol1e-13.npy")
    dX = np.diff(X, axis=0)
    u0 = np.log(np.maximum(np.linalg.norm(dX, axis=1), 1e-12) / 0.5)
    z = np.concatenate([X[1:-1].ravel(), u0])
    bnds = [(S1 + 1e-9, 3.6), (0.05, 3.6)] * (N - 1) + [(-40, 12)] * N
    for _ in range(3):
        z = minimize(mam.action, z, jac=True, method="L-BFGS-B", bounds=bnds,
                     options={"maxiter": 400000, "maxfun": 800000, "ftol": 1e-16, "gtol": 1e-13}).x
    X = np.vstack([mam.XA, z[:2 * (N - 1)].reshape(N - 1, 2), mam.XB])
    tau = np.exp(z[2 * (N - 1):])
    xm = 0.5 * (X[1:] + X[:-1]); v = np.diff(X, axis=0) / tau[:, None]
    l, m, _, _ = mam.lm(xm)
    p, _ = mam.pstar(v, l, m)
    t = np.concatenate([[0], np.cumsum(tau)])
    tm = 0.5 * (t[1:] + t[:-1])
    return tm, np.column_stack([xm, p])                     # (t, [x1, x2, p1, p2]) at segment midpoints


def trim(t, Y, r):
    """Keep the stretch of the guess farther than r from both end points (in x)."""
    dA = np.linalg.norm(Y[:, :2] - A[:2], axis=1); dB = np.linalg.norm(Y[:, :2] - B[:2], axis=1)
    k = (dA > r) & (dB > r)
    i0, i1 = np.argmax(k), len(k) - np.argmax(k[::-1])
    return t[i0:i1] - t[i0], Y[i0:i1]


DROP = int(os.environ.get("BVP_DROP", "0"))           # which of B's two unstable-mode rows energy makes redundant


def solve(t, Y, tol, epsA):
    """Fixed T; phase condition |x(0) - A| = epsA (removes the time-translation near-degeneracy); one of B's two
    projection rows is dropped (H = 0 makes it redundant) and its residual is returned as a check."""
    keep = [k for k in range(len(LB)) if k != DROP]

    def fun(_, y):
        d = f4(y[:4])
        return np.vstack([d, y[2] * d[0] + y[3] * d[1]])

    def bc(ya, yb):
        return np.concatenate([LA @ (ya[:4] - A), LB[keep] @ (yb[:4] - B),
                               [(ya[:2] - A[:2]) @ (ya[:2] - A[:2]) - epsA ** 2,
                                ya[4] - 0.5 * (ya[:2] - A[:2]) @ ya[2:4]]])

    S = np.concatenate([[0.0], np.cumsum(0.5 * (Y[1:, 2] + Y[:-1, 2]) * np.diff(Y[:, 0])
                                         + 0.5 * (Y[1:, 3] + Y[:-1, 3]) * np.diff(Y[:, 1]))])
    S += 0.5 * (Y[0, :2] - A[:2]) @ Y[0, 2:4]
    r = solve_bvp(fun, bc, t, np.vstack([Y.T, S]), tol=tol, max_nodes=400_000, bc_tol=tol * 1e-2)
    yT = r.sol(t[-1]); y0e = r.sol(0.0)
    W = yT[4] - 0.5 * (yT[:2] - B[:2]) @ yT[2:4]
    dropped = float(abs(LB[DROP] @ (yT[:4] - B))) / max(np.linalg.norm(yT[:4] - B), 1e-300)
    return r, float(W), float(np.linalg.norm(y0e[:2] - A[:2])), float(np.linalg.norm(yT[:2] - B[:2])), dropped


def _eig(J, which):
    w, V = np.linalg.eig(J)
    k = w.real > 0 if which == "unstable" else w.real < 0
    return w[k], V[:, k]


WA, VA = _eig(jac(A), "unstable")
WB, VB = _eig(jac(B), "stable")


def extend(r, dA, dB, n=200):
    """Continuation: extend the start by dA (backwards) along A's unstable modes and the end by dB along B's stable
    modes -- the linearised flow, exact to first order near each hyperbolic point."""
    t = np.linspace(0, r.x[-1], 3001)                       # resample: solve_bvp never removes nodes
    Y = r.sol(t)[:4].T
    from scipy.linalg import expm
    JA = jac(A)
    cB = np.linalg.lstsq(VB, (Y[-1] - B).astype(complex), rcond=None)[0]
    sa = np.linspace(-dA, 0, n + 1)[:-1] if dA > 0 else np.zeros(0)
    sb = np.linspace(0, dB, n + 1)[1:] if dB > 0 else np.zeros(0)
    Ya = np.array([A + expm(JA * v) @ (Y[0] - A) for v in sa]).reshape(-1, 4)   # linear flow, defective-safe
    Yb = (B[:, None] + VB @ (cB[:, None] * np.exp(WB[:, None] * sb[None, :]))).real.T
    tn = np.concatenate([sa, t, t[-1] + sb]) + dA
    Yn = np.vstack([Ya, Y, Yb])
    return tn, Yn, float(np.linalg.norm(Yn[0, :2] - A[:2]))


if __name__ == "__main__":
    tg, Yg = mam_guess()
    out = {"HILL_K": KH, "V1_s1": float(V1(S1)), "drop": DROP, "cells": []}
    t, Y = trim(tg, Yg, 3e-2)
    tol = float(os.environ.get("BVP_TOL", "1e-8"))
    dA, dB = float(os.environ.get("BVP_DA", "0.2")), float(os.environ.get("BVP_DB", "1.0"))
    epsA = float(np.linalg.norm(Y[0, :2] - A[:2]))
    for step in range(60):
        r, W, oa, ob, dr = solve(t, Y, tol, epsA)
        G = float(V1(S1) - W)
        out["cells"].append({"T": float(r.x[-1]), "tol": tol, "status": int(r.status), "nodes": len(r.x), "W": W,
                             "Gamma": G, "offset_A": oa, "offset_B": ob, "dropped_row_residual": dr})
        print(f"T={r.x[-1]:7.2f} st={r.status} nodes={len(r.x):6d} W={W:.13f} Gamma={G:.11e} "
              f"offA={oa:.1e} offB={ob:.1e} dropped-row={dr:.1e}", flush=True)
        if r.status != 0 or max(oa, ob) < 3e-7:
            break
        t, Y, epsA = extend(r, dA if oa > 3e-7 else 0.0, dB if ob > 3e-7 else 0.0)
    # the tolerance axis (rule 13), at the final T and epsA: re-solve the converged orbit, never re-continue
    t = np.linspace(0, r.x[-1], 6001); Y = r.sol(t)[:4].T
    for tt in (1e-7, 1e-8, 1e-9, 1e-10):
        r2, W2, oa2, ob2, dr2 = solve(t, Y, tt, epsA)
        out["cells"].append({"T": float(r2.x[-1]), "tol": tt, "status": int(r2.status), "nodes": len(r2.x), "W": W2,
                             "Gamma": float(V1(S1) - W2), "offset_A": oa2, "offset_B": ob2,
                             "dropped_row_residual": dr2, "axis": "tol"})
        print(f"  tol axis: tol={tt:g} st={r2.status} nodes={len(r2.x)} W={W2:.13f} Gamma={V1(S1) - W2:.11e} "
              f"offB={ob2:.1e}", flush=True)
    name = os.environ.get("BVP_NAME") or (("bvp_x" if KH == 1.0 else "bvp_xk") + f"_tol{tol:g}_drop{DROP}.json")
    (HERE.parents[1] / "results" / name).write_text(json.dumps(out, indent=1))

# OUTCOME (§134.8). Development: fixed-T projection BCs alone went near-singular as the end offsets shrank (time
# translation is only weakly broken), so a phase condition |x(0) - A| = epsA was added and one of B's two projection
# rows dropped (H = 0 makes one redundant). Dropping row 0 (B's deterministic +1.09 mode) converges; dropping row 1
# diverged. Continuation pads both ends with the linearised flow.
# B1 HELD (X): W* = 0.1540536705 (Gamma 1.68823963e-3), steady to ~1e-11 over offsets 1e-3 -> 1e-7 and tol 1e-7 -> 1e-10;
#    the rays' 1.688372e-3 is 1.3e-7 above, inside their stated 2e-6.
# B2 (XK): Gamma = 7.3381616e-3, steady to ~1e-11 on all three axes. It lies in [7.3369, 7.3411]e-3: "a really rises",
#    M7's "refuted" stands, and M6 FAILS (needs Gamma >= 7.3388e-3).
# My expectation (within 1e-6 of the rays' 7.3402e-3) was WRONG: the rays are 2.0e-6 high. Every real ray hit on Sigma
# has S >= 0.1484127596 (above the BVP's W, as an upper bound must be); only the rays' |p| -> 0 extrapolation
# (0.1484025122, linear in |p|^2 from |p| >= 1.9e-3) falls below it.
