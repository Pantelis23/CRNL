"""Independent cross-check of W*: direct minimum-action (MAM) with free segment times.
Path from the joint stable point (3,3) to the joint saddle (3, 0.606) [= inf over Sigma of W_c,
since W_c is non-increasing along the deterministic flow and Sigma flows into that saddle],
x1 >= s1 enforced by bounds. Action = sum_i tau_i L(x_mid, dx_i/tau_i), tau_i = exp(u_i).
Analytic gradient; L-BFGS-B.

Prediction (before running): agrees with ray value W* ~ 0.154054 from above, discretisation
error decreasing with N; path min x1 ~ 1.13, not at s1.
Usage: mam.py N [init_path.npy]
"""
import sys
import numpy as np
from scipy.optimize import minimize
from model import *

N = int(sys.argv[1])
XA = np.array([3.0, 3.0]); XB = np.array([3.0, R2_2])


def lm(x):
    x1, x2 = x[..., 0], x[..., 1]
    l = np.stack([lam1(x1), lam2(x1, x2)], -1)
    m = np.stack([mu1(x1), mu2(x2)], -1)
    # partials: dl[...,j,k] = d l_j / d x_k
    dl = np.zeros(x.shape + (2,)); dm = np.zeros(x.shape + (2,))
    dl[..., 0, 0] = dlam1(x1)
    dl[..., 1, 0] = dhill(x1) * A2 * x2 * x2
    dl[..., 1, 1] = hill(x1) * A2 * 2 * x2
    dm[..., 0, 0] = dmu1(x1)
    dm[..., 1, 1] = 3 * x2 * x2 + C2
    return l, m, dl, dm


def pstar(v, l, m):
    s = np.sqrt(v * v + 4 * l * m)
    return np.where(v >= 0, np.log((v + s) / (2 * l)), np.log(2 * m / (s - v))), s


def action(z, grad=True):
    X = np.vstack([XA, z[:2 * (N - 1)].reshape(N - 1, 2), XB])
    u = z[2 * (N - 1):]
    tau = np.exp(u)
    dX = np.diff(X, axis=0)
    xm = 0.5 * (X[1:] + X[:-1])
    v = dX / tau[:, None]
    l, m, dl, dm = lm(xm)
    p, s = pstar(v, l, m)
    L = (v * p - s + l + m).sum(1)
    S = (tau * L).sum()
    if not grad:
        return S
    Hm = (l * np.expm1(p) + m * np.expm1(-p)).sum(1)
    dLdl = -np.expm1(p); dLdm = -np.expm1(-p)
    dLdx = np.einsum('ij,ijk->ik', dLdl, dl) + np.einsum('ij,ijk->ik', dLdm, dm)
    gX = np.zeros_like(X)
    gX[1:] += p + 0.5 * tau[:, None] * dLdx
    gX[:-1] += -p + 0.5 * tau[:, None] * dLdx
    gu = tau * (L - (p * v).sum(1))   # = -tau * H
    return S, np.concatenate([gX[1:-1].ravel(), gu])


def init_from(path_file):
    P = np.load(path_file)            # rows: t, x1, x2, p1, p2, S
    xy = P[1:3].T
    xy = np.vstack([XA, xy, XB])
    # add the Sigma-hugging tail implicitly by straight segment to XB; reparametrise by
    # a monitor that clusters nodes near the fixed points
    seg = np.linalg.norm(np.diff(xy, axis=0), axis=1)
    sarc = np.concatenate([[0], np.cumsum(seg)]); sarc /= sarc[-1]
    q = 0.5 - 0.5 * np.cos(np.pi * np.linspace(0, 1, N + 1))
    X = np.column_stack([np.interp(q, sarc, xy[:, 0]), np.interp(q, sarc, xy[:, 1])])
    return X


if __name__ == '__main__':
    X = init_from(sys.argv[2])
    # initial times: deterministic-speed scale
    dX = np.diff(X, axis=0)
    u0 = np.log(np.maximum(np.linalg.norm(dX, axis=1), 1e-12) / 0.5)
    z0 = np.concatenate([X[1:-1].ravel(), u0])
    bnds = [(S1 + 1e-9, 3.6), (0.05, 3.6)] * (N - 1) + [(-40, 12)] * N
    best = None
    for rnd in range(6):
        r = minimize(action, z0, jac=True, method='L-BFGS-B', bounds=bnds,
                     options={'maxiter': 200000, 'maxfun': 400000, 'ftol': 1e-16, 'gtol': 1e-12})
        z0 = r.x
        print(f'  round {rnd}: S={r.fun:.10f} nit={r.nit} {r.message}', flush=True)
        if best is not None and abs(best - r.fun) < 1e-12:
            break
        best = r.fun
    Xo = np.vstack([XA, r.x[:2 * (N - 1)].reshape(N - 1, 2), XB])
    np.save(f'mam{TAG}_N{N}.npy', Xo)
    print(f'N={N} W*_MAM={r.fun:.10f} Gamma={V1(S1) - r.fun:.6e} min x1={Xo[:, 0].min():.5f}'
          f' (at x2={Xo[np.argmin(Xo[:, 0]), 1]:.4f}); total time={np.exp(r.x[2*(N-1):]).sum():.1f}')
