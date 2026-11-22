"""§134.20 (T-CASC-ba, candidate 2's kill test): is the saddle excess made by NONLINEARITY inside the core?

In the exact chain (XK, cap 1.1, Q-process law P), take the box |x1 - 3| <= 0.3, |x2 - 0.606| <= 0.05 around the joint
saddle H. Keep the exact P on the box boundary. Inside, replace every jump rate (the Doob-transformed stage-1 rates
lt1, mt1 and stage-2's lam2, mu2) by its LINEARISATION at H (value and gradient at H's lattice point), and solve the
forward equation for the interior. Same lattice, same jumps -- only the rates' nonlinearity inside the box is removed.

Predictions, written before running (Omega = 400, 800, 1600):
  N1  r(Omega) = P_lin(H) / P(H).
        |ln r| <= 0.05 at every Omega, with no trend          -> nonlinearity inside the core is NOT the cause:
                                                                 candidate 2 DEAD; the excess enters through the
                                                                 box boundary (candidate 1's territory).
        ln r falling with Omega at ~ -0.37 per ln Omega         -> candidate 2 HOLDS: the excess is made inside by the
                                                                 rates' nonlinearity.
        otherwise UNRESOLVED.
  N2  (instrument) with the EXACT rates inside, the box solve must return P(H) to <= 1e-8 relative (the same chain
      restricted, so it must reproduce itself).
  My expectation: N1's first branch -- the excess comes in with the boundary data.

Usage: CRNL_CAP_MULT=1.1 PYTHONPATH=. python experiments/inner_linear_box.py 400,800,1600 kr
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import experiments.episode_committor as ec
from experiments.sigma_flux import ARMS, element


def box_solve(P, rate_fn, i_lo, i_hi, j_lo, j_hi):
    """Forward equation inside the box, P fixed on its boundary. rate_fn(i, j) -> (lt1, mt1, lam2, mu2) arrays for
    (i, j) index arrays (rates OUT of (i, j))."""
    ii = np.arange(i_lo + 1, i_hi); jj = np.arange(j_lo + 1, j_hi)
    I, J = np.meshgrid(ii, jj, indexing="ij"); I, J = I.ravel(), J.ravel()
    n = len(I)
    pos = -np.ones(P.shape, dtype=np.int64)
    pos[I, J] = np.arange(n)
    rows, cols, vals = [], [], []
    rhs = np.zeros(n)
    lt, mt, l2, m2 = rate_fn(I, J)
    out = lt + mt + l2 + m2
    rows.append(np.arange(n)); cols.append(np.arange(n)); vals.append(-out)
    # inflow into (I, J) from its four neighbours, each with the neighbour's own out-rate toward (I, J)
    for di, dj, k in ((-1, 0, 0), (1, 0, 1), (0, -1, 2), (0, 1, 3)):
        Ni, Nj = I + di, J + dj
        r = rate_fn(Ni, Nj)[k]                                 # neighbour's rate toward (I, J): up1 / dn1 / up2 / dn2
        q = pos[Ni, Nj]
        inside = q >= 0
        rows.append(np.nonzero(inside)[0]); cols.append(q[inside]); vals.append(r[inside])
        rhs[~inside] -= r[~inside] * P[Ni[~inside], Nj[~inside]]
    M = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsc()
    x = spla.spsolve(M, rhs)
    Pb = P.copy()
    Pb[I, J] = x
    return Pb


if __name__ == "__main__":
    if len(sys.argv) > 2:
        ec.SOLVER = sys.argv[2]
    A, B = element(ARMS["XK"][0]), element(ARMS["XK"][1])
    out_p = Path("results/inner_linear_box.json")
    out = {}
    for om in [int(v) for v in sys.argv[1].split(",")]:
        m = A.cap(om) + 1
        sad1, sad2 = int(np.ceil(A.r2 * om)), int(np.ceil(B.r2 * om))
        ec.qprocess_J(A, B, om, m, sad1, sad2, lower="rail")
        L = ec.qprocess_J.last
        P = L["P"] / ec.PIN_SCALE
        n1s, base2, M1, M2 = L["n1s"], L["base2"], L["M1"], L["M2"]
        lt1 = L["lt1"].reshape(M1, M2); mt1 = L["mt1"].reshape(M1, M2)
        n2s = base2 + np.arange(M2)

        def exact(I, J):
            l2, m2 = B.rates(n2s[J].astype(float), n1s[I].astype(float), om, False)
            return lt1[I, J], mt1[I, J], l2, m2

        iH = int(round(3.0 * om)) - int(n1s[0]); jH = int(round(B.r2 * om)) - base2
        f0 = [np.asarray(v, float) for v in exact(np.array([iH]), np.array([jH]))]
        g1 = [(np.asarray(a, float) - np.asarray(b, float)) / 2 for a, b in
              zip(exact(np.array([iH + 1]), np.array([jH])), exact(np.array([iH - 1]), np.array([jH])))]
        g2 = [(np.asarray(a, float) - np.asarray(b, float)) / 2 for a, b in
              zip(exact(np.array([iH]), np.array([jH + 1])), exact(np.array([iH]), np.array([jH - 1])))]

        def linear(I, J):
            return tuple(f0[k][0] + g1[k][0] * (I - iH) + g2[k][0] * (J - jH) for k in range(4))

        a1, a2 = int(round(0.3 * om)), int(round(0.05 * om))
        box = (iH - a1, iH + a1, jH - a2, jH + a2)
        Pe = box_solve(P, exact, *box)
        Pl = box_solve(P, linear, *box)
        rec = {"N2_exact_reproduction": float(abs(Pe[iH, jH] / P[iH, jH] - 1)),
               "ln_r": float(math.log(Pl[iH, jH] / P[iH, jH])), "box": [int(v) for v in box]}
        out[str(om)] = rec
        print(om, rec, flush=True)
        out_p.write_text(json.dumps(out, indent=1))

# OUTCOME (§134.20): N2 HELD (exact rates reproduce P(H) to 1e-13). N1 UNRESOLVED as registered: ln r = +0.057 / +0.125 /
# +0.207 at Omega = 400 / 800 / 1600 -- LINEAR rates make P(H) LARGER, by ~+0.10-0.12 per ln Omega: neither "no change"
# nor candidate 2's predicted fall. Nonlinearity inside the box SUPPRESSES part of the excess rather than making it.
# Caveat (rule 13): the box is fixed physically (0.3 in x1), so its stretched size grows as sqrt(Omega) and the
# linearisation of the cubic stage-1 rates is crude at its edges; box size is a second axis not varied here.
# My expectation (first branch) was WRONG.
#
# REGISTERED NEXT (§134.21, before running): where does the excess build up? L at H + delta e_n, e_n the normal
# (x2) direction, delta = r (Omega Z22)^{-1/2}, Z22 = 0.6598, r = 1, 3, 10 (above and below H), Omega = 400-1600.
#   If dL/dlnOmega at fixed r equals H's (+0.30-0.37) for all r <= 10 -> the excess lives on scales >> the inner zone.
#   If it falls with r toward 0 by r = 10                              -> it is built within a few tube widths of H.
