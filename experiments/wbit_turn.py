"""§128 -- T-WBIT-a′'s kill test: do §126's local slopes turn back toward §127's Gamma below sigma = 0.25?

§127 found W-bit's exact joint-action exponents, Gamma_Z = +0.0125 and Gamma_N = -0.1066. §126's local slopes
of ln R against 1/sigma^2 are not approaching them: Z's (0.039 / 0.018 / 0.009) has fallen through Gamma_Z, and
N's (-0.158 / -0.195 / -0.217) moves away from Gamma_N. If Gamma is the asymptote, both must turn.

INSTRUMENT. §126's direct sparse LU cannot go past M = 160. The push-forward of a Gaussian through a logistic
covers about half the interval, so the 2-D kernel is about 25% dense (2.6e7 nonzeros at M = 100). The product
structure gives an O(M^3) matrix-vector product instead:
    (P v)[i, j] = sum_j' K2^(i)[j, j'] (K1 v)[i, j'].
The exit probabilities then come from GMRES on (I - P_aa) q = b, with the same dip / fall / alive sets and start
state as `wbit_cascade.exits`.

PREDICTIONS, WRITTEN BEFORE RUNNING.
  V1  INSTRUMENT. GMRES reproduces the stored direct-solve R of §126 (Z and N, sigma = 0.30 and 0.25, M = 160) to
      a relative 1e-8, with q_fall + q_dip = 1 to 1e-10.
  V2  GRID (rule 13). At each sigma, R's change between successive M (160 / 240 / 320) shrinks. A slope is read
      only if the grid change at the finest M is below a tenth of the change in ln R between the neighbouring
      sigma cells; otherwise that slope is written as unresolved.
  T   THE TURN (the verdict rule in THEORIES T-WBIT-a′). Every sigma (0.34, 0.30, 0.27, 0.25, 0.22, 0.20, 0.18) is
      run at the same M, so the slopes are paired (rule 18).
      - If the large-deviation form holds, N's slope sequence reverses: some slope is less negative than the one
        before it. Z's also reverses: some slope is larger than the one before it.
      - If N's slope passes -0.267, or Z's goes negative, W* is not the exponent of §126's R.
      - A sequence that neither turns nor crosses those bounds is UNRESOLVED at this sigma.
      My expectation, stated so that it can be wrong: N turns within the range, and Z's slope goes below zero
      before it turns. Z's R would then pass a maximum.

ADDED AFTER THE SIGMA = 0.34-0.18 RUNS (M = 160-640), BEFORE THE EXTENSION RAN. Z turned: its slopes run
0.0394 / 0.0175 / 0.0085 / 0.0066 / 0.0100 / 0.0147, already above Gamma_Z. N did not turn: -0.159 / -0.196 /
-0.217 / -0.233 / -0.243 / -0.246, but its increments shrink (-0.037 / -0.021 / -0.016 / -0.0097 / -0.003). The
extension to sigma = 0.16, 0.15, 0.14, 0.13, 0.12 runs at M = 480 / 640 / 960 on the laptop.
  T-ext-N  N's slope reverses within the extension: some slope at sigma <= 0.18 is less negative than the one
           before it. It does not cross -0.267. If it does neither, N is UNRESOLVED at sigma = 0.12 and is
           written as such.
  T-ext-Z  Z's slope stays above Gamma_Z = 0.0125 through sigma = 0.12. It has overshot, and a return toward
           Gamma_Z, which the large-deviation form requires eventually, is not reached in this range.

OUTCOMES (appended after running; predictions above left as written).
  V1 HELD: the stored direct-solve R is reproduced to 3e-11 relative (mass 1 to 4e-12, ~35 GMRES iterations), here
     and on the laptop to 1e-14.
  V2 FAILED AS WRITTEN: over M = 160 / 240 / 320 the change grows at some cells, because M = 240 is an outlier
     (Z, sigma = 0.25: +1.1e-3 then -1.4e-3). It shrinks over M = 320 / 480 / 640. The readability clause holds:
     every slope agrees to 1e-4 across M >= 320 (to 1e-3 across all M).
  T  Z TURNED (slopes 0.0394 / 0.0175 / 0.0085 / 0.0066 / 0.0100 / 0.0147). N did not turn by sigma = 0.18, so it
     was UNRESOLVED there. My stated expectation that Z goes negative first was WRONG.
  T-ext-N HELD. N's slope reverses: -0.2458 (0.20-0.18), -0.2430 (0.18-0.16), -0.2357 (0.16-0.15), [-0.2313
     (0.15-0.14), not readable]. It never crosses -0.267.
  T-ext-Z HELD. 0.0185 / 0.0201 / 0.0204 / 0.0201 / 0.0195 (0.18 -> 0.12): above Gamma_Z throughout, and it has
     turned a second time, peaking near sigma = 0.145 and descending toward Gamma_Z.
  INSTRUMENT LIMIT, found here, not pre-registered. q_fall + q_dip departs from 1 as sigma falls (-6e-8 at 0.15,
     -3e-7 at 0.14, -3e-6 at 0.13, -2e-5 at 0.12, at every M). GMRES's residual is amplified by the solve's
     conditioning, about 1 / (escape rate). A cell is read only if that deficit is below 1% of q_fall. This was
     set after seeing the data. It drops N at sigma <= 0.14: 0.13 gives q_fall = 1.9e-14 where the neighbouring
     slope implies ~1e-6, an impossible 1e8 drop. All Z cells pass. FINDINGS §128.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse.linalg as spla

from experiments.wbit_cascade import A, DT, H, HI1, PAIRS, T
from experiments.wbit_cell import fixed_points, kernel


def blocks(M, sigma, c, b2):
    K1, x = kernel(M, A, T, DT, sigma, 0.0)
    K2 = np.empty((M, M, M))
    for i in range(M):
        K2[i] = kernel(M, A, T, DT, sigma, b2 + c * (2 * x[i] - 1))[0].toarray()
    return K1.toarray(), K2, x


def exits_mf(M, sigma, c, b2, rtol=1e-13):
    K1, K2, x = blocks(M, sigma, c, b2)
    lo2 = fixed_points(A, T, DT, H)[0]
    P1, P2 = np.meshgrid(x, x, indexing="ij")
    dip = P1 < 0.5
    fall = (~dip) & (P2 < lo2 + 0.02)
    alive = ~(dip | fall)

    def Pmul(V):
        return np.einsum("ijk,ik->ij", K2, K1 @ V)

    def Aop(v):
        V = v.reshape(M, M) * alive
        return (V - alive * Pmul(V)).ravel()

    op = spla.LinearOperator((M * M, M * M), matvec=Aop, dtype=float)
    out, iters = {}, {}
    for name, target in (("fall", fall), ("dip", dip)):
        b = (alive * Pmul(target.astype(float))).ravel()
        n = [0]
        q, info = spla.gmres(op, b, rtol=rtol, atol=0.0, restart=300, maxiter=50,
                             callback=lambda r: n.__setitem__(0, n[0] + 1), callback_type="pr_norm")
        if info != 0:
            raise RuntimeError(f"gmres did not converge ({name}, info {info})")
        out[name], iters[name] = q.reshape(M, M), n[0]
    hi2 = fixed_points(A, T, DT, H)[-1]
    d = np.where(alive, (P1 - HI1) ** 2 + (P2 - hi2) ** 2, np.inf)
    k = np.unravel_index(np.argmin(d), d.shape)
    qf, qd = float(out["fall"][k]), float(out["dip"][k])
    return {"q_fall": qf, "q_dip": qd, "R": qf / qd, "mass_check": qf + qd, "gmres_iters": iters}


def main():
    out = Path("results/wbit_turn.json")
    rep = json.loads(out.read_text()) if out.exists() else {}
    mode = sys.argv[1] if len(sys.argv) > 1 else "turn"
    if mode == "v1":
        cells = [(n, s, 160) for n in "ZN" for s in (0.30, 0.25)]
    else:
        Ms = [int(m) for m in sys.argv[2].split(",")] if len(sys.argv) > 2 else [160, 240, 320]
        sig = [float(s) for s in sys.argv[3].split(",")] if len(sys.argv) > 3 else \
            [0.34, 0.30, 0.27, 0.25, 0.22, 0.20, 0.18]
        cells = [(n, s, M) for M in Ms for n in "ZN" for s in sig]
    for name, s, M in cells:
        key = f"{name}:{s}:{M}"
        if key in rep:
            continue
        t = time.time()
        r = exits_mf(M, s, *PAIRS[name])
        r["seconds"] = time.time() - t
        rep[key] = r
        print(key, json.dumps(r), flush=True)
        out.write_text(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main()
