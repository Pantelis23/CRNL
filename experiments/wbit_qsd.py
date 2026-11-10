"""§129 -- T-WBIT-a″: following N's slope below sigma = 0.15 with a solver that never subtracts.

§128's GMRES solve loses q_fall when it is smaller than the solve's error, about rtol / (escape rate). Here the
exit probabilities come from positive arithmetic only:
    mu_0 = delta(start),  mu_{t+1} = mu_t P_aa  (transpose product, all entries >= 0)
    F_T = sum_{t<T} mu_t . b_fall,  D_T = sum_{t<T} mu_t . b_dip,  m_T = |mu_T|.
Once mu_T has relaxed to the quasi-stationary distribution pi (its exit split r = pi.b_fall / pi.b_dip is
constant in t), the remaining mass m_T exits in that ratio, exactly:
    q_fall = F_T + m_T r / (1 + r),   q_dip = D_T + m_T / (1 + r).
The only approximation is mu_T proportional to pi, which is controlled by watching r converge. No step subtracts,
so tiny q_fall keeps its relative precision. Mass is conserved to rounding by construction, so it is no longer a
test of anything.

PREDICTIONS, WRITTEN BEFORE RUNNING.
  Q1  INSTRUMENT. This reproduces §128's GMRES R at its readable cells (Z and N, sigma = 0.25, 0.20, 0.16,
      M = 320) to 1e-7 relative.
  Q2  INSTRUMENT, INDEPENDENT. At sigma = 0.14 and 0.12 (N), R does not depend on the stopping rule: stopping at
      r-stability 1e-10 or 1e-13 changes R by < 1e-8 relative. Grid: R changes by < 1% from M = 640 to 960.
  T″  THE KILL TEST (THEORIES T-WBIT-a″). Over sigma = 0.16, 0.15, 0.14, 0.13, 0.12, 0.11, 0.10, N's local slope keeps
      rising: each slope is less negative than the one before, and stays below Gamma_N = -0.107. A second reversal
      (a slope more negative than its predecessor) refutes the suspect that N's slow approach is the
      saddle-parking path's prefactor. Z's slope keeps descending toward +0.0125 without a third reversal.
      Rule 20: the criterion is the direction of the sequence, not a distance to Gamma.

OUTCOMES (appended after running; predictions above left as written).
  Q1 HELD at the four cells where it can be judged: 4e-12 / 2e-12 (Z), 3e-11 / 7e-9 (N) at sigma = 0.25 / 0.20,
     M = 320. The 7e-9 is GMRES's error: N at 0.20 is where its mass deficit first shows. Q1 named sigma = 0.16 at
     M = 320, but §128 never ran that cell, so it could not be judged there. That was an error in writing the
     prediction.
  Q2 HELD. The stopping rule moves R by <= 7e-15. The grid moves it by <= 0.33% (M = 640 -> 960) and by <= 0.12%
     (960 -> 1280), shrinking. Mass is conserved to 2e-14 in every cell.
  T″ HELD. N's slopes over 0.18 -> 0.10 are -0.2431 / -0.2353 / -0.2248 / -0.2093 / -0.1898 / -0.1712 / -0.1576
     (M = 960; M = 640 and 1280 agree to 1e-4). They rise at every step and stay below -0.107. Z's are 0.0185 /
     0.0201 / 0.0204 / 0.0201 / 0.0195 / 0.0185 / 0.0173: past their maximum and descending, with no third
     reversal. FINDINGS §129.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

from experiments.wbit_cascade import A, DT, H, HI1, PAIRS, T
from experiments.wbit_cell import fixed_points
from experiments.wbit_turn import blocks


def exits_pos(M, sigma, c, b2, tol=1e-13, tmax=200000, check=50):
    K1, K2, x = blocks(M, sigma, c, b2)
    lo2 = fixed_points(A, T, DT, H)[0]
    P1, P2 = np.meshgrid(x, x, indexing="ij")
    dip = P1 < 0.5
    fall = (~dip) & (P2 < lo2 + 0.02)
    alive = ~(dip | fall)

    def Pmul(V):                       # (P v)[i, j]
        return np.einsum("ijk,ik->ij", K2, K1 @ V)

    def Pt(mu):                        # (mu P)[i', j']
        return K1.T @ np.einsum("ij,ijk->ik", mu, K2)

    b_fall = alive * Pmul(fall.astype(float))
    b_dip = alive * Pmul(dip.astype(float))
    hi2 = fixed_points(A, T, DT, H)[-1]
    d = np.where(alive, (P1 - HI1) ** 2 + (P2 - hi2) ** 2, np.inf)
    mu = np.zeros((M, M))
    mu[np.unravel_index(np.argmin(d), d.shape)] = 1.0
    F = D = 0.0
    r_prev, t = None, 0
    while t < tmax:
        f, dd = float((mu * b_fall).sum()), float((mu * b_dip).sum())
        F += f; D += dd
        mu = alive * Pt(mu)
        t += 1
        if t % check == 0:
            r = float((mu * b_fall).sum()) / float((mu * b_dip).sum())
            if r_prev is not None and abs(r / r_prev - 1) < tol:
                break
            r_prev = r
    m = float(mu.sum())
    qf, qd = F + m * r / (1 + r), D + m / (1 + r)
    return {"q_fall": qf, "q_dip": qd, "R": qf / qd, "r_qsd": r, "mass_left": m, "steps": t,
            "mass_check": F + D + m}


def main():
    out = Path("results/wbit_qsd.json")
    rep = json.loads(out.read_text()) if out.exists() else {}
    Ms = [int(m) for m in sys.argv[1].split(",")]
    sig = [float(s) for s in sys.argv[2].split(",")]
    names = sys.argv[3] if len(sys.argv) > 3 else "ZN"
    tol = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-13
    for M in Ms:
        for name in names:
            for s in sig:
                key = f"{name}:{s}:{M}" + ("" if tol == 1e-13 else f":tol{tol:g}")
                if key in rep:
                    continue
                t0 = time.time()
                r = exits_pos(M, s, *PAIRS[name], tol=tol)
                r["seconds"] = time.time() - t0
                rep[key] = r
                print(key, json.dumps(r), flush=True)
                out.write_text(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main()
