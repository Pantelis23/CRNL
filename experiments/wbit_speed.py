"""§130 -- what the edge-below-saddle design costs: error exponent against follow latency, across the window.

§126's pairs hold cell 2's drive at a healthy input fixed (h = 0.1) and move only the coupling c. The edge
e(c) = HI1 - (DC + h) / (2c) sits below cell 1's saddle when c < c_half = (DC + h) / (2 HI1 - 1) = 0.3582.
Cell 2 follows a genuine fall of cell 1 only if its drive with cell 1 at its low well,
d_low(c) = h - c k with k = 2 (HI1 - LO1) = 1.7171, is below -DC, i.e. c > c_min = (DC + h) / k = 0.1791.
So the design window is c in (0.1791, 0.3582): the edge runs from cell 1's low well up to its saddle. This
section asks what moving through the window buys (error exponent) and costs (latency). Cell 2 copies falls only:
at the healthy input it is bistable, so it does not follow cell 1 back up, and no c changes that.

PREDICTIONS, WRITTEN BEFORE RUNNING (constants computed beforehand from the map, no run of the cascade).
  P1  CLOSED FORM ACROSS THE WINDOW. For every edge below 1/2, the never-dipped fall parks cell 1 on its saddle
      (§127), so Gamma(c) = -S2(b2(c)), with b2(c) = h - c (2 HI1 - 1) the drive cell 2 sees there and S2 its 1-D
      escape action at that drive. Collocation W*(c) - S agrees with S2(b2(c)) to 1e-3 at c = 0.19, 0.22, 0.26,
      0.30, 0.34, 0.355. Gamma -> 0 continuously as c -> c_half, because b2 -> -DC and S2 -> 0. Above the window
      (c = 0.40, 0.50, 0.6707), Gamma > 0 and the closed form does not apply.
      KILL: any edge-below-1/2 value off by more than 1e-3 at converged L.
  P2  LATENCY DIVERGES AT THE WINDOW'S FLOOR, BY THE SADDLE-NODE NORMAL FORM. Deterministic latency tau(c): cell 1
      at its low well, cell 2 starting at its healthy-input high well 0.9469, steps until p2 < 1/2. Near c_min,
      cell 2 passes the ghost of its vanished high state. The normal form f(p) - p = -g eps - |beta| (p - p*)^2
      with eps = k (c - c_min), g = 1/6, |beta| = 0.8660 gives tau sqrt(c - c_min) -> pi / sqrt(k g |beta|) = 6.310.
      Tested as a limit (rule 20): over c - c_min = 1e-2, 1e-3, 1e-4, 1e-5, 1e-6 the product approaches 6.310, with
      |product - 6.310| shrinking at every step.
  P3  NOISE CAPS THE DIVERGENCE. At sigma = 0.15 the exact mean latency (1-D chain of cell 2 at frozen drive
      d_low(c), mean first passage from 0.9469 to p2 < 1/2) stays finite as c -> c_min. At c - c_min = 1e-4 it is
      below the deterministic tau. Far from the floor (c = 0.30) it is within 10% of tau.
  P4  THE TRADE-OFF IS MONOTONE (the design table). At sigma = 0.15, along c = 0.19 ... 0.355, R(c) (exact,
      positive-arithmetic solver, M = 320) falls monotonically as c decreases and the latency rises monotonically.
      Neither has an interior optimum.

OUTCOMES (appended after running; predictions above left as written).
  P1 HELD once converged in L. The first pass (L <= 80, pads 0 / 20 / 40) was not converged: W* was still falling
     by 1-2e-3 per 20 steps at c >= 0.30, and it printed 1.1e-3 misses at c = 0.34 and 0.355. Those are not read
     (rule 13). A fast collocation (analytic gradient; reproduces §127's Z to 3e-6) at L up to 428 gives
     Gamma - Gamma_closed = +3.0 / +2.9 / +2.8 / +2.8 / +2.6 / -2.2 e-6 at c = 0.19 / 0.22 / 0.26 / 0.30 / 0.34 /
     0.355. The near-constant +2.8e-6 is the grid error in the 1-D actions S and S2. The window's closed form is
     Gamma(c) = -S2(b2(c)): -0.191 at edge 0.12, rising to -0.0005 at edge 0.496. Just above the window,
     Gamma = +1.1e-5 (c = 0.40), +0.0014 (0.50) and +0.0125 (0.6707), so Gamma passes through 0 at edge = 1/2.
  P2 FAILED AS WRITTEN; the normal form HELD. Measured from §126's C_MIN, the products were 5.40 / 6.04 / 6.23 /
     6.37 / 7.35, and the last step moved away from 6.310. §126's DC is rounded to 6 places (0.207546 against
     0.2075464553), which puts C_MIN 2.65e-7 below the true floor; at the 1e-6 step that is a quarter of the
     distance. Measured from the exact floor, the products are 6.2200 / 6.2834 / 6.3020 / 6.3078 / 6.3095 at
     1e-4 ... 1e-8, approaching 6.3104 monotonically. The instrument's own constant was the axis (rule 13).
  P3 HELD. At sigma = 0.15 the noisy mean latency is 44.5 steps at 1e-4 above the floor (deterministic: 622),
     and 9.64 against 10 at c = 0.30 (3.6%). (P3 was run from the exact floor, the value P2's failure found.)
  P4 HELD. Edge 0.496 -> 0.120 (c 0.355 -> 0.19): R(0.15) = 1.65e-3 / 9.7e-4 / 2.2e-4 / 5.1e-5 / 1.8e-5 /
     1.18e-5, and the noisy latency 7.2 / 7.7 / 9.6 / 12.9 / 19.8 / 33.5 steps. Both are monotone. FINDINGS §130.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from experiments.wbit_action import S_DIP, escape_action, grid_path, joint_action, path_action_states, boundary
from experiments.wbit_cascade import A, DC, DT, H, HI1, T
from experiments.wbit_cell import fixed_points, kernel, step_map

LO1 = 1 - HI1
K = 2 * (HI1 - LO1)
C_MIN = (DC + H) / K
C_HALF = (DC + H) / (2 * HI1 - 1)
HI2 = fixed_points(A, T, DT, H)[-1]
NF_LIMIT = 6.3104138415954125          # pi / sqrt(k g |beta|), computed before the run
D_C_TRUE = 0.20754645532203025         # the high-state saddle-node drive to 1e-15 (§126's DC is rounded to 6 places)
C_MIN_TRUE = (D_C_TRUE + H) / K        # 2.65e-7 above C_MIN: found by P2's failure, see OUTCOMES


def b2_of(c):
    return H - c * (2 * HI1 - 1)


def edge_of(c):
    return HI1 - (DC + H) / (2 * c)


def det_latency(c, tmax=10_000_000):
    d = H - c * K
    p, t = HI2, 0
    while p >= 0.5 and t < tmax:
        p = step_map(p, d, A, T, DT)
        t += 1
    return t


def noisy_latency(c, sigma, M=1600):
    """Exact mean steps from p2 = HI2 to p2 < 1/2, cell 2 at frozen drive d_low(c)."""
    P, x = kernel(M, A, T, DT, sigma, H - c * K)
    alive = np.nonzero(x >= 0.5)[0]
    tau = spla.spsolve((sp.identity(len(alive)) - P[alive][:, alive]).tocsc(), np.ones(len(alive)))
    return float(np.interp(HI2, x[alive], tau))


def colloc_fast(c, b2, P0, g1, s, maxiter=5000):
    """Collocation as `wbit_action.path_action_states`, with an analytic gradient and constant constraint
    Jacobians, so that L in the hundreds is affordable. Validated against the original (tests)."""
    from scipy.optimize import minimize
    st = P0[0]
    L = len(P0) - 1

    def unpack(x):
        P = np.vstack([st, x.reshape(L, 2)])
        return P[:, 0], P[:, 1]

    def parts(x):
        p1, p2 = unpack(x)
        z1 = (p1[1:] - (1 - DT) * p1[:-1]) / DT
        z2 = (p2[1:] - (1 - DT) * p2[:-1]) / DT
        return p1, p2, z1, z2

    def cost_grad(x):
        p1, p2, z1, z2 = parts(x)
        if np.any(z1 <= 0) or np.any(z1 >= 1) or np.any(z2 <= 0) or np.any(z2 >= 1):
            return 1e6, np.zeros_like(x)
        e1 = T * np.log(z1 / (1 - z1)) - A * (2 * p1[:-1] - 1)
        e2 = T * np.log(z2 / (1 - z2)) - b2 - c * (2 * p1[:-1] - 1) - A * (2 * p2[:-1] - 1)
        w1, w2 = T / (z1 * (1 - z1)) / DT, T / (z2 * (1 - z2)) / DT
        g1_, g2_ = np.zeros(L + 1), np.zeros(L + 1)
        g1_[1:] += 0.5 * e1 * w1                                  # d/dp1_{t+1}
        g1_[:-1] += 0.5 * e1 * (-w1 * (1 - DT) - 2 * A) + 0.5 * e2 * (-2 * c)
        g2_[1:] += 0.5 * e2 * w2
        g2_[:-1] += 0.5 * e2 * (-w2 * (1 - DT) - 2 * A)
        grad = np.column_stack([g1_[1:], g2_[1:]]).ravel()
        return 0.25 * float(e1 @ e1 + e2 @ e2), grad

    # z constraints are linear in x: z = Jz x + z0
    n = 2 * L
    Jz = np.zeros((2 * L, n))
    for t in range(L):
        for k in range(2):
            row = k * L + t
            Jz[row, 2 * t + k] = 1 / DT
            if t > 0:
                Jz[row, 2 * (t - 1) + k] = -(1 - DT) / DT
    z0 = np.zeros(2 * L)
    z0[0], z0[L] = -(1 - DT) * st[0] / DT, -(1 - DT) * st[1] / DT

    def term(x):
        return np.array([np.interp(x[-2], g1, s) - x[-1]])

    def term_jac(x):
        j = np.zeros((1, n))
        h = 1e-7
        j[0, -2] = (np.interp(x[-2] + h, g1, s) - np.interp(x[-2] - h, g1, s)) / (2 * h)
        j[0, -1] = -1.0
        return j

    cons = [{"type": "ineq", "fun": lambda x: Jz @ x + z0 - 1e-7, "jac": lambda x: Jz},
            {"type": "ineq", "fun": lambda x: 1 - 1e-7 - (Jz @ x + z0), "jac": lambda x: -Jz},
            {"type": "ineq", "fun": term, "jac": term_jac}]
    r = minimize(cost_grad, P0[1:].ravel(), jac=True, constraints=cons, bounds=[(0.5, 0.9995), (0.0005, 0.9995)] * L,
                 method="SLSQP", options={"maxiter": maxiter, "ftol": 1e-15})
    zz = Jz @ r.x + z0
    feas = bool(np.all(zz > 0) and np.all(zz < 1) and term(r.x)[0] >= -1e-9 and np.all(r.x[0::2] >= 0.5 - 1e-12))
    return float(cost_grad(r.x)[0]), feas, r


def gamma_colloc_long(c, pads, N=81):
    b2 = b2_of(c)
    _, _, fld = joint_action(N, c, b2, field=True)
    P = grid_path(fld)
    g1, s = boundary(c, b2)
    rows = []
    for extra in pads:
        P0 = np.vstack([np.repeat(P[:1], extra, axis=0), P])
        W, feas, r = colloc_fast(c, b2, P0, g1, s)
        rows.append({"L": len(P0) - 1, "W_star": W, "feasible": feas, "success": bool(r.success)})
        print("  ", c, rows[-1], flush=True)
    return rows


def gamma_colloc(c, N=81, pads=(0, 20, 40)):
    b2 = b2_of(c)
    _, _, fld = joint_action(N, c, b2, field=True)
    P = grid_path(fld)
    g1, s = boundary(c, b2)
    rows = []
    for extra in pads:
        P0 = np.vstack([np.repeat(P[:1], extra, axis=0), P])
        W, feas, _ = path_action_states(c, b2, len(P0) - 1, P0, g1, s)
        rows.append({"L": len(P0) - 1, "W_star": W, "feasible": feas})
    best = min(r["W_star"] for r in rows if r["feasible"])
    return best, rows


def main():
    out = Path("results/wbit_speed.json")
    rep = json.loads(out.read_text()) if out.exists() else {}
    part = sys.argv[1] if len(sys.argv) > 1 else "all"
    if part in ("all", "p1"):
        for c in (0.19, 0.22, 0.26, 0.30, 0.34, 0.355, 0.40, 0.50, 0.6707):
            key = f"P1:{c}"
            if key in rep:
                continue
            W, rows = gamma_colloc(c)
            r = {"c": c, "edge": edge_of(c), "b2": b2_of(c), "W_star": W, "Gamma": S_DIP - W, "rows": rows}
            if c < C_HALF:
                r["S2_closed"] = escape_action(b2_of(c))
                r["Gamma_closed"] = -r["S2_closed"]
            rep[key] = r
            print(key, {k: v for k, v in r.items() if k != "rows"}, flush=True)
            out.write_text(json.dumps(rep, indent=2))
    if part == "p1long":                              # rule 13: P1's first pass was not converged in L
        for c, pads in ((0.19, (100, 200)), (0.22, (100, 200)), (0.26, (100, 200)), (0.30, (100, 200)),
                        (0.34, (160, 260)), (0.355, (250, 400)), (0.40, (100, 200)), (0.50, (100, 200))):
            key = f"P1long:{c}"
            if key in rep:
                continue
            rows = gamma_colloc_long(c, pads)
            best = min(r["W_star"] for r in rows if r["feasible"])
            r = {"c": c, "edge": edge_of(c), "W_star": best, "Gamma": S_DIP - best, "rows": rows}
            if c < C_HALF:
                r["Gamma_closed"] = -escape_action(b2_of(c))
                r["diff"] = r["Gamma"] - r["Gamma_closed"]
            rep[key] = r
            print(key, {k: v for k, v in r.items() if k != "rows"}, flush=True)
            out.write_text(json.dumps(rep, indent=2))
    if part in ("all", "p2"):
        rows = []
        for e in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
            tau = det_latency(C_MIN + e)
            rows.append({"c_minus_cmin": e, "tau": tau, "tau_sqrt": tau * e ** 0.5})
            print("P2", rows[-1], flush=True)
        rows_true = []
        for e in (1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
            tau = det_latency(C_MIN_TRUE + e)
            rows_true.append({"c_minus_cmin_true": e, "tau": tau, "tau_sqrt": tau * e ** 0.5})
            print("P2 true", rows_true[-1], flush=True)
        rep["P2"] = {"limit_predicted": NF_LIMIT, "rows": rows, "rows_true_cmin": rows_true}
        out.write_text(json.dumps(rep, indent=2))
    if part in ("all", "p3"):
        rows = []
        for c in (C_MIN_TRUE + 1e-4, C_MIN_TRUE + 1e-3, C_MIN_TRUE + 1e-2, 0.22, 0.26, 0.30, 0.34, 0.355):
            rows.append({"c": c, "tau_det": det_latency(c), "tau_noisy": noisy_latency(c, 0.15)})
            print("P3", rows[-1], flush=True)
        rep["P3"] = rows
        out.write_text(json.dumps(rep, indent=2))
    if part in ("all", "p4"):
        from experiments.wbit_qsd import exits_pos
        rows = []
        for c in (0.19, 0.22, 0.26, 0.30, 0.34, 0.355):
            r = exits_pos(320, 0.15, c, b2_of(c))
            rows.append({"c": c, "edge": edge_of(c), "R": r["R"], "mass_check": r["mass_check"],
                         "tau_det": det_latency(c), "tau_noisy": noisy_latency(c, 0.15)})
            print("P4", rows[-1], flush=True)
        rep["P4"] = rows
        out.write_text(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main()
