"""§127 -- W-bit's joint minimum action: the absolute test of §126's trend.

Discrete-time large deviations of the two-cell W-bit cascade (§126). A path costs sum_t (eta1^2 + eta2^2)/4 per
1/sigma^2 (each logit noise is N(0, 2 sigma^2)). Each cell's noise is fixed by its own step, so the one-step cost
separates: c1(p1 -> p1') + c2(p2 -> p2' | p1), with eta = T logit((p' - (1-dt) p)/dt) - drive.
  S   cell 1's dip action, high well -> p1 < 1/2 (1-D; §125: 0.32891).
  W*  the least action from both high wells to the joint basin of (cell 1 high, cell 2 low) -- the states whose
      deterministic orbit (eta = 0) ends with cell 2 low -- with p1 >= 1/2 at every step (never dipped).
  Gamma = S - W*. R ~ e^{Gamma / sigma^2} (§124's form), so the local slope of ln R against 1/sigma^2 must
  approach Gamma.

Value iteration on an N x N grid of p1 in [1/2, 1), p2 in (0, 1); successors restricted to the grid. Grid bias is
the approximation's own parameter (rule 13): run at N = 81, 121, 161 and extrapolate.

PREDICTIONS, WRITTEN BEFORE RUNNING.
  A1  SIGNS. Gamma_Z > 0 and Gamma_N < 0 at every N -- the sign of §126's trend in each pair.
  A2  ABSOLUTE (rule 16). The extrapolated Gamma lies between §126's last measured local slope and 0 for Z
      (slopes 0.039 / 0.018 / 0.009 falling: if they approach Gamma from above, Gamma_Z in [0, 0.009]), and for
      N within 0.1 of the last slope (-0.22, still steepening). A Gamma_Z above 0.02 would mean the measured
      slope is falling AWAY from it, and the large-deviation form does not describe §126's range.

ADDED AFTER THE GRID RUNS AT N = 41-161 WERE SEEN, BEFORE THE ROUTE BELOW WAS RUN. Grid W* converges slowly
(Z 0.3455 / 0.3272 / 0.3222 / 0.3205; N 0.4966 / 0.4625 / 0.4508 / 0.4458 at N = 41 / 81 / 121 / 161) and the
1/N, 1/N^2 and free-exponent extrapolations of N = 81-161 give W*_Z in [0.3155, 0.3186], W*_N in [0.4307, 0.4394].
S_grid does not converge smoothly (0.33098 / 0.32967 / 0.32957: the high well's nearest grid point moves), so
Gamma is taken against the exact S = 0.32891 of §125.
  A3  A SECOND ROUTE, GRID-FREE. Minimising sum (eta1^2 + eta2^2)/4 directly over L-step noise sequences
      (SLSQP; constraints p1_t >= 1/2 at every step and the final state on the fall side of the deterministic
      basin boundary p2 = s(p1)) gives W* inside the extrapolation band widened by 0.002 on each side
      (Z [0.3135, 0.3206], N [0.4287, 0.4414]), converged in L (L = 20 / 40 / 60 agree to 1e-4 at the best).
      Every feasible path is an upper bound on the true W*, so a value ABOVE the band means the optimiser
      found a local minimum (re-start it, do not read it); a value BELOW means the grid extrapolation is wrong.

OUTCOMES (appended after running; predictions above left as written).
  W*_Z = 0.31646 (collocation, L = 46 / 66 agree to 1e-6; shooting 0.31652 at L = 40-100; grid 1/N^2 extrapolation
  of N = 201 / 241: 0.31643) -> Gamma_Z = +0.01245. W*_N = 0.43552 (collocation, L = 70 / 90: 0.43573 / 0.43552;
  grid extrapolations of N = 161-241: 0.4289-0.4364) -> Gamma_N = -0.10661; independently, cell 2's own escape
  action at the drive it sees with cell 1 parked on its saddle is 0.10659, so W*_N = S + S2(b2) to 2e-5.
  A1 FAILED AS WRITTEN at N = 41 (Gamma_Z = -0.013: W*'s grid bias, +0.029, exceeds Gamma_Z itself); held at
     every N >= 81 and in the converged values.
  A2 FAILED in both pairs. Gamma_Z = 0.0125 is outside [0, 0.009] (not above 0.02); Gamma_N = -0.107 is 0.11 from
     -0.217. §126's local slopes do not approach Gamma: Z's (0.039 / 0.018 / 0.009) has already fallen through
     it, N's (-0.158 / -0.195 / -0.217) moves away from it. A2 assumed a monotone approach; if Gamma is the
     asymptote, both slopes must turn back below sigma = 0.25.
  A3 HELD in both pairs (inside the bands). Its L clause failed for shooting at L = 20 (Z 0.31755, 1e-3 high);
     shooting never found N's minimum (best feasible 0.4888 at L = 20, 0.5892 at L = 40; S2(h) = 0.5957);
     collocation over states was needed and is the reported route. FINDINGS §127.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from experiments.wbit_cascade import A, DT, PAIRS, T
from experiments.wbit_cell import fixed_points

S_DIP = 0.32890998229862456


def eta_cost(p, q, drive):
    """(eta^2)/4 for a step p -> q under drive (arrays broadcast); inf if q is unreachable."""
    z = (q - (1 - DT) * p) / DT
    ok = (z > 0) & (z < 1)
    zz = np.where(ok, z, 0.5)
    e = T * np.log(zz / (1 - zz)) - drive - A * (2 * p - 1)
    return np.where(ok, e * e / 4, np.inf)


def det_target(g1, g2, c, b2, steps=400):
    """States (p1 >= 1/2) whose eta = 0 orbit ends with cell 2 low (cell 1 always relaxes high)."""
    from scipy.special import expit
    P1, P2 = np.meshgrid(g1, g2, indexing="ij")
    for _ in range(steps):
        P1, P2 = ((1 - DT) * P1 + DT * expit((A * (2 * P1 - 1)) / T),
                  (1 - DT) * P2 + DT * expit((b2 + c * (2 * P1 - 1) + A * (2 * P2 - 1)) / T))
    return P2 < 0.397


def joint_action(N, c, b2, sweeps=400, field=False):
    g1 = np.linspace(0.5, 0.9995, N)
    g2 = np.linspace(0.0005, 0.9995, N)
    Wv = np.where(det_target(g1, g2, c, b2), 0.0, np.inf)
    C1 = eta_cost(g1[:, None], g1[None, :], 0.0)                    # (p1, p1')
    C2 = np.stack([eta_cost(g2[:, None], g2[None, :], b2 + c * (2 * p1 - 1)) for p1 in g1])  # (p1, p2, p2')
    for it in range(sweeps):
        Wn = Wv.copy()
        for i in range(N):
            # min over p2' of C2[i, p2, p2'] + W[p1', p2'] -> (p2, p1')
            B = np.min(C2[i][:, None, :] + Wv[None, :, :], axis=2)
            Wn[i] = np.minimum(Wv[i], np.min(C1[i][None, :] + B, axis=1))
        done = np.array_equal(np.isfinite(Wn), np.isfinite(Wv)) and np.nanmax(
            np.abs(np.where(np.isfinite(Wn), Wn - Wv, 0.0))) < 1e-13
        Wv = Wn
        if done:
            break
    hi1 = fixed_points(A, T, DT)[-1]
    hi2 = fixed_points(A, T, DT, 0.1)[-1]
    i1, i2 = int(np.argmin(abs(g1 - hi1))), int(np.argmin(abs(g2 - hi2)))
    if field:
        return float(Wv[i1, i2]), it, (g1, g2, Wv, C1, C2, i1, i2)
    return float(Wv[i1, i2]), it


def grid_path(fld, maxlen=400):
    """Greedy descent of the grid value function from the start: the grid's optimal path (p1_t, p2_t)."""
    g1, g2, Wv, C1, C2, i, j = fld
    path = [(g1[i], g2[j])]
    while Wv[i, j] > 0 and len(path) < maxlen:
        tot = C1[i][:, None] + C2[i][j][None, :] + Wv
        i, j = np.unravel_index(np.argmin(tot), tot.shape)
        path.append((g1[i], g2[j]))
    return np.array(path)


def path_to_noise(path, c, b2):
    """The noise sequence that produces a given path exactly (inverse logistic per cell)."""
    p1, p2 = path[:-1, 0], path[:-1, 1]
    q1, q2 = path[1:, 0], path[1:, 1]
    z1 = np.clip((q1 - (1 - DT) * p1) / DT, 1e-12, 1 - 1e-12)
    z2 = np.clip((q2 - (1 - DT) * p2) / DT, 1e-12, 1 - 1e-12)
    e1 = T * np.log(z1 / (1 - z1)) - A * (2 * p1 - 1)
    e2 = T * np.log(z2 / (1 - z2)) - b2 - c * (2 * p1 - 1) - A * (2 * p2 - 1)
    return e1, e2


def dip_action_grid(N):
    """S on the same p1 grid (plus a sink just below 1/2), so its grid bias matches W*'s."""
    g1 = np.linspace(0.5, 0.9995, N)
    q = np.concatenate([[0.4999], g1])
    C = eta_cost(g1[:, None], q[None, :], 0.0)
    Av = np.full(N, np.inf)
    for _ in range(2000):
        An = np.min(C + np.concatenate([[0.0], Av])[None, :], axis=1)
        if np.array_equal(np.isfinite(An), np.isfinite(Av)) and np.max(np.abs(np.where(np.isfinite(An), An - Av, 0))) < 1e-14:
            Av = An
            break
        Av = An
    return float(Av[int(np.argmin(abs(g1 - fixed_points(A, T, DT)[-1])))])


def escape_action(drive, n=4000):
    """One cell's minimum action at constant drive, from its high well to just below its saddle (1-D value
    iteration). escape_action(0) is S; for N, W* = S + escape_action(b2) (cell 1 parked on its saddle)."""
    lo, sad, hi = fixed_points(A, T, DT, drive)
    g = np.linspace(sad, 0.9995, n)
    q = np.concatenate([[sad - 1e-4], g])
    C = eta_cost(g[:, None], q[None, :], drive)
    Av = np.full(n, np.inf)
    for _ in range(4000):
        An = np.min(C + np.concatenate([[0.0], Av])[None, :], axis=1)
        done = np.array_equal(np.isfinite(An), np.isfinite(Av)) and np.max(
            np.abs(np.where(np.isfinite(An), An - Av, 0.0))) < 1e-14
        Av = An
        if done:
            break
    return float(np.interp(hi, g, Av))


def boundary(c, b2, n=801, steps=400):
    """s(p1): the p2 below which the eta = 0 orbit from (p1, p2) ends with cell 2 low (1.0 if all of it does)."""
    from scipy.special import expit
    g1 = np.linspace(0.5, 0.9995, n)
    lo, hi = np.zeros(n), np.ones(n)
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        P1, P2 = g1.copy(), mid.copy()
        for _ in range(steps):
            P1, P2 = ((1 - DT) * P1 + DT * expit((A * (2 * P1 - 1)) / T),
                      (1 - DT) * P2 + DT * expit((b2 + c * (2 * P1 - 1) + A * (2 * P2 - 1)) / T))
        falls = P2 < 0.397
        lo, hi = np.where(falls, mid, lo), np.where(falls, hi, mid)
    return g1, 0.5 * (lo + hi)


def path_action(c, b2, L, x0, g1=None, s=None):
    """Grid-free route (A3): minimise sum (eta1^2 + eta2^2)/4 over L-step noise sequences from both high wells,
    with p1_t >= 1/2 at every step and the end state on the fall side of the boundary s(p1). SLSQP."""
    from scipy.optimize import minimize
    from scipy.special import expit
    if g1 is None:
        g1, s = boundary(c, b2)
    p10, p20 = fixed_points(A, T, DT)[-1], fixed_points(A, T, DT, 0.1)[-1]

    def orbit(x):
        e1, e2 = x[:L], x[L:]
        p1, p2, out = p10, p20, []
        for t in range(L):
            p1, p2 = ((1 - DT) * p1 + DT * expit((A * (2 * p1 - 1) + e1[t]) / T),
                      (1 - DT) * p2 + DT * expit((b2 + c * (2 * p1 - 1) + A * (2 * p2 - 1) + e2[t]) / T))
            out.append((p1, p2))
        return np.array(out)

    cons = [{"type": "ineq", "fun": lambda x: orbit(x)[:, 0] - 0.5},
            {"type": "ineq", "fun": lambda x: np.interp(orbit(x)[-1, 0], g1, s) - orbit(x)[-1, 1]}]
    r = minimize(lambda x: 0.25 * np.dot(x, x), x0, jac=lambda x: 0.5 * x, constraints=cons,
                 method="SLSQP", options={"maxiter": 2000, "ftol": 1e-13})
    feas = bool(np.all(orbit(r.x)[:, 0] >= 0.5 - 1e-9)) and \
        bool(np.interp(orbit(r.x)[-1, 0], g1, s) - orbit(r.x)[-1, 1] >= -1e-9)
    return float(r.fun), feas, r


def path_action_states(c, b2, L, P0=None, g1=None, s=None):
    """A3 by collocation: optimise over the states (p1_t, p2_t), t = 1..L, with the noise an explicit function
    of consecutive states (path_to_noise). Shooting over the noise (path_action) is ill-conditioned when the
    path hugs cell 1's saddle, which amplifies errors by 1.25 per step. P0: (L+1) x 2 start guess."""
    from scipy.optimize import minimize
    if g1 is None:
        g1, s = boundary(c, b2)
    st = np.array([fixed_points(A, T, DT)[-1], fixed_points(A, T, DT, 0.1)[-1]])

    def full(x):
        return np.vstack([st, x.reshape(L, 2)])

    def zs(x):
        P = full(x)
        return np.concatenate([(P[1:, 0] - (1 - DT) * P[:-1, 0]) / DT, (P[1:, 1] - (1 - DT) * P[:-1, 1]) / DT])

    def cost(x):
        P = full(x)
        e1, e2 = path_to_noise(P, c, b2)
        return 0.25 * float(e1 @ e1 + e2 @ e2) if np.all((zs(x) > 0) & (zs(x) < 1)) else 1e6

    cons = [{"type": "ineq", "fun": lambda x: zs(x) - 1e-7},
            {"type": "ineq", "fun": lambda x: 1 - 1e-7 - zs(x)},
            {"type": "ineq", "fun": lambda x: np.interp(x[-2], g1, s) - x[-1]}]
    if P0 is None:
        P0 = np.column_stack([np.linspace(st[0], 0.55, L + 1), np.linspace(st[1], 0.5, L + 1)])
    x0 = P0[1:].ravel()
    bounds = [(0.5, 0.9995), (0.0005, 0.9995)] * L
    r = minimize(cost, x0, constraints=cons, bounds=bounds, method="SLSQP",
                 options={"maxiter": 3000, "ftol": 1e-14})
    feas = bool(np.all(zs(r.x) > 0) and np.all(zs(r.x) < 1) and np.interp(r.x[-2], g1, s) - r.x[-1] >= -1e-9
                and np.all(r.x[0::2] >= 0.5 - 1e-12))
    return float(cost(r.x)), feas, r


def main():
    out = Path("results/wbit_action.json")
    rep = json.loads(out.read_text()) if out.exists() else {}
    Ns = [int(n) for n in sys.argv[1].split(",")] if len(sys.argv) > 1 else [81, 121, 161]
    for name, (c, b2) in PAIRS.items():
        for N in Ns:
            key = f"{name}:{N}"
            if key in rep:
                continue
            W, it = joint_action(N, c, b2)
            Sg = dip_action_grid(N)
            rep[key] = {"W_star": W, "S_grid": Sg, "Gamma_grid": Sg - W, "Gamma_vs_exact_S": S_DIP - W, "sweeps": it}
            print(key, rep[key], flush=True)
            out.write_text(json.dumps(rep, indent=2))
    # A3, the grid-free route: collocation warm-started from the N = 121 grid path, padded at the start
    # (sitting at the fixed point is free) to test convergence in L.
    for name, (c, b2) in PAIRS.items():
        key = f"colloc:{name}"
        if key in rep:
            continue
        _, _, fld = joint_action(121, c, b2, field=True)
        P = grid_path(fld)
        g1, s = boundary(c, b2)
        rows = []
        for extra in (0, 20, 40):
            P0 = np.vstack([np.repeat(P[:1], extra, axis=0), P])
            W, feas, r = path_action_states(c, b2, len(P0) - 1, P0, g1, s)
            rows.append({"L": len(P0) - 1, "W_star": W, "feasible": feas, "Gamma": S_DIP - W})
            print(key, rows[-1], flush=True)
        rep[key] = rows
        out.write_text(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main()
