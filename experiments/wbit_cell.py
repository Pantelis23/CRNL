"""§125 -- W-bit's analog cell as a CRNL element: an exact transfer operator, validated against W-bit.

W-bit (../W-bit, read-only here) relaxes radix-R cells by `AnalogWbitNetwork.step`:
    s_i <- s_i + dt * (softmax((theta_i + sum_j Theta_ij s_j + xi_i) / T) - s_i),   xi_i ~ N(0, sigma^2 I_R)
For R = 2 the state is one number, p = s_i[1], and the softmax of two inputs is a logistic of their
difference, whose noise is N(0, 2 sigma^2). With a self-coupling Theta_ii = a I (and drive b = theta_i[1]
- theta_i[0], input c (2 p_up - 1) from an upstream cell through Theta_ij = c I):
    p' = (1 - dt) p + dt * logistic((b + a (2p - 1) + c (2 p_up - 1) + eta) / T),   eta ~ N(0, 2 sigma^2).
Bistable when a > 2T at b = 0. This module gives that map an exact discretisation: the next-state law
is the push-forward of a Gaussian through a logistic, so bin masses follow from the normal CDF.

PREDICTIONS, WRITTEN BEFORE RUNNING.
  W1  FAITHFULNESS (rule 10: the harness must not do work the substrate does not). Simulating W-bit's
      own AnalogWbitNetwork (one self-coupled R = 2 cell) and this module's reduced map with the same
      parameters, the one-step law of p' from a fixed p agrees: mean and sd within 3 standard errors
      over 20000 draws, at three starting points (low well, saddle, high well).
  W2  EXACT vs MONTE CARLO. The flip rate of the discretised chain (1/MFPT from the high well to the
      saddle's far side) converges in the grid size M (M = 400 / 800 / 1600: relative change shrinking)
      and agrees with direct simulation of the W-bit map within 2 standard errors at a noise level where
      flips are observable (>= 200 flips).
  W3  ARRHENIUS IN 1/sigma^2. ln(rate) is linear in 1/sigma^2 over the resolvable range, and its slope
      equals minus the discrete-time minimum action (min sum eta_t^2 / 4 over noise paths from the high
      well to the saddle, solved by direct optimisation) within 10% -- the CME's "Omega" is 1/sigma^2 here.

OUTCOMES (appended after running; predictions above left as written). W1 HELD: W-bit's own class and
the reduced map agree (mean 0.2-1.5 se, sd within 1.7%). W2 HELD: 1.574971 / 1.574898 / 1.574880 e-2
at M = 400 / 800 / 1600; Monte Carlo (126045 flips) 1.5756e-2 +- 4.4e-5. W3 HELD: the local Arrhenius
slope approaches the minimum action S = 0.32891 from below, 0.910 -> 0.950 (sigma 0.40 -> 0.23) ->
0.995 S at sigma = 0.13 (rate 7.1e-10). FINDINGS §125.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import expit, ndtr

WBIT_SRC = Path(__file__).resolve().parents[2] / "W-bit" / "src"


def step_map(p, drive, a, T, dt):
    """Deterministic part of the reduced map at noise eta = 0 (drive = b + c (2 p_up - 1))."""
    return (1 - dt) * p + dt * expit((drive + a * (2 * p - 1)) / T)


def kernel(M, a, T, dt, sigma, drive=0.0):
    """Row-stochastic M x M transfer matrix on bin centres of [0, 1]. Exact bin masses: P(p' <= y) =
    Phi((T logit(z) - b - a(2p-1)) / (sigma sqrt 2)), z = (y - (1-dt) p) / dt, clipped to (0, 1)."""
    edges = np.linspace(0.0, 1.0, M + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    rows, cols, vals = [], [], []
    s = sigma * math.sqrt(2.0)
    for i, p in enumerate(centres):
        z = (edges - (1 - dt) * p) / dt
        z = np.clip(z, 1e-300, 1 - 1e-16)
        cdf = np.where(edges <= (1 - dt) * p, 0.0,
                       np.where(edges >= (1 - dt) * p + dt, 1.0,
                                ndtr((T * np.log(z / (1 - z)) - drive - a * (2 * p - 1)) / s)))
        m = np.diff(cdf)
        nz = np.nonzero(m > 1e-300)[0]
        rows += [i] * len(nz); cols += list(nz); vals += list(m[nz] / m[nz].sum())
    return sp.csr_matrix((vals, (rows, cols)), shape=(M, M)), centres


def flip_rate(M, a, T, dt, sigma, drive=0.0):
    """Exact per-step escape rate from the high well: 1 / MFPT from the high fixed point to p < 1/2."""
    P, x = kernel(M, a, T, dt, sigma, drive)
    alive = np.nonzero(x >= 0.5)[0]
    A = (sp.identity(len(alive)) - P[alive][:, alive]).tocsc()
    tau = spla.spsolve(A, np.ones(len(alive)))
    hi = fixed_points(a, T, dt, drive)[-1]
    return 1.0 / float(np.interp(hi, x[alive], tau))


def fixed_points(a, T, dt, drive=0.0):
    g = np.linspace(1e-6, 1 - 1e-6, 200001)
    f = step_map(g, drive, a, T, dt) - g
    idx = np.nonzero(np.sign(f[:-1]) != np.sign(f[1:]))[0]
    return [float(g[i]) for i in idx]


def simulate_flips(a, T, dt, sigma, steps, seed=0, drive=0.0, chains=2000):
    """Vectorised direct simulation of the reduced W-bit map; counts high -> below-1/2 escapes."""
    rng = np.random.default_rng(seed)
    hi = fixed_points(a, T, dt, drive)[-1]
    p = np.full(chains, hi)
    flips, t_alive = 0, 0
    for _ in range(steps):
        eta = rng.normal(0.0, sigma * math.sqrt(2.0), chains)
        p = (1 - dt) * p + dt * expit((drive + a * (2 * p - 1) + eta) / T)
        t_alive += chains
        f = p < 0.5
        flips += int(f.sum())
        p[f] = hi                                   # restart at the high well after each escape
    return flips, t_alive


def wbit_one_step(a, T, dt, sigma, p0, draws, seed=0):
    """W-bit's own AnalogWbitNetwork: one self-coupled R = 2 cell, one step from s = (1-p0, p0)."""
    sys.path.insert(0, str(WBIT_SRC))
    from wbit.analog_network import AnalogWbitNetwork
    import random
    random.seed(seed)
    out = []
    net = AnalogWbitNetwork(1, 2)
    net.set_interaction_weights(0, 0, [[a, 0.0], [0.0, a]])
    for _ in range(draws):
        net.state = [[1 - p0, p0]]
        net.step(temperature=T, noise_level=sigma, dt=dt)
        out.append(net.state[0][1])
    return np.array(out)


def min_action(a, T, dt, n=3000):
    """Discrete-time minimum action from the high well to p < 1/2: A(p) = min_q [eta(p,q)^2/4 + A(q)],
    eta the noise that maps p to q (inverse logistic); value iteration on a grid of [1/2, 1)."""
    g = np.linspace(0.5, 0.999, n)
    A = np.full(n, np.inf); A[0] = 0.0
    P, Q = np.meshgrid(g, np.concatenate([[0.4999], g]), indexing="ij")
    z = (Q - (1 - dt) * P) / dt
    ok = (z > 0) & (z < 1)
    zz = np.where(ok, z, 0.5)
    C = np.where(ok, (T * np.log(zz / (1 - zz)) - a * (2 * P - 1)) ** 2 / 4, np.inf)
    for _ in range(1000):
        An = np.min(C + np.concatenate([[0.0], A])[None, :], axis=1)
        if np.max(np.abs(np.where(np.isfinite(An) & np.isfinite(A), An - A, 0.0))) < 1e-14:
            break
        A = An
    return float(np.interp(fixed_points(a, T, dt)[-1], g, A))


def main():
    import json
    a, T, dt = 1.5, 0.5, 0.5
    rep = {"params": {"a": a, "T": T, "dt": dt}, "fixed_points": fixed_points(a, T, dt)}
    rng = np.random.default_rng(1)
    w1 = []
    for p0 in (rep["fixed_points"][0], 0.5, rep["fixed_points"][-1]):
        w = wbit_one_step(a, T, dt, 0.4, p0, 20000, seed=3)
        r = (1 - dt) * p0 + dt * expit((a * (2 * p0 - 1) + rng.normal(0, 0.4 * math.sqrt(2), 20000)) / T)
        se = math.hypot(w.std(), r.std()) / math.sqrt(20000)
        w1.append({"p0": p0, "dmean_over_se": (w.mean() - r.mean()) / se, "dsd_rel": w.std() / r.std() - 1})
    rep["W1"] = w1
    rep["W2"] = {"rates_M": {str(M): flip_rate(M, a, T, dt, 0.4) for M in (400, 800, 1600)}}
    f, n = simulate_flips(a, T, dt, 0.4, 4000, seed=7)
    rep["W2"]["mc"] = {"flips": f, "steps": n, "rate": f / n, "se": math.sqrt(f) / n}
    S = min_action(a, T, dt)
    w3 = []
    for sig, M in ((0.40, 1600), (0.35, 1600), (0.30, 1600), (0.27, 3200), (0.25, 3200), (0.23, 3200),
                   (0.20, 4000), (0.17, 6400), (0.15, 8000), (0.13, 12000)):
        w3.append({"sigma": sig, "M": M, "rate": flip_rate(M, a, T, dt, sig)})
        print(w3[-1], flush=True)
    for i in range(1, len(w3)):
        w3[i]["local_slope_over_S"] = (math.log(w3[i - 1]["rate"]) - math.log(w3[i]["rate"])) / (
            1 / w3[i]["sigma"] ** 2 - 1 / w3[i - 1]["sigma"] ** 2) / S
    rep["W3"] = {"S": S, "rows": w3}
    Path("results/wbit_cell.json").write_text(json.dumps(rep, indent=2, default=float))
    print("wrote results/wbit_cell.json")


if __name__ == "__main__":
    main()
