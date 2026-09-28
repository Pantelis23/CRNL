"""§120 -- can the near-miss constant be predicted from single-element data alone?

§119 found, out of sample, that a restoring stage's failure probability given correct input grows
as phi * p_d * k1 t, with phi = w_d/(k1 t) a stage-1 quantity and p_d -- stage 2's low probability
on dipped-and-recovered trajectories -- NOT set by zone width alone (Q1: equal-width pairs a and Z
differ ~5x). Rule 16: until the constant is computed in absolute terms from something other than
the joint chain it is a fitted shape, not a tested model.

THE REDUCTION. Stage 2 as a two-state element whose high -> low hazard is kappa2(x) = its own decay
rate absorbed at its saddle with its input FROZEN at x (the §92/§96 frozen limit, pointwise). Stage
1's exact 1-D chain, with §118's dip tag, carries a stage-2-alive density that is killed at rate
kappa2(n1/Omega). Stage 2's low -> high recovery is ignored (it is exponentially slow for every
downstream element here except S's). Then

    p_d_red = 1 - alive(d = 1, stage 1 high) / all(d = 1, stage 1 high)   at t = 16,

and p0_red the same on d = 0. No joint chain enters. Stage 1's seed = its QSD; stage 2 starts alive.

PREDICTIONS, WRITTEN BEFORE RUNNING.

  R1  ORDER. p_d_red orders all eight arms of §119 the way the exact p_d does, at both Omega --
      in particular a < Z, the pair that zone width got wrong. (Criterion kill test: a reduction
      that only sees width would put a ~ Z and fail here.)
  R2  ABSOLUTE, rule 16. p_d_red / p_d lies within [0.5, 2] in every arm with w > 0 at every
      Omega. The frozen limit was 2-7x off for §102's rates, but those were barrier crossings under
      a fluctuating input; here most of the damage is a barrier-free fall during a slow crawl, which
      is where freezing the input should be closest to right. If R2 fails in a way that is constant
      across arms, report the constant; if it is not constant, the reduction is not the mechanism.
  R3  NEVER-DIPPED TOO. p0_red / p0 in [0.5, 2] wherever p0 > 1e-9: §119 found p0 on stage 1's
      clock in arm X, and the same frozen hazard along never-dipped excursions should give it.

OUTCOMES, appended after running (every prediction block left as written). R1, R2, R3 FAILED: the
frozen hazard over-predicts p_d 1.1-117x and p0 2x-3.5e4x, and swaps a and c. R4 FAILED by its letter
(a at 0.34-0.39) with every positive-zone ratio in 0.34-0.67 over a 90x range of p_d. R5 held at
Omega = 40; a and d tie within noise at 55. R6 FAILED by its letter in X (0.66-0.70), and every
other cell sits in [0.7, 1.06], none overshooting; N becomes exact within error. FINDINGS §120.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from experiments.weak_upstream import ARMS, MARKS, dip_rate, element

T_END = MARKS[-1]


def hazard(B, om, cap):
    """kappa2 at every stage-1 count 0..cap: stage 2's absorbed decay rate, input frozen there."""
    return np.array([dip_rate(B, om, n / om, False) for n in range(cap + 1)])


def reduced(A, B, om):
    Q, m, _ = A._gen(om, 0.0, True)                    # stage 1 alone, 1-D, reflecting box
    cap = m - 1
    sad1 = int(np.ceil(A.r2 * om))
    kap = hazard(B, om, cap)
    C = Q.tocoo()
    r, c, v = C.row, C.col, C.data
    dip = (r != c) & (r >= sad1) & (c < sad1)
    # tagged 1-D chain (d = 0 block, d = 1 block), same construction as §118
    T = sp.coo_matrix((np.concatenate([v[~dip], v[dip], v]),
                       (np.concatenate([r[~dip], r[dip], r + m]),
                        np.concatenate([c[~dip], c[dip] + m, c + m]))), shape=(2 * m, 2 * m)).tocsr()
    K = sp.diags(np.concatenate([kap, kap]))
    keep, w = A.qsd(om, 0.0, first=True)
    p = np.zeros(2 * m); p[keep] = w
    all_ = p.copy()                                    # stage-1 law, no killing
    alive = p.copy()                                   # stage-1 law x P(stage 2 still high)
    Ta = ((T - K).T * 2.0).tocsr()
    Tt = (T.T * 2.0).tocsr()
    t = 0.0
    while t < T_END - 1e-9:
        all_ = spla.expm_multiply(Tt, all_)
        alive = spla.expm_multiply(Ta, alive)
        t += 2.0
    hi = np.arange(m) >= sad1
    out = {}
    for d, sl in ((0, slice(0, m)), (1, slice(m, 2 * m))):
        tot = all_[sl][hi].sum()
        out[f"p{d}_red"] = float(1.0 - alive[sl][hi].sum() / tot)
    return out


MEMORY = """
§120b -- written after R1-R3 FAILED (frozen hazard over-predicts p_d 1.1x-117x, not a constant, the
ratio growing with Omega; a and c swapped; p0 over-predicted 2x-35000x). Reading: stage 2 is not a
memoryless two-state switch. A fall is a continuous drift that takes time and REVERSES if stage 1
returns before stage 2 has passed its separatrix, so a frozen hazard over-counts every short
excursion -- more so at large Omega, where excursions below edge_2 are shallower and briefer.

THE SECOND REDUCTION, still single-element data only. Stage 1 stochastic from its own 1-D chain
(tau-leaped, many paths), started one step below its saddle at the moment of the dip; stage 2
DETERMINISTIC, its own rate equation, driven by x1(t) along each path, started at its high fixed
point with input at stage 1's operating point. Dip times uniform over the window, so each path runs
for T_r ~ U(0, 16). p_d_ode = P(x2 < saddle_2 at T_r | stage 1 back above its saddle at T_r).

  R4  ABSOLUTE. p_d_ode / p_d in [0.5, 2] for every arm with w > 0, at every Omega. It omits the
      ~24% that §118 P11b put on stage 2's pre-dip sag, and stage 2's own noise, so it should run
      LOW, not high -- the opposite side from R2's failure. A ratio above 2 anywhere means the
      memory is still not the missing piece.
  R5  ORDER. p_d_ode orders the eight arms as the exact p_d does, including a < c and a < Z.
"""

TAU = 0.002
PATHS = 20000


def ode_rate(B, x2, x1):
    k1a, k1r, k2b, k2r = B.c
    return B.hill(x1) * k1a * x2 ** 2 + k2b - k1r * x2 ** 3 - k2r * x2


def high_fixed_point(B, x_in):
    k1a, k1r, k2b, k2r = B.c
    r = np.roots([-k1r, B.hill(x_in) * k1a, -k2r, k2b])
    return float(max(z.real for z in r if abs(z.imag) < 1e-9))


def memory_reduction(A, B, om, seed=0):
    rng = np.random.default_rng(seed)
    sad1 = int(np.ceil(A.r2 * om))
    mu1 = A.operating_points(om, 1)[0][0]
    n1 = np.full(PATHS, float(sad1 - 1))
    x2 = np.full(PATHS, high_fixed_point(B, mu1))
    t_r = rng.uniform(0.0, T_END, PATHS)
    done = np.zeros(PATHS, bool)
    fin_n1, fin_x2 = np.zeros(PATHS), np.zeros(PATHS)
    t = 0.0
    while not done.all():
        lam, mu = A.rates(n1, np.zeros(PATHS), om, True)
        n1 = np.maximum(n1 + rng.poisson(lam * TAU) - rng.poisson(mu * TAU), 0.0)
        x2 = x2 + TAU * ode_rate(B, x2, n1 / om)
        t += TAU
        hit = (~done) & (t >= t_r)
        fin_n1[hit], fin_x2[hit] = n1[hit], x2[hit]
        done |= hit
    back = fin_n1 >= sad1
    low = fin_x2 < B.r2
    return {"p_d_ode": float(low[back].mean()) if back.any() else float("nan"),
            "back": int(back.sum()), "paths": PATHS}


NOISE = """
§120c -- the kill test named in §120's summary, first half, written BEFORE running. The memory
reduction runs 0.34-0.67 of the exact p_d. Suspects: stage 2's own noise, and the pre-dip sag
(§118 P11b: ~24% of the effect). Here stage 2 is made STOCHASTIC -- tau-leaped from its own birth-
death rates, driven by the same stage-1 paths, started from its own quasi-stationary law with input
at stage 1's operating point -- and nothing else changes. This isolates the noise; the sag stays out.

  R6  Every positive-zone ratio rises from 0.34-0.67 into [0.7, 1.05]: noise supplies most of the
      missing factor and the sag (~0.76-0.8 of the effect) the rest, so it must NOT overshoot 1
      by more than Monte Carlo error. N's ratio becomes non-zero. If ratios stay near 0.5, noise is
      not the missing piece; if they overshoot well past 1, the dip-time/conditioning shortcuts in
      the reduction are carrying errors of their own and the residual cannot be read as the sag.
"""


def noisy_reduction(A, B, om, seed=0):
    rng = np.random.default_rng(seed)
    sad1, sad2 = int(np.ceil(A.r2 * om)), int(np.ceil(B.r2 * om))
    mu1 = A.operating_points(om, 1)[0][0]
    keep, w = B.qsd(om, mu1, first=False)
    n1 = np.full(PATHS, float(sad1 - 1))
    n2 = rng.choice(keep, size=PATHS, p=w).astype(float)
    t_r = rng.uniform(0.0, T_END, PATHS)
    done = np.zeros(PATHS, bool)
    fin_n1, fin_n2 = np.zeros(PATHS), np.zeros(PATHS)
    t = 0.0
    while not done.all():
        lam, mu = A.rates(n1, np.zeros(PATHS), om, True)
        l2, u2 = B.rates(n2, n1, om, False)
        n1 = np.maximum(n1 + rng.poisson(lam * TAU) - rng.poisson(mu * TAU), 0.0)
        n2 = np.maximum(n2 + rng.poisson(l2 * TAU) - rng.poisson(u2 * TAU), 0.0)
        t += TAU
        hit = (~done) & (t >= t_r)
        fin_n1[hit], fin_n2[hit] = n1[hit], n2[hit]
        done |= hit
    back = fin_n1 >= sad1
    return {"p_d_noisy": float((fin_n2[back] < sad2).mean()), "back": int(back.sum())}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("results/frozen_downstream.json"))
    args = ap.parse_args()
    ref = json.loads(pathlib.Path("results/weak_upstream.json").read_text())
    rep = {}
    for key in sorted(k for k in ref if ":" in k):
        arm, om = key.split(":")
        A, B = element(ARMS[arm][0]), element(ARMS[arm][1])
        r = reduced(A, B, int(om))
        s = ref[key]["series"][-1]
        r.update(p_d=s["p_d"], p0=s["p0"], ratio_d=r["p1_red"] / s["p_d"],
                 ratio_0=r["p0_red"] / s["p0"] if s["p0"] > 1e-9 else float("nan"))
        m = memory_reduction(A, B, int(om))
        m["se"] = float(np.sqrt(max(m["p_d_ode"] * (1 - m["p_d_ode"]), 1e-12) / m["back"]))
        r.update(p_d_ode=m["p_d_ode"], se_ode=m["se"], back=m["back"], ratio_ode=m["p_d_ode"] / s["p_d"])
        z = noisy_reduction(A, B, int(om))
        r.update(p_d_noisy=z["p_d_noisy"], back_noisy=z["back"], ratio_noisy=z["p_d_noisy"] / s["p_d"],
                 se_noisy=float(np.sqrt(max(z["p_d_noisy"] * (1 - z["p_d_noisy"]), 1e-12) / z["back"])))
        rep[key] = r
        print(f"{key:>5}: p_d exact {s['p_d']:.4g} | frozen {r['p1_red']:.4g} ({r['ratio_d']:.3f}x)"
              f" | memory {m['p_d_ode']:.4g} +- {m['se']:.2g} ({r['ratio_ode']:.3f}x)"
              f"   p0 frozen/exact {r['ratio_0']:.3f}", flush=True)
    global TAU
    TAU = TAU / 2                                       # rule 13: the reduction's own step
    A, B = element("E-070"), element("E-084")
    rep["X:40/tau_half"] = memory_reduction(A, B, 40)
    TAU = TAU * 2
    args.out.write_text(json.dumps(rep, indent=2, default=float))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
