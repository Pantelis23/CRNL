"""§123 -- the exact p_d as a joint-chain committor: a route past the window to Omega > Omega_c.

§122's universal-floor suspect (T-CASC-ap) needs an EXACT p_d in a narrow-zone pair beyond its
crossover Omega_c = lam(s1)/(ls w^2) -- 227 for Z -- where evolving §119's whole window would need
~10^7 tagged states. A dip is an episode: it starts when stage 1 steps below its saddle, and ends
when stage 1 commits -- to its high rail or its low one. So, exactly on the joint chain,

    p_d_ep = sum_n2 J(n2) h(sad1-1, n2) / sum_n2 J(n2) c(sad1-1, n2),

  J(n2)    the stage-2 law AT the dip: the crossing flux sad1 -> sad1-1 out of the quasi-stationary
           law of the joint chain with stage 1 alive (n1 >= sad1) -- so the pre-dip sag is in it;
  c        P(stage 1 reaches its high rail before its low rail), from that state, on the joint chain;
  h        P(the same, AND stage 2 below its saddle at that moment).

c and h are each one sparse linear solve. Nothing is evolved, and no window enters.

PREDICTIONS, WRITTEN BEFORE RUNNING.

  E1  VALIDATION. At Omega = 40 and 55, p_d_ep reproduces §119's stored windowed exact p_d for X
      and Z within 10%. (They differ in definition: the window counts trajectories high at t = 16,
      which may dip again; the episode ends at commitment. If they disagree by more, the episode is
      not the window's quantity and E2 cannot be read against §122.)
  E2  THE KILL TEST (T-CASC-ap), conditional on E1. Z at Omega = 230 and 400 (past Omega_c = 227):
      the exact p_d_ep RISES from its Omega = 55 value, turning the decline that §119 measured
      (0.066 -> 0.054). If it keeps falling past Omega_c, the universal-floor suspect is refuted in
      the exact chain, whatever the reductions said.
  E3  phi_ep = sum J c / sum J -- the episode's upstream factor -- equals the stage-1 committor to
      within 2%: stage 2 does not feed back on stage 1, so it must.

  E4  Written after Z:800 returned p_d_ep = 0.44 with lam0 = -5e-15 and J_low_frac = 0.30 (against
      3e-8 at 400) -- the alive chain has two near-absorbing pockets (stage 2 high / low), and once
      both decay rates are below double precision the sigma = 0 eigenvector is an arbitrary mixture.
      Replacement with no eigenproblem: J from the OCCUPATION before the first dip, starting at the
      healthy product seed, x = (-Q_alive^T)^{-1} p_seed -- exactly stage 2's law at the first dip.
      Predict it reproduces the eigenvector-J p_d_ep within 1% at Omega <= 230, and returns
      J_low_frac at 800 to the ~1e-8 scale. Z:800 is read for E2 only from this J.

  E5  Written after E4's occupation solve was itself found broken at Omega >= 200 (relative
      residual ~1e2, total occupation FALLING from 200 to 400 where 1/k1 grows, J's shape jumping):
      stage 1 is autonomous, so the alive chain has TWO eigenvalues within a whisker of zero (stage 1
      dying, and that plus stage 2's own slow switching), and any solve near zero mixes them. All
      cells with lam0 <~ 1e-12 -- X at Omega >= 200, Z at Omega >= 400 -- are withdrawn, including
      E2's "still falling at 800" and X's non-monotone 0.596 / 0.538 / 0.674.
      Fix: make stage-2-low absorbing too. The chain alive AND stage 2 high has ONE tiny eigenvalue,
      O(1)-separated from the rest, so shift-invert is well-conditioned. It omits the fraction of
      stage 2 already low at the dip (<= 1.5% wherever it was reliably measured), carried as an
      error bar. Predict: matches the validated values within 1.5% at X <= 100, Z <= 230; and its
      own residual is <= 1e-8 at every Omega it is run at.

OUTCOMES (appended after running; every prediction above left as written). E1 held (window vs
episode within 6%). E3 held (phi_ep = stage-1 committor to 4 digits). E4's first half held at
Omega >= 55 (1.3% at 40); its second half was moot -- the occupation route itself broke at
lam0 <~ 1e-12. E5 held (within 1.2% of every validated cell; residuals <= 1.4e-10 at every Omega).
E2 FAILED as written: exact Z keeps falling through Omega = 400 = 1.76 Omega_c (0.065 -> 0.032);
Z:800 was not computed (stopped for system memory). Exact X rises monotonically to 0.807 at
Omega = 800 and phi * p_d to 0.380. FINDINGS §123.

  E6  Written before Z:800 ran (on the user's 64 GB laptop, after it reproduced Z:230/hi to 1.1e-10
      relative). T-CASC-ap needs exact Z to turn eventually; the reduction turned near Omega ~ 800
      and exact/reduction is falling. Predict Z:800 p_d >= 0.0306 (no more than 5% below Z:400's
      0.0322): flattening or rising. A fall of more than 10% (below 0.0290) weakens T-CASC-ap for
      narrow zones at accessible Omega; it would not refute a limit, but it would say Omega_c is
      not where the exact chain turns.
      -> HELD, and then some: Z:800 = 0.186 (from 0.032 at 400), eig_res 1.1e-10, on the laptop
         after it reproduced Z:230/hi to 1.1e-10. A 5.8x step is itself suspect (exact/reduction
         went 2.35 -> ~14), hence E7.

  E7  Written before running. If the turn is real it is smooth in Omega: Z:600 lies strictly
      between Z:400 and Z:800 (0.032 < p_d < 0.186), and Z:1000 lies above Z:800. Z:600 outside
      that band, or Z:1000 below Z:800, marks the 800 value as an isolated spike -- withdrawn
      until an independent route reproduces it.
      -> HELD: Z:600 = 0.0354, Z:1000 = 0.562 (laptop; eig_res <= 1.5e-10, |lambda_1| ~ 2.3).
"""


from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from experiments.weak_upstream import ARMS, element, pair_generator


def qsd_alive(Q, m, sad1):
    """Left principal eigenvector of the joint chain restricted to stage 1 alive (n1 >= sad1)."""
    N = Q.shape[0]
    alive = np.nonzero(np.arange(N) // m >= sad1)[0]
    Qa = Q[alive][:, alive].T.tocsc()
    vals, vecs = spla.eigs(Qa, k=1, sigma=0.0, which="LM")
    p = np.abs(np.real(vecs[:, 0]))
    full = np.zeros(N)
    full[alive] = p / p.sum()
    return full, float(-np.real(vals[0]))


def occupation_alive(Q, m, sad1, seed_vec):
    """Expected occupation of each alive state before stage 1 first steps below its saddle."""
    N = Q.shape[0]
    alive = np.nonzero(np.arange(N) // m >= sad1)[0]
    Qa = Q[alive][:, alive].T.tocsc()
    x = spla.spsolve(-Qa, seed_vec[alive])
    full = np.zeros(N)
    full[alive] = np.abs(x)
    return full / full.sum(), float("nan")


def qsd_alive_high(Q, m, sad1, sad2):
    """QSD of the chain with stage 1 alive AND stage 2 high: one tiny eigenvalue, well separated."""
    N = Q.shape[0]
    idx = np.arange(N)
    keep = np.nonzero((idx // m >= sad1) & (idx % m >= sad2))[0]
    Qk = Q[keep][:, keep].T.tocsc()
    vals, vecs = spla.eigs(Qk, k=2, sigma=0.0, which="LM")
    order = np.argsort(np.abs(vals))
    v = np.real(vecs[:, order[0]])
    v = v * np.sign(v.sum())
    res = float(np.abs(Qk @ v - np.real(vals[order[0]]) * v).max() / np.abs(v).max())
    full = np.zeros(N)
    full[keep] = np.clip(v, 0.0, None)
    gap = float(np.abs(vals[order[1]]) / max(np.abs(vals[order[0]]), 1e-300))
    return full / full.sum(), float(-np.real(vals[order[0]])), res, gap


def episode(A, B, om, method="qsd"):
    from experiments.weak_upstream import seed
    Q, m = pair_generator(A, B, om)
    N = m * m
    sad1, sad2 = int(np.ceil(A.r2 * om)), int(np.ceil(B.r2 * om))
    lo1, hi1 = int(round(A.r1 * om)), int(round(A.r3 * om))
    t0 = time.time()
    if method == "qsd":
        pi, lam0 = qsd_alive(Q, m, sad1)
    elif method == "occ":
        pi, lam0 = occupation_alive(Q, m, sad1, seed(A, B, om, m)[0])
    else:
        pi, lam0, eig_res, eig_gap = qsd_alive_high(Q, m, sad1, sad2)
    n2 = np.arange(m)
    _, mu1 = A.rates(np.full(m, float(sad1)), np.zeros(m), om, True)
    J = pi.reshape(m, m)[sad1, :] * mu1                     # flux sad1 -> sad1-1, by n2
    J /= J.sum()
    # backward problems on interior n1 in (lo1, hi1); boundaries n1 = lo1 (0) and n1 = hi1 (target)
    idx = np.arange(N)
    n1 = idx // m
    inner = np.nonzero((n1 > lo1) & (n1 < hi1))[0]
    Qi = Q[inner][:, inner].tocsc()
    top = np.nonzero(n1 == hi1)[0]
    Qt = Q[inner][:, top]
    rhs_c = -np.asarray(Qt.sum(axis=1)).ravel()
    low2 = (top % m) < sad2
    rhs_h = -np.asarray(Qt[:, low2].sum(axis=1)).ravel()
    lu = spla.splu(Qi)
    c, h = lu.solve(rhs_c), lu.solve(rhs_h)
    pos = {s: i for i, s in enumerate(inner)}
    start = [pos[(sad1 - 1) * m + j] for j in range(m)]
    cs, hs = c[start], h[start]
    return {"omega": om, "states": int(len(inner)), "secs": time.time() - t0,
            "phi_ep": float((J * cs).sum()), "p_d_ep": float((J * hs).sum() / (J * cs).sum()),
            "J_low_frac": float(J[:sad2].sum()), "lam0": lam0,
            **({"eig_res": eig_res, "eig_gap": eig_gap} if method == "hi" else {})}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("results/episode_committor.json"))
    ap.add_argument("--cells", default="X:40,X:55,Z:40,Z:55")
    ap.add_argument("--method", default="hi", choices=["qsd", "occ", "hi"],
                    help="hi (default, E5): the only route valid at every Omega; qsd/occ kept as the record of E4's failure")
    args = ap.parse_args()
    rep = json.loads(args.out.read_text()) if args.out.exists() else {}
    for cell in args.cells.split(","):
        key = cell if args.method == "qsd" else f"{cell}/{args.method}"
        if key in rep:
            continue
        arm, om = cell.split(":")
        A, B = element(ARMS[arm][0]), element(ARMS[arm][1])
        r = episode(A, B, int(om), args.method)
        rep[key] = r
        print(key, json.dumps(r), flush=True)
        args.out.write_text(json.dumps(rep, indent=2, default=float))


if __name__ == "__main__":
    main()
