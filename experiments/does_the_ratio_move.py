"""§110 (T-CASC-x) -- does the fast/frozen POSITION track the timescale ratio it is named for?

§109 left the drift sharply posed. With matched seeding, stage 2's effective escape rate runs
from AT the fast limit at Omega = 14 (position ~0.0) to PAST the frozen limit at Omega = 70
(position ~1.16): the system traverses the entire fast-to-frozen bracket as Omega grows, and
exits it.

The framing inherited from §92/§93 says that position is set by a TIMESCALE RATIO -- how fast
the upstream fluctuates against how long the downstream's escape takes to happen. If that is
right, the ratio must move the way the position moves. **The naive version of the argument
predicts the opposite**: the upstream's correlation time is a macroscopic, Omega-independent
quantity, while the downstream's escape time grows like exp(A*Omega), so the upstream should
become RELATIVELY FASTER as Omega grows and the position should move toward the FAST end. It
moves the other way.

That naive version compares against the wrong downstream time. What the upstream has to be
fast against is not the WAITING time between escapes but the TRAVERSAL time of one escape --
how long the downstream spends crossing, during which a fluctuating input can help or hinder.
That is §105's conditional first-passage time, and it is not obviously Omega-independent.

Both quantities are 1-D and neither needs a joint solve.

PREDICTIONS, WRITTEN BEFORE RUNNING.

  P1  THE UPSTREAM CLOCK. tau_up is the intra-basin relaxation time of stage 1 -- the
      sub-dominant relaxation of its generator restricted above the saddle, NOT its escape
      rate. It is a macroscopic quantity, so it must be nearly Omega-INDEPENDENT. §100 found
      exactly this for the walled stage (gap 1.13x over Omega = 14-55) and it must hold here
      too. If tau_up instead moves with Omega by more than ~20%, the "upstream clock" is not a
      single number and the whole framing needs restating before anything else is read.

  P2  THE DOWNSTREAM CLOCK. tau_cross is the CONDITIONAL traversal time from the rail to the
      saddle, given the trajectory gets there rather than falling back (§105's h-transform).
      I predict it GROWS with Omega -- slowly, roughly like log(Omega) -- because a sharper
      instanton spends longer near the saddle where the drift vanishes. Reported as a curve.

  P3  THE DISCRIMINATOR, and it can print either verdict. The measured position moves from
      ~0.0 to ~1.16 as Omega goes 14 -> 70, i.e. toward FROZEN. Resonant activation says frozen
      means the upstream is SLOW compared to the crossing, so tau_up / tau_cross must RISE
      with Omega by something like the same factor. **If instead the ratio falls, or is flat,
      the position is not set by this ratio and the §92/§93 framing does not carry over to the
      escape rate at all** -- which would retire the last inherited explanation for the drift,
      the way §108/§109 retired the averaging family.

  P4  IF THE RATIO DOES MOVE, HOW WELL? Report position against ratio directly rather than
      asserting a correlation from two endpoints (rule 18 -- read the cells). Five Omegas is
      few; a monotone co-movement across all five is suggestive, three out of five is not.

WHAT THIS CANNOT SETTLE. It tests one named candidate. A refutation does not say what DOES set
the position, and P3's wording is deliberately about direction rather than magnitude -- the
bracket's own width changes with Omega (2.39 -> 1.83), so a position and a ratio cannot be
compared as numbers, only as trends.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import experiments.chemical_cascade as cc
from experiments.chain_without_a_joint_solve import chain_operating_points
from experiments.depth_compounding import C, R2, R3

OMEGAS = (14, 20, 30, 55, 70)


def _stage_generator(om, x_up, first=False, cap_mult=1.25):
    cap = int(np.ceil(cap_mult * R3 * om))
    m = cap + 1
    rows, cols, vals = [], [], []
    for s in range(m):
        tot = 0.0
        l, u = cc.rates_stage(float(s), x_up * om, om, C, R3, first, "hill")
        if s < cap and l > 0:
            rows.append(s); cols.append(s + 1); vals.append(l); tot += l
        if s > 0 and u > 0:
            rows.append(s); cols.append(s - 1); vals.append(u); tot += u
        rows.append(s); cols.append(s); vals.append(-tot)
    return sp.csr_matrix((vals, (rows, cols)), shape=(m, m)), m, cap


def upstream_relaxation(om):
    """Intra-basin relaxation time of stage 1: 1/(lambda_2 - lambda_1) above its saddle.

    lambda_1 is the QSD decay (the escape). What a downstream stage sees fluctuating is the
    NEXT mode -- how fast the input forgets where it was inside the high basin.
    """
    Q, m, _ = _stage_generator(om, 0.0, first=True)
    keep = np.where(np.arange(m) > R2 * om)[0]
    A = Q[keep][:, keep].toarray()
    ev = np.sort(-np.real(np.linalg.eigvals(A)))
    return 1.0 / float(ev[1] - ev[0]), float(ev[0]), float(ev[1])


def downstream_crossing(om, cap_mult=1.25):
    """Conditional traversal time of STAGE 2, rail -> saddle, given it gets there.

    §105's `conditional_traversal` hardwires stage-1 rates (no upstream); stage 2's landscape
    depends on its input, so the h-transform is redone here with the upstream pinned at stage
    1's operating point. The start is ONE SITE below the rail: the rail is itself the other
    absorbing boundary, so starting exactly on it is degenerate.
    """
    mus, _ = chain_operating_points(om, 2)
    x_up = mus[0]
    a = int(np.ceil(R2 * om))                 # target: the saddle
    b = int(round(R3 * om))                   # the rail -- the other absorbing end
    n0 = b - 1
    idx = np.arange(a, b + 1)
    m = len(idx)
    L = np.zeros((m, m))
    for i, s in enumerate(idx):
        if s in (a, b):
            L[i, i] = 1.0
            continue
        lam, mu = cc.rates_stage(float(s), x_up * om, om, C, R3, False, "hill")
        L[i, i] = -(lam + mu)
        L[i, i + 1] = lam
        L[i, i - 1] = mu
    rhs_h = np.zeros(m); rhs_h[0] = 1.0
    h = np.linalg.solve(L, rhs_h)
    rhs_v = -h.copy(); rhs_v[0] = 0.0; rhs_v[-1] = 0.0
    v = np.linalg.solve(L, rhs_v)
    j = list(idx).index(n0)
    return float(v[j] / h[j]), float(h[j]), x_up


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/does_the_ratio_move.json"))
    args = ap.parse_args()

    # §109's measured positions, matched seeding (results/it_was_the_seed.json p3)
    src = json.loads(pathlib.Path("results/it_was_the_seed.json").read_text())
    pos = {}
    for r in src["p3"]:
        top = r["bracket_top"]
        pos[r["omega"]] = float(np.log(r["qsd"]) / np.log(top))

    print("P1/P2 -- the two clocks, both 1-D, neither needing a joint solve")
    print(f"{'Omega':>7}{'tau_up':>11}{'lam1 (escape)':>15}{'lam2':>10}"
          f"{'tau_cross':>12}{'ratio':>10}{'position':>11}")
    rows = []
    for om in OMEGAS:
        tau_up, l1, l2 = upstream_relaxation(om)
        tau_x, h, mu = downstream_crossing(om)
        rows.append({"omega": om, "tau_up": tau_up, "lam1": l1, "lam2": l2,
                     "tau_cross": tau_x, "split_prob": h,
                     "ratio": tau_up / tau_x, "position": pos.get(om, float("nan"))})
        print(f"{om:>7}{tau_up:>11.4f}{l1:>15.4e}{l2:>10.4f}{tau_x:>12.4f}"
              f"{tau_up / tau_x:>10.4f}{pos.get(om, float('nan')):>11.4f}", flush=True)
        args.out.write_text(json.dumps(rows, indent=2, default=float))

    ups = [r["tau_up"] for r in rows]
    xs = [r["tau_cross"] for r in rows]
    rt = [r["ratio"] for r in rows]
    ps = [r["position"] for r in rows]

    print(f"\nP1: tau_up spans {max(ups)/min(ups):.4f}x over Omega = {OMEGAS[0]}-{OMEGAS[-1]}"
          f"   -> {'HOLDS' if max(ups)/min(ups) < 1.2 else 'FAILS'} (macroscopic)")
    print(f"P2: tau_cross spans {max(xs)/min(xs):.4f}x, "
          f"{'rising' if xs == sorted(xs) else 'not monotone'}")

    print("\nP3 -- does the ratio move the way the position moves?")
    print(f"  position {ps[0]:.4f} -> {ps[-1]:.4f}   (toward FROZEN)")
    print(f"  ratio    {rt[0]:.4f} -> {rt[-1]:.4f}   "
          f"({'rising' if rt[-1] > rt[0] else 'FALLING'})")
    co = (rt == sorted(rt)) == (ps == sorted(ps))
    print(f"  both monotone in the same direction: {co}"
          f"   -> P3 {'consistent with resonant activation' if co else 'REFUTES it'}")

    print("\nP4 -- the cells, not the endpoints")
    for r in rows:
        print(f"  Omega {r['omega']:>3}: ratio {r['ratio']:.4f}   position {r['position']:.4f}")

    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
