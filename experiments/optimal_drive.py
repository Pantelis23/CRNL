"""An OPTIMAL DRIVE for restoration — located, and tested against the trap that killed §9.2

`cost_floor` found that R = Sigma/L, the entropy produced per nat of reliability, has an
INTERIOR MINIMUM in the drive: 26.01 k_B/nat at gamma = 0.0025, falling to 16.90 at
gamma = 0.04, rising to 73.99 at gamma = 0.32. Both limits diverge and the mechanism is
explicit in the components:

    s / ln(1/gamma) = 0.77..0.80, essentially CONSTANT   (entropy per molecule tracks
                                                          the affinity A/3 = ln(1/gamma))
    c saturates at ~0.19 as gamma -> 0, and collapses to 0 as gamma -> gamma_c

so **R ~ 0.79 * ln(1/gamma) / c(gamma)**, divergent at both ends: drive too hard and each
cycle dissipates ln(1/gamma) while the barrier has already saturated; drive too softly
and the landscape shallows faster than the saving. **The minimum is a design principle --
the drive at which a restoring switch is cheapest per nat of reliability.**

**THIS PROJECT HAS CLAIMED AN OPTIMAL DRIVE BEFORE AND WITHDRAWN IT.** THEORIES §4:
"Dissipation has a minimum near gamma ~ 0.3 -- a clean U-shaped curve", killed because
the decision threshold was held fixed while delta*(gamma) shrank, so above gamma ~ 0.42
the threshold sat OUTSIDE the landscape. That withdrawal stands: the minimum found here
is at gamma ~ 0.04, nearly an order of magnitude away, and §9.2's curve was an artifact.
But the coincidence of shape is exactly why this needs the harder test before it is
believed at all.

**THE TEST THAT DECIDES IT.** A real optimum is a property of the chemistry. A protocol
artifact moves when the protocol moves. So gamma* and R* are measured across a grid of
eps (where the decision starts) and theta (where it is declared finished). If the optimum
is real, both should shift it only weakly; if either drags it, this is §9.2 again with a
different threshold and must be reported as such.

PREDICTIONS, written before running:

  P1  A finer gamma grid with proper Omega extrapolation (R_inf = s/c from separate
      linear fits) locates gamma* and R* to within the grid spacing.
  P2  THE DECIDING TEST. gamma* moves by less than a factor of ~1.5 across
      eps in {0.25, 0.35, 0.50} and theta in {0.70, 0.80, 0.90}, and R* by less than
      ~15%. **If gamma* tracks theta the way §9.2's did, this is the same artifact and
      the result is withdrawn on the spot.**
  P3  The mechanism constant k = s / ln(1/gamma) stays ~0.79 across gamma AND across
      protocol. k is where the affinity enters, so if k is protocol-dependent the
      decomposition is not clean and R's minimum has no simple reading.
  P4  With R* in hand, the founding claim becomes a number: the free energy needed for
      one decision at transistor-grade reliability, P(error) = 1e-15, is
      R* x ln(1e15) k_B T. Reported as an arithmetic consequence, not a new measurement.
  P5  If gamma* sits close to a recognisable value it is reported as a number and NOT
      decorated. §28.2's power law and §35.1's -1/2 both came from reading structure into
      a fitted quantity, and both were withdrawn.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from experiments.cost_of_reliability import cell


def r_inf(gamma, eps, theta, omegas):
    L, S, om = [], [], []
    for o in omegas:
        try:
            r = cell(gamma, o, eps, theta)
        except Exception:
            continue
        if np.isfinite(r["L"]) and np.isfinite(r["Sigma"]) and r["L"] > 0:
            L.append(r["L"]); S.append(r["Sigma"]); om.append(float(o))
    if len(om) < 4:
        return None
    om, L, S = np.array(om), np.array(L), np.array(S)
    cL, cS = np.polyfit(om, L, 1), np.polyfit(om, S, 1)
    r2 = min(1 - np.var(L - np.polyval(cL, om)) / np.var(L),
             1 - np.var(S - np.polyval(cS, om)) / np.var(S))
    return {"c": float(cL[0]), "s": float(cS[0]),
            "R": float(cS[0] / cL[0]), "r2": float(r2)}


def locate(gammas, eps, theta, omegas):
    pts = []
    for g in gammas:
        v = r_inf(g, eps, theta, omegas)
        if v:
            pts.append((g, v))
    if len(pts) < 3:
        return None, None, pts
    Rs = np.array([v["R"] for _, v in pts])
    i = int(np.argmin(Rs))
    if 0 < i < len(pts) - 1:                      # parabolic refinement in ln gamma
        x = np.log([pts[j][0] for j in (i - 1, i, i + 1)])
        y = Rs[[i - 1, i, i + 1]]
        a, b, _ = np.polyfit(x, y, 2)
        gstar = float(np.exp(-b / (2 * a))) if a > 0 else pts[i][0]
        rstar = float(np.polyval([a, b, _], np.log(gstar)))
    else:
        gstar, rstar = pts[i][0], float(Rs[i])
    return gstar, rstar, pts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.015, 0.025, 0.035, 0.05, 0.07, 0.10, 0.14])
    ap.add_argument("--omegas", type=int, nargs="+",
                    default=[120, 200, 280, 360, 440, 520])
    ap.add_argument("--eps-fracs", type=float, nargs="+", default=[0.25, 0.35, 0.50])
    ap.add_argument("--thetas", type=float, nargs="+", default=[0.70, 0.80, 0.90])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/optimal_drive.json"))
    args = ap.parse_args()

    t0 = time.time()
    print("=== P1: locate the optimum at the reference protocol (eps=0.35, theta=0.80)")
    print(f"{'gamma':>8}{'c':>12}{'s':>12}{'R_inf':>10}{'k=s/ln(1/g)':>14}{'R^2':>10}")
    g0, r0, pts = locate(args.gammas, 0.35, 0.80, args.omegas)
    for g, v in pts:
        print(f"{g:>8.4f}{v['c']:>12.6f}{v['s']:>12.6f}{v['R']:>10.4f}"
              f"{v['s']/np.log(1/g):>14.6f}{v['r2']:>10.6f}")
    print(f"\n  gamma* = {g0:.5f}   R* = {r0:.4f} k_B per nat"
          f"   (A* = {-3*np.log(g0):.4f})")

    print(f"\n=== P2 THE DECIDING TEST: does the optimum move with the protocol?")
    print(f"{'eps':>7}{'theta':>7}{'gamma*':>10}{'R*':>10}{'worst R^2':>11}")
    grid = []
    for eps in args.eps_fracs:
        for th in args.thetas:
            gs, rs, p = locate(args.gammas, eps, th, args.omegas)
            if gs is None:
                print(f"{eps:>7.2f}{th:>7.2f}   SKIPPED")
                continue
            w = min(v["r2"] for _, v in p)
            grid.append({"eps": eps, "theta": th, "gamma_star": gs, "R_star": rs,
                         "worst_r2": w})
            print(f"{eps:>7.2f}{th:>7.2f}{gs:>10.5f}{rs:>10.4f}{w:>11.6f}")

    if grid:
        gsv = np.array([x["gamma_star"] for x in grid])
        rsv = np.array([x["R_star"] for x in grid])
        fg = gsv.max() / gsv.min()
        dr = 100 * (rsv.max() - rsv.min()) / rsv.mean()
        print(f"\n  gamma* spans {gsv.min():.5f}..{gsv.max():.5f}  "
              f"-> a factor of {fg:.2f}")
        print(f"  R*     spans {rsv.min():.4f}..{rsv.max():.4f}  -> {dr:.1f}%")
        ok = fg < 1.5 and dr < 15.0
        print(f"  -> P2 {'HOLDS: the optimum is a property of the chemistry' if ok else 'FAILS: the optimum tracks the protocol -- this is §9.2 again'}")

        # does gamma* track theta specifically? that was §9.2's failure mode
        for name, key in (("theta", "theta"), ("eps", "eps")):
            vals = sorted({x[key] for x in grid})
            med = [np.median([x["gamma_star"] for x in grid if x[key] == v])
                   for v in vals]
            print(f"  gamma* vs {name}: "
                  + "  ".join(f"{v}:{m:.5f}" for v, m in zip(vals, med))
                  + f"   spread {100*(max(med)-min(med))/np.mean(med):.1f}%")

    print(f"\n=== P3: is the mechanism constant k = s/ln(1/gamma) protocol-independent?")
    ks = [v["s"] / np.log(1 / g) for g, v in pts]
    print(f"  at the reference protocol, k over gamma: {min(ks):.4f}..{max(ks):.4f}"
          f"  spread {100*(max(ks)-min(ks))/np.mean(ks):.1f}%")

    print(f"\n=== P4: the founding claim, as a number")
    for target, label in ((1e-15, "transistor-grade 1e-15"), (1e-9, "1e-9"),
                          (1e-3, "1e-3")):
        nats = -np.log(target)
        print(f"  P(error) = {label:>22}: {nats:>6.2f} nats x {r0:.2f} = "
              f"{nats*r0:>8.1f} k_B T per decision")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"reference": {"gamma_star": g0, "R_star": r0,
                                                  "points": [(g, v) for g, v in pts]},
                                    "protocol_grid": grid}, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
