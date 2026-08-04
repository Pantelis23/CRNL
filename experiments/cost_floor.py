"""Is there a FLOOR on the cost of reliability? — the gamma -> 0 limit

`cost_of_reliability` measured R = Sigma/L, the entropy produced per nat of reliability
bought, and found it monotone in gamma -- confirming §9.2's withdrawal, no interior
minimum, diverging to 892 k_B/nat as gamma -> gamma_c. **But at the small-gamma end it
FLATTENS: 16.156 at gamma = 0.05 against 16.329 at gamma = 0.10.**

That flattening is the whole question, because it should not happen for free. Two things
move in opposite directions as the drive rises:

  * the entropy per dissipative jump grows like ln(1/gamma) = A/3, since the reverse
    propensity that appears in ln(a_fwd/a_rev) is vanishing;
  * the barrier `c` grows too, because a smaller gamma deepens the landscape.

R ~ (jumps x ln(1/gamma)) / (c Omega), so a finite limit requires those to CANCEL. If
they do, **that limit is the minimum thermodynamic cost of reliability for this network**
-- a floor that no amount of drive can beat, which is exactly the quantity the founding
question is about. If instead R diverges logarithmically, there is no floor and
reliability is bought by paying unbounded affinity per cycle.

TWO INSTRUMENT FIXES, both needed before the limit means anything:

  (1) R AT THE LARGEST Omega IS NOT R. Both sides carry Omega-independent offsets,
      L = c*Omega + L0 and Sigma = s*Omega + Sigma0, so R(Omega) = (s Omega + Sigma0) /
      (c Omega + L0) and the asymptotic cost is **R_inf = s/c from SEPARATE linear fits**,
      not the ratio at one Omega. The previous run's drifts were 1-19%, so this matters.
  (2) The gamma range must reach far enough below 0.05 to tell a flattening from a slow
      logarithm. ln(1/gamma) only doubles between gamma = 0.05 and 0.0025.

PREDICTIONS, written before running:

  P1  GATE. L and Sigma are each LINEAR in Omega with high R^2 -- that is what makes
      R_inf = s/c well defined. Curvature in either means the offsets are not constant
      and the extrapolation is not licensed.
  P2  R_inf differs from R at the largest Omega by roughly the drifts already seen
      (1-19%), and is the number that should have been quoted.
  P3  THE TEST. R_inf approaches a FINITE LIMIT as gamma -> 0. **I expect this**, on the
      strength of 16.156 against 16.329, and if it holds the limit is the minimum cost of
      a nat of reliability in AM.
  P4  IF R_inf instead grows like ln(1/gamma) at small gamma, there is no floor. The
      apparent flattening over 0.05-0.10 would then be the crossover of two competing
      logarithms, which is precisely the trap THEORIES §4 records as "one spare parameter
      can eat a logarithm" -- and the range 0.05..0.10 is far too short to tell them
      apart, which is why this run goes to 0.0025.
  P5  The cancellation is reported COMPONENTWISE -- s (entropy per molecule) and c
      (nats per molecule) separately against gamma -- so that a finite limit is visibly a
      cancellation between two diverging quantities rather than an asserted one.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from experiments.cost_of_reliability import cell


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.0025, 0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32])
    ap.add_argument("--omegas", type=int, nargs="+",
                    default=[120, 200, 280, 360, 440, 520])
    ap.add_argument("--eps-frac", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/cost_floor.json"))
    args = ap.parse_args()

    t0 = time.time()
    print("R_inf = s/c from SEPARATE linear fits of Sigma(Omega) and L(Omega)")
    print(f"{'gamma':>8}{'A/3=ln(1/g)':>13}{'c (nats/mol)':>14}{'s (k_B/mol)':>13}"
          f"{'R_inf':>10}{'R at max Om':>13}{'R^2 L':>9}{'R^2 S':>9}")
    rows = []
    for g in args.gammas:
        L, S, om = [], [], []
        for o in args.omegas:
            try:
                r = cell(g, o, args.eps_frac, args.theta)
            except Exception:
                continue
            if np.isfinite(r["L"]) and np.isfinite(r["Sigma"]) and r["L"] > 0:
                L.append(r["L"]); S.append(r["Sigma"]); om.append(float(o))
        if len(om) < 4:
            print(f"{g:>8.4f}   SKIPPED ({len(om)} usable cells)")
            continue
        om, L, S = np.array(om), np.array(L), np.array(S)
        cL = np.polyfit(om, L, 1); cS = np.polyfit(om, S, 1)
        r2L = 1 - np.var(L - np.polyval(cL, om)) / np.var(L)
        r2S = 1 - np.var(S - np.polyval(cS, om)) / np.var(S)
        Rinf = cS[0] / cL[0]
        rows.append({"gamma": g, "c": float(cL[0]), "s": float(cS[0]),
                     "R_inf": float(Rinf), "R_last": float(S[-1] / L[-1]),
                     "r2_L": float(r2L), "r2_S": float(r2S),
                     "L0": float(cL[1]), "S0": float(cS[1])})
        print(f"{g:>8.4f}{np.log(1/g):>13.4f}{cL[0]:>14.6f}{cS[0]:>13.6f}"
              f"{Rinf:>10.4f}{S[-1]/L[-1]:>13.4f}{r2L:>9.6f}{r2S:>9.6f}")

    print(f"\n=== P1 gate: are L and Sigma linear in Omega?")
    worst = min(min(r["r2_L"], r["r2_S"]) for r in rows)
    print(f"  worst R^2 across all cells: {worst:.6f}"
          f"   -> P1 {'HOLDS' if worst > 0.999 else 'FAILS'}")

    print(f"\n=== P3/P4: does R_inf converge as gamma -> 0?")
    rows.sort(key=lambda r: r["gamma"])
    print(f"{'gamma':>8}{'R_inf':>10}{'change':>10}")
    for i, r in enumerate(rows):
        ch = "" if i == 0 else f"{100*(r['R_inf']-rows[i-1]['R_inf'])/rows[i-1]['R_inf']:+.2f}%"
        print(f"{r['gamma']:>8.4f}{r['R_inf']:>10.4f}{ch:>10}")
    small = [r for r in rows if r["gamma"] <= 0.04]
    if len(small) >= 3:
        x = np.log([1 / r["gamma"] for r in small])
        y = np.array([r["R_inf"] for r in small])
        p = np.polyfit(x, y, 1)
        r2 = 1 - np.var(y - np.polyval(p, x)) / np.var(y)
        print(f"\n  R_inf vs ln(1/gamma) over the smallest gammas: "
              f"slope {p[0]:+.4f}  R^2 {r2:.4f}")
        print(f"    a slope near 0 means a FLOOR (P3); a clear positive slope means "
              f"logarithmic growth (P4)")
        spread = 100 * (y.max() - y.min()) / y.mean()
        print(f"    R_inf spread over that range: {spread:.2f}% "
              f"while ln(1/gamma) changes by {100*(x.max()-x.min())/x.mean():.0f}%")
        print(f"  -> {'P3: a FLOOR at R_inf ~ ' + f'{y.mean():.3f} k_B per nat' if abs(p[0]) < 0.15 else 'P4: no floor, R_inf grows logarithmically'}")

    print(f"\n=== P5: the cancellation, componentwise")
    print(f"{'gamma':>8}{'ln(1/g)':>10}{'c':>12}{'s':>12}{'s/ln(1/g)':>12}")
    for r in rows:
        print(f"{r['gamma']:>8.4f}{np.log(1/r['gamma']):>10.4f}{r['c']:>12.6f}"
              f"{r['s']:>12.6f}{r['s']/np.log(1/r['gamma']):>12.6f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
