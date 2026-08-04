"""The cost of restoration is per E-FOLD OF GAIN, and that one is preparation-free

§37 measured R = Sigma/L, entropy per nat of reliability, and found it varies 203% with
protocol. **The reason is structural, not a sensitivity.** Measured directly at fixed
gamma and Omega, as the input margin eps rises:

    Sigma FALLS   1094 -> 205 k_B      (a start nearer the threshold needs fewer
    L     RISES   15.0 -> 101.7 nats    reactions, and is more reliable)

so R = Sigma/L collapses by a factor of **36** for a trivial reason. **"Cost per nat of
reliability" is not a quantity.** Reliability is bought with input margin, which is free;
dissipation buys something else. §37's withdrawal of the universal cost was right, and
this is why.

**WHAT DISSIPATION ACTUALLY BUYS IS GAIN.** A restoring switch takes an input margin
eps*delta* and delivers theta*delta* -- it AMPLIFIES. In the same table,

    Sigma / ln(theta/eps)   varies only 10.9% (gamma=0.07) and 16.4% (gamma=0.25)

against R's 3600%. Entropy scales with the LOGARITHM of the gain, not the distance
travelled -- which is what an exponential amplifier should do, since it traverses margin
multiplicatively. So the candidate invariant is

    **G = Sigma / ( Omega * ln(theta/eps) )   [k_B per molecule per e-fold of gain]**

and it is preparation-free BY CONSTRUCTION: eps enters only through the gain it defines.

**THE RESIDUAL 11-16% IS EXPECTED AND IS THE THING TO REMOVE.** Sigma carries an
Omega-independent offset just as L and Sigma did in §37, so Sigma = G*Omega*ln(gain) +
Sigma_0 and dividing rather than fitting leaves the offset drifting through ln(gain).
The honest estimator is the SLOPE of Sigma against Omega*ln(theta/eps), fitted jointly
across eps AND Omega, which is what this measures.

PREDICTIONS, written before running:

  P1  GATE. Sigma is LINEAR in the single product `Omega * ln(theta/eps)` jointly over a
      grid of eps and Omega -- one slope, one intercept, high R^2. **That is what makes G
      a quantity at all.** If the data need separate eps and Omega terms, restoration has
      no single cost and the whole framing joins R in the bin.
  P2  G is margin-independent and Omega-independent by construction if P1 holds, and the
      residual scatter should fall well below the 11-16% seen when dividing instead of
      fitting.
  P3  G(gamma) is then the founding question's number: **the free energy a restoring
      switch spends, per molecule, to amplify a signal margin by a factor e.** Reported
      as measured, and NOT decorated against recognisable constants -- §28.2's power law
      and §35.1's -1/2 both came from reading structure into a fitted quantity.
  P4  G has a minimum in gamma. **If its location differs from §37's gamma* ~ 0.07, then
      §37 optimised the wrong quantity** and the design principle moves. That is the
      outcome that costs a result published an hour ago, and it is the reason to run this
      at the same gammas.
  P5  If P1 FAILS -- if Sigma needs eps and Omega separately -- then the cost of
      restoration is irreducibly preparation-dependent, and the honest conclusion is that
      the founding question has no single-number answer. That would be a real finding and
      it is reported as one, not as a failure to measure.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from experiments.cost_of_reliability import cell


def grid(gamma, epss, omegas, theta):
    X, Y, meta = [], [], []
    for eps in epss:
        for om in omegas:
            try:
                r = cell(gamma, om, eps, theta)
            except Exception:
                continue
            if not np.isfinite(r["Sigma"]):
                continue
            X.append(om * np.log(theta / eps))
            Y.append(r["Sigma"])
            meta.append({"eps": eps, "omega": om, "Sigma": r["Sigma"],
                         "L": r["L"], "gain": float(np.log(theta / eps))})
    return np.array(X), np.array(Y), meta


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.02, 0.04, 0.07, 0.12, 0.20, 0.30, 0.40])
    ap.add_argument("--epss", type=float, nargs="+",
                    default=[0.20, 0.30, 0.40, 0.50, 0.60])
    ap.add_argument("--omegas", type=int, nargs="+", default=[150, 250, 350, 450])
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/cost_per_efold.json"))
    args = ap.parse_args()

    t0 = time.time()
    print("G = slope of Sigma against Omega*ln(theta/eps), fitted jointly over eps and Omega")
    print(f"{'gamma':>7}{'cells':>7}{'G (k_B/mol/e-fold)':>21}{'intercept':>12}"
          f"{'R^2':>10}{'naive spread':>14}")
    rows = []
    for g in args.gammas:
        X, Y, meta = grid(g, args.epss, args.omegas, args.theta)
        if len(X) < 8:
            print(f"{g:>7.2f}   SKIPPED ({len(X)} cells)")
            continue
        p = np.polyfit(X, Y, 1)
        res = Y - np.polyval(p, X)
        r2 = 1 - res.var() / Y.var()
        naive = np.array([m["Sigma"] / (m["omega"] * m["gain"]) for m in meta])
        spread = 100 * (naive.max() - naive.min()) / naive.mean()
        rows.append({"gamma": g, "G": float(p[0]), "intercept": float(p[1]),
                     "r2": float(r2), "n": len(X), "naive_spread": float(spread),
                     "rms": float(np.sqrt((res ** 2).mean()))})
        print(f"{g:>7.2f}{len(X):>7}{p[0]:>21.5f}{p[1]:>12.2f}{r2:>10.6f}"
              f"{spread:>13.1f}%")

    print(f"\n=== P1 gate: is Sigma linear in the single product?")
    worst = min(r["r2"] for r in rows)
    print(f"  worst R^2 across gammas: {worst:.6f}"
          f"   -> P1 {'HOLDS -- G is a quantity' if worst > 0.995 else 'FAILS -- the cost is irreducibly preparation-dependent (P5)'}")

    print(f"\n=== P2: does fitting beat dividing?")
    for r in rows:
        print(f"  gamma={r['gamma']:.2f}: naive Sigma/(Om ln gain) spreads "
              f"{r['naive_spread']:.1f}%, fitted slope R^2 = {r['r2']:.6f}")

    print(f"\n=== P3/P4: the cost of one e-fold of restoration, against the drive")
    print(f"{'gamma':>7}{'affinity A':>12}{'G':>12}")
    for r in rows:
        print(f"{r['gamma']:>7.2f}{-3*np.log(r['gamma']):>12.4f}{r['G']:>12.5f}")
    Gs = np.array([r["G"] for r in rows]); gs = np.array([r["gamma"] for r in rows])
    i = int(np.argmin(Gs))
    interior = 0 < i < len(Gs) - 1
    print(f"\n  minimum G = {Gs[i]:.5f} k_B per molecule per e-fold at gamma = {gs[i]:.2f}")
    print(f"  interior minimum? {'YES' if interior else 'NO -- monotone'}")
    print(f"  §37's optimum for R was gamma* ~ 0.07"
          f"   -> {'SAME location' if abs(gs[i]-0.07) < 0.05 else 'DIFFERENT -- §37 optimised the wrong quantity (P4)'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
