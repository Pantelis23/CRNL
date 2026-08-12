"""T15-l: is §63's width exponent -1/2 the PITCHFORK, or a coincidence?

§63 measured the width of the restoration threshold as w ~ Omega^(-0.497), from an
instrument that mentions no reference at all. That exponent is not free: if the escape
action vanishes at the boundary as

    A(gamma) ~ (gamma_c - gamma)^nu

then the gamma at which Omega*A = k is gamma_c - (k/Omega)^(1/nu), so the width between
two fixed metastability levels goes as **Omega^(-1/nu)**. §63's -0.497 therefore says
**nu = 2**, and nothing in §63 tested that independently -- it is the same measurement
read a second way.

**BUT nu = 2 IS PREDICTED, AND BY A SECTION THAT PREDATES ALL OF THIS.** §9.1 established
that the AM landscape loses its two rails at gamma_c = 1/2 through a PITCHFORK, measuring
delta* proportional to sqrt(gamma_c - gamma) and counting exactly 3 fixed points below and 1
above from 1830 simplex starts. A pitchfork normal form U = a x^2/2 + b x^4/4 with
a proportional to (gamma_c - gamma) has barrier a^2/(4b), so

    **A(gamma) proportional to (gamma_c - gamma)^2,  i.e. nu = 2,  i.e. w ~ Omega^(-1/2)**

follows from §9.1 alone, with no reference to §63. So §63's exponent is a PREDICTION of the
deterministic bifurcation type, and this measures whether the prediction is met.

PREDICTIONS, written before running.

  P1  GATE, and rule 13 first. A(gamma) is extracted as -ln|lambda_A|/Omega, which is an
      Omega -> infinity limit. **Convergence WITHIN Omega must be checked before comparing
      ACROSS gamma.** Report A at successive Omega per gamma; any gamma whose A has not
      settled to 1% is excluded and counted, not fitted.

      SECOND PASS, before re-running. **The gate fired and refused every gamma**: the ratio
      estimator -ln|lambda_A|/Omega drifted 8-39% across Omega = 120..500, so no nu was
      reported. The cause is not slow convergence of A but the ESTIMATOR: ln|lambda_A| =
      -Omega*A + b*ln(Omega) + c, so dividing by Omega leaves the whole prefactor as an
      O(ln Omega / Omega) contamination. The repair is the LOCAL SLOPE,

          A_eff(W1, W2) = -[ln|lambda_A(W2)| - ln|lambda_A(W1)|] / (W2 - W1)

      which cancels c exactly and suppresses b. **A three-parameter fit is deliberately NOT
      used**: §35.3 proved Omega and ln(Omega) are collinear over any bounded window here,
      and §35.1's b values were withdrawn for exactly that. The local slope asks only for
      the leading term, which is the one that survives collinearity. Convergence is then
      demanded of A_eff itself, and the Omega window is chosen PER GAMMA, since deep in the
      metastable phase lambda_A falls under the precision floor (P6) and near gamma_c it
      needs the largest Omega available.
  P2  **THE TEST, absolute (rule 16). nu = 2.000.** Fit ln A against ln(gamma_c - gamma)
      with gamma_c = 1/2 held FIXED at the value §62 proves exactly -- gamma_c is not a
      free parameter here, which is what makes this a test rather than a two-parameter fit.
  P3  **VERDICT RULE (rule 19).** nu is continuous, so "close to 2" is not a criterion and a
      single fit is not either: near-critical fits drift with the window. The criterion is
      **convergence of nu as the window shrinks toward gamma_c**. Report nu for a nested
      sequence of windows. Data that would make it print the other answer: nu drifting AWAY
      from 2, or settling on a different constant. A drift toward 2 that has not arrived is
      reported as such, not rounded.

      SECOND VERSION, and the first was wrong in the way rule 19 names. It compared
      |nu - 2| at the widest window against the narrowest and declared "nu -> 2" if the
      narrowest came in under 0.05. The measured values were 1.9507, 1.9465, 1.9498, 1.9517
      -- **FLAT, non-monotone, scatter 0.005** -- and it printed "P2 HOLDS: nu -> 2" because
      0.0483 < 0.05. A rule that reads a 0.001 wobble as a trend cannot tell a drift from a
      constant, which is the only thing it was there to do. Replaced by a rule that compares
      the TREND against the SCATTER: nu is drifting only if the change across nested windows
      exceeds the spread among them; otherwise nu is flat, and the question is whether 2 lies
      inside the scatter. Data that would now print "-> 2": nested values marching toward 2
      by more than their own spread.
  P4  If nu = 2 holds, §63's Omega^(-1/2) is explained rather than merely measured, and the
      two are one fact -- which must then be SAID, because reporting them as independent
      confirmations would be double-counting one measurement.
  P5  The amplitude is reported but not predicted: the normal-form coefficient b is not
      derived here, so A/(gamma_c-gamma)^2 is a measured constant. **It is named as measured**
      -- §49's precedent of declining pi, and §61's of declining a prefactor that was not
      universal.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.networks.am_reversible import am_reversible
from experiments.threshold_sharpness import lambda_A

GAMMA_C = 0.5


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omegas", type=int, nargs="+", default=[120, 200, 320, 500])
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.20, 0.25, 0.30, 0.34, 0.38, 0.41, 0.43, 0.45, 0.46, 0.47])
    ap.add_argument("--tol", type=float, default=0.01)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/pitchfork_action.json"))
    args = ap.parse_args()
    t0 = time.time()

    print("=== P1 GATE (rule 13): has the LOCAL SLOPE A_eff converged IN Omega?")
    print(f"{'gamma':>7}{'Omega used':>22}{'A_eff series':>34}{'drift%':>9}{'used':>6}")
    A, rows = {}, []
    for g in args.gammas:
        oms, lns = [], []
        for om in args.omegas:
            lam = lambda_A(am_reversible(float(g)), om)
            if lam is None or lam >= 0 or abs(lam) < 1e-12:
                continue                       # P6 floor: deep metastable, unresolvable
            oms.append(om); lns.append(-np.log(abs(lam)))
        if len(oms) < 3:
            print(f"{g:>7.3f}{'too few resolvable':>22}{'':>34}{'--':>9}{'no':>6}")
            continue
        slopes = [(lns[i + 1] - lns[i]) / (oms[i + 1] - oms[i])
                  for i in range(len(oms) - 1)]
        drift = abs(slopes[-1] - slopes[-2]) / max(abs(slopes[-1]), 1e-30)
        used = drift <= args.tol
        A[g] = float(slopes[-1])
        rows.append({"gamma": g, "A": A[g], "drift": float(drift), "used": bool(used),
                     "omegas": oms, "slopes": [float(s) for s in slopes]})
        win = f"{min(oms)}..{max(oms)} ({len(oms)})"
        ser = " ".join(f"{s:.6f}" for s in slopes[-3:])
        print(f"{g:>7.3f}{win:>22}{ser:>34}{100*drift:>9.2f}{'yes' if used else 'NO':>6}")
    good = [r for r in rows if r["used"]]
    print(f"  -> P1: {len(good)} of {len(rows)} gamma converged to {100*args.tol:g}%;"
          f" {len(rows)-len(good)} excluded and counted")

    print(f"\n=== P2/P3: nu from ln A vs ln(1/2 - gamma), gamma_c held FIXED at 1/2")
    print(f"{'window':>22}{'n':>4}{'nu':>10}{'amplitude':>12}")
    nus = []
    for lo in (0.20, 0.30, 0.38, 0.41, 0.43):
        sel = [r for r in good if r["gamma"] >= lo]
        if len(sel) < 3:
            print(f"{f'gamma >= {lo}':>22}{len(sel):>4}   too few")
            continue
        x = np.log(GAMMA_C - np.array([r["gamma"] for r in sel]))
        y = np.log(np.array([r["A"] for r in sel]))
        nu, c = np.polyfit(x, y, 1)
        nus.append((lo, float(nu)))
        hi = max(r["gamma"] for r in sel)
        label = f"gamma in [{lo:.2f}, {hi:.2f}]"
        print(f"{label:>22}{len(sel):>4}{nu:>10.4f}{np.exp(c):>12.4f}")
    if len(nus) >= 3:
        v = np.array([n for _, n in nus])
        print(f"\n  nu across nested windows: " + ", ".join(f"{n:.4f}" for _, n in nus))
        scatter = float(v.max() - v.min())
        trend = float(abs(v[0] - 2.0) - abs(v[-1] - 2.0))     # >0 means moving toward 2
        gap = float(abs(v.mean() - 2.0))
        print(f"  scatter across windows {scatter:.4f};"
              f" net movement toward 2 {trend:+.4f};"
              f" distance of the mean from 2 {gap:.4f} ({gap/max(scatter,1e-9):.1f}x scatter)")
        if trend > scatter:
            print(f"  -> P2 nu is DRIFTING toward 2 (movement exceeds scatter);"
                  f" the limit is not reached at {v[-1]:.4f} and is not rounded to it")
        elif gap <= scatter:
            print(f"  -> P2 HOLDS: nu = 2 within the scatter,"
                  f" so §63's Omega^(-1/2) IS §9.1's pitchfork")
        else:
            print(f"  -> P2 nu is FLAT at {v.mean():.4f} +- {scatter/2:.4f} and 2 is"
                  f" {gap/max(scatter,1e-9):.1f}x the scatter away:"
                  f" **2 is EXCLUDED over this window**")
            print(f"     1/nu = {1/v.mean():.4f}, to be compared with §63's width exponent")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"rows": rows, "nu": nus}, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
