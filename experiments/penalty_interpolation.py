"""T-CASC-k: is the position between the two limits a function of the barrier depth alone?

§97 left the composition penalty BRACKETED but not closed: the frozen-upstream average and the
fast-upstream limit are both computable, and the measurement sits between them. Where between is
the open question, and there are two data points:

    calibration element, A*Omega = 5.7 :  fast 1.139, measured 4.442, frozen 4.845
    §97's new element,   A*Omega = 13.2:  fast 1.010, measured 4.474, frozen 8.925

In log space the measurement sits **94%** of the way to the frozen limit on the first and **68%**
on the second. Deeper barrier -> further from the frozen limit, which is what an
excursion-duration argument would give: a deeper barrier needs a deeper and therefore RARER and
SHORTER upstream excursion, so fewer excursions last long enough for the escape to complete, and
the rate-average over-counts.

**Two points cannot establish a law and fitting them would be §91's mistake exactly.** Omega is
the clean axis: it changes A*Omega without touching the landscape, the transfer map, the margin in
concentration units, or the coupling. So the curve is traced on ONE element and the OTHER element's
point is then checked against it, with nothing fitted.

PREDICTIONS, written before running.

  P1  GATE. At Omega = 30 this must reproduce §92's stored frozen 4.845, fast 1.139 and measured
      4.442. Different numbers mean a different instrument.
  P2  **THE SWEEP.** Position vs A*Omega across Omega. **Predicted: falls monotonically** -- deeper
      barriers sit further from the frozen limit.
  P3  **THE TEST, OUT OF SAMPLE AND NOTHING FITTED.** Does §97's element, at A*Omega = 13.2 and
      position 0.683, land on the curve traced by the FIRST element's Omega sweep? **If it does,
      the penalty closes**: frozen and fast are computable and the position is a function of the
      barrier depth. **If it does not, A*Omega is not the interpolation variable**, and the honest
      output is the bracket.
  P4  **RULE 15, and it binds here.** If P3 fails, do NOT fit a two-parameter interpolation to
      two elements. Report the bracket and stop. §91's slope was exactly such a fit and §97 showed
      what it was worth.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import scipy.sparse.linalg as spla

import experiments.chemical_cascade as cc
from experiments.margin_law import (
    R1, R2, R3, build_reflect, predict, stage1_stationary, upstream_qsd,
)
from experiments.timescale_ratio import pinned_reference

# §92 and §97's stored values (rule 16).
CAL_OM, CAL_FROZEN, CAL_FAST, CAL_MEAS = 30, 4.845, 1.139, 4.442
NEW_AOM, NEW_FROZEN, NEW_FAST, NEW_MEAS = 13.2, 8.925, 1.010, 4.474
A_CAL = 0.190241


def position(frozen, fast, meas):
    """Where the measurement sits between the limits, in LOG space: 0 = fast, 1 = frozen."""
    lf, la, lm = np.log(fast), np.log(frozen), np.log(meas)
    return float((lm - lf) / (la - lf))


def measure(om, t=2.0):
    cc.HILL_N, cc.HILL_K = 4.0, 1.0
    xs, px = upstream_qsd(om)
    frozen = predict("hill", 4.0, 1.0, om, xs, px)
    fast = predict("hill", 4.0, 1.0, om, xs, px, quenched=True)
    _, pi1 = stage1_stationary(om)
    Q, up, m2, cap = build_reflect(om, 1.0)
    p = np.zeros(len(up) * m2)
    for a, w in enumerate(pi1):
        p[a * m2 + int(round(R3 * om))] = w
    p = spla.expm_multiply(Q.T * t, p)
    lo = (np.arange(len(up) * m2) % m2) < R2 * om
    meas = float(p[lo].sum()) / pinned_reference(om, 1.0, t)
    return frozen, fast, meas


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omegas", type=int, nargs="+", default=[14, 20, 30, 40, 55, 70, 85])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/penalty_interpolation.json"))
    args = ap.parse_args()
    out = {}

    print("=== P1 GATE: reproduce §92's stored numbers at Omega = 30")
    f0, a0, m0 = measure(CAL_OM)
    print(f"  frozen {f0:.3f} (§92: {CAL_FROZEN});  fast {a0:.3f} ({CAL_FAST});"
          f"  measured {m0:.3f} ({CAL_MEAS})")
    ok = (abs(f0 / CAL_FROZEN - 1) < 0.02 and abs(a0 / CAL_FAST - 1) < 0.02
          and abs(m0 / CAL_MEAS - 1) < 0.02)
    print(f"  -> P1 {'HOLDS: same instrument' if ok else 'FAILS'}")
    assert ok

    print("\n=== P2: sweep Omega on ONE element. A*Omega moves; the landscape does not.")
    print(f"{'Omega':>7}{'A*Omega':>10}{'fast':>9}{'measured':>11}{'frozen':>9}{'position':>11}")
    rows = []
    for om in args.omegas:
        f, a, m = measure(om)
        pos = position(f, a, m)
        rows.append({"omega": om, "AOm": A_CAL * om, "frozen": f, "fast": a,
                     "meas": m, "pos": pos})
        print(f"{om:>7}{A_CAL*om:>10.2f}{a:>9.3f}{m:>11.3f}{f:>9.3f}{pos:>11.4f}")
    out["sweep"] = rows
    poss = [r["pos"] for r in rows]
    mono = all(poss[i + 1] < poss[i] for i in range(len(poss) - 1))
    print(f"  position runs {poss[0]:.4f} -> {poss[-1]:.4f}"
          f"  ({'monotone falling' if mono else 'NOT monotone'})")
    print(f"  -> P2 {'HOLDS: deeper barriers sit further from the frozen limit' if mono else 'FAILS: the position is not monotone in the barrier depth'}")

    print("\n=== P3: does §97's OTHER element land on this curve? Nothing fitted.")
    new_pos = position(NEW_FROZEN, NEW_FAST, NEW_MEAS)
    aoms = np.array([r["AOm"] for r in rows])
    ps = np.array(poss)
    inside = aoms.min() <= NEW_AOM <= aoms.max()
    # np.interp REQUIRES increasing x. `aoms` already increases with Omega; the first version
    # passed aoms[::-1], which decreases, and np.interp then returned the first element (1.3089,
    # the shallowest cell) and printed FAILS off it. Same family as §86.1(3) and rule 19: the
    # tell was that the "prediction" was exactly a row of the table.
    assert np.all(np.diff(aoms) > 0), "A*Omega must increase for np.interp"
    pred = float(np.interp(NEW_AOM, aoms, ps)) if inside else np.nan
    print(f"  §97's element: A*Omega = {NEW_AOM}, measured position = {new_pos:.4f}")
    print(f"  the sweep spans A*Omega = {aoms.min():.2f}..{aoms.max():.2f}"
          f"   -- {'INSIDE, so this is interpolation' if inside else 'OUTSIDE, so this would be extrapolation and is not run (rule 19)'}")
    out["new"] = {"AOm": NEW_AOM, "pos": new_pos, "pred": pred}
    if inside:
        print(f"  curve at that depth predicts position = {pred:.4f}"
              f"   -> off by {100*(pred-new_pos)/new_pos:+.1f}%")
        ok3 = abs(pred - new_pos) < 0.08
        print(f"  -> P3 {'HOLDS: a SECOND element lands on the first elements curve, so the position is a function of the barrier depth and the penalty CLOSES -- frozen, fast and position are all computable' if ok3 else 'FAILS: the second element does not land on the curve, so A*Omega is not the interpolation variable'}")
        if not ok3:
            print("  -> P4: the bracket is the output. NO two-parameter interpolation is fitted")
            print("     to two elements -- §91's slope was exactly that and §97 showed its worth.")
    else:
        print("  -> P3 UNDECIDED without extending the sweep to cover the second element's depth.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
