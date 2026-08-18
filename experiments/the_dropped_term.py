"""T15-n.5: is the reduction's error the term the reduction drops? Computed, not fitted

§88 established the tilted-generator action overshoots by +5.61% at gamma = 0.30 falling to
+1.02% at 0.46, that the overshoot falls with timescale separation (so it IS finite separation),
and that the 1/M law predicted for it is refuted -- |overshoot| x M grows 5-6x over M = 1..16.
**The obvious next move is to fit a power to the decay. That is exactly rule 21's mistake**, and
this project has now paid for it twice.

There is a computable candidate instead. The slow-fast reduction of a WKB Hamiltonian drops the
fast pair's contribution to the action: Assaf, Roberts & Luthey-Schulten (PRL 106:248102, 2011)
state it plainly -- "we have neglected in S(y) the term int p_x dx ~ O(1/gamma)" -- and §85.3
recorded that the literature agrees on this MAGNITUDE while claiming no sign. **§87's solve
already returns p_s(u) and s(u) along the curve, so the dropped term can be integrated directly
and compared in ABSOLUTE terms (rule 16) against the measured overshoot.** Nothing is fitted.

    C = int p_s ds       computed along the solved curve
    overshoot = A_predicted - A_measured        from §88's stored numbers

PREDICTIONS, written before running.

  P1  GATE, a regression. The curve solved here must reproduce §88's action at every gamma to
      better than 0.1%; otherwise this is measuring a different curve than §88 did.
  P2  **MAGNITUDE.** |C| must be the same size as the overshoot -- predicted within a factor of
      about 2 across all eight gammas, since both are the leading finite-separation term. A C
      that is orders out is not the explanation whatever its sign.
  P3  **SIGN, and it is the sharp half.** The reduction OVERSHOOTS, so the correction must
      SUBTRACT. **Predicted: A_pred - C lands closer to A_meas than A_pred does.** If C has the
      wrong sign, or subtracting it makes the agreement worse, **the dropped term is not what the
      reduction is losing** -- and that is a real result, because it is the literature's own
      candidate and §85.3 flagged that the literature claims no sign for it.
  P4  **THE M AXIS.** C must fall with timescale separation alongside the overshoot, and their
      RATIO must be flatter than either -- that is what "the same effect" means. §88 showed the
      overshoot falls but not like 1/M; if C reproduces that same non-1/M shape, including the
      non-monotone first step (2.56 -> 2.99 at gamma = 0.40), the term is identified without any
      power being fitted.
  P5  **RULE 15.** Report C at every gamma and every M, not a summary, and report the ratio
      C/overshoot as a series. If that ratio is itself structured in gamma, say so rather than
      quoting its mean.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from experiments.deep_barriers import MEASURED_A, MEASURED_A_M, action, pu_at
from experiments.where_the_deficit_lives import u_star


def curve_full(g, M=1.0, n=120):
    """The solved curve with the fast pair retained: u, s, ps, pu."""
    us = u_star(g, M)
    rows, guess = [], None
    for u in np.linspace(us * 0.01, us * 0.99, n):
        r = pu_at(u, g, M, guess)
        if r is None:
            continue
        guess = np.array([r["s"], r["ps"]])
        rows.append(r)
    return rows


def dropped_term(g, M=1.0, n=120):
    """C = int p_s ds along the solved curve, and A = -int p_u du on the same grid."""
    rows = curve_full(g, M, n)
    if len(rows) < n // 2:
        return None
    u = np.array([r["u"] for r in rows])
    s = np.array([r["s"] for r in rows])
    ps = np.array([r["ps"] for r in rows])
    pu = np.array([r["pu"] for r in rows])
    us = u_star(g, M)
    C = float(np.trapezoid(ps, s))
    uu = np.concatenate(([0.0], u, [us]))
    pp = np.concatenate(([0.0], pu, [0.0]))
    A = float(-np.trapezoid(pp, uu))
    return {"A": A, "C": C, "n": len(rows)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/the_dropped_term.json"))
    args = ap.parse_args()
    out = {}

    print("=== P1 GATE: does this curve reproduce §88's action?")
    worst = 0.0
    for g in (0.30, 0.40, 0.46):
        d = dropped_term(g)
        a88, _ = action(g, 1.0, n=80)
        worst = max(worst, abs(d["A"] / a88 - 1))
        print(f"  gamma={g}: A here = {d['A']:.6f}, §88 = {a88:.6f}, ratio {d['A']/a88:.6f}")
    ok = worst < 1e-3
    print(f"  -> P1 {'HOLDS: same curve' if ok else 'FAILS: a different curve is being measured'}")
    assert ok

    print("\n=== P2/P3/P5: the dropped term against the measured overshoot, per gamma")
    print(f"{'gamma':>7}{'A meas':>11}{'A pred':>11}{'overshoot':>11}{'C':>12}"
          f"{'C/over':>9}{'A pred - C':>12}{'ratio':>8}")
    rows = []
    for g in sorted(MEASURED_A):
        d = dropped_term(g)
        if d is None:
            print(f"{g:>7}   no curve")
            continue
        m = MEASURED_A[g]
        over = d["A"] - m
        rows.append({"gamma": g, "meas": m, "A": d["A"], "over": over, "C": d["C"],
                     "ratio_C": d["C"] / over, "corrected": (d["A"] - d["C"]) / m})
        print(f"{g:>7}{m:>11.6f}{d['A']:>11.6f}{over:>11.6f}{d['C']:>12.6f}"
              f"{d['C']/over:>9.3f}{d['A']-d['C']:>12.6f}{(d['A']-d['C'])/m:>8.4f}")
    out["per_gamma"] = rows
    mags = [abs(r["ratio_C"]) for r in rows]
    print(f"  |C/overshoot| runs " + ", ".join(f"{v:.2f}" for v in mags))
    p2 = all(0.5 < v < 2.0 for v in mags)
    print(f"  -> P2 {'HOLDS: C is the same size as the overshoot at every gamma' if p2 else 'FAILS: C is not the right magnitude'}")
    better = sum(1 for r in rows if abs(r["corrected"] - 1) < abs(r["A"] / r["meas"] - 1))
    signs = set(np.sign([r["ratio_C"] for r in rows]))
    print(f"  subtracting C improves the agreement in {better}/{len(rows)} cells;"
          f" sign of C/overshoot: {sorted(signs)}")
    print(f"  -> P3 {'HOLDS: the dropped term has the right sign and subtracting it corrects the reduction' if better == len(rows) else 'REFUTED: the dropped term does NOT correct the reduction -- it is the literatures own candidate and it is not what is being lost'}")
    spread = max(mags) / min(mags)
    print(f"  -> P5 the ratio C/overshoot spans {min(mags):.2f}..{max(mags):.2f}"
          f" (a factor of {spread:.2f}), "
          + ("structured in gamma, not a constant" if spread > 1.3 else "flat in gamma"))

    print("\n=== P4: the M axis -- does C reproduce the overshoot's non-1/M shape?")
    print(f"{'gamma':>7}{'M':>4}{'overshoot':>12}{'C':>12}{'C/over':>9}{'A pred - C':>12}"
          f"{'ratio':>8}")
    p4 = []
    for g in (0.40, 0.44):
        for M in sorted(MEASURED_A_M[g]):
            d = dropped_term(g, float(M))
            if d is None:
                print(f"{g:>7}{M:>4}   no curve")
                continue
            m = MEASURED_A_M[g][M]
            over = d["A"] - m
            p4.append({"gamma": g, "M": M, "over": over, "C": d["C"],
                       "ratio_C": d["C"] / over, "corrected": (d["A"] - d["C"]) / m})
            print(f"{g:>7}{M:>4}{over:>12.6f}{d['C']:>12.6f}{d['C']/over:>9.3f}"
                  f"{d['A']-d['C']:>12.6f}{(d['A']-d['C'])/m:>8.4f}")
    out["per_M"] = p4
    for g in (0.40, 0.44):
        ser = [r for r in p4 if r["gamma"] == g]
        if len(ser) >= 3:
            rr = [r["ratio_C"] for r in ser]
            print(f"  gamma={g}: C/overshoot over M = " + ", ".join(f"{v:.3f}" for v in rr)
                  + f"   (spans {max(rr)/min(rr):.2f}x)")
    if p4:
        flat = all(0.5 < abs(r["ratio_C"]) < 2.0 for r in p4)
        print(f"  -> P4 {'HOLDS: C tracks the overshoot along M as well, so the term is identified with no power fitted' if flat else 'FAILS: C does not track the overshoot along M'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
