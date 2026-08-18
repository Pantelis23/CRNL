"""T15-n.6: is §89's flat 0.1% residual the instrument? Re-extract A without a two-point slope

§89 corrected the tilted-generator action with the term the reduction drops and landed at
**1.0007 +- 0.0004 in every cell, gamma and M alike** -- the correction removed a 7x spread in the
raw overshoot and left a CONSTANT. §89.1 named that constant a suspect, not a finding, on the
grounds that a constant offset surviving a correction which fixed everything else is more likely
to sit on the MEASUREMENT side.

**The measurement side is specific and checkable.** §84's A is the last two-point local slope
d ln T/dOmega, admitted by a 2% convergence gate. The WKB form is

    ln T = A*Omega + b*ln Omega + c + O(1/Omega)

so a two-point slope over [Om1, Om2] returns A + b*ln(Om2/Om1)/(Om2-Om1) = A + O(1/Omega), a
BIAS OF FIXED SIGN that a convergence gate does not remove -- the gate only checks that successive
slopes stop moving quickly, not that they have stopped moving.

**T15-n's own history warns against the alternative.** Fitting A, b and c together is the
collinearity §35.3 proved unresolvable over a bounded window, and §64's action route -- the one
that gave the 2.19 outlier -- was the heaviest extrapolation of the three. So this does NOT fit
the three-parameter form. It takes the local slopes on a ladder of Omega and extrapolates them
against 1/Omega, which is one parameter and is exactly Richardson.

PREDICTIONS, written before running.

  P1  GATE. ln T must be positive, increasing in Omega, and below 28 at every point (§84's
      ceiling: the double-precision limit bends the slope upward near ln T ~ 35 and turned one
      cell NEGATIVE at Omega = 650 in scouting).
  P2  **THE BIAS MUST BE VISIBLE BEFORE IT IS REMOVED.** The local slopes must drift
      systematically with Omega and be roughly linear in 1/Omega. If they scatter instead, there
      is no 1/Omega bias to extrapolate and the whole premise is wrong.
  P3  **THE TEST, and its direction is fixed in advance.** §89's corrected ratio is 1.0007, so
      closing it requires the measured A to be **LARGER** by about 0.07%. **Predicted: the
      Richardson-extrapolated A exceeds §84's two-point value by roughly that, at all three
      gammas.** A correction of the wrong sign, or one an order too small, refutes the instrument
      account.
  P4  **THE CONSEQUENCE.** Recompute §89's corrected ratio against the re-extracted A.
      **Predicted: it moves toward 1.000.** Per rule 20 the criterion is that it MOVES TOWARD,
      not that it lands inside any tolerance -- the extrapolation is itself approximate.
  P5  **IF IT DOES NOT CLOSE.** Then there is a second physical term at the 0.1% level and the
      right response is to say so and name it, **not to fit it** -- §84's P3 became a claim about
      a phenomenon that was its own Taylor term precisely that way.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from experiments.deep_barriers import MEASURED_A
from experiments.ode_does_not_determine_it import am_ln_mfpt
from experiments.the_dropped_term import dropped_term

# §84's two-point values, quoted as stored numbers (rule 16).
LADDER = {0.35: [100, 150, 200, 250, 300, 350],
          0.40: [200, 300, 400, 500, 650, 800],
          0.44: [500, 700, 900, 1100, 1400, 1700]}


def ln_T_ladder(g, omegas):
    out = []
    for om in omegas:
        v = am_ln_mfpt(g, 0.0, om)
        if v is not None and 0.0 < v < 28.0:
            out.append((om, v))
    return out


def slopes(pts):
    """Successive local slopes and the harmonic midpoint they belong to."""
    return [((p2[1] - p1[1]) / (p2[0] - p1[0]), 2.0 / (1.0 / p1[0] + 1.0 / p2[0]))
            for p1, p2 in zip(pts, pts[1:])]


def richardson(pts):
    """Extrapolate the local slopes to 1/Omega -> 0. One parameter, not three (§35.3)."""
    sl = slopes(pts)
    if len(sl) < 3:
        return None, sl
    x = np.array([1.0 / m for _, m in sl])
    y = np.array([s for s, _ in sl])
    a, b = np.polyfit(x, y, 1)
    return float(b), sl          # intercept at 1/Omega = 0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/richardson_A.json"))
    args = ap.parse_args()
    out, rows = {}, []

    print("=== P1/P2: the ladder, and is the bias visible before it is removed?")
    for g, oms in LADDER.items():
        pts = ln_T_ladder(g, oms)
        ok = all(0 < v < 28 for _, v in pts) and all(
            pts[i + 1][1] > pts[i][1] for i in range(len(pts) - 1))
        sl = slopes(pts)
        print(f"\n  gamma={g}  ln T at Omega " + ", ".join(f"{o}:{v:.3f}" for o, v in pts))
        print(f"    local slopes: " + ", ".join(f"{s:.6f}" for s, _ in sl))
        x = np.array([1.0 / m for _, m in sl])
        y = np.array([s for s, _ in sl])
        r = np.corrcoef(x, y)[0, 1] if len(x) > 2 else np.nan
        A_ext, _ = richardson(pts)
        rows.append({"gamma": g, "pts": pts, "slopes": [s for s, _ in sl],
                     "A_rich": A_ext, "corr": float(r), "clean": bool(ok)})
        print(f"    P1 {'clean' if ok else 'FAILS: ln T not positive/increasing/below 28'};"
              f"  slope vs 1/Omega correlation {r:+.4f}")
    out["ladders"] = rows
    p1 = all(r["clean"] for r in rows)
    p2 = all(abs(r["corr"]) > 0.9 for r in rows)
    print(f"\n  -> P1 {'HOLDS' if p1 else 'FAILS'};  P2 "
          f"{'HOLDS: the slopes are linear in 1/Omega, so there IS a bias to remove' if p2 else 'FAILS: the slopes scatter -- no 1/Omega bias, and the premise is wrong'}")

    print("\n=== P3/P4: the re-extracted A, and what it does to §89's residual")
    print(f"{'gamma':>7}{'§84 two-point':>15}{'Richardson':>13}{'shift %':>10}"
          f"{'§89 ratio':>12}{'new ratio':>11}")
    p34 = []
    for r in rows:
        g = r["gamma"]
        if r["A_rich"] is None:
            continue
        old = MEASURED_A[g]
        d = dropped_term(g, n=80)
        corrected = d["A"] - d["C"]
        p34.append({"gamma": g, "A_old": old, "A_new": r["A_rich"],
                    "shift": 100 * (r["A_rich"] / old - 1),
                    "ratio_old": corrected / old, "ratio_new": corrected / r["A_rich"]})
        print(f"{g:>7}{old:>15.6f}{r['A_rich']:>13.6f}"
              f"{100*(r['A_rich']/old-1):>10.3f}{corrected/old:>12.4f}"
              f"{corrected/r['A_rich']:>11.4f}")
    out["p34"] = p34

    print("\n=== RULE 15: every candidate extrapolation, not the one that closes the residual")
    print(f"{'gamma':>7}{'last slope':>12}{'lin all':>10}{'lin last3':>11}{'quad':>9}"
          f"{'lin last4':>11}{'spread %':>10}")
    p15 = []
    for r in rows:
        g = r["gamma"]
        pts = r["pts"]
        sl = slopes(pts)
        x = np.array([1.0 / m for _, m in sl])
        y = np.array([s for s, _ in sl])
        cands = {"last": float(y[-1]),
                 "lin all": float(np.polyfit(x, y, 1)[-1]),
                 "lin last3": float(np.polyfit(x[-3:], y[-3:], 1)[-1]),
                 "quad": float(np.polyfit(x, y, 2)[-1]) if len(x) >= 4 else np.nan,
                 "lin last4": float(np.polyfit(x[-4:], y[-4:], 1)[-1]) if len(x) >= 4 else np.nan}
        v = [c for c in cands.values() if np.isfinite(c)]
        d = dropped_term(g, n=80)
        corr = d["A"] - d["C"]
        ratios = {k: corr / c for k, c in cands.items() if np.isfinite(c)}
        p15.append({"gamma": g, "cands": cands, "ratios": ratios,
                    "spread": 100 * (max(v) / min(v) - 1)})
        print(f"{g:>7}{cands['last']:>12.6f}{cands['lin all']:>10.6f}"
              f"{cands['lin last3']:>11.6f}{cands['quad']:>9.6f}{cands['lin last4']:>11.6f}"
              f"{100*(max(v)/min(v)-1):>10.3f}")
        print("        §89 corrected/A for each: "
              + ", ".join(f"{k}={vv:.4f}" for k, vv in ratios.items()))
    out["p15"] = p15
    straddle = all(min(r["ratios"].values()) < 1.0 < max(r["ratios"].values()) for r in p15)
    worst = max(r["spread"] for r in p15)
    print(f"  ansatz spread in the measured A: up to {worst:.3f}%, against §89's residual of")
    print(f"  +0.07..+0.18%. **The extrapolations STRADDLE 1 at every gamma.**"
          if straddle else "  the extrapolations agree on a side")
    print(f"  -> T15-n.6 {'RESOLVES AS UNRESOLVABLE HERE: the measured A is not determined to better than its own ansatz spread, which EXCEEDS the residual it was being used to judge. §89s residual is below the resolution of the instrument that would test it -- attributable to neither side' if straddle else 'the residual survives the ansatz sweep'}")
    print("  And §89.1's flatness argument is weaker than it looked: every cell used the SAME")
    print("  two-point estimator, so its bias is COMMON MODE and a flat residual is what that")
    print("  produces. Flatness across cells sharing an estimator is not evidence about physics.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
