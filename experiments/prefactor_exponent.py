"""T14-d: is the prefactor exponent -1/2, and is its gamma-dependence real?

§35 fitted `ln P = c*Omega + b*ln(Omega) + a` to the exact collapse and found b =
-0.4484 / -0.4394 / -0.4089 / -0.3964 at gamma = 0.20 / 0.25 / 0.30 / 0.35 -- near
WKB's -1/2, drifting, and correlated with BOTH gamma and the number of decades the fit
spanned (29.21 / 21.20 / 14.30 / 8.79).

**THE CONFOUND IS TIGHTER THAN "two things move together".** §35 used the SAME Omega
grid (150..1500) at every gamma, so the decade count is a deterministic function of
gamma: decades = |slope| * dOmega / ln 10. Varying gamma at a fixed Omega window cannot
separate them even in principle -- they are the same variable wearing two names. The
only way out is to give each gamma a DIFFERENT Omega window chosen so the decade counts
match, which is §30.2's lesson applied deliberately: break a confound with a second
sweep whose structure is opposite, not with a better fit to the first.

  SWEEP A -- decades held FIXED, gamma varies. Each gamma gets an Omega window sized to
             span the same number of decades. If b still tracks gamma here, it is
             physics. If the spread collapses, b was never a function of gamma.
  SWEEP B -- gamma held FIXED, decades vary. Sub-windows of one gamma's data at 9, 14,
             21, 29 decades. If b tracks decades here, the three-term ansatz is
             incomplete and the fit is absorbing higher-order curvature.

**Note what the exact arithmetic rules out.** The residuals are 1e-3 in ln P on data
good to 1e-13 relative, so this is not sampling noise and not an ill-conditioned solve.
If b moves with the window while the data are exact, the MODEL is missing a term --
that is the only remaining explanation, and it is testable directly (P4).

PREDICTIONS, written before running:

  P1  THE TEST. At matched decade count the exponents across gamma agree to within a
      few percent, and the residual gamma-dependence largely vanishes. **I expect this
      outcome**, because the alternative requires a physical mechanism that makes an
      algebraic prefactor depend on the drive, and none of §15's closed forms contains
      one.
  P2  In sweep B, b moves systematically with window length at FIXED gamma -- toward
      -1/2 as the window lengthens. That is the same effect seen from the other side,
      and P1 and P2 must agree with each other or neither reading is safe.
  P3  If P1 and P2 both hold, extrapolating b against 1/decades to an infinite lever arm
      gives a gamma-independent limit. Reported against -1/2 as an absolute
      comparison -- WKB predicts exactly -1/2 for this class and nothing here was fitted
      to make it so.
  P4  RULE 15, and the direct test of "the model is missing a term". A four-term fit
      adding `d/Omega` should give a b that is MORE window-independent than the
      three-term fit. If it does, the missing-term reading is confirmed outright; if b
      still drifts under the richer model, something else is going on and P1's
      interpretation is not safe either. All fits reported, never only the flattering
      one.
  P5  If instead b tracks gamma at matched decades AND is window-independent at fixed
      gamma, the gamma-dependence is real, -1/2 is wrong for this network, and that is a
      genuine finding needing a mechanism -- which will NOT be proposed in the same
      breath as the measurement (rule 17).
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from experiments.deep_tail import solve_direct

TARGETS = (9.0, 14.0, 21.0, 29.0)


def collapse_curve(gamma, eps, theta, omegas):
    out = []
    for om in omegas:
        r = solve_direct(gamma, int(om), eps, theta)
        if np.isfinite(r["p"]) and 0.0 < r["p"] < 1.0:
            out.append(r)
    return out


def fit_window(rows, terms="log"):
    om = np.array([r["omega"] for r in rows], float)
    lp = np.log([r["p"] for r in rows])
    er = np.array([r["eps_realised"] for r in rows])
    ec = er - er.mean()
    cols = [om, np.log(om), om * ec, np.ones_like(om)]
    if terms == "log+inv":
        cols.insert(2, 1.0 / om)
    A = np.vstack(cols).T
    c, *_ = np.linalg.lstsq(A, lp, rcond=None)
    res = lp - A @ c
    return {"c": float(c[0]), "b": float(c[1]),
            "rms": float(np.sqrt((res ** 2).mean())),
            "decades": float((lp.max() - lp.min()) / np.log(10)),
            "o_lo": float(om.min()), "o_hi": float(om.max()), "n": len(om)}


def sub_window(rows, target_decades):
    """Longest prefix of `rows` (from the shallow end) spanning ~target decades."""
    lp = np.log([r["p"] for r in rows])
    span = (lp[0] - lp) / np.log(10)
    idx = np.where(span <= target_decades)[0]
    return rows[: max(idx[-1] + 1, 6)] if len(idx) else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.20, 0.25, 0.30, 0.35])
    ap.add_argument("--eps-frac", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--omega-lo", type=int, default=150)
    ap.add_argument("--omega-cap", type=int, default=2000)
    ap.add_argument("--points", type=int, default=14)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/prefactor_exponent.json"))
    args = ap.parse_args()

    t0 = time.time()
    print("T14-d: matched-decade sweeps to separate gamma from lever arm")
    curves = {}
    for g in args.gammas:
        oms = [int(round(o)) for o in
               np.linspace(args.omega_lo, args.omega_cap, args.points)]
        rows = collapse_curve(g, args.eps_frac, args.theta, oms)
        curves[g] = rows
        dec = (np.log(rows[0]["p"]) - np.log(rows[-1]["p"])) / np.log(10)
        print(f"  gamma={g}: {len(rows)} cells, Omega {oms[0]}..{oms[-1]}, "
              f"{dec:.2f} decades, P down to {rows[-1]['p']:.3e}"
              f"   ({time.time()-t0:.0f}s)")

    print(f"\n=== SWEEP A: decades held FIXED, gamma varies")
    A_rows = []
    for target in TARGETS:
        line = []
        for g in args.gammas:
            sub = sub_window(curves[g], target)
            if sub is None or len(sub) < 6:
                continue
            f = fit_window(sub)
            if abs(f["decades"] - target) > 0.35 * target:
                continue
            line.append((g, f))
        if len(line) < 2:
            continue
        print(f"\n  target {target:.0f} decades:")
        print(f"  {'gamma':>7}{'achieved dec':>14}{'Omega window':>16}{'b':>10}"
              f"{'c':>12}{'rms':>9}")
        for g, f in line:
            win = f"{int(f['o_lo'])}-{int(f['o_hi'])}"
            print(f"  {g:>7.2f}{f['decades']:>14.2f}{win:>16}"
                  f"{f['b']:>10.4f}{f['c']:>12.6f}{f['rms']:>9.5f}")
            A_rows.append({"target": target, "gamma": g, **f})
        bs = np.array([f["b"] for _, f in line])
        print(f"    spread in b across gamma at matched decades: "
              f"{100*(bs.max()-bs.min())/abs(bs.mean()):.2f}%"
              f"   (§35's unmatched spread was 11.6%)")

    print(f"\n=== SWEEP B: gamma held FIXED, decades vary")
    B_rows = []
    for g in args.gammas:
        line = []
        for target in TARGETS:
            sub = sub_window(curves[g], target)
            if sub is None or len(sub) < 6:
                continue
            f = fit_window(sub)
            if abs(f["decades"] - target) > 0.35 * target:
                continue
            line.append(f)
        if len(line) < 3:
            continue
        print(f"\n  gamma = {g}:")
        print(f"  {'decades':>9}{'Omega window':>16}{'b (3-term)':>12}"
              f"{'b (4-term)':>12}{'rms3':>9}{'rms4':>9}")
        for f in line:
            sub = sub_window(curves[g], f["decades"] + 0.01)
            f4 = fit_window(sub, terms="log+inv")
            win = f"{int(f['o_lo'])}-{int(f['o_hi'])}"
            print(f"  {f['decades']:>9.2f}{win:>16}"
                  f"{f['b']:>12.4f}{f4['b']:>12.4f}{f['rms']:>9.5f}{f4['rms']:>9.5f}")
            B_rows.append({"gamma": g, **f, "b4": f4["b"], "rms4": f4["rms"]})
        bs = np.array([f["b"] for f in line])
        print(f"    b moves {100*(bs.max()-bs.min())/abs(bs.mean()):.2f}% "
              f"across window length at FIXED gamma")

    print(f"\n=== P3: extrapolate b to an infinite lever arm")
    print(f"  {'gamma':>7}{'b at 1/dec -> 0':>18}{'R^2':>9}{'points':>8}")
    lims = []
    for g in args.gammas:
        rs = [r for r in B_rows if r["gamma"] == g]
        if len(rs) < 3:
            continue
        x = np.array([1.0 / r["decades"] for r in rs])
        y = np.array([r["b"] for r in rs])
        p = np.polyfit(x, y, 1)
        r2 = 1 - np.var(y - np.polyval(p, x)) / np.var(y)
        lims.append(p[1])
        print(f"  {g:>7.2f}{p[1]:>18.4f}{r2:>9.5f}{len(rs):>8}")
    if lims:
        lims = np.array(lims)
        print(f"  limits across gamma: {lims.min():.4f}..{lims.max():.4f}, "
              f"mean {lims.mean():.4f}   against WKB's -0.5 -> "
              f"{abs(lims.mean()/-0.5):.4f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"sweepA": A_rows, "sweepB": B_rows},
                                   indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
