"""T14-c-iv: does the excess really cross zero near gamma ~ 0.357? -- reopened by §35.

§28.3 fitted `excess = 0.2240 - 0.6276*gamma` at R^2 = 0.9905, crossing zero near
gamma ~ 0.357, and offered it explicitly as a DESCRIPTION rather than a law. It could
not test the crossing, because at large gamma the collapse is so shallow that four
decades of P needs Omega ~ 2000, and `OMEGA_CAP` was 900. **That cap existed because of
the probability floor, not because of Omega**: `p_cme` computes the error as
`1 - split`, which dies to cancellation near 1e-12, so a shallow-gamma cell could not be
pushed deep enough to fit a slope at all. §35 removes that floor by solving for the
wrong outcome directly, validated componentwise to 1e-13 at P ~ 1e-33. The cells §28.3
reported as skipped are now reachable.

**§35 ALSO CHANGES WHAT "THE SLOPE" MEANS, and this experiment has to respect that.**
The local slope drifts with Omega, so a fitted slope is a property of its window. §28.3's
numbers were measured over P = 1e-2 -> 1e-6; any comparison against them must use the
SAME window or it is comparing different quantities (rule 11's spirit -- a control must
share a clock with its arm). So every gamma here is measured twice:

  * `matched`  -- the identical P = 1e-2 -> 1e-6 window §28.3 used, so the new gammas
                  extend that table on its own terms;
  * `deep`     -- P = 1e-2 -> 1e-20, only possible after §35, to show how much of the
                  "excess" is window and how much is gamma.

PREDICTIONS, written before running:

  P1  All three of §28.3's unreachable gammas (0.38, 0.41, 0.44) now yield a fitted
      slope on the matched window. If any still fails it is an Omega/state-space limit,
      not a probability limit, and is reported as such.
  P2  THE TEST, on the matched window so it is commensurable with §28.3. The line
      predicts excess = -0.0145, -0.0333, -0.0521 at gamma = 0.38, 0.41, 0.44. Either
      (a) the excess really does go negative and roughly follows the line -- the
      description survives into new territory; or (b) it flattens near zero -- the line
      was local and the truth is saturation, meaning §28.3's crossing point was an
      artifact of extrapolating a straight line off the end of its data. **I expect (b)**,
      because a negative excess means §15's closed form is too SHALLOW there, and §22.4
      and §28 both found it consistently too steep. A sign change in a discrepancy that
      has one sign everywhere else needs more than a two-parameter fit to be believed.
  P3  The `deep` window gives a SMALLER |slope| than the `matched` window at every
      gamma, following §35's drift, so the deep-window excess is systematically larger.
      If the two windows give the same excess the drift does not affect this observable
      and §28.3's numbers need no window caveat.
  P4  eps-controlled fits throughout, and the achieved P window reported per cell so
      matching can be checked rather than assumed -- §28.1's entire error was an
      unmatched grid nobody had checked.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.networks.am_reversible import GAMMA_C, am_reversible, delta_star
from experiments.collapse_slope_absolute import V_exact
from experiments.deep_tail import solve_direct

WINDOWS = {"matched": (1e-2, 1e-6), "deep": (1e-2, 1e-20)}


def p_at(gamma, omega, eps, theta):
    try:
        r = solve_direct(gamma, int(omega), eps, theta)
    except Exception:
        return None
    p = r["p"]
    return r if (np.isfinite(p) and 0.0 < p < 1.0) else None


def find_omega(gamma, eps, theta, target, cap):
    """Smallest Omega whose P drops to `target`. None if beyond the cap."""
    lo, hi = 20, cap
    r = p_at(gamma, hi, eps, theta)
    if r is None or r["p"] > target:
        return None
    while hi - lo > 8:
        mid = (lo + hi) // 2
        rm = p_at(gamma, mid, eps, theta)
        if rm is None or rm["p"] > target:
            lo = mid
        else:
            hi = mid
    return int(hi)


def fit_cell(gamma, eps, theta, window, cap, n_cells):
    p_hi, p_lo = WINDOWS[window]
    o_hi = find_omega(gamma, eps, theta, p_hi, cap)
    o_lo = find_omega(gamma, eps, theta, p_lo, cap)
    if o_hi is None or o_lo is None or o_lo <= o_hi:
        return None
    grid = sorted({int(round(o)) for o in np.linspace(o_hi, o_lo, n_cells)})
    om, lp, er = [], [], []
    for o in grid:
        r = p_at(gamma, o, eps, theta)
        if r is None:
            continue
        om.append(float(o)); lp.append(np.log(r["p"]))
        er.append(r["eps_realised"])
    om, lp, er = np.array(om), np.array(lp), np.array(er)
    if len(om) < 6:
        return None
    A = np.vstack([om, om * (er - er.mean()), np.ones_like(om)]).T
    c, *_ = np.linalg.lstsq(A, lp, rcond=None)
    res = lp - A @ c
    return {"slope": float(c[0]), "r2": float(1 - res.var() / lp.var()),
            "n": len(om), "o_lo": int(om.min()), "o_hi": int(om.max()),
            "decades": float((lp.max() - lp.min()) / np.log(10)),
            "p_hi": float(np.exp(lp.max())), "p_lo": float(np.exp(lp.min()))}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.25, 0.30, 0.35, 0.38, 0.41, 0.44])
    ap.add_argument("--eps-frac", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--cap", type=int, default=2600)
    ap.add_argument("--cells", type=int, default=10)
    ap.add_argument("--windows", type=str, nargs="+", default=["matched", "deep"])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/collapse_slope_deep.json"))
    args = ap.parse_args()

    t0 = time.time()
    print(f"direct solve (no subtraction), Omega cap {args.cap}, "
          f"eps = {args.eps_frac}, {args.cells} cells per fit")
    print(f"§28.3's line: excess = 0.2240 - 0.6276*gamma  (crossing at gamma ~ 0.357)")
    out = {}
    for window in args.windows:
        lo, hi = WINDOWS[window]
        print(f"\n=== window: P = {lo:g} -> {hi:g}")
        print(f"{'gamma':>7}{'Omega':>14}{'dec':>6}{'measured':>12}{'predicted':>12}"
              f"{'ratio':>8}{'excess':>9}{'line says':>11}{'R^2':>9}")
        rows = []
        for g in args.gammas:
            if g >= GAMMA_C:
                print(f"{g:>7.2f}   SKIPPED (gamma >= gamma_c)"); continue
            f = fit_cell(g, args.eps_frac, args.theta, window, args.cap, args.cells)
            if f is None:
                print(f"{g:>7.2f}   SKIPPED (needs Omega > {args.cap})"); continue
            net = am_reversible(g)
            pred = -2.0 * V_exact(net, args.eps_frac * delta_star(g))
            ratio = pred / f["slope"]
            line = 0.2240 - 0.6276 * g
            rows.append({"gamma": g, "window": window, "predicted": pred,
                         "ratio": ratio, "excess": ratio - 1.0, "line": line, **f})
            print(f"{g:>7.2f}{f['o_lo']:>8}-{f['o_hi']:<5}{f['decades']:>6.1f}"
                  f"{f['slope']:>12.6f}{pred:>12.6f}{ratio:>8.4f}{ratio-1:>9.4f}"
                  f"{line:>11.4f}{f['r2']:>9.6f}")
        out[window] = rows

    print(f"\n=== P2: does the excess cross zero, or flatten?")
    m = out.get("matched", [])
    new = [r for r in m if r["gamma"] > 0.36]
    if new:
        print(f"{'gamma':>7}{'measured excess':>18}{'line predicts':>15}"
              f"{'meas - line':>13}")
        for r in new:
            print(f"{r['gamma']:>7.2f}{r['excess']:>18.4f}{r['line']:>15.4f}"
                  f"{r['excess']-r['line']:>13.4f}")
        neg = [r for r in new if r["excess"] < -0.005]
        print(f"  cells with a genuinely NEGATIVE excess: {len(neg)}/{len(new)}")
        print(f"  -> {'the crossing is real, the line survives past its data'
                     if len(neg) == len(new) else
                     'the excess FLATTENS rather than crossing -- §28.3 extrapolated a '
                     'straight line off the end of its data'}")
    else:
        print("  no cells past gamma = 0.36 were reachable")

    if len(args.windows) > 1 and out.get("deep"):
        print(f"\n=== P3: how much of the excess is the WINDOW rather than gamma?")
        print(f"{'gamma':>7}{'matched slope':>15}{'deep slope':>13}"
              f"{'matched excess':>16}{'deep excess':>13}")
        dd = {r["gamma"]: r for r in out["deep"]}
        for r in m:
            d = dd.get(r["gamma"])
            if d:
                print(f"{r['gamma']:>7.2f}{r['slope']:>15.6f}{d['slope']:>13.6f}"
                      f"{r['excess']:>16.4f}{d['excess']:>13.4f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
