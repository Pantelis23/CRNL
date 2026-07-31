"""T14-c-i: does §28's discrepancy collapse against the timescale separation?

§28 predicted the collapse slope from closed forms with no free parameter and got
4.1% at gamma = 0.30 but 24-36% at gamma = 0.15, and showed the drift is NOT
finite-Omega (refitting on the upper Omega half moves it further from 1). The named
suspect is the prediction's 1-D reduction, which slaves the pool to its nullcline:
§24.1a measured that reduction's error shrinking monotonically with the timescale
separation `sep = 3(1+2g)/(1-2g)`, which is 5.6 at gamma=0.15 against 12.0 at 0.30 --
worse reduction exactly where the prediction is worse.

TWO POINTS FIT ANYTHING. On those two, `excess = ratio - 1` times `sep` gives 1.34 and
0.49 (not constant) while `excess * sep^2` gives 7.4 and 5.9 (closer) -- which is
suggestive of `excess ~ sep^-2` and worth exactly nothing as evidence. This measures
five gammas so the suggestion becomes a curve or dies.

PRE-COMMITTED CRITERIA, fixed before any ratio is computed, because §24.3 and §28
both turned on cells that should never have been in the fit:

  * DECADES: a gamma counts only if its exact sweep spans >= 2.0 decades. Below that
    the slope is under-determined -- this is why §28 excluded gamma = 0.45 at 0.40
    decades.
  * PRECISION FLOOR: only cells with P >= 1e-12 are used. `p_cme` returns `1 - split`,
    so a tail at 1e-17 is pure cancellation noise in double precision; at 1e-12 about
    four significant digits survive. Checked empirically -- gamma = 0.15 stays
    monotone to 2.0e-13 at Omega = 400 -- but the floor is set above where it was
    checked, not at it.
  * eps CONTROL: every fit is eps-controlled as in §27, because the integer lattice
    makes realised eps wobble ~10% and that alone bounced §27's raw local slopes 40%.

PREDICTIONS, written before running:

  P1  The excess collapses against `sep` on some simple power -- most likely between
      sep^-1 and sep^-2 given the two-point hint. If it collapses, the 1-D slaved
      reduction is the culprit, §15's closed form is intact, and §28's gamma-drift is
      an artifact of the reduction rather than of the physics.
  P2  THE OUTCOME THAT MATTERS MORE. If the excess does NOT collapse -- if it varies
      with gamma in a way the separation cannot account for -- then the
      gamma-dependence lives in `kappa*delta*^2` itself. **That would be the most
      consequential withdrawal available in this project**, because §15 corrected §12
      and §16-§28 all lean on it. I am not expecting this, which is exactly why the
      criteria above are fixed in advance.
  P3  The ratio should exceed 1 at every gamma (prediction too steep), continuing
      §28's pattern and §22.4's finding that `kappa delta^2` is stiffer than the exact
      barrier. A ratio crossing below 1 somewhere would mean the sign of the reduction
      error flips with gamma, which no current account predicts.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from crnl.information import wall_coefficient
from crnl.networks.am_reversible import am_reversible, delta_star
from experiments.approximation_hierarchy import _setup, p_cme
from experiments.collapse_slope_absolute import V_exact

MIN_DECADES = 2.0
P_FLOOR = 1e-12

OMEGAS = {
    0.15: [40, 80, 120, 160, 200, 240, 280, 320, 360],
    0.20: [40, 80, 120, 160, 200, 260, 320, 380, 440],
    0.25: [40, 100, 160, 220, 280, 340, 400, 460, 520],
    0.30: [40, 120, 200, 280, 360, 440, 520, 600, 680],
    0.35: [40, 140, 240, 340, 440, 540, 640, 740, 840],
}


def sweep(gamma, eps, theta, omegas):
    ds = delta_star(gamma)
    om, lp, er, ps = [], [], [], []
    for o in omegas:
        try:
            p = p_cme(gamma, o, eps, theta)
        except Exception:
            continue
        if not np.isfinite(p) or p < P_FLOOR or p >= 1.0:
            continue
        n0, _t, _ = _setup(gamma, o, eps, theta)
        om.append(float(o)); lp.append(np.log(p)); ps.append(p)
        er.append(int(n0[0] - n0[1]) / (ds * o))
    om, lp, er = np.array(om), np.array(lp), np.array(er)
    if len(om) < 4:
        return None
    dec = float(np.log10(max(ps) / min(ps)))
    A = np.vstack([om, om * (er - er.mean()), np.ones_like(om)]).T
    c, *_ = np.linalg.lstsq(A, lp, rcond=None)
    r = lp - A @ c
    return {"slope": float(c[0]), "r2": float(1 - r.var() / lp.var()),
            "decades": dec, "n_cells": len(om),
            "p_hi": float(max(ps)), "p_lo": float(min(ps))}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eps-frac", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/collapse_slope_collapse.json"))
    args = ap.parse_args()

    print(f"criteria fixed in advance: >= {MIN_DECADES} decades, P >= {P_FLOOR:g}, "
          f"eps-controlled fits\n")
    print(f"{'gamma':>6}{'sep':>7}{'cells':>6}{'decades':>8}{'P range':>22}"
          f"{'measured':>11}{'predicted':>11}{'ratio':>8}{'R^2':>9}")
    rows = []
    for g in sorted(OMEGAS):
        s = sweep(g, args.eps_frac, args.theta, OMEGAS[g])
        if s is None:
            print(f"{g:>6.2f}   no usable cells"); continue
        net = am_reversible(g); ds = delta_star(g); x0 = args.eps_frac * ds
        pred = -2.0 * V_exact(net, x0)
        sep = 3 * (1 + 2 * g) / (1 - 2 * g)
        ratio = pred / s["slope"]
        keep = s["decades"] >= MIN_DECADES
        rows.append({"gamma": g, "sep": sep, "predicted": pred, "ratio": ratio,
                     "kept": keep, "kappa": wall_coefficient(g), **s})
        flag = "" if keep else "   EXCLUDED (< min decades)"
        print(f"{g:>6.2f}{sep:>7.2f}{s['n_cells']:>6}{s['decades']:>8.2f}"
              f"{s['p_hi']:>10.2e}{s['p_lo']:>12.2e}"
              f"{s['slope']:>11.6f}{pred:>11.6f}{ratio:>8.3f}{s['r2']:>9.6f}{flag}")

    use = [r for r in rows if r["kept"]]
    if len(use) >= 3:
        sep = np.array([r["sep"] for r in use])
        exc = np.array([r["ratio"] for r in use]) - 1.0
        print(f"\n=== P1: does the excess collapse against the separation?")
        print(f"{'gamma':>6}{'sep':>8}{'excess':>10}{'exc*sep':>10}"
              f"{'exc*sep^2':>12}")
        for r, s_, e_ in zip(use, sep, exc):
            print(f"{r['gamma']:>6.2f}{s_:>8.2f}{e_:>10.4f}{e_*s_:>10.3f}"
                  f"{e_*s_**2:>12.2f}")
        pos = exc > 0
        if pos.sum() >= 3:
            p = np.polyfit(np.log(sep[pos]), np.log(exc[pos]), 1)
            res = np.log(exc[pos]) - np.polyval(p, np.log(sep[pos]))
            r2 = 1 - res.var() / np.log(exc[pos]).var()
            print(f"\n  excess = {np.exp(p[1]):.3f} * sep^({p[0]:.3f})   R^2 = {r2:.4f}")
            print(f"  (a clean power law supports P1: the 1-D reduction is the culprit"
                  f" and §15's closed form is intact)")
        else:
            print(f"\n  P3 VIOLATED: {int((~pos).sum())} of {len(use)} gammas have "
                  f"ratio < 1, so the sign of the discrepancy flips. No current "
                  f"account predicts that.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
