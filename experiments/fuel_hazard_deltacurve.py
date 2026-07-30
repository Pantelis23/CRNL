"""T10b-iii-d, part 2: does a theta-dependent imposed separation close the residue?

§23.6 seeded the past-`gamma_c` hazard cells at a CONSTANT imposed separation
`delta_past` and closed 39.6% of the integral's exponent error, leaving +0.1445.
§23.7 then measured what the cascade actually carries there and found the constant
is wrong in both directions: the real separation DECAYS along the past-`gamma_c`
sequence -- 0.30 -> 0.13 pooled -- so a constant 0.20 overstates the hazard early
and understates it late.

§23.7 also found the decay is a function of `theta`, not of survival history:
`delta ~ theta` alone gives weighted R^2 = 0.9177 across three budgets and eight
survival indices, and adding `ln Phi` lifts it only to 0.9675 with a coefficient of
-0.028 (about 23% of the delta range). So the repair keeps the model MARKOVIAN in
`theta`, which is the opposite of what §23.6 expected when it wrote that any fix
would make the object non-Markovian.

This imposes `delta_past(theta)` interpolated from §23.7's measured medians instead
of a constant, re-measures `q(theta)` cell by cell, and re-integrates. Still a
one-parameter-per-cell empirical input, so still NOT the absolute test §23.4 ran --
it is now a measured-curve model and is reported as one.

PREDICTIONS, written before running. §23.6's residue is not one error but two, in
two different regions, and this run can only touch one of them:

  P1  Small tanks gain depth. At Phi/Omega = 25 the integral predicts 5.47 against
      a measured 7, and most past-`gamma_c` stages are at LOW k where the real
      separation (0.30) exceeds the imposed 0.20 -- 2465 trials at k=1 against 1048
      at k=3 -- so correcting it lowers the early hazard. Predicted depth rises to
      roughly **6-6.5**.
  P2  Large tanks barely move. At Phi/Omega = 400 only 9.5% of trials reach
      `gamma_eff >= gamma_c` at all and the median count is 0 stages (§23.4), so its
      predicted 48.23 is set almost entirely in the PRE-`gamma_c` region and
      `delta_past` cannot touch it.
  P3  Hence the exponent falls from 0.7919 toward **~0.75**, closing perhaps another
      quarter of the residue and leaving ~0.10.
  P4  AND THE RESIDUE THAT SURVIVES IS THE LARGE-TANK OVER-PREDICTION, which is a
      pre-`gamma_c` error: 48.23 predicted against 44 measured, +10%, in a region
      where the rail exists, the state is quasi-static, and §23.5 already showed
      re-seeding onto that rail changes nothing. If the exponent after this fix is
      still above ~0.72, the remaining discrepancy is NOT about the dying landscape
      at all and the whole past-`gamma_c` line of attack has been chasing the
      smaller of the two errors.
  P5  WHAT WOULD FALSIFY THE MARKOVIAN CLAIM. If imposing `delta_past(theta)` leaves
      the exponent at 0.79, then `theta` does not determine the hazard even with the
      right separation, the R^2 = 0.9177 collapse is not doing the work it appears
      to, and §23.7's Markovian reading must be withdrawn.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.networks.am_reversible import delta_star
from experiments.cascade_fuel_vs_noise import stage_time
from experiments.fuel_hazard import GAMMA0, GAMMA_INF, predict_depth
from experiments.fuel_hazard_pastgc import (
    FUEL_CONCS, depth_continuous, fit, hazard_at,
)


def delta_curve(path: pathlib.Path, min_n: int = 50):
    """delta_past(theta) from §23.7's measured medians, pooled over budgets."""
    d = json.load(open(path))
    pts = [(b["theta_median"], b["delta_median"])
           for r in d for b in r["by_k"] if b["entered"] >= min_n]
    pts.sort()
    th = np.array([p[0] for p in pts])
    dl = np.array([p[1] for p in pts])
    # collapse duplicate thetas so np.interp gets a monotone x
    uniq, inv = np.unique(np.round(th, 4), return_inverse=True)
    med = np.array([np.median(dl[inv == i]) for i in range(len(uniq))])
    lo, hi = float(med[0]), float(med[-1])
    return uniq, med, (lambda t: float(np.interp(t, uniq, med, left=lo, right=hi)))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=40)
    ap.add_argument("--sigma-frac", type=float, default=0.03)
    ap.add_argument("--n-theta", type=int, default=14)
    ap.add_argument("--trials", type=int, default=800)
    ap.add_argument("--seed", type=int, default=20260730)
    ap.add_argument("--bias", type=pathlib.Path,
                    default=pathlib.Path("results/fuel_survivor_bias.json"))
    ap.add_argument("--measured", type=pathlib.Path,
                    default=pathlib.Path("results/fuel_depth_scaling.json"))
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/fuel_hazard_deltacurve.json"))
    args = ap.parse_args()

    om, sf = args.omega, args.sigma_frac
    sigma_counts = sf * delta_star(GAMMA0) * om
    theta0 = GAMMA0 / (1.0 + GAMMA0)
    thetas = np.linspace(theta0, 0.48, args.n_theta)
    phis = np.array([int(round(fc * om)) for fc in FUEL_CONCS], float)

    xs, ys, dfun = delta_curve(args.bias)
    print(f"delta_past(theta) from §23.7, {len(xs)} pooled points:")
    print("  " + "  ".join(f"{a:.3f}->{b:.3f}" for a, b in zip(xs, ys)))

    meas = json.load(open(args.measured)) if args.measured.exists() else []
    dmeas = []
    for phi in phis:
        g = [r for r in meas if r.get("arm") == "fueled" and r.get("phi") == int(phi)
             and abs(r.get("sigma_frac", -1) - sf) < 1e-9]
        dmeas.append(g[0]["depth_median"] if g else float("nan"))

    t0 = time.time()
    print(f"\nOmega={om} sigma/delta*={sf} {args.trials} trials/cell")
    print(f"  reference: measured 0.6474+-0.022 | constant delta_past=0.20 gave "
          f"0.7919+-0.022 | hard stop 0.8867+-0.016")
    print(f"  P3 predicts ~0.75; P1 predicts small-tank depth 6-6.5 (was 5.47)\n")

    rows, conts, ints = [], [], []
    for fc, phi in zip(FUEL_CONCS, phis):
        ts = stage_time(GAMMA0, fc, 2.0)
        curve = []
        for i, t in enumerate(thetas):
            dp = dfun(float(t))
            c = hazard_at(float(t), om, int(phi), sigma_counts, ts,
                          args.trials, args.seed + i * 31 + int(fc), dp)
            c["delta_past_used"] = dp
            curve.append(c)
        conts.append(depth_continuous(curve, int(phi), theta0))
        ints.append(predict_depth(curve, int(phi), theta0)["depth_pred"])
        rows.append({"phi": int(phi), "sigma_frac": sf, "omega": om,
                     "depth_continuous": conts[-1], "depth_pred": ints[-1],
                     "curve": curve})
        dm = dmeas[len(conts) - 1]
        print(f"  Phi={int(phi):>6}  D_pred={conts[-1]:>6.2f}  D_meas={dm:>5.1f}"
              f"  ratio={conts[-1] / dm:>5.3f}")

    p, se = fit(phis, np.array(conts))
    print(f"\n  exponent {p:.4f} +- {se:.3f}    "
          f"[measured 0.6474, constant-delta 0.7919, hard stop 0.8867]")
    closed = (0.7919 - p) / (0.7919 - 0.6474) if p == p else float("nan")
    print(f"  closes a further {closed:.1%} of the residue left by §23.6")
    print(f"  depths: " + " ".join(f"{v:.2f}" for v in conts))
    print(f"  measured: " + " ".join(f"{v:.0f}" for v in dmeas))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(
        {"exponent": p, "exponent_se": se, "depths": conts,
         "depths_integer": ints, "measured": dmeas,
         "delta_curve": {"theta": xs.tolist(), "delta": ys.tolist()},
         "cells": rows}, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
