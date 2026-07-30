"""T10b-iii-f: is the integral's burn rate a conditional-mean error? (rule 12, inside the instrument)

Raised by an independent analysis of §23's stored data, and the code confirms it.
`fuel_hazard.hazard_at` accumulates `waste.append(...)` for EVERY trial in the cell,
including the ones that lost the bit, and reports `c = mean(waste)`. The integral
then advances the SURVIVORS' burn fraction with a mean pooled over survivors and
losers together.

That is not a harmless average, because the network's own stoichiometry couples burn
to the thing being lost. `am_fueled`'s dominant consuming reaction is

    f1:  X + Y + F -> 2B + W        propensity proportional to n_X * n_Y * n_F

which is MAXIMAL at delta = 0. A trial that loses the bit rides delta ~ 0 through the
stage and burns hard; a survivor holds delta large and barely fires f1. In the killing
zone the measured `q` reaches 0.29-0.52, so the pooled mean is heavily weighted by
loser burn -- and the integral spends it on survivors. **This is rule 12 -- a
conditional quantity summarised unconditionally -- sitting inside the instrument, and
every sensitivity sweep in §23.6 and §23.8 swept `delta`'s effect on `q` while nobody
swept its effect on `c`.**

It also has the right shape. The correction is proportional to the fraction of life
spent where `q` is large: large at Phi/Omega = 25 (which dies at theta ~ 0.42 where
q ~ 0.4-0.5) and smaller at Phi/Omega = 400 (which dies at theta ~ 0.30 where
q ~ 0.10-0.24). A shape error across the budget is exactly what §23.7 said the residue
is and what no parametric repair has moved.

This measures `c_surv(theta)` and `c_lost(theta)` separately -- the ONLY change from
§23.9's instrument -- and re-integrates with the survivors' burn.

PREDICTIONS, written before running:

  P1  `c_lost / c_surv > 1`, and the excess is concentrated where `q` is large. If
      the ratio is 1 within a few percent everywhere, the stoichiometric argument is
      wrong and this candidate dies immediately.
  P2  Re-integrating with `c_surv` RAISES depth at every budget, because the tank
      lasts longer once survivors are charged only their own burn.
  P3  It raises the SMALL tanks proportionally more, because they spend a larger
      FRACTION of their life in the high-`q` region, so the exponent falls. But not
      to 0.6474 -- the large tanks also see q ~ 0.24 near the end, so they gain too,
      and the two partly cancel the way every correction in this thread has. I
      predict the exponent lands between **0.68 and 0.74**, i.e. this is worth more
      than theta-dispersion's 19.2% but does not close the residue alone.
  P4  The independent analysis predicts specifically that this reproduces ARM A
      (§23.5's rail-reseeded diagnostic: depths 7, 10, 16, 29, 52, exponent 0.7077),
      leaving inter-stage memory to carry the rest down to the measured 0.6474. That
      is a sharper prediction than mine and is recorded so it can be scored: if the
      depths land near 7/10/16/28/52 the decomposition is settled.
  P5  WHAT KILLS IT. If `c_surv` moves the exponent by less than the ~0.028 that
      theta-dispersion bought, then the burn side is no better than the hazard side
      and the residue is genuinely structural in the independence assumption alone.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.networks.am_fueled import am_fueled, gamma_effective, initial_counts
from crnl.networks.am_reversible import GAMMA_C, delta_star
from crnl.vectorized import compile_network, gillespie_fast
from experiments.cascade_fuel_vs_noise import _kick, stage_time
from experiments.fuel_hazard import GAMMA0, GAMMA_INF
from experiments.fuel_hazard_pastgc import (
    FUEL_CONCS, depth_continuous, fit, seed_past_gc,
)

DELTA_PAST = 0.20     # §23.8: the integral is insensitive to this


def cell(theta: float, omega: int, phi: int, sigma_counts: float, t_stage: float,
         trials: int, seed: int) -> dict:
    """q(theta) plus burn split by OUTCOME -- the one change from §23.9."""
    w = int(round(theta * phi))
    g = float(gamma_effective(phi - w, w, GAMMA_INF))
    n_start = (seed_past_gc(omega, phi, w, DELTA_PAST) if g >= GAMMA_C
               else initial_counts(omega, phi, w, gamma_inf=GAMMA_INF))
    if int(n_start[0]) <= int(n_start[1]):
        return {"theta": theta, "q": 1.0, "c": float("nan"),
                "c_surv": float("nan"), "c_lost": float("nan"),
                "n_surv": 0, "gamma_eff": g}
    comp = compile_network(am_fueled(GAMMA_INF), float(omega))
    rng = np.random.default_rng(seed)
    surv, lost = [], []
    for _ in range(trials):
        n = _kick(n_start.copy(), sigma_counts, rng)
        nf = gillespie_fast(comp, n, rng, t_max=t_stage).n_final
        burn = int(nf[4]) - int(n_start[4])
        (lost if int(nf[0]) <= int(nf[1]) else surv).append(burn)
    allb = np.array(surv + lost, float)
    return {
        "theta": theta, "gamma_eff": g,
        "q": len(lost) / trials,
        "c": float(allb.mean()),
        "c_surv": float(np.mean(surv)) if surv else float("nan"),
        "c_lost": float(np.mean(lost)) if lost else float("nan"),
        "n_surv": len(surv), "n_lost": len(lost),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=40)
    ap.add_argument("--sigma-frac", type=float, default=0.03)
    ap.add_argument("--n-theta", type=int, default=14)
    ap.add_argument("--trials", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260730)
    ap.add_argument("--measured", type=pathlib.Path,
                    default=pathlib.Path("results/fuel_depth_scaling.json"))
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/fuel_burn_conditioning.json"))
    args = ap.parse_args()

    om, sf = args.omega, args.sigma_frac
    sigma_counts = sf * delta_star(GAMMA0) * om
    theta0 = GAMMA0 / (1.0 + GAMMA0)
    thetas = np.linspace(theta0, 0.48, args.n_theta)
    phis = np.array([int(round(fc * om)) for fc in FUEL_CONCS], float)

    meas = json.load(open(args.measured)) if args.measured.exists() else []
    dmeas = []
    for phi in phis:
        g = [r for r in meas if r.get("arm") == "fueled" and r.get("phi") == int(phi)
             and abs(r.get("sigma_frac", -1) - sf) < 1e-9]
        dmeas.append(g[0]["depth_median"] if g else float("nan"))

    t0 = time.time()
    print(f"Omega={om} sigma/delta*={sf} {args.trials} trials/cell")
    print(f"  reference: measured 0.6474 | pooled-c integral 0.7919 | "
          f"ensemble 0.7641 | Arm A 0.7077")
    print(f"  P3 predicts 0.68-0.74; the independent analysis predicts ~0.71 "
          f"with depths 7/10/16/28/52\n")

    rows, d_all, d_surv = [], [], []
    for fc, phi in zip(FUEL_CONCS, phis):
        ts = stage_time(GAMMA0, fc, 2.0)
        curve = [cell(float(t), om, int(phi), sigma_counts, ts, args.trials,
                      args.seed + i * 31 + int(fc)) for i, t in enumerate(thetas)]
        c_s = [dict(c, c=c["c_surv"] if np.isfinite(c.get("c_surv", np.nan))
                    else c["c"]) for c in curve]
        d_all.append(depth_continuous(curve, int(phi), theta0))
        d_surv.append(depth_continuous(c_s, int(phi), theta0))
        rows.append({"phi": int(phi), "curve": curve,
                     "depth_pooled": d_all[-1], "depth_surv": d_surv[-1],
                     "depth_meas": dmeas[len(d_all) - 1]})
        print(f"  Phi={int(phi):>6}  pooled-c {d_all[-1]:>6.2f}  "
              f"surv-c {d_surv[-1]:>6.2f}  measured {dmeas[len(d_all)-1]:>5.1f}")

    print(f"\n=== P1: c_lost / c_surv, and where q is large (Phi=1000)")
    print(f"{'theta':>7} {'q':>7} {'c_surv':>8} {'c_lost':>8} {'ratio':>7} {'n_surv':>7}")
    for c in rows[0]["curve"]:
        if not np.isfinite(c.get("c_surv", np.nan)) or not np.isfinite(c.get("c_lost", np.nan)):
            continue
        print(f"{c['theta']:>7.4f} {c['q']:>7.4f} {c['c_surv']:>8.2f} "
              f"{c['c_lost']:>8.2f} {c['c_lost']/c['c_surv']:>7.3f} {c['n_surv']:>7}")

    pa, sa = fit(phis, np.array(d_all))
    ps, ss = fit(phis, np.array(d_surv))
    print(f"\n  pooled-c  exponent {pa:.4f} +- {sa:.3f}")
    print(f"  surv-c    exponent {ps:.4f} +- {ss:.3f}")
    print(f"  measured           0.6474 +- 0.022")
    print(f"  closes {(pa-ps)/(pa-0.6474):.1%} of the pooled integral's residue "
          f"(theta-dispersion bought 19.2%)")
    print(f"  depths surv-c: " + " ".join(f"{v:.2f}" for v in d_surv))
    print(f"  measured:      " + " ".join(f"{v:.0f}" for v in dmeas))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(
        {"exponent_pooled": pa, "exponent_pooled_se": sa,
         "exponent_surv": ps, "exponent_surv_se": ss,
         "measured_exponent": 0.6474, "cells": rows}, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
