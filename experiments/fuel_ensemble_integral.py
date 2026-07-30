"""T10b-iii-e: is the integral's error its FORM rather than its numbers?

Three repairs to §23.4's hazard integral have now been tested and all changed only
the values it is fed: imposing quasi-staticity on the state (§23.5, 1.6 sigma from
doing nothing), removing the hard stop past `gamma_c` (§23.6, 39.6%), and correcting
the separation imposed there (§23.8, 2.7%). +0.14 of exponent survives, as -22% at
the smallest tank and +10% at the largest -- a SHAPE error across `Phi` that no
parametric correction has moved.

So look at the form instead. The integral makes exactly two structural approximations:

    theta_{d+1} = theta_d + c(theta_d)/Phi      <- deterministic burn
    S_{d+1}     = S_d * (1 - q(theta_d))        <- independent stages

and the two surviving errors have the signs and the `Phi`-scalings those two
approximations would produce.

  (a) DETERMINISTIC THETA. Real burn is stochastic, so at stage `d` there is a
      DISTRIBUTION of `theta`, not one value. Trials that burn fast reach high
      `theta` and die, so the surviving population is biased toward slow burners
      sitting BELOW the mean. `q(theta)` rises steeply, so evaluating it at the mean
      overstates the hazard among survivors and the integral under-predicts depth.
      Relative spread of accumulated burn goes like `1/sqrt(Phi)`, so this is
      strongest at small tanks -- the sign and the scaling of the -22%.
  (b) INDEPENDENT STAGES. `delta` carries memory across stages: a bit knocked down
      stays down. Correlated trials die sooner than an independent product predicts,
      so the integral over-predicts -- the +10%, present everywhere and dominant at
      large tanks where nothing else acts.

This tests (a) ONLY, and does so without touching a single hazard value. The same
measured `q(theta)` and the same measured burn are used; the only change is that
`theta` is propagated as an ensemble of stochastic trajectories rather than one
deterministic path, with the per-stage burn drawn from its MEASURED mean and spread
at that `theta`. If (a) is the small-tank error, the ensemble median depth rises at
small `Phi` and barely moves at large `Phi`. (b) is not addressed here: the ensemble
still applies `q(theta)` independently each stage.

PREDICTIONS, written before running:

  P1  Direction and scaling: ensemble depth >= deterministic depth at every budget,
      with the gain FALLING as `Phi` rises. This is close to forced by Jensen given a
      convex `q`, so it is not much of a test on its own -- P2 is the test.
  P2  MAGNITUDE, and I am genuinely unsure. A crude estimate says this may be too
      small: per-stage burn is ~25 molecules, so if it were POISSON the per-stage sd
      would be ~5, and over the ~7 stages a small tank lasts the accumulated spread
      is `5*sqrt(7) = 13` molecules against ~175 burned, i.e. **7.5%** -- against
      **2.9%** at the largest tank. A 2.6x ratio in a 7.5% dispersion producing a 22%
      depth error requires the hazard to be very steep in `theta` where it matters,
      or the burn to be much burstier than Poisson. **I predict the ensemble closes
      LESS THAN HALF of the -22% at Phi/Omega = 25 unless the measured burn sd is
      well above the Poisson value**, and the measured sd is reported alongside so
      the two can be told apart.
  P3  The exponent falls from 0.79 but NOT to 0.6474, because (b) is untouched and
      acts in the opposite direction at the large-tank end. Landing at ~0.72-0.75
      would support the two-error picture; landing at 0.647 would mean (a) is the
      whole story and (b) is not needed, which would be a cleaner result than I
      expect.
  P4  WHAT KILLS THE READING. If the ensemble changes the exponent by less than its
      fit error, then theta-dispersion is not the small-tank error either, the
      structural reading in §23.8's closing paragraph is wrong, and the residue has
      no surviving suspect of any kind -- parametric or structural.
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
from experiments.fuel_hazard_pastgc import FUEL_CONCS, depth_continuous, fit, seed_past_gc

DELTA_PAST = 0.20     # §23.8: the integral is insensitive to this, so it is fixed


def cell(theta: float, omega: int, phi: int, sigma_counts: float, t_stage: float,
         trials: int, seed: int) -> dict:
    """q(theta), and the burn distribution -- mean AND sd, which is the new part."""
    w = int(round(theta * phi))
    g = float(gamma_effective(phi - w, w, GAMMA_INF))
    n_start = (seed_past_gc(omega, phi, w, DELTA_PAST) if g >= GAMMA_C
               else initial_counts(omega, phi, w, gamma_inf=GAMMA_INF))
    if int(n_start[0]) <= int(n_start[1]):
        return {"theta": theta, "q": 1.0, "c": float("nan"), "c_sd": float("nan"),
                "gamma_eff": g}
    comp = compile_network(am_fueled(GAMMA_INF), float(omega))
    rng = np.random.default_rng(seed)
    lost, waste = 0, []
    for _ in range(trials):
        n = _kick(n_start.copy(), sigma_counts, rng)
        nf = gillespie_fast(comp, n, rng, t_max=t_stage).n_final
        if int(nf[0]) <= int(nf[1]):
            lost += 1
        waste.append(int(nf[4]) - int(n_start[4]))
    w_arr = np.array(waste, float)
    return {"theta": theta, "q": lost / trials, "c": float(w_arr.mean()),
            "c_sd": float(w_arr.std(ddof=1)), "gamma_eff": g,
            "poisson_sd": float(np.sqrt(max(w_arr.mean(), 0.0)))}


def ensemble_depth(curve: list[dict], phi: int, theta0: float, n_traj: int,
                   rng: np.random.Generator, max_depth: int = 4000) -> float:
    """Median depth over stochastic theta trajectories. Same q, same burn.

    The crossing of S = 1/2 is INTERPOLATED, matching `depth_continuous`. Taking
    `np.median` of integer depths instead quantises the answer and inflated this
    experiment's headline from 19.2% to 27.3% on the first pass -- the same fault
    §23.6 caught in `predict_depth`, reintroduced here two sections later.
    """
    ths = np.array([c["theta"] for c in curve])
    qs = np.array([c["q"] for c in curve])
    cs = np.array([c["c"] for c in curve])
    sds = np.array([c["c_sd"] for c in curve])
    ok = np.isfinite(cs)
    if not ok.any():
        return float("nan")
    th_max = ths[ok][-1]

    theta = np.full(n_traj, theta0)
    alive = np.ones(n_traj, bool)
    surv = [1.0]
    for _ in range(max_depth):
        if not alive.any():
            break
        idx = np.flatnonzero(alive)
        t = theta[idx]
        died = rng.random(t.size) < np.interp(t, ths, qs)
        # burn drawn from the measured mean and spread at the current theta
        c = np.interp(t, ths[ok], cs[ok])
        s = np.interp(t, ths[ok], sds[ok])
        theta[idx] = t + np.maximum(rng.normal(c, s), 0.0) / phi
        alive[idx[died]] = False
        alive[idx[theta[idx] > th_max]] = False     # rail gone, bit cannot be held
        surv.append(alive.sum() / n_traj)
    surv = np.array(surv)
    for i in range(len(surv) - 1):
        if surv[i] >= 0.5 >= surv[i + 1]:
            return i + (surv[i] - 0.5) / (surv[i] - surv[i + 1])
    return float(len(surv) - 1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=40)
    ap.add_argument("--sigma-frac", type=float, default=0.03)
    ap.add_argument("--n-theta", type=int, default=14)
    ap.add_argument("--trials", type=int, default=800)
    ap.add_argument("--n-traj", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=20260730)
    ap.add_argument("--measured", type=pathlib.Path,
                    default=pathlib.Path("results/fuel_depth_scaling.json"))
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/fuel_ensemble_integral.json"))
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
    rng = np.random.default_rng(args.seed + 991)
    print(f"Omega={om} sigma/delta*={sf} {args.trials} trials/cell, "
          f"{args.n_traj} ensemble trajectories")
    print(f"  reference: measured 0.6474+-0.022 | deterministic integral ~0.79\n")
    print(f"{'Phi':>7} {'determ':>8} {'ensemble':>9} {'measured':>9} "
          f"{'det/meas':>9} {'ens/meas':>9} {'burn sd/Poisson':>17}")

    rows, det, ens = [], [], []
    for fc, phi in zip(FUEL_CONCS, phis):
        ts = stage_time(GAMMA0, fc, 2.0)
        curve = [cell(float(t), om, int(phi), sigma_counts, ts, args.trials,
                      args.seed + i * 31 + int(fc)) for i, t in enumerate(thetas)]
        d = depth_continuous(curve, int(phi), theta0)
        e = ensemble_depth(curve, int(phi), theta0, args.n_traj, rng)
        det.append(d)
        ens.append(e)
        rat = [c["c_sd"] / c["poisson_sd"] for c in curve
               if np.isfinite(c.get("c_sd", np.nan)) and c.get("poisson_sd", 0) > 0]
        dm = dmeas[len(det) - 1]
        print(f"{int(phi):>7} {d:>8.2f} {e:>9.2f} {dm:>9.1f} "
              f"{d/dm:>9.3f} {e/dm:>9.3f} {np.median(rat):>17.2f}")
        rows.append({"phi": int(phi), "depth_det": d, "depth_ens": e,
                     "depth_meas": dm, "curve": curve})

    pd_, sd_ = fit(phis, np.array(det))
    pe, se = fit(phis, np.array(ens))
    print(f"\n  deterministic exponent {pd_:.4f} +- {sd_:.3f}")
    print(f"  ensemble      exponent {pe:.4f} +- {se:.3f}")
    print(f"  measured               0.6474 +- 0.022")
    gap = pd_ - 0.6474
    print(f"  closes {(pd_-pe)/gap:.1%} of the deterministic integral's residue"
          if gap else "")
    small = (ens[0] - det[0]) / (dmeas[0] - det[0]) if dmeas[0] != det[0] else float("nan")
    print(f"  P2: closes {small:.1%} of the -22% at the smallest tank "
          f"(predicted < 50% unless burn is super-Poisson)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(
        {"exponent_det": pd_, "exponent_det_se": sd_, "exponent_ens": pe,
         "exponent_ens_se": se, "measured_exponent": 0.6474,
         "cells": rows}, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
