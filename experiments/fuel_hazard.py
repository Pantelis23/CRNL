"""T10b-iii-a: does a hazard integral PREDICT the fuel-limited depth? — §23's mechanism, tested absolutely

§23 measured that the budget exponent drifts with the channel and the population,
so `D_fuel` is not a property of the tank. §23.3 offered a mechanism -- extra
stages each carry their own loss probability, so a bigger tank loses the bit at a
shallower `gamma_eff` -- and THEORIES T10b-iii-a proposed to test it by collapsing
the exponent onto a curve in the per-stage loss probability.

THAT TEST WAS BADLY POSED, and the reason matters. §23.3's wording says the extra
stages are "spent in the healthy part of the tank". But §23's own control arm
survives a mean of 377 stages at sigma_ch/delta* = 0.03, so the healthy-tank
per-stage loss probability there is at most ~0.0027, and over the 44 stages the
largest tank actually lasts that accumulates to ~12%. **It cannot be what kills the
big-tank cascade at a quiet channel.** Whatever the exponent is a function of, it is
not an accumulation of HEALTHY-tank hazard. §23.3 is corrected by whatever this
measures.

The repair keeps the structure and fixes the rate. The thing that accumulates is
the per-stage loss probability *at the tank's current state*, `q(theta)`, which
climbs as the burn fraction `theta` rises and `gamma_eff = gamma_inf*w/f` walks the
landscape toward death. Then the cascade's median depth is fixed with NO FREE
PARAMETER by a survival product:

    theta_{d+1} = theta_d + c(theta_d)/Phi ,     S_{d+1} = S_d * (1 - q(theta_d))

with `theta_0 = gamma_0/(1+gamma_0)` the seeded waste, and `q(theta)`, `c(theta)`
both MEASURED here in single stages from the attractor of that theta. Predicted
median depth is where `S` crosses 1/2. Rule 16: this is an absolute prediction
against a quantity already measured five times over in §23, not a fitted slope.

PREDICTIONS, written before running:

  P1  Healthy-tank hazard is negligible at the quiet channel -- `q(theta_0) <
      1/500` at sigma_ch/delta* = 0.03 -- confirming that §23.3's "healthy part of
      the tank" phrasing is wrong and has to be corrected rather than defended.
  P2  `q(theta)` rises steeply and is the whole story: the hazard integral
      reproduces §23's measured median depths within 1.5x at the LARGE budgets,
      with no fitted parameter.
  P3  IT SHOULD FAIL AT THE SMALL BUDGETS, and for a stated reason. The hazard
      integral assumes the state sits on the rail `delta*(gamma_eff)` of its
      current tank -- quasi-static. A small tank drains in a handful of stages, and
      `stage_time` is set by the relaxation time at `gamma_0`, not at the degraded
      `gamma_eff`, so the state LAGS the collapsing rail. §23 already shows the
      symptom: the smallest tank dies at `theta = 0.43`, i.e. `gamma_eff = 0.75`,
      well past `gamma_c = 0.5` where the rail does not exist at all. A bit cannot
      sit on a rail that is gone, so it is surviving on kinetics, and a
      quasi-static hazard integral must UNDER-predict its depth there.
  P4  Hence the exponent's drift needs no new physics: `q` depends on sigma and
      Omega, the tank sets how many stages are spent at each `theta`, and the two
      compose. If the integral reproduces the exponents at both sigma to within
      their fit error at the large budgets, the drift is explained.

WHAT WOULD KILL THE WHOLE PICTURE: if the integral misses the large-budget depths
by more than ~2x, or misses in the WRONG DIRECTION (over-predicting depth, when the
rail-lag argument says it can only under-predict), then accumulated hazard along the
tank's decline is not the mechanism and §23.3's paragraph should be withdrawn
outright rather than corrected.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.networks.am_fueled import am_fueled, death_waste_fraction, initial_counts
from crnl.networks.am_reversible import delta_star, lambda_antisym
from crnl.vectorized import compile_network, gillespie_fast
from experiments.cascade_fuel_vs_noise import _kick, stage_time

GAMMA0 = 0.3
GAMMA_INF = 1.0


def hazard_at(theta: float, omega: int, phi: int, sigma_counts: float,
              t_stage: float, trials: int, seed: int) -> dict:
    """One stage from the attractor of burn fraction `theta`: P(lose), waste made.

    The start state is `initial_counts` with `waste0 = theta*phi`, which places
    (X, Y, B) on the rail `delta*(gamma_eff(theta))` -- the quasi-static
    assumption P3 says will fail for small tanks. Reported `lost_at_start` is the
    fraction of cells where that rail has already collapsed to zero separation,
    because there the hazard is 1 by construction and carries no information.
    """
    w0 = int(round(theta * phi))
    n_start = initial_counts(omega, phi, w0, gamma_inf=GAMMA_INF)
    if int(n_start[0]) <= int(n_start[1]):
        return {"theta": theta, "q": 1.0, "c": float("nan"),
                "rail_gone": True, "trials": 0}

    comp = compile_network(am_fueled(GAMMA_INF), float(omega))
    rng = np.random.default_rng(seed)
    lost, waste = 0, []
    for _ in range(trials):
        n = _kick(n_start.copy(), sigma_counts, rng)
        nf = gillespie_fast(comp, n, rng, t_max=t_stage).n_final
        if int(nf[0]) <= int(nf[1]):
            lost += 1
        waste.append(int(nf[4]) - int(n_start[4]))
    return {"theta": theta, "q": lost / trials, "c": float(np.mean(waste)),
            "rail_gone": False, "trials": trials}


def predict_depth(curve: list[dict], phi: int, theta0: float,
                  max_depth: int = 4000) -> dict:
    """Median depth from the survival product. No fitted parameter."""
    ths = np.array([r["theta"] for r in curve])
    qs = np.array([r["q"] for r in curve])
    cs = np.array([r["c"] for r in curve])
    ok = np.isfinite(cs)
    theta, S, d = theta0, 1.0, 0
    while S > 0.5 and d < max_depth:
        q = float(np.interp(theta, ths, qs))
        c = float(np.interp(theta, ths[ok], cs[ok])) if ok.any() else 0.0
        if c <= 0:
            break
        S *= (1.0 - q)
        theta += c / phi
        d += 1
        if theta >= ths[-1]:            # ran off the measured curve
            return {"depth_pred": float(d), "ran_off": True, "S": S}
    return {"depth_pred": float(d), "ran_off": False, "S": S}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=40)
    ap.add_argument("--sigma-fracs", type=float, nargs="+", default=[0.03, 0.15])
    ap.add_argument("--fuel-concs", type=float, nargs="+",
                    default=[25.0, 50.0, 100.0, 200.0, 400.0])
    ap.add_argument("--n-theta", type=int, default=16)
    ap.add_argument("--trials", type=int, default=1500)
    ap.add_argument("--seed", type=int, default=20260729)
    ap.add_argument("--measured", type=pathlib.Path,
                    default=pathlib.Path("results/fuel_depth_scaling.json"))
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/fuel_hazard.json"))
    args = ap.parse_args()

    om = args.omega
    ds0 = delta_star(GAMMA0)
    theta0 = GAMMA0 / (1.0 + GAMMA0)
    wd = death_waste_fraction(GAMMA_INF)
    thetas = np.linspace(theta0, 0.48, args.n_theta)
    t0 = time.time()

    print(f"Omega={om} gamma0={GAMMA0} theta_0={theta0:.4f} "
          f"landscape death at theta={wd:.4f}  ({args.trials} trials/cell)")

    measured = json.load(open(args.measured)) if args.measured.exists() else []
    rows = []
    for sf in args.sigma_fracs:
        sigma_counts = sf * ds0 * om
        print(f"\n=== sigma/delta* = {sf}   (sigma = {sigma_counts:.2f} molecules)")
        # the hazard curve is measured once per (sigma, Phi) because c and the
        # lattice both depend on Phi; the largest tank sets the reference curve
        for fc in args.fuel_concs:
            phi = int(round(fc * om))
            ts = stage_time(GAMMA0, fc, 2.0)
            curve = [hazard_at(float(t), om, phi, sigma_counts, ts,
                               args.trials, args.seed + i * 31 + int(fc))
                     for i, t in enumerate(thetas)]
            pred = predict_depth(curve, phi, theta0)
            got = [r for r in measured
                   if r.get("arm") == "fueled" and r.get("phi") == phi
                   and abs(r.get("sigma_frac", -1) - sf) < 1e-9]
            dm = got[0]["depth_median"] if got else float("nan")
            ratio = pred["depth_pred"] / dm if dm == dm and dm else float("nan")
            q0 = curve[0]["q"]
            print(f"  Phi={phi:>6}  q(theta_0)={q0:<8.5f} "
                  f"D_pred={pred['depth_pred']:>7.1f}  D_meas={dm:>6.1f}  "
                  f"ratio={ratio:>6.3f}" + ("   [ran off curve]" if pred["ran_off"] else ""))
            rows.append({"omega": om, "sigma_frac": sf, "phi": phi, "fuel_conc": fc,
                         "t_stage": ts, "q_theta0": q0, "theta0": theta0,
                         "depth_pred": pred["depth_pred"], "depth_meas": dm,
                         "ratio": ratio, "ran_off": pred["ran_off"],
                         "curve": curve})

    print(f"\n=== P1: is healthy-tank hazard negligible at the quiet channel?")
    for sf in args.sigma_fracs:
        qs = [r["q_theta0"] for r in rows if abs(r["sigma_frac"] - sf) < 1e-9]
        print(f"  sigma/d*={sf}: q(theta_0) = {min(qs):.5f}-{max(qs):.5f}"
              f"   (1/500 = 0.00200)")

    print(f"\n=== P2/P3: does the hazard integral predict the measured depth?")
    for sf in args.sigma_fracs:
        rs = [r for r in rows if abs(r["sigma_frac"] - sf) < 1e-9
              and r["ratio"] == r["ratio"]]
        if not rs:
            continue
        print(f"  sigma/d*={sf}: " + "  ".join(
            f"Phi={r['phi']}:{r['ratio']:.2f}x" for r in rs))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
