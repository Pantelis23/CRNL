"""Does the bit outrun its own landscape? — the correction §23.3 needs

§23.3 said the fuel-limited depth is sublinear in the budget because the extra
stages a bigger tank buys each carry their own per-stage loss probability, so the
bit is lost while the drive is still strong. `fuel_hazard.py` tested that as an
absolute hazard integral and it is mostly WRONG:

  - accumulated hazard cuts the predicted depth by 0.0% at four of five budgets
    (7.7% at the largest, quiet channel), so it cannot be carrying the effect;
  - with hazard switched off entirely, "burn until the rail is gone" gives
    exponent 0.916 -- close to the linear value §23's P5 argued for;
  - the integral reproduces the sigma-DRIFT of the exponent (-0.121 predicted vs
    -0.107 measured) but over-predicts its LEVEL by 0.25 (0.893 vs 0.647);
  - and the entire remaining gap is that the integral UNDER-predicts the small
    tanks by up to 1.75x (ratio 0.571 at Phi = 1000) while matching the large ones.

So the sublinearity is mostly not big tanks under-performing. It is **small tanks
over-performing**, exactly where `fuel_hazard.py`'s P3 said a quasi-static
prediction must fail. §23 already recorded the symptom without reading it: the
smallest tank loses the bit at `theta = 0.43`, i.e. `gamma_eff = 0.75`, past
`gamma_c = 0.5` -- where `delta*` does not exist and the landscape is MONOSTABLE.
A bit cannot sit on a rail that is gone, so it is being carried by kinetics: the
stage time is two relaxation times at `gamma_0`, but the relaxation rate
`lambda(gamma_eff)` vanishes at `gamma_c`, so as the tank empties the state needs
ever longer to follow the collapsing landscape and simply stops following it.

This measures the lag directly instead of inferring it from a residual.

PREDICTIONS, written before running:

  P1  The carried separation exceeds the rail, and by more for smaller tanks:
      `delta_meas / delta*(gamma_eff)` climbs above 1 as the tank drains, and at
      matched `theta` is LARGER at small Phi.
  P2  Small tanks spend a substantial number of stages at `gamma_eff > gamma_c`,
      holding a bit in a monostable landscape; large tanks spend ~none. This is
      the whole of the small-tank over-performance.
  P3  The large tank is quasi-static: at Phi/Omega = 400 the loss happens at
      `gamma_eff` where the rail still exists and `delta_meas ~ delta*`, which is
      why `fuel_hazard.py` matched it to 1.09x.
  P4  THE ONE THAT COULD KILL THE READING. If instead `delta_meas/delta*` is ~1
      for every budget right up to loss, then there is no lag, the bit is not
      outrunning its landscape, and the small-tank excess has some other source --
      in which case §23.3 should be withdrawn outright with no replacement rather
      than corrected.
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

GAMMA0 = 0.3
GAMMA_INF = 1.0


def trace(omega: int, fuel_conc: float, sigma_frac: float, trials: int,
          seed: int, max_depth: int = 600) -> dict:
    """Per-stage (theta, delta, rail) until the bit is lost."""
    ds0 = delta_star(GAMMA0)
    sigma_counts = sigma_frac * ds0 * omega
    phi = int(round(fuel_conc * omega))
    w0 = int(round(phi * GAMMA0 / (1.0 + GAMMA0)))
    ts = stage_time(GAMMA0, fuel_conc, 2.0)
    comp = compile_network(am_fueled(GAMMA_INF), float(omega))
    n_start = initial_counts(omega, phi, w0, gamma_inf=GAMMA_INF)
    rng = np.random.default_rng(seed)

    depths, past_gc, g_loss, lag_loss, curves = [], [], [], [], []
    for t in range(trials):
        n = n_start.copy()
        stages = []
        for d in range(1, max_depth + 1):
            n = _kick(n, sigma_counts, rng)
            n = gillespie_fast(comp, n, rng, t_max=ts).n_final
            w = int(n[4])
            g = float(gamma_effective(phi - w, w, GAMMA_INF))
            delta = (int(n[0]) - int(n[1])) / omega
            rail = delta_star(g) if g < GAMMA_C else 0.0
            stages.append((w / phi, g, delta, rail))
            if int(n[0]) <= int(n[1]):
                break
        depths.append(len(stages))
        past_gc.append(sum(1 for s in stages if s[1] >= GAMMA_C))
        # state one stage BEFORE loss: the last stage still holding the bit
        hold = stages[-2] if len(stages) >= 2 else stages[-1]
        g_loss.append(hold[1])
        lag_loss.append(hold[2] / hold[3] if hold[3] > 1e-9 else np.inf)
        if t < 12:
            curves.append(stages)
    fin = [v for v in lag_loss if np.isfinite(v)]
    return {
        "omega": omega, "fuel_conc": fuel_conc, "phi": phi,
        "sigma_frac": sigma_frac, "t_stage": ts, "trials": trials,
        "depth_median": float(np.median(depths)),
        "stages_past_gc_median": float(np.median(past_gc)),
        "frac_trials_past_gc": float(np.mean([p > 0 for p in past_gc])),
        "gamma_eff_at_hold_median": float(np.median(g_loss)),
        "frac_hold_beyond_rail": float(np.mean([not np.isfinite(v) for v in lag_loss])),
        "lag_at_hold_median": float(np.median(fin)) if fin else float("inf"),
        "example_traces": curves,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=40)
    ap.add_argument("--sigma-fracs", type=float, nargs="+", default=[0.03, 0.15])
    ap.add_argument("--fuel-concs", type=float, nargs="+",
                    default=[25.0, 50.0, 100.0, 200.0, 400.0])
    ap.add_argument("--trials", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20260729)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/fuel_rail_lag.json"))
    args = ap.parse_args()

    t0 = time.time()
    print(f"Omega={args.omega} gamma_c={GAMMA_C}  (last stage still holding the bit)")
    rows = []
    for sf in args.sigma_fracs:
        print(f"\n=== sigma/delta* = {sf}")
        print(f"{'Phi':>7} {'D_med':>7} {'g_eff@hold':>11} {'delta/rail':>11} "
              f"{'held past rail':>15} {'stages past g_c':>16} {'any past g_c':>13}")
        for fc in args.fuel_concs:
            r = trace(args.omega, fc, sf, args.trials, args.seed + int(fc))
            lag = r["lag_at_hold_median"]
            print(f"{r['phi']:>7} {r['depth_median']:>7.1f} "
                  f"{r['gamma_eff_at_hold_median']:>11.4f} "
                  f"{lag:>11.3f} {r['frac_hold_beyond_rail']:>15.1%} "
                  f"{r['stages_past_gc_median']:>16.1f} "
                  f"{r['frac_trials_past_gc']:>13.1%}")
            rows.append(r)

    print("\n=== P2: are the small tanks' extra stages spent past gamma_c?")
    for sf in args.sigma_fracs:
        rs = [r for r in rows if abs(r["sigma_frac"] - sf) < 1e-9]
        print(f"  sigma/d*={sf}: " + "  ".join(
            f"Phi={r['phi']}:{r['stages_past_gc_median']:.0f}/"
            f"{r['depth_median']:.0f}" for r in rs))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
