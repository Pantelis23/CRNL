"""T10b-iii-b: is the small-tank excess really a kinetic lag? — §23.4's attribution, tested

§23.4 attributed the 0.25 exponent gap -- quasi-static hazard integral 0.8925 vs
measured 0.6474 at sigma_ch/delta* = 0.03 -- to the bit being held past `gamma_c` on
kinetics, because the stage is two relaxation times at `gamma_0` while
`lambda(gamma_eff) = (1-2 gamma_eff)/3` vanishes at `gamma_c`. That was an
attribution, not a calculation. Two arms test it.

ARM A -- RAIL, the diagnostic. At the start of every stage the (X, Y, B) state is
re-seeded onto the attractor of its CURRENT `gamma_eff`, keeping the tank (F, W)
exactly as the chemistry left it. This imposes the quasi-static assumption by hand.

  **THIS ARM CONTAINS A FREE RESTORING ELEMENT AND IS NOT A PHYSICAL RESULT.**
  Re-seeding onto the rail is precisely the "harness doing work the chemistry
  cannot" that cost this project three withdrawn results, and it is used here only
  as an instrument: it makes the simulation obey the same assumption the hazard
  integral makes, so the two can be compared. No claim about a real network may be
  read off Arm A. It re-seeds (X, Y, B) only -- the tank is untouched and
  `n_X+n_Y+n_B = Omega` is preserved, so it is not also a fresh blank pool.

ARM B -- ADAPTIVE, the physical one. The stage time is set from the CURRENT state,
`t_stage = relax / (lambda(gamma_eff) * f)`, so every stage is two of the tank's own
present relaxation times rather than two of `gamma_0`'s. The state then has time to
follow the collapsing landscape. This diverges as `gamma_eff -> gamma_c`, which is
not a nuisance but the physics, so it is capped at `cap` times the baseline stage and
the cap is SWEPT -- rule 13, an approximation's own parameter is a second axis.

Both arms keep everything else identical to §23: same channel, same
`am_fueled(gamma_inf=1)`, same loss criterion `n_X <= n_Y`, same budgets. Arm B's
longer stages burn more fuel per stage, so its absolute depths are NOT comparable to
§23's; only the EXPONENT in `Phi` is, which is the quantity under test.

PREDICTIONS, written before running:

  P1  Arm A reproduces the hazard integral: exponent 0.89 +- a few percent, against
      §23's measured 0.6474. If Arm A does not land near 0.89, then the difference
      between the integral and the measurement is NOT the quasi-static assumption
      and §23.4's attribution fails immediately -- this is the cheap decisive test,
      and it also independently checks the integrator itself.
  P2  Arm B's exponent rises MONOTONICALLY with the cap, from ~0.65 at cap = 1
      (which is §23's own configuration, so cap = 1 must reproduce 0.647 -- a
      built-in control on the whole harness) toward Arm A's value.
  P3  Arm B never reaches Arm A. Full quasi-staticity needs `t_stage -> inf` at
      `gamma_c`, so at any finite cap some lag survives. If Arm B instead SATURATES
      well below Arm A -- say at 0.75 -- then time-scale lag is only part of the
      excess and something else holds the bit past `gamma_c`.
  P4  THE KILLER. If Arm B stays flat near 0.65 across a 30x cap sweep while Arm A
      sits at 0.89, the excess is not a time-scale effect at all, and §23.4's
      kinetic attribution must be withdrawn -- leaving the 0.25 gap named but
      unexplained, which is the honest state to be in.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.networks.am_fueled import am_fueled, gamma_effective, initial_counts
from crnl.networks.am_reversible import GAMMA_C, delta_star, lambda_antisym
from crnl.vectorized import compile_network, gillespie_fast
from experiments.cascade_fuel_vs_noise import _kick, stage_time

GAMMA0 = 0.3
GAMMA_INF = 1.0
RELAX = 2.0


def run_arm(arm: str, omega: int, fuel_conc: float, sigma_frac: float,
            trials: int, seed: int, *, cap: float = 1.0,
            max_depth: int = 600) -> dict:
    ds0 = delta_star(GAMMA0)
    sigma_counts = sigma_frac * ds0 * omega
    phi = int(round(fuel_conc * omega))
    w0 = int(round(phi * GAMMA0 / (1.0 + GAMMA0)))
    ts_base = stage_time(GAMMA0, fuel_conc, RELAX)
    comp = compile_network(am_fueled(GAMMA_INF), float(omega))
    n_start = initial_counts(omega, phi, w0, gamma_inf=GAMMA_INF)
    rng = np.random.default_rng(seed)

    depths, reseeds = [], []
    for _ in range(trials):
        n = n_start.copy()
        lost, nres = None, 0
        for d in range(1, max_depth + 1):
            w = int(n[4])
            g = float(gamma_effective(phi - w, w, GAMMA_INF))

            if arm == "rail" and g < GAMMA_C:
                # DIAGNOSTIC ONLY: impose the quasi-static assumption by hand.
                # (X, Y, B) go to the attractor of the current gamma_eff; the tank
                # is left exactly as the chemistry left it.
                seeded = initial_counts(omega, phi, w, gamma_inf=GAMMA_INF)
                if int(seeded[0]) > int(seeded[1]):
                    n = np.concatenate([seeded[:3], n[3:]])
                    nres += 1

            if arm == "adaptive":
                lam = lambda_antisym(g)
                f = max(int(n[3]), 1) / omega
                ts = RELAX / (lam * f) if lam > 1e-9 else np.inf
                ts = float(min(ts, cap * ts_base))
            else:
                ts = ts_base

            n = _kick(n, sigma_counts, rng)
            n = gillespie_fast(comp, n, rng, t_max=ts).n_final
            if int(n[0]) <= int(n[1]):
                lost = d
                break
        depths.append(lost if lost is not None else max_depth)
        reseeds.append(nres)
    return {"arm": arm, "cap": cap, "omega": omega, "fuel_conc": fuel_conc,
            "phi": phi, "sigma_frac": sigma_frac, "trials": trials,
            "t_stage_base": ts_base,
            "depth_median": float(np.median(depths)),
            "depth_mean": float(np.mean(depths)),
            "censored": float(np.mean([d >= max_depth for d in depths])),
            "reseeds_median": float(np.median(reseeds))}


def fit(x, y) -> tuple[float, float]:
    lx, ly = np.log(x), np.log(y)
    p, c = np.polyfit(lx, ly, 1)
    r = ly - (p * lx + c)
    se = np.sqrt(max(r @ r, 1e-30) / (len(x) - 2) / np.sum((lx - lx.mean()) ** 2))
    return float(p), float(se)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=40)
    ap.add_argument("--sigma-frac", type=float, default=0.03)
    ap.add_argument("--fuel-concs", type=float, nargs="+",
                    default=[25.0, 50.0, 100.0, 200.0, 400.0])
    ap.add_argument("--caps", type=float, nargs="+", default=[1.0, 3.0, 10.0, 30.0])
    ap.add_argument("--trials", type=int, default=80)
    ap.add_argument("--seed", type=int, default=20260729)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/fuel_quasistatic.json"))
    args = ap.parse_args()

    t0 = time.time()
    phis = np.array([int(round(fc * args.omega)) for fc in args.fuel_concs], float)
    print(f"Omega={args.omega} sigma/delta*={args.sigma_frac} "
          f"{args.trials} trials   [§23 measured exponent 0.6474 +- 0.022, "
          f"hazard integral 0.8925 +- 0.017]")
    rows = []

    print("\n=== ARM A: rail-reseeded (DIAGNOSTIC, contains a free restoring element)")
    ra = [run_arm("rail", args.omega, fc, args.sigma_frac, args.trials,
                  args.seed + int(fc)) for fc in args.fuel_concs]
    rows += ra
    print("  " + "  ".join(f"Phi={int(r['phi'])}:{r['depth_median']:.1f}"
                           f"(re{r['reseeds_median']:.0f})" for r in ra))
    pa, sa = fit(phis, np.array([r["depth_mean"] for r in ra]))
    print(f"  exponent {pa:.4f} +- {sa:.3f}   [P1 wants ~0.89]")

    print("\n=== ARM B: adaptive stage time, cap swept")
    print(f"{'cap':>6} {'depths (median)':>44} {'exponent':>16}")
    for cap in args.caps:
        rb = [run_arm("adaptive", args.omega, fc, args.sigma_frac, args.trials,
                      args.seed + int(fc), cap=cap) for fc in args.fuel_concs]
        rows += rb
        pb, sb = fit(phis, np.array([r["depth_mean"] for r in rb]))
        ds = " ".join(f"{r['depth_median']:.0f}" for r in rb)
        print(f"{cap:>6.0f} {ds:>44} {pb:>9.4f}+-{sb:<5.3f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
