"""Which ceiling binds for a CASCADE -- exhaustion or noise? (THEORIES T10b-ii)

FINDINGS 20 measured a fuel-limited LIFETIME for one held bit; 12.1 measured a
noise-limited DEPTH for a cascade. They are different units and were never put
side by side. This runs both through ONE harness so the comparison is not an
artifact of two conventions.

TWO ARMS, identical in every respect except whether the drive can run out:

  fueled   `am_fueled` with a finite tank Phi; gamma_eff = gamma_inf*w/f rises as
           it empties.
  control  `am_reversible(gamma0, k=f0)` -- the same drive held forever, with the
           rate scaled by the initial fuel concentration f0 so both arms share a
           clock. FINDINGS 20.1 explains why that scaling is not optional: the
           fueled network is third order, and comparing unscaled lifetimes is the
           mismatched-control error behind 10.3's withdrawn result.

THE CHANNEL. Between stages the decision coordinate is kicked by a Gaussian of sd
`sigma_frac * delta*`, applied by moving molecules between X and Y ONLY, so
n_X + n_Y, n_B and the tank are all untouched. FINDINGS 12's exact cascade
re-seeds onto an input lattice at a canonical blank level; not resetting n_B here
is deliberate -- a reset would be the harness handing the chemistry a fresh blank
pool it did not earn. Both arms use the identical channel, which is what makes the
comparison meaningful even though the convention differs slightly from 12's.

PREDICTIONS, written before running:

  P1  The fueled depth D_fuel is roughly INDEPENDENT of the channel noise: the
      tank is drained by the restoring chemistry, not by the channel.
  P2  The control depth D_noise grows like exp(delta*^2 / 2 sigma^2) as the
      channel quietens -- 12.1's ceiling -- so it is unbounded as sigma -> 0
      while D_fuel is not.
  P3  Hence a CROSSOVER: noise binds at large sigma, exhaustion binds at small
      sigma. Below the crossover a cascade dies of exhaustion, which is a failure
      mode this project has never seen.
  P4  THE SHARP ONE. The crossover sits at sigma_c^2 = delta*^2 / (2 ln D_fuel),
      and D_fuel is linear in Phi, so sigma_c moves only as 1/sqrt(ln Phi).
      **Fuel is a logarithmically weak lever**: tripling the budget should move
      the crossover noise by ~10%, not by 3x.

P1 was REFUTED (§20.3: D_fuel falls 2.2x-3.8x across sigma -- the ceilings
compound). P4's premise -- `D_fuel` linear in `Phi` -- was never tested; the two
budgets on hand gave `~Phi^0.67`. That is T10b-iii, and the second prediction
block below is for it.

PREDICTIONS for T10b-iii, written before reading any output of the multi-budget
run:

  P5  THE EXPONENT IS 1, not 0.67. Mechanism: the cascade does not die when the
      tank hits zero, it dies when `gamma_eff = gamma_inf*w/f` climbs far enough
      that `kappa(gamma_eff)*delta*(gamma_eff)^2` no longer beats the channel.
      That is a condition on the BURN FRACTION `w/(w+f)`, which is scale-free in
      `Phi`. Fuel burned per stage is set by how many molecules must be flipped,
      i.e. by Omega, not by the tank. Fixed usable fraction / fixed per-stage
      cost = linear. The measured 0.67 rests on two points and I expect it to be
      an artifact of exactly the compounding that killed P1: at the smaller
      budget the tank is a larger fraction spent by any given depth, so the
      channel bites earlier and drags the low-Phi point DOWN -- which bends the
      two-point chord UP, not down. So 0.67 is not even the sign I would predict
      for that contamination, and if the five-budget fit lands near 0.67 the
      burn-fraction picture is wrong.
  P6  THE LOAD-BEARING ONE. Whatever the exponent is, it is **the same at both
      sigma** within fit error. `D_fuel` is a property of the budget; sigma sets
      the prefactor (P1's refutation) but not the scaling. If the exponents
      differ between sigma/delta* = 0.03 and 0.08 by more than their uncertainty,
      then there is no such thing as "the fuel-limited depth" as a budget
      property, and 20.3's two-ceiling framing has to be reworded rather than
      refined.

OUTCOME of P5/P6 (five budgets, sigma/delta* = 0.03 and 0.08): P5 REFUTED, the
exponent is 0.664 +- 0.017 and 0.607 +- 0.050, with linear 20 and 8 sigma away.
P6 survived at 1.09 sigma -- but the test was weak, see P7.

WHERE THE SUBLINEARITY LIVES. Decomposing `D = (theta - theta_0) Phi / c` with
`theta` the waste fraction at loss and `theta_0 = gamma_0/(1+gamma_0)` the seeded
waste: fuel burned per stage `c` is FLAT (25-32 molecules over a 16x budget), as
the a-priori argument requires -- burn rate is third order so it goes like `f`,
and `stage_time` goes like `1/f`. The whole effect is that the USABLE burn
fraction collapses, 0.198 -> 0.069. Equivalently the cascade dies at a shallower
`gamma_eff`: 0.75 at Phi=1000 (past gamma_c, the landscape was already dead) but
only 0.43 at Phi=16000 (landscape still alive).

  MECHANISM: a bigger tank does not buy proportionally more depth because the
  extra stages each carry their own per-stage loss probability. More stages
  spent in the healthy part of the tank means more accumulated chances to lose
  the bit, so the bit goes while the drive is still strong. This is 20.3's
  compounding, seen from the budget side.

  P7, THE FORK, written before running the wider sweep. That mechanism makes the
  exponent a function of the CHANNEL, not of the budget: raise sigma and each
  stage is likelier to lose the bit, so the accumulation bites sooner and the
  exponent must FALL. **I predict the exponent drops monotonically across
  sigma/delta* = 0.03, 0.15, 0.30, refuting P6 once the channel is wide enough
  to matter.** P6 passed only because 0.03 and 0.08 sit on the SAME side of
  20.3's crossover -- both exhaustion-bound, depths differing by 1.1x -- so it
  never asked the question. If instead the exponent holds near 0.66 out to 0.30,
  the accumulation mechanism above is wrong and D_fuel really is a budget
  property, which is the answer T10b-iii asked for but not the one I expect.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.networks.am_fueled import (
    am_fueled, death_waste_fraction, gamma_effective, initial_counts,
)
from crnl.networks.am_reversible import (
    am_reversible, delta_star, lambda_antisym,
)
from crnl.vectorized import compile_network, gillespie_fast


def _kick(n: np.ndarray, sigma_counts: float, rng) -> np.ndarray:
    """Channel: move molecules between X and Y, conserving n_X + n_Y."""
    d = int(round(rng.normal(0.0, sigma_counts)))
    if d == 0:
        return n
    m = n.copy()
    move = min(abs(d), int(m[1]) if d > 0 else int(m[0]))
    if d > 0:
        m[0] += move
        m[1] -= move
    else:
        m[0] -= move
        m[1] += move
    return m


def stage_time(gamma0: float, fuel_conc: float, relax: float) -> float:
    """A stage of `relax` relaxation times.

    Fixing t_stage in ABSOLUTE time is what made the first run useless: at fuel
    concentration 10 the whole tank lasts ~7.5 time units (FINDINGS 20) while a
    stage was 8, so the cascade died at depth 1. The meaningful depth is the
    lifetime measured in RESTORATION times, and the relaxation time is
    1/(lambda * f0), which itself shrinks as the tank gets richer. Holding the
    stage at a fixed multiple of it is what makes depths comparable across fuel
    concentrations -- and it is why D_fuel grows with the budget even though
    FINDINGS 20's absolute lifetime shrinks.
    """
    f0 = fuel_conc * (1.0 - gamma0 / (1.0 + gamma0))
    return float(relax / (lambda_antisym(gamma0) * f0))


def run_arm(arm: str, omega: int, sigma_frac: float, t_stage: float,
            max_depth: int, trials: int, seed: int, *,
            fuel_conc: float, gamma0: float, gamma_inf: float) -> dict:
    ds = delta_star(gamma0)
    sigma_counts = sigma_frac * ds * omega
    f0 = fuel_conc * (1.0 - gamma0 / (1.0 + gamma0))
    phi = int(round(fuel_conc * omega))
    w0 = int(round(phi * gamma0 / (1.0 + gamma0)))

    if arm == "fueled":
        net, n_start = am_fueled(gamma_inf), initial_counts(
            omega, phi, w0, gamma_inf=gamma_inf)
    else:
        net = am_reversible(gamma0, k=f0)
        base = initial_counts(omega, phi, w0, gamma_inf=gamma_inf)
        n_start = base[:3].copy()
    comp = compile_network(net, float(omega))
    rng = np.random.default_rng(seed)

    depths, wfracs, survived = [], [], 0
    for _ in range(trials):
        n = n_start.copy()
        lost = None
        for d in range(1, max_depth + 1):
            n = _kick(n, sigma_counts, rng)
            n = gillespie_fast(comp, n, rng, t_max=t_stage).n_final
            if int(n[0]) <= int(n[1]):
                lost = d
                break
        if lost is None:
            survived += 1
            depths.append(max_depth)
        else:
            depths.append(lost)
        if arm == "fueled":
            wfracs.append(float(n[4]) / phi)
    out = {
        "arm": arm, "omega": omega, "sigma_frac": sigma_frac,
        "t_stage": t_stage, "trials": trials, "max_depth": max_depth,
        "depth_median": float(np.median(depths)),
        "depth_mean": float(np.mean(depths)),
        "censored": survived / trials,
    }
    if arm == "fueled":
        out["fuel_conc"] = fuel_conc
        out["phi"] = phi
        out["waste_frac_at_loss"] = float(np.median(wfracs))
        out["waste_frac_at_death"] = death_waste_fraction(gamma_inf)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=40)
    ap.add_argument("--sigma-fracs", type=float, nargs="+",
                    default=[0.15, 0.22, 0.30, 0.40, 0.50])
    ap.add_argument("--fuel-concs", type=float, nargs="+", default=[10.0, 30.0])
    ap.add_argument("--gamma0", type=float, default=0.30)
    ap.add_argument("--gamma-inf", type=float, default=1.0)
    ap.add_argument("--stage-relax", type=float, default=2.0,
                    help="stage length in relaxation times (see stage_time)")
    ap.add_argument("--max-depth", type=int, default=250)
    ap.add_argument("--trials", type=int, default=40)
    ap.add_argument("--seed", type=int, default=20260728)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/cascade_fuel_vs_noise.json"))
    args = ap.parse_args()

    t0 = time.time()
    ds = delta_star(args.gamma0)
    print(f"Omega={args.omega} gamma0={args.gamma0} delta*={ds:.4f} "
          f"stage={args.stage_relax} relaxation times, max_depth={args.max_depth}")
    for fc in args.fuel_concs:
        print(f"  fuel conc {fc}: t_stage = "
              f"{stage_time(args.gamma0, fc, args.stage_relax):.4f}")
    rows = []

    print(f"\n{'sigma/d*':>9} {'CONTROL depth':>14} {'cens':>6}", end="")
    for fc in args.fuel_concs:
        print(f" | {'FUEL ' + str(fc):>12} {'cens':>6} {'w/wd':>6}", end="")
    print()
    for sf in args.sigma_fracs:
        ts0 = stage_time(args.gamma0, args.fuel_concs[0], args.stage_relax)
        c = run_arm("control", args.omega, sf, ts0, args.max_depth,
                    args.trials, args.seed, fuel_conc=args.fuel_concs[0],
                    gamma0=args.gamma0, gamma_inf=args.gamma_inf)
        rows.append(c)
        line = f"{sf:>9.2f} {c['depth_median']:>14.1f} {c['censored']:>6.2f}"
        for fc in args.fuel_concs:
            f = run_arm("fueled", args.omega, sf,
                        stage_time(args.gamma0, fc, args.stage_relax),
                        args.max_depth, args.trials, args.seed + 7, fuel_conc=fc,
                        gamma0=args.gamma0, gamma_inf=args.gamma_inf)
            rows.append(f)
            line += (f" | {f['depth_median']:>12.1f} {f['censored']:>6.2f} "
                     f"{f['waste_frac_at_loss']/f['waste_frac_at_death']:>6.2f}")
        print(line)

    print("\n=== P1: is the fueled depth flat in the channel noise?")
    for fc in args.fuel_concs:
        fs = [r for r in rows if r["arm"] == "fueled" and r["fuel_conc"] == fc
              and r["censored"] < 0.2]
        if len(fs) < 2:
            print(f"  fuel {fc}: fewer than 2 uncensored cells; nothing quotable")
            continue
        d = [r["depth_median"] for r in fs]
        print(f"  fuel {fc}: depth {min(d):.0f}-{max(d):.0f} "
              f"({max(d)/max(min(d),1e-9):.2f}x) over sigma/d* "
              f"{min(r['sigma_frac'] for r in fs):.2f}-"
              f"{max(r['sigma_frac'] for r in fs):.2f}")
    print("\n=== P3/P4: where does the binding ceiling change hands?")
    for fc in args.fuel_concs:
        prev = None
        for sf in args.sigma_fracs:
            c = next(r for r in rows if r["arm"] == "control"
                     and r["sigma_frac"] == sf)
            f = next(r for r in rows if r["arm"] == "fueled"
                     and r["sigma_frac"] == sf and r["fuel_conc"] == fc)
            binds = "fuel" if f["depth_median"] < 0.8 * c["depth_median"] else "noise"
            if prev is not None and binds != prev:
                print(f"  fuel {fc}: hands over between sigma/d* "
                      f"{prev_sf:.2f} and {sf:.2f}")
            prev, prev_sf = binds, sf
        if prev is not None:
            print(f"  fuel {fc}: at the quietest channel tested "
                  f"({min(args.sigma_fracs):.2f}) the binding ceiling is "
                  f"{'fuel' if True else ''}"
                  f"{'' if True else ''}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
