"""A second ceiling: restoration that dies of exhaustion, not of noise.

FINDINGS 12.1's depth ceiling is purely noise-limited, because `gamma` has always
been an infinite reservoir -- set it once and it holds forever at no cost.
FINDINGS 9 therefore measures what restoration *dissipates* while nothing ever
runs down. `networks/am_fueled` makes the drive a finite fuel species, so the
effective drive is gamma_eff = gamma_inf * w/f and it RISES as the tank empties --
the mirror of FINDINGS 19's cooling, which drove gamma down.

PREDICTIONS, written before running (they also appear in the network module):

  P1  Restoration dies at waste fraction 1/(1 + 2 gamma_inf), = 1/3 at
      gamma_inf = 1, where gamma_eff crosses gamma_c = 1/2. Parameter-free.
  P2  Fuel consumed = net forward firings, the quantity FINDINGS 9 already
      measures. Its "cost of remembering" becomes a LIFETIME rather than a rate.
  P3  THE HEADLINE. At fixed fuel CONCENTRATION the budget and the burn rate are
      both extensive, so the fuel-limited lifetime is Omega-INDEPENDENT, while the
      noise-limited lifetime grows like exp(Omega). Beyond a crossover Omega
      restoration is fuel-limited and MORE MOLECULES BUY NOTHING -- the exact
      mirror of FINDINGS 1's wall, where molecules bought exponential reliability.
  P4  The bit is lost BEFORE the formal death point, because gamma_eff rises
      continuously and the barrier degrades all the way up to it. The loss
      fraction w_loss/w_death should rise toward 1 as Omega grows.

WHAT THE FIRST PROBE GOT WRONG, kept because it is the same class of error this
file keeps making. Runs were let go to completion, and they do not stop at the
death point -- the chemistry sails through it and relaxes to equilibrium
(gamma_eff -> 1, waste -> 0.49, committed populations equal). The death point is
where RESTORATION ends, not where reactions end, so the measurement has to stop at
the loss event. Reading the end state instead would have reported "waste 0.49" and
buried P1.

Loss is defined as the first state with n_X <= n_Y from a start committed to X.
That is a readout, never an intervention -- the run is not restarted, nudged, or
re-railed, which is the failure mode (harness doing the chemistry) behind three
withdrawn results in this project.
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
from crnl.networks.am_reversible import am_reversible, reverse_pairing
from crnl.thermo import gillespie_instrumented
from crnl.vectorized import compile_network

MAX_STEPS = 40_000_000


def fueled_lifetime(omega: int, fuel_conc: float, gamma0: float,
                    gamma_inf: float, trials: int, seed: int) -> dict:
    """Time until the bit is lost, with a finite tank."""
    phi = int(round(fuel_conc * omega))
    w0 = int(round(phi * gamma0 / (1.0 + gamma0)))
    n0 = initial_counts(omega, phi, w0, gamma_inf=gamma_inf)
    net = am_fueled(gamma_inf)
    comp = compile_network(net, float(omega))
    pair = reverse_pairing(net)
    rng = np.random.default_rng(seed)

    times, wfracs, hit_budget = [], [], 0
    for _ in range(trials):
        r = gillespie_instrumented(comp, n0, rng, pair,
                                   stop=lambda n: n[0] <= n[1],
                                   max_steps=MAX_STEPS)
        if r.steps >= MAX_STEPS:
            hit_budget += 1
            continue
        times.append(r.t_final)
        wfracs.append(float(r.n_final[4]) / phi)
    if not times:
        return {}
    return {
        "arm": "fueled", "omega": omega, "fuel_conc": fuel_conc, "phi": phi,
        "gamma0": gamma0, "gamma_inf": gamma_inf, "trials": trials,
        "n_ok": len(times), "budget_hits": hit_budget,
        "lifetime_median": float(np.median(times)),
        "lifetime_mean": float(np.mean(times)),
        "lifetime_sem": float(np.std(times) / np.sqrt(len(times))),
        "waste_frac_at_loss": float(np.median(wfracs)),
        "waste_frac_at_death": death_waste_fraction(gamma_inf),
        "loss_over_death": float(np.median(wfracs)) / death_waste_fraction(gamma_inf),
    }


def noise_lifetime(omega: int, gamma: float, trials: int, seed: int,
                   rate_scale: float = 1.0) -> dict:
    """Control: the SAME drive held forever. Loss can only be by noise.

    `rate_scale` MUST be the initial fuel concentration f0, and this is not a
    detail. The fueled network is third order, so its forward rate carries a
    factor f -- at fuel concentration 10 with gamma_eff(0) = 0.3 that is f0 = 7.7,
    making the fueled chemistry 7.7x faster in absolute time than an unscaled
    `am_reversible`. Comparing raw lifetimes across the two arms without this is
    comparing different clocks, which is exactly the mismatched-control error that
    produced FINDINGS 10.3's withdrawn result. With the scale applied, the two
    arms share chemistry, initial drive AND initial rate, and differ only in
    whether the drive depletes.
    """
    net = am_reversible(gamma, k=rate_scale)
    comp = compile_network(net, float(omega))
    pair = reverse_pairing(net)
    from crnl.networks.am_reversible import delta_star
    nb = int(round(omega * gamma / (1.0 + gamma)))
    rest = omega - nb
    sep = int(round(delta_star(gamma) * omega))
    if (rest - sep) % 2:
        sep -= 1
    n0 = np.array([(rest + sep) // 2, rest - (rest + sep) // 2, nb],
                  dtype=np.int64)
    rng = np.random.default_rng(seed)
    times, hit_budget = [], 0
    for _ in range(trials):
        r = gillespie_instrumented(comp, n0, rng, pair,
                                   stop=lambda n: n[0] <= n[1],
                                   max_steps=MAX_STEPS)
        if r.steps >= MAX_STEPS:
            hit_budget += 1
            continue
        times.append(r.t_final)
    if not times:
        return {"arm": "noise", "omega": omega, "gamma": gamma,
                "n_ok": 0, "budget_hits": hit_budget}
    return {
        "arm": "noise", "omega": omega, "gamma": gamma, "trials": trials,
        "n_ok": len(times), "budget_hits": hit_budget,
        "lifetime_median": float(np.median(times)),
        "lifetime_mean": float(np.mean(times)),
        "lifetime_sem": float(np.std(times) / np.sqrt(len(times))),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omegas", type=int, nargs="+", default=[30, 45, 60, 90, 120])
    # The noise control costs exp(Omega) steps and is measurable only at the small
    # end. That asymmetry IS the result, not a sampling shortfall.
    ap.add_argument("--noise-omegas", type=int, nargs="+", default=[30, 36, 42, 48])
    ap.add_argument("--fuel-concs", type=float, nargs="+", default=[10.0, 30.0])
    ap.add_argument("--gamma0", type=float, default=0.30)
    ap.add_argument("--gamma-inf", type=float, default=1.0)
    ap.add_argument("--trials", type=int, default=120)
    ap.add_argument("--seed", type=int, default=20260728)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/fuel_ceiling.json"))
    args = ap.parse_args()

    t0 = time.time()
    print(f"gamma_inf={args.gamma_inf}  gamma_eff(0)={args.gamma0}  "
          f"death at waste/phi={death_waste_fraction(args.gamma_inf):.4f}")

    rows = []
    for fc in args.fuel_concs:
        print(f"\n--- FUELED, fuel concentration {fc} (P3: should be FLAT in Omega)")
        print(f"{'Omega':>6} {'phi':>7} {'lifetime':>11} {'+-sem':>10} "
              f"{'w/phi@loss':>11} {'/death':>8} {'budget':>7}")
        for om in args.omegas:
            r = fueled_lifetime(om, fc, args.gamma0, args.gamma_inf,
                                args.trials, args.seed + om)
            if not r:
                print(f"{om:>6}   all trials hit the step budget")
                continue
            rows.append(r)
            print(f"{om:>6} {r['phi']:>7} {r['lifetime_median']:>11.3f} "
                  f"{r['lifetime_sem']:>10.3f} {r['waste_frac_at_loss']:>11.4f} "
                  f"{r['loss_over_death']:>8.3f} {r['budget_hits']:>7}")

    # One control PER fuel concentration: the arms must share a clock, and the
    # fueled arm's clock is set by f0 (see noise_lifetime's docstring).
    for fc in args.fuel_concs:
        f0 = fc * (1.0 - args.gamma0 / (1.0 + args.gamma0))
        print(f"\n--- NOISE control for fuel conc {fc}: gamma={args.gamma0} held "
              f"forever, rate scaled by f0={f0:.3f} so both arms share a clock")
        print(f"{'Omega':>6} {'lifetime':>11} {'+-sem':>10} {'budget':>7}")
        for om in args.noise_omegas:
            r = noise_lifetime(om, args.gamma0, args.trials,
                               args.seed + 500 + om, rate_scale=f0)
            r["fuel_conc"] = fc
            rows.append(r)
            if not r.get("n_ok"):
                print(f"{om:>6}   exceeds the step budget in every trial "
                      "(lifetime unmeasurably long)")
                continue
            print(f"{om:>6} {r['lifetime_median']:>11.3f} "
                  f"{r['lifetime_sem']:>10.3f} {r['budget_hits']:>7}")

    print("\n=== the two ceilings, per fuel concentration (same clock in each row)")
    for fc in args.fuel_concs:
        ns = [r for r in rows if r["arm"] == "noise" and r.get("fuel_conc") == fc
              and r.get("n_ok", 0) > 5]
        fs = [r for r in rows if r["arm"] == "fueled" and r["fuel_conc"] == fc]
        if len(ns) < 3 or len(fs) < 2:
            continue
        x = np.array([r["omega"] for r in ns], dtype=float)
        y = np.log(np.array([r["lifetime_median"] for r in ns]))
        sl, ic = np.polyfit(x, y, 1)
        res = y - (sl * x + ic)
        lv = float(np.median([r["lifetime_median"] for r in fs]))
        spread = (max(r["lifetime_median"] for r in fs)
                  / min(r["lifetime_median"] for r in fs))
        cross = (np.log(lv) - ic) / sl if sl > 0 else float("nan")
        om_max = max(r["omega"] for r in fs)
        gap = np.exp(sl * om_max + ic) / lv
        print(f"\nfuel conc {fc}:")
        print(f"  noise ceiling  ln(lifetime) = {sl:.4f}*Omega + {ic:.3f}  "
              f"R^2 = {1 - res.var()/y.var():.4f}   (exponential)")
        print(f"  fuel  ceiling  plateau {lv:.2f}, spread {spread:.2f}x over a "
              f"{om_max/min(r['omega'] for r in fs):.0f}x Omega range   (flat)")
        print(f"  they cross at Omega ~ {cross:.0f}; by Omega = {om_max} the noise "
              f"ceiling is {gap:.3g}x further away")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
