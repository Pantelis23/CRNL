"""How does the restoration barrier scale with alphabet size? c(n) and Omega_req(n).

`radix_wall.py` measures c at one alphabet size on a fixed Omega grid. That is not
enough to get the *scaling*, because the fit band moves: as n grows the barrier
shrinks, so the error-rich window (where the wall is both visible and cleanly
exponential) marches to much larger Omega. A single shared grid saturates at high n.

This experiment therefore adapts per n:

  phase 1  probe a geometric Omega ladder at low trials to LOCATE the band where
           champion-loss sits in ~[0.03, 0.42] -- error-rich but not saturated
  phase 2  run a refined 7-point Omega grid inside that band at high trials
  fit      log(loss) = intercept - c(n)*Omega on the clean band (crossover-aware
           fitter reused from radix_wall)

Outputs c(n) and the population cost Omega_required(n) to hold champion-loss at a
fixed target (default 5%), and writes results/radix_cn.json.

Headline result (committed figure, delta=0.10): c(n) falls ~7x from n=2 to n=32
but SATURATES at a floor ~0.0022 -- confirmed out to n=64. Under a fixed pairwise
margin the radix penalty on the margin is bounded; the unbounded cost is Omega.

    python -m experiments.radix_scaling --quick
    python -m experiments.radix_scaling --ns 2 3 4 6 8 12 16 24 32 --trials 6000
"""

from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from experiments.radix_wall import run_sweep, fit_c


def losses_of(points):
    """(omega, champion-loss fraction) for each swept point."""
    return [(p["omega"], (p["trials"] - p["champion_wins"]) / p["trials"])
            for p in points]


def locate_band(probe_rows, lo_frac=0.03, hi_frac=0.42):
    """Omega values whose champion-loss is error-rich but not saturated.

    Falls back to progressively looser windows so a stubborn n still yields a
    grid rather than crashing; returns the full ladder as a last resort.
    """
    ls = losses_of(probe_rows)
    band = [om for om, l in ls if lo_frac < l < hi_frac]
    if len(band) < 3:
        band = [om for om, l in ls if 0.008 < l < 0.5]
    if len(band) < 2:
        band = [om for om, _ in ls]
    return min(band), max(band)


def run_n(n, delta, ladder, probe_trials, run_trials, base_seed, jobs, grid_points=7):
    """Two-phase measurement of c(n): locate the band, then fit inside it."""
    probe = run_sweep(n, ladder, delta, probe_trials, base_seed, jobs)
    lo, hi = locate_band(probe)
    grid = sorted(set(int(round(x)) for x in np.linspace(lo, hi, grid_points)))
    run = run_sweep(n, grid, delta, run_trials, base_seed + 100, jobs)
    loss_counts = [p["trials"] - p["champion_wins"] for p in run]
    fit = fit_c(grid, loss_counts, run_trials)
    return {
        "n": n,
        "probe_losses": losses_of(probe),
        "grid": grid,
        "run_losses": losses_of(run),
        "fit": fit,
    }


def omega_required(fit, target=0.05):
    """Omega needed to hold champion-loss at `target`, from the fitted wall."""
    if not fit:
        return float("nan")
    return (fit["intercept"] - np.log(target)) / fit["c"]


def make_figure(results, delta, target, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ns = sorted(int(k) for k, v in results.items() if v["fit"])
    cs = [results[str(n)]["fit"]["c"] for n in ns]
    om = [omega_required(results[str(n)]["fit"], target) for n in ns]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.loglog(ns, cs, "o-", color="#1f77b4")
    if len(ns) >= 3:
        p = np.polyfit(np.log(ns), np.log(cs), 1)
        ax1.plot(ns, np.exp(p[1]) * np.array(ns, float) ** p[0], "--",
                 color="#d62728", label=f"c(n) ~ n^{p[0]:.2f}")
        ax1.legend(fontsize=9)
    ax1.set_xlabel("alphabet size n")
    ax1.set_ylabel("barrier c(n)")
    ax1.set_title(f"Radix wall: barrier vs alphabet size (delta={delta})")
    ax1.grid(True, which="both", alpha=0.3)

    ax2.loglog(ns, om, "s-", color="#d62728")
    ax2.set_xlabel("alphabet size n")
    ax2.set_ylabel(f"Omega to hold champion-loss <= {target:.0%}")
    ax2.set_title("Population cost of radix")
    ax2.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote figure -> {out_path}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ns", type=int, nargs="+",
                   default=[2, 3, 4, 6, 8, 12, 16, 24, 32])
    p.add_argument("--delta", type=float, default=0.10,
                   help="fixed pairwise margin (0.10 = 55/45 at n=2)")
    p.add_argument("--omega-cap", type=int, default=1100,
                   help="largest Omega to probe; bounds per-point cost")
    p.add_argument("--probe-trials", type=int, default=1200)
    p.add_argument("--trials", type=int, default=6000)
    p.add_argument("--target", type=float, default=0.05,
                   help="champion-loss target for Omega_required(n)")
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--out", default=os.path.join(
        os.path.dirname(__file__), "radix_scaling.png"))
    p.add_argument("--data", default=os.path.join(
        os.path.dirname(__file__), os.pardir, "results", "radix_cn.json"))
    args = p.parse_args()

    if args.quick:
        args.ns = [2, 4, 8]
        args.probe_trials, args.trials, args.omega_cap = 400, 1500, 600

    ladder = [w for w in
              (20, 30, 45, 65, 95, 140, 200, 290, 420, 600, 850, 1100)
              if w <= args.omega_cap]

    print(f"radix scaling  (delta={args.delta}, trials={args.trials}, "
          f"jobs={args.jobs})")
    t0 = time.time()
    results = {}
    for n in args.ns:
        tn = time.time()
        res = run_n(n, args.delta, ladder, args.probe_trials, args.trials,
                    args.seed, args.jobs)
        results[str(n)] = res
        fit = res["fit"]
        if fit:
            print(f"  n={n:>3}  c={fit['c']:.5g}  R2={fit['r2']:.3f}  "
                  f"band=[{fit['omega_lo']},{fit['omega_hi']}]  "
                  f"Om_req={omega_required(fit, args.target):.0f}  "
                  f"({time.time()-tn:.0f}s)", flush=True)
        else:
            print(f"  n={n:>3}  (no clean band)  ({time.time()-tn:.0f}s)",
                  flush=True)

    os.makedirs(os.path.dirname(os.path.abspath(args.data)), exist_ok=True)
    with open(args.data, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"wrote data -> {args.data}   (total {time.time()-t0:.0f}s)")
    make_figure(results, args.delta, args.target, args.out)


if __name__ == "__main__":
    main()
