"""Why restoration matters: signal survival across a deep cascade.

The founding claim of the whole project is about keeping a signal distinguishable
"across a deep cascade." Here we finally build one. A one-bit signal (sign = the
bit, magnitude = confidence) is passed through D stages. Each stage injects
channel noise; the two cascade types differ only in what the stage then does:

  * restoring (AM): feed the noisy signal in as an AM initial bias, run to
    consensus at finite Omega, and emit the winning rail at FULL magnitude -- the
    analog drift is snapped back to +-1 every stage. Errors are corrected, not
    accumulated.
  * non-restoring: pass the noisy analog value straight through. Channel noise
    accumulates as a random walk; the sign eventually flips.

Measure P(final bit correct) vs depth D. Restoration should hold it near 1 for
many stages while the non-restoring cascade decays toward a coin flip -- the
concrete reason binary/restoring logic enables deep computation.

    python -m experiments.cascade            # sweep + figure
    python -m experiments.cascade --quick
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from crnl.networks import approximate_majority
from crnl.vectorized import compile_network, gillespie_fast
from crnl.stochastic import seed_for


def _am_restore(s, omega, comp, names, rng):
    """One restoring stage: bias AM by signal s in [-1,1], return +-1 (winner)."""
    s = float(np.clip(s, -1.0, 1.0))
    fx = (1.0 + s) / 2.0
    nx = int(round(fx * omega))
    n0 = np.array([nx, omega - nx, 0])
    r = gillespie_fast(comp, n0, rng, species=names)
    return 1.0 if r.n_final[0] > 0 else -1.0     # full-magnitude clean rail


def run_cascade(mode, depth, sigma, omega, trials, base_seed, s_init=0.3):
    """P(correct) at each depth 1..depth. Intended bit is +1 (sign of s_init)."""
    net = approximate_majority()
    names = list(net.species)
    comp = compile_network(net, omega)
    correct = np.zeros(depth, dtype=np.int64)
    for t in range(trials):
        rng = seed_for(omega, t, base=base_seed)
        s = s_init
        for d in range(depth):
            s = s + rng.normal(0.0, sigma)          # channel noise (same for both modes)
            if mode == "restoring":
                s = _am_restore(s, omega, comp, names, rng)
            else:                                    # non-restoring: analog passthrough
                s = float(np.clip(s, -1.0, 1.0))
            correct[d] += (np.sign(s) == 1.0)
    return correct / trials


def _worker(args):
    mode, depth, sigma, omega, trials, seed = args
    return (mode, omega, run_cascade(mode, depth, sigma, omega, trials, seed))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--depth", type=int, default=40)
    p.add_argument("--sigma", type=float, default=0.35, help="per-stage channel noise")
    p.add_argument("--omegas", type=int, nargs="+", default=[20, 40, 80])
    p.add_argument("--trials", type=int, default=4000)
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "cascade.png"))
    args = p.parse_args()

    if args.quick:
        args.depth, args.omegas, args.trials = 25, [20, 60], 1500

    print(f"cascade  (depth={args.depth}, sigma={args.sigma}, trials={args.trials})")
    tasks = [("nonrestoring", args.depth, args.sigma, args.omegas[-1], args.trials, args.seed)]
    tasks += [("restoring", args.depth, args.sigma, w, args.trials, args.seed)
              for w in args.omegas]
    import multiprocessing as mp
    with mp.Pool(args.jobs) as pool:
        results = pool.map(_worker, tasks)

    depths = np.arange(1, args.depth + 1)
    nonrest = next(r for m, w, r in results if m == "nonrestoring")
    rest = {w: r for m, w, r in results if m == "restoring"}

    print(f"{'depth':>6} {'non-rest':>9} " + " ".join(f"AM Ω={w:<4}" for w in args.omegas))
    for i in range(0, args.depth, max(1, args.depth // 12)):
        row = " ".join(f"{rest[w][i]:8.3f}" for w in args.omegas)
        print(f"{depths[i]:>6} {nonrest[i]:>9.3f}  {row}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(depths, nonrest, "o-", color="#7f7f7f", lw=2,
            label="non-restoring (analog passthrough)")
    colors = plt.cm.viridis(np.linspace(0.15, 0.8, len(args.omegas)))
    for c, w in zip(colors, args.omegas):
        ax.plot(depths, rest[w], "o-", color=c, lw=2, label=f"restoring AM (Ω={w})")
    ax.axhline(0.5, color="red", lw=0.8, ls=":", label="coin flip")
    ax.set_xlabel("cascade depth (stages)")
    ax.set_ylabel("P(final bit correct)")
    ax.set_title(f"Signal survival across a deep cascade  (per-stage noise σ={args.sigma})\n"
                 "restoration keeps error bounded; analog drift accumulates")
    ax.set_ylim(0.4, 1.02)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print(f"wrote figure -> {args.out}")


if __name__ == "__main__":
    main()
