"""Chemical freeze-out: does expansion prevent the decision? (crnl/expanding.py)

Run Approximate Majority from a balanced 50/50 start in an exponentially
expanding volume Omega(t) = Omega0 * exp(H t), sweeping the expansion rate H.
Below a critical H the consensus completes (a clean rail); above it the reaction
freezes mid-decision, locking in a relic minority abundance -- the chemical
analogue of cosmological freeze-out.

Observables vs H:
  * order parameter  <|n_X - n_Y| / (n_X + n_Y)>   (1 = decided, 0 = undecided)
  * relic abundance  <min(n_X, n_Y) / Omega>       (the frozen-in minority)
  * P(frozen), P(fully resolved)
Run at several Omega0 to expose the finite-size sharpening of the transition.

    python -m experiments.expansion            # default sweep + figure
    python -m experiments.expansion --quick
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from crnl.networks import approximate_majority
from crnl.vectorized import compile_network
from crnl.expanding import gillespie_expanding
from crnl.stochastic import seed_for


def run_point(omega: int, hubble: float, trials: int, base_seed: int) -> dict:
    net = approximate_majority()
    names = list(net.species)
    comp = compile_network(net, omega)
    n0 = np.array([omega // 2, omega - omega // 2, 0])
    orders = np.empty(trials)
    relics = np.empty(trials)
    xfrac = np.empty(trials)
    frozen = 0
    resolved = 0  # minority exactly 0 (clean decision)
    key = int(round(hubble * 1000))
    for t in range(trials):
        r = gillespie_expanding(comp, n0, seed_for(omega, t, base=base_seed + key),
                                hubble=hubble, species=names)
        x, y = int(r.n_final[0]), int(r.n_final[1])
        tot = x + y
        orders[t] = abs(x - y) / tot if tot else 1.0
        relics[t] = min(x, y) / omega
        xfrac[t] = x / tot if tot else 0.5
        frozen += r.frozen
        resolved += (min(x, y) == 0)
    return {
        "omega": int(omega), "hubble": float(hubble), "trials": trials,
        "order": float(orders.mean()), "order_sd": float(orders.std()),
        "relic": float(relics.mean()),
        "p_frozen": frozen / trials, "p_resolved": resolved / trials,
        "xfrac": xfrac,
    }


def _worker(args):
    return run_point(*args)


def run_sweep(omegas, hubbles, trials, base_seed, jobs):
    tasks = [(int(w), float(h), trials, base_seed) for w in omegas for h in hubbles]
    if jobs and jobs > 1:
        import multiprocessing as mp
        with mp.Pool(jobs) as pool:
            res = pool.map(_worker, tasks)
    else:
        res = [_worker(t) for t in tasks]
    by_omega = {}
    for r in res:
        by_omega.setdefault(r["omega"], []).append(r)
    for w in by_omega:
        by_omega[w].sort(key=lambda r: r["hubble"])
    return by_omega


def make_figure(by_omega, hist_omega, hist_hubbles, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(by_omega)))

    # panel 1: order parameter vs H (the phase transition)
    ax = axes[0]
    for c, (omega, rows) in zip(colors, sorted(by_omega.items())):
        H = [r["hubble"] for r in rows]
        o = [r["order"] for r in rows]
        ax.semilogx(H, o, "o-", color=c, label=f"Ω={omega}")
    ax.set_xlabel("expansion rate H  (units of k)")
    ax.set_ylabel("order parameter  ⟨|X−Y|/(X+Y)⟩")
    ax.set_title("Freeze-out transition\n(decided → frozen undecided)")
    ax.axhline(0.5, color="grey", lw=0.8, ls=":")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.25)

    # panel 2: relic minority abundance vs H
    ax = axes[1]
    for c, (omega, rows) in zip(colors, sorted(by_omega.items())):
        H = [r["hubble"] for r in rows]
        rel = [r["relic"] for r in rows]
        ax.semilogx(H, rel, "s-", color=c, label=f"Ω={omega}")
    ax.set_xlabel("expansion rate H  (units of k)")
    ax.set_ylabel("relic minority  ⟨min(X,Y)/Ω⟩")
    ax.set_title("Frozen-in relic abundance\n(0 = clean decision, 0.5 = fully undecided)")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.25)

    # panel 3: frozen X-fraction distribution at 3 H (one Omega): bimodal -> central
    ax = axes[2]
    rows = {r["hubble"]: r for r in by_omega[hist_omega]}
    bins = np.linspace(0, 1, 26)
    for h in hist_hubbles:
        # nearest available hubble
        hh = min(rows, key=lambda x: abs(x - h))
        ax.hist(rows[hh]["xfrac"], bins=bins, histtype="step", lw=2,
                density=True, label=f"H={hh:g}")
    ax.set_xlabel("frozen X-fraction  X/(X+Y)")
    ax.set_ylabel("density")
    ax.set_title(f"Frozen composition (Ω={hist_omega})\nslow: rails; fast: central relic")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote figure -> {out_path}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--omegas", type=int, nargs="+", default=[40, 80, 200, 400])
    p.add_argument("--hubbles", type=float, nargs="+", default=None)
    p.add_argument("--trials", type=int, default=4000)
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--out", default=os.path.join(
        os.path.dirname(__file__), "expansion.png"))
    args = p.parse_args()

    if args.hubbles is None:
        args.hubbles = [0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.22,
                        0.32, 0.5, 0.75, 1.0, 2.0, 5.0]
    if args.quick:
        args.omegas = [40, 160]
        args.hubbles = [0.01, 0.03, 0.07, 0.1, 0.15, 0.3, 0.7, 2.0]
        args.trials = 1500

    print(f"expansion freeze-out  (trials={args.trials}, jobs={args.jobs})")
    by_omega = run_sweep(args.omegas, args.hubbles, args.trials, args.seed, args.jobs)

    for omega, rows in sorted(by_omega.items()):
        print(f"-- Omega={omega} --")
        for r in rows:
            print(f"   H={r['hubble']:6.3f}  order={r['order']:.3f}  "
                  f"relic={r['relic']:.4f}  P(frozen)={r['p_frozen']:.2f}  "
                  f"P(resolved)={r['p_resolved']:.3f}")

    hist_omega = args.omegas[-1]
    make_figure(by_omega, hist_omega, [0.03, 0.15, 1.0], args.out)


if __name__ == "__main__":
    main()
