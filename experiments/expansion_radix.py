"""Does a bigger alphabet freeze more easily? (n-winner freeze-out)

Run n-winner Approximate Majority from a balanced start (all n committed species
equal, B=0) in an exponentially expanding volume, sweeping the expansion rate H
for several alphabet sizes n. Ties two threads together: the discovery sweep
showed consensus *time* grows with n, so we expect a bigger alphabet to freeze at
a *lower* critical H (longer to decide -> easier to freeze half-made).

Order parameter (winner dominance, normalized so it is comparable across n):

    D = (max_i x_i - 1/n) / (1 - 1/n),   x_i = committed fraction of species i

    D = 1  -> one species won (decided);  D = 0  -> all equal (frozen undecided).

Also reports the mean number of surviving committed species vs H.

    python -m experiments.expansion_radix            # default sweep + figure
    python -m experiments.expansion_radix --quick
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from crnl.networks import n_winner
from crnl.vectorized import compile_network
from crnl.expanding import gillespie_expanding
from crnl.stochastic import seed_for


def run_point(n: int, omega: int, hubble: float, trials: int, base_seed: int) -> dict:
    net = n_winner(n)
    names = list(net.species)
    comp = compile_network(net, omega)
    # balanced start: committed as equal as integers allow, blank = 0
    base = omega // n
    counts = np.full(n, base, dtype=np.int64)
    counts[: omega - base * n] += 1        # spread the remainder (fixed, small bias ~n/omega)
    n0 = np.concatenate([counts, [0]])

    D = np.empty(trials)
    survivors = np.empty(trials)
    resolved = 0
    key = int(round(hubble * 1000))
    for t in range(trials):
        r = gillespie_expanding(comp, n0, seed_for(omega, t, base=base_seed + key),
                                hubble=hubble, species=names)
        committed = r.n_final[:-1]
        tot = int(committed.sum())
        s = int((committed > 0).sum())
        survivors[t] = s
        if tot == 0:
            D[t] = 0.0                      # all-blank: fully undecided
        else:
            maxfrac = committed.max() / tot
            D[t] = (maxfrac - 1.0 / n) / (1.0 - 1.0 / n)
        resolved += (s == 1)
    return {
        "n": n, "omega": int(omega), "hubble": float(hubble), "trials": trials,
        "D": float(D.mean()), "survivors": float(survivors.mean()),
        "p_resolved": resolved / trials,
    }


def _worker(args):
    return run_point(*args)


def run_sweep(ns, omega, hubbles, trials, base_seed, jobs):
    tasks = [(int(n), int(omega), float(h), trials, base_seed)
             for n in ns for h in hubbles]
    if jobs and jobs > 1:
        import multiprocessing as mp
        with mp.Pool(jobs) as pool:
            res = pool.map(_worker, tasks)
    else:
        res = [_worker(t) for t in tasks]
    by_n = {}
    for r in res:
        by_n.setdefault(r["n"], []).append(r)
    for n in by_n:
        by_n[n].sort(key=lambda r: r["hubble"])
    return by_n


def critical_h(rows, level=0.5):
    """H at which the order parameter D crosses `level` (linear interp in log H)."""
    H = np.array([r["hubble"] for r in rows])
    D = np.array([r["D"] for r in rows])
    for i in range(len(H) - 1):
        if D[i] >= level >= D[i + 1]:
            f = (D[i] - level) / (D[i] - D[i + 1]) if D[i] != D[i + 1] else 0.0
            return float(np.exp(np.log(H[i]) + f * (np.log(H[i + 1]) - np.log(H[i]))))
    return float("nan")


def make_figure(by_n, omega, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    ns = sorted(by_n)
    colors = plt.cm.plasma(np.linspace(0, 0.85, len(ns)))

    ax = axes[0]
    for c, n in zip(colors, ns):
        rows = by_n[n]
        ax.semilogx([r["hubble"] for r in rows], [r["D"] for r in rows],
                    "o-", color=c, label=f"n={n}")
    ax.axhline(0.5, color="grey", lw=0.8, ls=":")
    ax.set_xlabel("expansion rate H  (units of k)")
    ax.set_ylabel("winner dominance  D")
    ax.set_title(f"Freeze-out vs alphabet size (Ω={omega})\nD=1 decided, D=0 undecided")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[1]
    for c, n in zip(colors, ns):
        rows = by_n[n]
        ax.semilogx([r["hubble"] for r in rows], [r["survivors"] for r in rows],
                    "s-", color=c, label=f"n={n}")
    ax.set_xlabel("expansion rate H  (units of k)")
    ax.set_ylabel("mean surviving species")
    ax.set_title("Frozen survivor count\n(1 = decided, n = balanced relic)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.25)

    # critical H vs n
    ax = axes[2]
    hstar = [critical_h(by_n[n]) for n in ns]
    ax.loglog(ns, hstar, "D-", color="#d62728")
    ax.set_xlabel("alphabet size n")
    ax.set_ylabel("critical expansion rate H*  (D=0.5)")
    ax.set_title("Bigger alphabet freezes easier?\nH* vs n")
    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote figure -> {out_path}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ns", type=int, nargs="+", default=[2, 3, 4, 6, 8, 12, 16])
    p.add_argument("--omega", type=int, default=160)
    p.add_argument("--hubbles", type=float, nargs="+", default=None)
    p.add_argument("--trials", type=int, default=3000)
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--out", default=os.path.join(
        os.path.dirname(__file__), "expansion_radix.png"))
    args = p.parse_args()

    if args.hubbles is None:
        args.hubbles = [0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.2, 0.35, 0.6, 1.0, 2.0]
    if args.quick:
        args.ns = [2, 4, 8, 16]
        args.hubbles = [0.01, 0.03, 0.06, 0.1, 0.2, 0.4, 1.0]
        args.trials = 1200

    print(f"n-winner freeze-out  (Ω={args.omega}, trials={args.trials}, jobs={args.jobs})")
    by_n = run_sweep(args.ns, args.omega, args.hubbles, args.trials, args.seed, args.jobs)
    for n in sorted(by_n):
        hs = critical_h(by_n[n])
        print(f"-- n={n:>2}  H*(D=0.5)={hs:.3f} --")
        for r in by_n[n]:
            print(f"     H={r['hubble']:6.3f}  D={r['D']:.3f}  "
                  f"survivors={r['survivors']:.2f}  P(resolved)={r['p_resolved']:.3f}")
    make_figure(by_n, args.omega, args.out)


if __name__ == "__main__":
    main()
