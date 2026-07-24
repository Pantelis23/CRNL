"""Symmetric-start discovery: characterize outcome distribution + dynamics.

spec 2026-07-24 §4.2. All n committed species as equal as integers allow (B=0);
the C mod n remainder is placed on RANDOM distinct species per trial so no single
label is systematically favored (eps_eff ~ n/C -> negligible for C >> n). No
exponent is fit (symmetry admits no barrier); we report the outcome distribution
{single winner / blank / coexist / undecided} and consensus time, under two
resource conventions:
  * fixed total Omega ("same beaker, more symbol types"), and
  * fixed per-species density Omega = n*m ("fair resources per symbol").

    python -m experiments.radix_discovery --quick
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from crnl.networks import n_winner
from crnl.vectorized import compile_network, gillespie_fast
from crnl.classify import classify_winner
from crnl.stochastic import seed_for


def symmetric_counts(n: int, omega: int, rng: np.random.Generator) -> np.ndarray:
    base = omega // n
    counts = np.full(n, base, dtype=np.int64)
    rem = int(omega - base * n)
    if rem:
        extra = rng.choice(n, size=rem, replace=False)
        counts[extra] += 1
    return np.concatenate([counts, [0]])


def run_point(n, omega, trials, base_seed, max_steps):
    net = n_winner(n)
    names = list(net.species)
    compiled = compile_network(net, omega)
    bins = {"single": 0, "blank": 0, "coexist": 0, "undecided": 0}
    times = []
    for t in range(trials):
        rng = seed_for(omega, t, base=base_seed + n)
        n0 = symmetric_counts(n, omega, rng)
        res = gillespie_fast(compiled, n0, rng, max_steps=max_steps, species=names)
        w = classify_winner(res, blank="B")
        if w == "blank":
            bins["blank"] += 1
        elif w == "coexist":
            bins["coexist"] += 1
        elif w == "undecided":
            bins["undecided"] += 1
        else:
            bins["single"] += 1
        if res.absorbed:
            times.append(res.t_final)
    med = float(np.median(times)) if times else float("nan")
    return {"n": n, "omega": int(omega), "trials": trials,
            **bins, "median_consensus_time": med}


def _worker(args):
    return run_point(*args)


def run_sweep(points, trials, base_seed, jobs, max_steps):
    tasks = [(n, int(w), trials, base_seed, max_steps) for (n, w) in points]
    if jobs and jobs > 1:
        import multiprocessing as mp
        with mp.Pool(jobs) as pool:
            return pool.map(_worker, tasks)
    return [_worker(t) for t in tasks]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ns", type=int, nargs="+",
                   default=[2, 3, 4, 6, 8, 12, 16, 24, 32])
    p.add_argument("--omega-total", type=int, default=120,
                   help="fixed total Omega for the 'same beaker' sweep")
    p.add_argument("--density", type=int, default=20,
                   help="per-species m for the fixed-density sweep (Omega = n*m)")
    p.add_argument("--trials", type=int, default=4000)
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--max-steps", type=int, default=5_000_000)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--out", default=os.path.join(
        os.path.dirname(__file__), "radix_discovery.png"))
    args = p.parse_args()

    ns = [2, 3, 4, 6, 8] if args.quick else args.ns
    trials = 1000 if args.quick else args.trials

    print("== fixed total Omega =", args.omega_total, "==")
    fixed = run_sweep([(n, args.omega_total) for n in ns],
                      trials, args.seed, args.jobs, args.max_steps)
    for r in fixed:
        _print_row(r)

    print("== fixed density Omega = n *", args.density, "==")
    dens = run_sweep([(n, n * args.density) for n in ns],
                     trials, args.seed + 1000, args.jobs, args.max_steps)
    for r in dens:
        _print_row(r)

    _make_figure(ns, fixed, dens, args.out)


def _print_row(r):
    tr = r["trials"]
    print(f"  n={r['n']:>3} Om={r['omega']:>4}  "
          f"single={r['single'] / tr:.3f} blank={r['blank'] / tr:.3f} "
          f"coexist={r['coexist'] / tr:.3f} undecided={r['undecided'] / tr:.3f}  "
          f"t_med={r['median_consensus_time']:.2f}")


def _make_figure(ns, fixed, dens, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def frac(rows, key):
        return [r[key] / r["trials"] for r in rows]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for rows, label, ax in ((fixed, "fixed total Omega", axes[0]),
                            (dens, "fixed density", axes[1])):
        ax.plot(ns, frac(rows, "single"), "o-", label="single winner")
        ax.plot(ns, frac(rows, "blank"), "D-", label="all-blank")
        ax.plot(ns, frac(rows, "coexist"), "s-", label="coexist")
        ax.plot(ns, frac(rows, "undecided"), "^-", label="undecided")
        ax.set_xlabel("alphabet size n"); ax.set_ylabel("outcome fraction")
        ax.set_title(label); ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
    axes[2].plot(ns, [r["median_consensus_time"] for r in fixed], "o-",
                 label="fixed total")
    axes[2].plot(ns, [r["median_consensus_time"] for r in dens], "s-",
                 label="fixed density")
    axes[2].set_xlabel("alphabet size n")
    axes[2].set_ylabel("median consensus time")
    axes[2].set_title("Consensus time vs alphabet size")
    axes[2].grid(True, alpha=0.3); axes[2].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out_path, dpi=130)
    print(f"wrote figure -> {out_path}")


if __name__ == "__main__":
    main()
