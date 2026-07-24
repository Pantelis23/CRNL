"""Champion-vs-field: measure the barrier c(n) and population cost Omega_req(n).

spec 2026-07-24 §4.1. One champion (X1) leads EACH rival by a fixed pairwise
margin delta (anchored to 55/45 at n=2, delta=0.10). Error = champion does not
win. As n grows at fixed delta, basins crowd and c(n) falls -- the cost of radix.

    python -m experiments.radix_wall --quick
    python -m experiments.radix_wall --trials 20000 --jobs 8
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from crnl.networks import n_winner
from crnl.vectorized import compile_network, gillespie_fast
from crnl.classify import classify_winner
from crnl.stochastic import seed_for


def champion_counts(n: int, omega: int, delta: float) -> np.ndarray:
    """Integer counts (X1..Xn, B=0) with champion X1 strictly leading each rival.

    Fixed pairwise margin: f_champ = 1/n + delta*(n-1)/n, f_other = 1/n - delta/n.
    Remainder from rounding is handed to the largest fractional parts. When
    delta*omega is small the fractional remainder alone can leave the champion
    merely TIED with a rival (e.g. n=4, omega=11, delta=0.05 -> [3,3,3,2]); a
    final guard moves one molecule from a leading rival to the champion so the
    champion strictly leads at t=0 for any delta>0 (as long as omega has a
    molecule to move). On the default experiment grid no tie ever occurs and the
    guard is inert.
    """
    f = np.full(n, 1.0 / n - delta / n)
    f[0] = 1.0 / n + delta * (n - 1) / n
    exact = f * omega
    counts = np.floor(exact).astype(np.int64)
    rem = int(omega - counts.sum())
    frac = exact - counts
    order = np.argsort(-frac)
    for i in range(rem):
        counts[order[i]] += 1
    # strict-lead guard (fable-5 review): floor is monotonic so counts[0] is the
    # max floor and rivals can at most TIE it; if so, move one unit to the champion.
    rivals_max = int(counts[1:].max())
    if rivals_max >= 1 and counts[0] <= rivals_max:
        donor = 1 + int(np.argmax(counts[1:]))
        counts[donor] -= 1
        counts[0] += 1
    return np.concatenate([counts, [0]])


def run_point(n, omega, delta, trials, base_seed):
    net = n_winner(n)
    names = list(net.species)
    compiled = compile_network(net, omega)
    n0 = champion_counts(n, omega, delta)
    champ = "X1"
    champion_wins = other_wins = blank = coexist = 0
    for t in range(trials):
        res = gillespie_fast(compiled, n0, seed_for(omega, t, base=base_seed + n),
                             species=names)
        w = classify_winner(res, blank="B")
        if w == champ:
            champion_wins += 1
        elif w == "blank":
            blank += 1
        elif w == "coexist" or w == "undecided":
            coexist += 1
        else:
            other_wins += 1
    return {"n": n, "omega": int(omega), "trials": trials,
            "champion_wins": champion_wins, "other_wins": other_wins,
            "blank": blank, "coexist": coexist}


def _champ_loss_counts(point):
    return point["trials"] - point["champion_wins"]


def fit_c(omegas, champ_loss_counts, trials, min_events=15, max_frac=0.4):
    """Log-linear fit log(loss_frac) = intercept - c*Omega on the clean band.

    Keeps points that are error-rich enough (>= min_events losses) and not
    saturated (loss fraction < max_frac). Drops are the crossover/flat regions.
    """
    omegas = np.asarray(omegas, dtype=float)
    counts = np.asarray(champ_loss_counts, dtype=float)
    fr = counts / trials
    keep = (counts >= min_events) & (fr > 0) & (fr < max_frac)
    if keep.sum() < 2:
        return None
    xs = omegas[keep]
    ys = np.log(fr[keep])
    slope, intercept = np.polyfit(xs, ys, 1)
    pred = slope * xs + intercept
    ss_res = np.sum((ys - pred) ** 2)
    ss_tot = np.sum((ys - ys.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {"c": -slope, "intercept": intercept, "r2": r2,
            "omega_lo": int(xs.min()), "omega_hi": int(xs.max()),
            "n_points": int(keep.sum())}


def _worker(args):
    return run_point(*args)


def run_sweep(n, omegas, delta, trials, base_seed, jobs):
    tasks = [(n, int(w), delta, trials, base_seed) for w in omegas]
    if jobs and jobs > 1:
        import multiprocessing as mp
        with mp.Pool(jobs) as pool:
            return pool.map(_worker, tasks)
    return [_worker(t) for t in tasks]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--delta", type=float, default=0.10)
    p.add_argument("--ns", type=int, nargs="+",
                   default=[2, 3, 4, 6, 8, 12, 16, 24, 32])
    p.add_argument("--omegas", type=int, nargs="+",
                   default=[20, 40, 60, 80, 100, 140, 180, 240, 300])
    p.add_argument("--trials", type=int, default=8000)
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--out", default=os.path.join(
        os.path.dirname(__file__), "radix_wall.png"))
    args = p.parse_args()

    ns = [2, 3, 4, 6, 8] if args.quick else args.ns
    omegas = [20, 40, 60, 80, 100] if args.quick else args.omegas
    trials = 2000 if args.quick else args.trials

    print(f"radix wall (delta={args.delta}, trials={trials}, jobs={args.jobs})")
    c_of_n = {}
    for n in ns:
        pts = run_sweep(n, omegas, args.delta, trials, args.seed, args.jobs)
        loss = [_champ_loss_counts(pt) for pt in pts]
        fit = fit_c(omegas, loss, trials)
        cval = fit["c"] if fit else float("nan")
        c_of_n[n] = (cval, fit)
        print(f"  n={n:>3}  c(n)={cval:.4g}  "
              + ("R2=%.3f band=[%d,%d]" % (fit["r2"], fit["omega_lo"], fit["omega_hi"])
                 if fit else "(no clean band)"))
    _make_figure(c_of_n, args.delta, args.out)


def _make_figure(c_of_n, delta, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ns = sorted(n for n, (c, _) in c_of_n.items() if np.isfinite(c))
    cs = [c_of_n[n][0] for n in ns]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.loglog(ns, cs, "o-", color="#1f77b4")
    ax1.set_xlabel("alphabet size n"); ax1.set_ylabel("barrier c(n)")
    ax1.set_title(f"Radix wall: barrier vs alphabet size (delta={delta})")
    ax1.grid(True, which="both", alpha=0.3)
    # Omega to hold champion-loss at 5%: Omega* = (intercept - log 0.05)/c
    target = 0.05
    omreq = []
    for n in ns:
        _, fit = c_of_n[n]
        omreq.append((fit["intercept"] - np.log(target)) / fit["c"])
    ax2.loglog(ns, omreq, "s-", color="#d62728")
    ax2.set_xlabel("alphabet size n")
    ax2.set_ylabel("Omega to hold champion-loss <= 5%")
    ax2.set_title("Population cost of radix")
    ax2.grid(True, which="both", alpha=0.3)
    fig.tight_layout(); fig.savefig(out_path, dpi=130)
    print(f"wrote figure -> {out_path}")


if __name__ == "__main__":
    main()
