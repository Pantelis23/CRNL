"""Is chemical freeze-out a genuine phase transition? Finite-size scaling.

`expansion.py` shows the freeze-out transition sharpening as Omega grows. That is
suggestive but not proof: a crossover also sharpens. The test that distinguishes a
real transition is a *data collapse* -- curves at different system sizes must fall
onto one master curve under the scaling ansatz

    D(H, Omega) = F( (H - Hc) * Omega^a )

with a single critical rate Hc and width exponent a. If that works across a wide
range of Omega, the sharpening is finite-size scaling around a true critical point,
and the finite-Omega crossings obey H*(Omega) = Hc + const * Omega^(-a).

Method: sweep D(H, Omega) on a fine H grid for a geometric ladder of Omega, then
find (Hc, a) minimizing the Bhattacharjee-Seno collapse residual -- each point's
distance from the linear interpolation of its neighbours in the scaled variable.

Headline result (committed figure): Hc ~ 0.055, a ~ 0.38, with six Omega spanning
x32 collapsing cleanly onto one curve. NOTE the exponent is a two-parameter fit
without error bars and sits between 1/3 and 2/5 -- the universality class is NOT
pinned down; see FINDINGS.md.

    python -m experiments.freezeout_scaling --quick
    python -m experiments.freezeout_scaling --omegas 40 80 160 320 640 1280
"""

from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from experiments.expansion import run_point


def _worker(args):
    omega, h, trials, seed = args
    r = run_point(omega, h, trials, seed)
    return {"omega": int(omega), "hubble": float(h), "order": r["order"],
            "p_frozen": r["p_frozen"], "relic": r["relic"]}


def sweep(omegas, hubbles, trials, seed, jobs):
    tasks = [(w, h, trials, seed) for w in omegas for h in hubbles]
    if jobs and jobs > 1:
        import multiprocessing as mp
        with mp.Pool(jobs) as pool:
            res = pool.map(_worker, tasks)
    else:
        res = [_worker(t) for t in tasks]
    by = {}
    for r in res:
        by.setdefault(r["omega"], []).append(r)
    for w in by:
        by[w].sort(key=lambda r: r["hubble"])
    return by


def collapse_residual(by_omega, Hc, a):
    """Bhattacharjee-Seno residual for the ansatz D = F((H-Hc)*Omega^a).

    Pool all (scaled H, D) points, sort by the scaled variable, and score each
    interior point against the linear interpolation of its two neighbours. A good
    collapse means every size's points lie on the same curve, so the residual is
    small. Lower is better.
    """
    us, ds = [], []
    for w, rows in by_omega.items():
        H = np.array([r["hubble"] for r in rows])
        D = np.array([r["order"] for r in rows])
        us.append((H - Hc) * float(w) ** a)
        ds.append(D)
    u = np.concatenate(us)
    d = np.concatenate(ds)
    order = np.argsort(u)
    u, d = u[order], d[order]
    res, n = 0.0, 0
    for i in range(1, len(u) - 1):
        if u[i + 1] == u[i - 1]:
            continue
        w_lo = (u[i] - u[i - 1]) / (u[i + 1] - u[i - 1])
        pred = d[i - 1] + (d[i + 1] - d[i - 1]) * w_lo
        res += (d[i] - pred) ** 2
        n += 1
    return res / max(n, 1)


def fit_collapse(by_omega, hc_range=(0.02, 0.16), a_range=(0.0, 0.9), grid=41):
    """Coarse-then-refined grid search for (Hc, a)."""
    best = None
    for Hc in np.linspace(*hc_range, grid):
        for a in np.linspace(*a_range, grid):
            r = collapse_residual(by_omega, Hc, a)
            if best is None or r < best[0]:
                best = (r, Hc, a)
    r0, Hc0, a0 = best
    for Hc in np.linspace(Hc0 - 0.01, Hc0 + 0.01, grid):
        for a in np.linspace(max(0.0, a0 - 0.1), a0 + 0.1, grid):
            r = collapse_residual(by_omega, Hc, a)
            if r < best[0]:
                best = (r, Hc, a)
    return {"residual": best[0], "Hc": best[1], "a": best[2]}


def crossing(rows, level=0.5):
    """H where the order parameter crosses `level` (log-interpolated)."""
    H = np.array([r["hubble"] for r in rows])
    D = np.array([r["order"] for r in rows])
    for i in range(len(H) - 1):
        if D[i] >= level >= D[i + 1]:
            if D[i] == D[i + 1]:
                return float(H[i])
            f = (D[i] - level) / (D[i] - D[i + 1])
            return float(np.exp(np.log(H[i]) + f * (np.log(H[i + 1]) - np.log(H[i]))))
    return float("nan")


def make_figure(by_omega, fit, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    omegas = sorted(by_omega)
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(omegas)))
    Hc, a = fit["Hc"], fit["a"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    for c, w in zip(colors, omegas):
        rows = by_omega[w]
        ax1.semilogx([r["hubble"] for r in rows], [r["order"] for r in rows],
                     "o-", color=c, ms=4, label=f"Ω={w}")
    ax1.axhline(0.5, color="grey", lw=0.8, ls=":")
    ax1.set_xlabel("expansion rate H")
    ax1.set_ylabel("order parameter D")
    ax1.set_title("Raw freeze-out curves")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.25)

    for c, w in zip(colors, omegas):
        rows = by_omega[w]
        H = np.array([r["hubble"] for r in rows])
        D = np.array([r["order"] for r in rows])
        ax2.plot((H - Hc) * float(w) ** a, D, "o", color=c, ms=4, label=f"Ω={w}")
    ax2.set_xlabel(f"(H − Hc)·Ω^a   [Hc={Hc:.3f}, a={a:.2f}]")
    ax2.set_ylabel("D")
    ax2.set_title("Finite-size-scaling collapse")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote figure -> {out_path}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--omegas", type=int, nargs="+",
                   default=[40, 80, 160, 320, 640, 1280])
    p.add_argument("--h-range", type=float, nargs=2, default=[0.03, 0.4])
    p.add_argument("--h-points", type=int, default=22)
    p.add_argument("--trials", type=int, default=6000)
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--out", default=os.path.join(
        os.path.dirname(__file__), "freezeout_scaling.png"))
    p.add_argument("--data", default=os.path.join(
        os.path.dirname(__file__), os.pardir, "results", "freezeout_fss.json"))
    args = p.parse_args()

    if args.quick:
        args.omegas, args.h_points, args.trials = [40, 80, 160], 10, 800

    hubbles = [round(h, 4) for h in
               np.geomspace(args.h_range[0], args.h_range[1], args.h_points)]
    print(f"freeze-out FSS  ({len(args.omegas)}x{len(hubbles)} points, "
          f"trials={args.trials}, jobs={args.jobs})")
    t0 = time.time()
    by_omega = sweep(args.omegas, hubbles, args.trials, args.seed, args.jobs)

    print("H\\Omega  " + "  ".join(f"{w:>6}" for w in sorted(by_omega)))
    for i, h in enumerate(hubbles):
        row = "  ".join(f"{by_omega[w][i]['order']:6.3f}" for w in sorted(by_omega))
        print(f"{h:7.4f}  {row}")

    fit = fit_collapse(by_omega)
    print(f"\ncollapse: Hc={fit['Hc']:.4f}  a={fit['a']:.3f}  "
          f"residual={fit['residual']:.3e}")
    print("Omega  H*(D=0.5)")
    for w in sorted(by_omega):
        print(f"{w:>5}  {crossing(by_omega[w]):.4f}")

    os.makedirs(os.path.dirname(os.path.abspath(args.data)), exist_ok=True)
    with open(args.data, "w") as fh:
        json.dump({str(w): rows for w, rows in by_omega.items()}, fh, indent=2)
    print(f"wrote data -> {args.data}   (total {time.time()-t0:.0f}s)")
    make_figure(by_omega, fit, args.out)


if __name__ == "__main__":
    main()
