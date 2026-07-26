"""Is there a critical expansion rate? Extending FINDINGS Sec.5 by 128x in Omega.

Sec.5 fitted `D(H,Omega) = F((H-Hc)*Omega^a)` over Omega = 40..1280 and reported
Hc ~ 0.055, a ~ 0.38 -- "a genuine transition". `crnl/freezeout.py` shows the
expanding SSA is *exactly* ordinary SSA stopped at internal time tau = 1/H, so
D(H,Omega) is nothing but the AM consensus-time distribution read at 1/H. From a
symmetric start that time diverges logarithmically:

    d = x - y obeys d' = b*d EXACTLY, b -> 1/3, and d(0) ~ Omega^{-1/2} (shot
    noise), so  tau*(Omega) = (1/lambda) * (1/2) * ln Omega + O(1)
                            = (3/2) ln Omega + O(1),   lambda = 1/3.

Prediction, with no fitted parameter: **1/H*(Omega) = (3/2) ln Omega + B**, hence
**Hc = 0** -- no critical rate at all, and in the Omega -> infinity limit ANY
expansion freezes the decision. Sec.5's Hc > 0 should be a degeneracy of its 32x
range.

Three things are measured here, all from the same runs:

  1. tau*(Omega) at four crossing levels, out to Omega = 163840, with error bars
     from independent replicates. Log law (2 params) vs Hc + C Omega^-a (3).
  2. Curvature of tau* in ln Omega -- a positive Hc *must* bend it over toward
     the ceiling 1/Hc, so the curvature bounds Hc from above.
  3. A zero-parameter collapse D = F(tau - (3/2) ln Omega) scored against
     Sec.5's two-parameter FSS collapse by the same Bhattacharjee-Seno residual.

Bonus check with a different predicted number: the median *absorption* time
should scale as (5/2) ln Omega, because clearing the last molecules off the rail
adds its own ln Omega at unit rate on top of the (3/2) ln Omega of symmetry
breaking.

    python -m experiments.freezeout_law --quick
    python -m experiments.freezeout_law                 # ~1 h
"""

from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from crnl.freezeout import (
    am_observables,
    am_order_exact,
    crossing_tau,
    deterministic_times,
    internal_clock_sweep,
    n_winner_observables,
)
from crnl.networks import approximate_majority, n_winner
from crnl.vectorized import compile_network

LEVELS = (0.25, 0.5, 0.75, 0.9)


def _start(omega: int, n_sym: int, bias: float = 0.0) -> list:
    """Start state. bias = 0 is symmetric; bias > 0 gives species 1 a head start.

    The distinction is the whole point of the biased control. From an EXACTLY
    symmetric start the seed of the decision is shot noise, size Omega^{-1/2}, so
    the decision time grows like ln Omega and H* -> 0. From a start with an
    Omega-INDEPENDENT margin the seed does not shrink, so the decision time is
    Omega-independent and H* tends to a positive constant. Same instrument, same
    network, opposite conclusions -- which is how you tell a critical point from a
    shrinking initial condition.

    `bias` is the pairwise margin (x1 - x_j)/(sum of committed) at n_sym = 2.
    """
    net = n_winner(n_sym) if n_sym > 2 else approximate_majority()
    n0 = np.zeros(len(net.species), dtype=np.int64)
    if not 0.0 <= bias < 1.0:
        raise ValueError(f"bias must be in [0, 1), got {bias}")
    lead = int(round(omega * (1.0 / n_sym + bias * (n_sym - 1) / n_sym)))
    lead = min(max(lead, 1), omega - (n_sym - 1))
    base, rem = divmod(omega - lead, n_sym - 1)
    n0[0] = lead
    n0[1:n_sym] = base
    n0[1:1 + rem] += 1
    return net, n0.tolist()


def run_rep(omega: int, taus, trials: int, rep: int, seed: int,
            n_sym: int = 2, bias: float = 0.0) -> dict:
    """One independent replicate at one Omega. The unit of parallelism."""
    net, n0 = _start(omega, n_sym, bias)
    comp = compile_network(net, omega)
    obs = am_observables if n_sym == 2 else n_winner_observables(n_sym)
    rng = np.random.default_rng([seed, omega, rep, n_sym, int(bias * 1e6)])
    out = internal_clock_sweep(comp, n0, taus, trials, rng, obs)
    return {"omega": int(omega), "rep": rep,
            "order": out["means"][0], "second": out["means"][1],
            "p_absorbed": out["p_absorbed"]}


def _worker(a):
    return run_rep(*a)


def combine(omega: int, reps_out: list, taus, trials: int, n_sym: int) -> dict:
    """Pool replicates; the scatter ACROSS replicates is the error bar.

    Replicate scatter rather than a within-run SEM, because tau* is a nonlinear
    functional of the whole curve and its grid points are correlated.
    """
    curves = np.array([r["order"] for r in reps_out])
    second = np.array([r["second"] for r in reps_out])
    absorbed = np.array([r["p_absorbed"] for r in reps_out])
    reps = len(reps_out)
    taustars = np.array([[crossing_tau(taus, c, L) for L in LEVELS]
                         for c in curves])
    tmeds = np.array([crossing_tau(taus, a, 0.5) for a in absorbed])
    return {
        "omega": int(omega), "n_sym": int(n_sym), "trials": trials, "reps": reps,
        "order": curves.mean(axis=0).tolist(),
        "order_sem": (curves.std(axis=0, ddof=1) / np.sqrt(reps)).tolist(),
        "p_absorbed": absorbed.mean(axis=0).tolist(),
        # AM: relic min(X,Y)/Omega.  n-winner: surviving committed species count,
        # which is what distinguishes a collective decision (all n still alive at
        # the crossing) from an extinction-driven one (several already gone).
        "second": second.mean(axis=0).tolist(),
        "tau_star": np.nanmean(taustars, axis=0).tolist(),
        "tau_star_sem": (np.nanstd(taustars, axis=0, ddof=1)
                         / np.sqrt(reps)).tolist(),
        "tau_absorb_median": float(np.nanmean(tmeds)),
        "tau_absorb_median_sem": float(np.nanstd(tmeds, ddof=1) / np.sqrt(reps)),
    }


# --------------------------------------------------------------------------- #
# analysis


def fit_log_law(omegas, taus, sems):
    """tau* = A ln Omega + B, weighted; A is the prediction under test (3/2)."""
    x = np.log(np.asarray(omegas, dtype=float))
    y = np.asarray(taus, dtype=float)
    w = 1.0 / np.maximum(np.asarray(sems, dtype=float), 1e-9) ** 2
    X = np.vstack([x, np.ones_like(x)]).T
    W = np.diag(w)
    cov = np.linalg.inv(X.T @ W @ X)
    beta = cov @ (X.T @ W @ y)
    resid = y - X @ beta
    chi2 = float(resid @ (w * resid))
    return {"A": float(beta[0]), "B": float(beta[1]),
            "A_sd": float(np.sqrt(cov[0, 0])), "B_sd": float(np.sqrt(cov[1, 1])),
            "chi2": chi2, "dof": len(x) - 2,
            "rms": float(np.sqrt(np.mean(resid ** 2)))}


def fit_power_law(omegas, taus, sems):
    """H* = Hc + C Omega^-a, fitted on tau* = 1/H* with the same weights."""
    from scipy.optimize import least_squares

    w = np.asarray(omegas, dtype=float)
    y = np.asarray(taus, dtype=float)
    s = np.maximum(np.asarray(sems, dtype=float), 1e-9)

    def resid(p):
        Hc, C, a = p
        return (1.0 / (Hc + C * w ** (-a)) - y) / s

    best = None
    for a0 in (0.15, 0.3, 0.45, 0.7):
        for hc0 in (1e-4, 0.01, 0.04, 0.06):
            try:
                r = least_squares(resid, [hc0, 0.5, a0],
                                  bounds=([0, 1e-6, 0.01], [1, 1e3, 4]))
            except Exception:
                continue
            if best is None or r.cost < best.cost:
                best = r
    Hc, C, a = best.x
    return {"Hc": float(Hc), "C": float(C), "a": float(a),
            "chi2": float(2 * best.cost), "dof": len(w) - 3,
            "ceiling_tau": float(1.0 / Hc) if Hc > 0 else float("inf")}


def curvature_test(omegas, taus, sems):
    """Quadratic term in ln Omega. A positive Hc forces tau* to bend DOWN."""
    x = np.log(np.asarray(omegas, dtype=float))
    x = x - x.mean()
    y = np.asarray(taus, dtype=float)
    w = 1.0 / np.maximum(np.asarray(sems, dtype=float), 1e-9) ** 2
    X = np.vstack([x ** 2, x, np.ones_like(x)]).T
    W = np.diag(w)
    cov = np.linalg.inv(X.T @ W @ X)
    beta = cov @ (X.T @ W @ y)
    return {"quad": float(beta[0]), "quad_sd": float(np.sqrt(cov[0, 0])),
            "sigma": float(beta[0] / np.sqrt(cov[0, 0]))}


def exact_route(omegas, tau_max: float, n_grid: int = 401) -> dict:
    """Second, independent route: solve the CME instead of sampling it.

    dp/dtau = Q^T p on the conserved simplex, so tau*(Omega) carries NO sampling
    error at all. Costs ~Omega^3 (states x matvecs), which caps it well below the
    SSA's reach -- the point is to confirm the SSA's tau* and the slope where both
    instruments overlap, not to extend the range.
    """
    out = []
    for w in omegas:
        t0 = time.time()
        # cost is ~Omega^3 * tau_max, so integrate only as far as the crossings
        # actually reach: 1.5 ln Omega is the prediction, +10 covers D = 0.9.
        span = min(tau_max, 1.5 * np.log(w) + 10.0)
        taus, curve = am_order_exact(int(w), span, n_grid)
        ts = [crossing_tau(taus[1:], curve[1:], L) for L in LEVELS]
        out.append({"omega": int(w), "tau_star": ts,
                    "seconds": round(time.time() - t0, 2)})
        print(f"  exact Ω={w:6d}  τ*(0.5) = {ts[LEVELS.index(0.5)]:7.4f}"
              f"   ({out[-1]['seconds']:.1f}s)")
    ws = [r["omega"] for r in out]
    fits = {}
    for i, L in enumerate(LEVELS):
        y = [r["tau_star"][i] for r in out]
        fits[str(L)] = fit_log_law(ws, y, [1e-6] * len(ws))
        print(f"  exact D={L}: A = {fits[str(L)]['A']:.4f}"
              f"  rms {fits[str(L)]['rms']:.4f}")
    return {"rows": out, "fits": fits}


#: FINDINGS Sec.6's measured H*(n) at Omega=160, D=0.5, from the EXPANDING sweep.
SEC6_HSTAR = {2: 0.121, 3: 0.110, 4: 0.101, 6: 0.089, 8: 0.079, 12: 0.074,
              16: 0.071}


def sec6_check(omega: int = 160, trials: int = 2000, reps: int = 4,
               tau_max: float = 60.0, seed: int = 7) -> dict:
    """Does the reduction reproduce Sec.6's own H*(n), with no expansion at all?

    Sec.6 measured H*(n) by running the *expanding* SSA at many H. If the time
    change is right, `1/H*(n)` must simply be the ORDINARY n-winner consensus
    time at the same Omega and level -- so measuring that and inverting it is a
    check of the reduction against a published table it never saw.
    """
    taus = np.geomspace(1.0, tau_max, 300)
    out = []
    for n, h6 in SEC6_HSTAR.items():
        net, n0 = _start(omega, n)
        comp = compile_network(net, omega)
        obs = am_observables if n == 2 else n_winner_observables(n)
        ts = [crossing_tau(taus, internal_clock_sweep(
                  comp, n0, taus, trials,
                  np.random.default_rng([seed, n, rep]), obs)["means"][0], 0.5)
              for rep in range(reps)]
        t = float(np.mean(ts))
        e = float(np.std(ts, ddof=1) / np.sqrt(reps))
        out.append({"n": n, "sec6_hstar": h6, "tau_star": t, "tau_star_sem": e,
                    "hstar": 1.0 / t, "ratio": t * h6})
    return {"omega": omega, "rows": out}


def extinction_check(n_sym: int, spec=None, trials: int = 2000,
                     seed: int = 31) -> dict:
    """Why the n-winner slope falls short of (2n-1)/2 at small Omega.

    The collective route needs all n committed species alive while the slow
    symmetry-breaking mode (rate lambda = 1/(2n-1)) grows. A species that
    fluctuates to zero is gone for good -- X_i is only regenerated by
    B + X_i -> 2X_i, which needs X_i > 0 -- and losing species decides the contest
    faster than lambda allows. So: how many are still alive when D crosses 0.5?
    Collective => ~n. Extinction-driven => well below n.
    """
    if spec is None:
        spec = [(60, 60), (120, 60), (240, 70), (480, 80), (960, 90),
                (3840, 110), (15360, 130), (46080, 150)]
    out = []
    for omega, tau_max in spec:
        net, n0 = _start(omega, n_sym)
        comp = compile_network(net, omega)
        taus = np.geomspace(1.0, tau_max, 300)
        r = internal_clock_sweep(comp, n0, taus, trials,
                                 np.random.default_rng([seed, n_sym, omega]),
                                 n_winner_observables(n_sym))
        ts = crossing_tau(taus, r["means"][0], 0.5)
        out.append({"omega": omega, "per_species": omega / (2 * n_sym - 1),
                    "tau_star": ts,
                    "alive_at_crossing": float(np.interp(ts, taus,
                                                         r["means"][1]))})
    return {"n_sym": n_sym, "trials": trials, "rows": out}


def fixed_H_table(rows, taus, hs=(0.1219, 0.0879, 0.055, 0.04)) -> dict:
    """The model-free half of the argument: D at FIXED H, as Omega grows.

    No fit is involved. If Hc > 0 then for H < Hc the system must still decide as
    Omega -> infinity, so D(H, Omega) must stay near 1. But at fixed internal time
    tau = 1/H the chain converges (Kurtz) to the mass-action ODE, which from an
    exactly symmetric start stays exactly symmetric -- so D(H, Omega) -> 0 for
    EVERY H > 0. Watching D fall monotonically toward 0 at a fixed H below Sec.5's
    claimed Hc is a direct refutation that needs no functional form at all.
    """
    out = {}
    for h in hs:
        out[str(h)] = [float(np.interp(1.0 / h, taus, r["order"])) for r in rows]
    return out


def collapse_residual(coords, values):
    """Bhattacharjee-Seno residual, same estimator as freezeout_scaling.py."""
    u = np.concatenate(coords)
    d = np.concatenate(values)
    o = np.argsort(u)
    u, d = u[o], d[o]
    res, n = 0.0, 0
    for i in range(1, len(u) - 1):
        if u[i + 1] == u[i - 1]:
            continue
        f = (u[i] - u[i - 1]) / (u[i + 1] - u[i - 1])
        res += (d[i] - (d[i - 1] + (d[i + 1] - d[i - 1]) * f)) ** 2
        n += 1
    return res / max(n, 1)


def score_collapses(rows, taus, h_grid):
    """Zero-parameter shift collapse vs Sec.5's two-parameter FSS collapse.

    Both are scored on the SAME points -- D sampled at tau = 1/H over a common H
    grid -- with the same residual functional, so the comparison is like for like.
    """
    from scipy.optimize import minimize

    sample = []
    for r in rows:
        curve = np.interp(1.0 / h_grid, taus, r["order"])
        sample.append((float(r["omega"]), curve))

    def shift(A):
        return collapse_residual([(1.0 / h_grid) - A * np.log(w)
                                  for w, _ in sample], [c for _, c in sample])

    def fss(p):
        Hc, a = p
        if not (-1.0 < Hc < 1.0 and 0.0 <= a < 2.0):
            return 1e9            # keep Nelder-Mead out of the overflow corners
        return collapse_residual([(h_grid - Hc) * w ** a for w, _ in sample],
                                 [c for _, c in sample])

    best_fss = min((minimize(fss, x0, method="Nelder-Mead",
                             options={"xatol": 1e-5, "fatol": 1e-12})
                    for x0 in ([0.055, 0.38], [0.001, 0.15], [0.02, 0.25])),
                   key=lambda r: r.fun)
    best_shift = minimize(lambda p: shift(p[0]), [1.5], method="Nelder-Mead",
                          options={"xatol": 1e-5, "fatol": 1e-12})
    return {
        "shift_fixed_1.5": float(shift(1.5)),
        "shift_best": {"A": float(best_shift.x[0]), "residual": float(best_shift.fun)},
        "fss_best": {"Hc": float(best_fss.x[0]), "a": float(best_fss.x[1]),
                     "residual": float(best_fss.fun)},
    }


# --------------------------------------------------------------------------- #


def make_figure(rows, taus, fits, collapse, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    omegas = np.array([r["omega"] for r in rows], dtype=float)
    ts = np.array([r["tau_star"][LEVELS.index(0.5)] for r in rows])
    es = np.array([r["tau_star_sem"][LEVELS.index(0.5)] for r in rows])

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5))
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(rows)))

    ax = axes[0]
    for c, r in zip(colors, rows):
        ax.plot(taus, r["order"], color=c, lw=1.2, label=f"Ω={r['omega']}")
    ax.axhline(0.5, color="grey", ls=":", lw=0.8)
    ax.set_xlabel("internal time τ = 1/H")
    ax.set_ylabel("order parameter D")
    ax.set_title("One SSA pass per Ω gives every H\n(D(H,Ω) = D at τ = 1/H)")
    ax.legend(fontsize=6, ncol=2)
    ax.grid(alpha=0.25)

    ax = axes[1]
    x = np.log(omegas)
    ax.errorbar(x, ts, yerr=es, fmt="o", ms=5, color="k", capsize=2,
                label="measured τ*(D=0.5)")
    xx = np.linspace(x.min(), x.max(), 100)
    ax.plot(xx, fits["log"]["A"] * xx + fits["log"]["B"], "-", color="C0",
            label=f"(A={fits['log']['A']:.3f}±{fits['log']['A_sd']:.3f}) lnΩ + B")
    ax.plot(xx, 1.5 * xx + fits["log"]["B"], "--", color="C2",
            label="predicted slope 3/2")
    p = fits["pow"]
    ax.plot(xx, 1.0 / (p["Hc"] + p["C"] * np.exp(-p["a"] * xx)), ":", color="C3",
            lw=2, label=f"Hc+CΩ⁻ᵃ (Hc={p['Hc']:.4f})")
    old = fits["pow_old_range"]
    ax.plot(xx, 1.0 / (old["Hc"] + old["C"] * np.exp(-old["a"] * xx)), "-.",
            color="C1", lw=1.5,
            label=f"§5 fit on Ω≤1280 (Hc={old['Hc']:.4f}), extrapolated")
    ax.set_xlabel("ln Ω")
    ax.set_ylabel("τ* = 1/H*")
    ax.set_title("No critical rate: τ* is linear in lnΩ")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)

    ax = axes[2]
    for c, r in zip(colors, rows):
        ax.plot(taus - 1.5 * np.log(r["omega"]), r["order"], color=c, lw=1.2)
    ax.set_xlim(-8, 12)
    ax.set_xlabel("τ − (3/2) lnΩ")
    ax.set_ylabel("D")
    ax.set_title("Zero-parameter collapse\n"
                 f"residual {collapse['shift_fixed_1.5']:.2e} "
                 f"vs FSS {collapse['fss_best']['residual']:.2e}")
    ax.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote figure -> {out_path}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--omegas", type=int, nargs="+", default=None)
    p.add_argument("--trials", type=int, default=2500)
    p.add_argument("--reps", type=int, default=8)
    p.add_argument("--tau-max", type=float, default=70.0)
    p.add_argument("--tau-points", type=int, default=260)
    p.add_argument("--n-sym", type=int, default=2)
    p.add_argument("--bias", type=float, default=0.0,
                   help="Omega-independent initial pairwise margin (control: "
                        "with bias > 0, H* must NOT drift with Omega)")
    p.add_argument("--extinction-check", action="store_true",
                   help="how many committed species survive to the crossing "
                        "(n-winner only; explains the small-Omega shortfall)")
    p.add_argument("--sec6-check", action="store_true",
                   help="reproduce FINDINGS §6's H*(n) with a NON-expanding SSA")
    p.add_argument("--exact-omegas", type=int, nargs="*", default=None,
                   help="also solve the CME exactly at these Omega (cost ~Omega^3)")
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--from-data", action="store_true",
                   help="re-run the analysis on the saved sweep, no simulation")
    p.add_argument("--out", default=None)
    p.add_argument("--data", default=None)
    args = p.parse_args()

    tag = "" if args.n_sym == 2 else f"_n{args.n_sym}"
    if args.bias:
        tag += f"_bias{args.bias:g}"
    if args.out is None:
        args.out = os.path.join(os.path.dirname(__file__), f"freezeout_law{tag}.png")
    if args.data is None:
        args.data = os.path.join(os.path.dirname(__file__), os.pardir, "results",
                                 f"freezeout_law{tag}.json")
    if args.omegas is None:
        args.omegas = [40 * 2 ** k for k in range(13)]
    if args.quick:
        args.omegas, args.trials, args.reps = [40, 160, 640, 2560], 800, 4

    t0 = time.time()
    if args.from_data:
        with open(args.data) as fh:
            saved = json.load(fh)
        taus = np.array(saved["taus"])
        rows = saved["rows"]
        args.trials, args.reps = saved["trials"], saved["reps"]
        args.n_sym = rows[0].get("n_sym", args.n_sym)
        print(f"re-analysing {args.data}: {len(rows)} Ω x {args.reps} reps x "
              f"{args.trials} trials, n_sym={args.n_sym}")
    else:
        taus = np.geomspace(1.0, args.tau_max, args.tau_points)
        print(f"freeze-out law: {len(args.omegas)} Ω x {args.reps} reps x "
              f"{args.trials} trials, n_sym={args.n_sym}, bias={args.bias}, "
              f"jobs={args.jobs}")
        # biggest Omega first so the long pole starts immediately
        tasks = [(w, taus, args.trials, rep, args.seed, args.n_sym, args.bias)
                 for w in sorted(args.omegas, reverse=True)
                 for rep in range(args.reps)]
        if args.jobs > 1:
            import multiprocessing as mp
            with mp.Pool(min(args.jobs, len(tasks))) as pool:
                out = pool.map(_worker, tasks, chunksize=1)
        else:
            out = [_worker(t) for t in tasks]
        by = {}
        for r in out:
            by.setdefault(r["omega"], []).append(r)
        rows = [combine(w, by[w], taus, args.trials, args.n_sym)
                for w in sorted(by)]
        print(f"  sweep done in {time.time()-t0:.0f}s")

    omegas = [r["omega"] for r in rows]
    fits = {}
    print("\n   Ω      lnΩ   " + "  ".join(f"τ*({L})" for L in LEVELS)
          + "   τ_abs(med)")
    for r in rows:
        cells = "  ".join(f"{t:7.3f}" for t in r["tau_star"])
        print(f"{r['omega']:7d} {np.log(r['omega']):7.3f}  {cells}   "
              f"{r['tau_absorb_median']:7.3f}")

    if args.n_sym == 2 and not args.bias:
        print("\n-- LOCAL slope dτ*/dlnΩ, against the ODE route (both -> 3/2) --")
        i50 = LEVELS.index(0.5)
        # The ODE route cannot predict the seed AMPLITUDE, only its Omega^-1/2
        # scaling; t_level(s) = t_level(1) + 3 ln s, so one point calibrates it
        # and every slope stays a prediction.
        base = deterministic_times(rows[-1]["omega"])["t_level"]
        seed = float(np.exp((rows[-1]["tau_star"][i50] - base) / 3.0))
        print(f"   (ODE seed amplitude calibrated at Ω={rows[-1]['omega']}: "
              f"δ₀ = {seed:.3f}·Ω^-1/2; slopes are unaffected by it. The ODE "
              f"route\n    pins the LIMIT 3/2, not the size of the finite-Ω "
              f"excess -- a quenched\n    seed cannot reproduce continuous noise "
              f"injection at small Ω.)")
        print("     Ω pair        measured           ODE       Hc>0 needs -> 0")
        det = {r["omega"]: deterministic_times(r["omega"], seed_scale=seed)
               for r in rows}
        loc = []
        for a, b in zip(rows, rows[1:]):
            dl = np.log(b["omega"] / a["omega"])
            m = (b["tau_star"][i50] - a["tau_star"][i50]) / dl
            e = np.hypot(a["tau_star_sem"][i50], b["tau_star_sem"][i50]) / dl
            o = (det[b["omega"]]["t_level"] - det[a["omega"]]["t_level"]) / dl
            loc.append({"lo": a["omega"], "hi": b["omega"], "slope": float(m),
                        "sem": float(e), "ode": float(o)})
            print(f"  {a['omega']:6d}→{b['omega']:<7d}  {m:6.3f} ± {e:.3f}"
                  f"       {o:6.3f}")
        fits["local_slopes"] = loc
        fits["ode_seed_scale"] = seed
        print("   ODE offsets after calibration: "
              + "  ".join(f"{det[r['omega']]['t_level'] - r['tau_star'][i50]:+.2f}"
                          for r in rows))

    print("\n-- tau* = A lnΩ + B  (prediction A = 3/2) --")
    for i, L in enumerate(LEVELS):
        f = fit_log_law(omegas, [r["tau_star"][i] for r in rows],
                        [r["tau_star_sem"][i] for r in rows])
        fits[f"log_{L}"] = f
        print(f"  D={L}:  A = {f['A']:.4f} ± {f['A_sd']:.4f}   B = {f['B']:.3f}"
              f"   rms {f['rms']:.4f}   ({(f['A']-1.5)/f['A_sd']:+.1f}σ from 3/2)")
    fits["log"] = fits[f"log_{0.5}"]

    fa = fit_log_law(omegas, [r["tau_absorb_median"] for r in rows],
                     [max(r["tau_absorb_median_sem"], 1e-4) for r in rows])
    fits["log_absorb"] = fa
    print(f"  absorption median: A = {fa['A']:.4f} ± {fa['A_sd']:.4f}"
          f"   (prediction 5/2)")

    i50 = LEVELS.index(0.5)
    t50 = [r["tau_star"][i50] for r in rows]
    s50 = [r["tau_star_sem"][i50] for r in rows]
    fits["pow"] = fit_power_law(omegas, t50, s50)
    old = [j for j, w in enumerate(omegas) if w <= 1280]
    fits["pow_old_range"] = fit_power_law([omegas[j] for j in old],
                                          [t50[j] for j in old],
                                          [s50[j] for j in old])
    fits["curvature"] = curvature_test(omegas, t50, s50)
    print(f"\n-- Hc + C Ω^-a on the FULL range: Hc = {fits['pow']['Hc']:.5f}, "
          f"a = {fits['pow']['a']:.3f}, chi2/dof = "
          f"{fits['pow']['chi2']/max(fits['pow']['dof'],1):.2f}"
          f"   (log law chi2/dof = "
          f"{fits['log']['chi2']/max(fits['log']['dof'],1):.2f})")
    print(f"-- same fit restricted to Ω ≤ 1280: Hc = "
          f"{fits['pow_old_range']['Hc']:.5f}, a = {fits['pow_old_range']['a']:.3f}"
          f"  -> predicts τ*(Ω={omegas[-1]}) = "
          f"{1.0/(fits['pow_old_range']['Hc'] + fits['pow_old_range']['C']*omegas[-1]**-fits['pow_old_range']['a']):.2f}"
          f", measured {t50[-1]:.2f}")
    c = fits["curvature"]
    print(f"-- curvature of τ* in lnΩ: {c['quad']:+.4f} ± {c['quad_sd']:.4f} "
          f"({c['sigma']:+.1f}σ); a positive Hc requires it to be negative")

    fixed = fixed_H_table(rows, taus)
    print("\n-- D at FIXED H as Ω grows (no fit; §5 claims Hc ≈ 0.055) --")
    print("      Ω   " + "  ".join(f"H={h:<7}" for h in fixed))
    for i, r in enumerate(rows):
        print(f"  {r['omega']:7d}  "
              + "  ".join(f"{fixed[h][i]:7.4f}  " for h in fixed))

    h_grid = 1.0 / np.geomspace(4.0, 40.0, 60)[::-1]
    collapse = score_collapses(rows, taus, h_grid)
    print(f"\n-- collapse residuals (lower is better, same estimator) --")
    print(f"  shift, A fixed at 3/2 (0 params): {collapse['shift_fixed_1.5']:.4e}")
    print(f"  shift, A fitted   (1 param):      "
          f"{collapse['shift_best']['residual']:.4e}  (A={collapse['shift_best']['A']:.3f})")
    print(f"  §5 FSS            (2 params):     "
          f"{collapse['fss_best']['residual']:.4e}  "
          f"(Hc={collapse['fss_best']['Hc']:.4f}, a={collapse['fss_best']['a']:.3f})")

    extinct = None
    if args.extinction_check:
        if args.n_sym < 3:
            raise SystemExit("--extinction-check needs --n-sym >= 3")
        print(f"\n-- committed species alive when D crosses 0.5 (n={args.n_sym}) --")
        extinct = extinction_check(args.n_sym)
        print("      Ω   Ω/(2n−1)   τ*(0.5)   alive/n")
        for r in extinct["rows"]:
            print(f"  {r['omega']:6d}   {r['per_species']:8.1f}   "
                  f"{r['tau_star']:7.2f}   "
                  f"{r['alive_at_crossing']/args.n_sym:6.3f}"
                  f"   ({r['alive_at_crossing']:.2f} of {args.n_sym})")

    sec6 = None
    if args.sec6_check:
        print("\n-- §6's H*(n) from an ORDINARY (non-expanding) SSA at Ω=160 --")
        sec6 = sec6_check()
        print("   n   §6 H*    1/τ* here    τ* measured        ratio τ*·H*_§6")
        for r in sec6["rows"]:
            print(f" {r['n']:3d}   {r['sec6_hstar']:.4f}   {r['hstar']:.4f}     "
                  f"{r['tau_star']:6.3f} ± {r['tau_star_sem']:.3f}     "
                  f"{r['ratio']:.3f}")

    exact = None
    if args.exact_omegas:
        print("\n-- exact CME route (no sampling error) --")
        exact = exact_route(args.exact_omegas, args.tau_max)
        ssa = {r["omega"]: r for r in rows}
        print("     Ω     exact τ*    SSA τ*            diff")
        for r in exact["rows"]:
            s = ssa.get(r["omega"])
            if s is None:
                continue
            i = LEVELS.index(0.5)
            print(f"  {r['omega']:6d}  {r['tau_star'][i]:9.4f}  "
                  f"{s['tau_star'][i]:7.4f} ± {s['tau_star_sem'][i]:.4f}   "
                  f"{r['tau_star'][i]-s['tau_star'][i]:+.4f}")

    os.makedirs(os.path.dirname(os.path.abspath(args.data)), exist_ok=True)
    with open(args.data, "w") as fh:
        json.dump({"taus": taus.tolist(), "levels": list(LEVELS),
                   "rows": rows, "fits": fits, "collapse": collapse,
                   "exact": exact, "fixed_H": fixed, "bias": args.bias,
                   "sec6": sec6, "extinction": extinct,
                   "trials": args.trials, "reps": args.reps}, fh)
    print(f"\nwrote data -> {args.data}   (total {time.time()-t0:.0f}s)")
    make_figure(rows, taus, fits, collapse, args.out)


if __name__ == "__main__":
    main()
