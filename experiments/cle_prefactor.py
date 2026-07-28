"""Is the CLE's small excess a prefactor, or just discretisation? (T11b-i)

FINDINGS 21.3a left one residue standing. Across all eight n = 2 cells at 40,000
trials the CLE sat ABOVE the exact CME error probability -- ratios 1.018 to 1.045,
every z positive, combining to roughly 3 sigma for a ~+3% uniform excess -- while
its error EXPONENT matched exact (0.980 and 1.001). A uniform factor on p is a
prefactor effect and leaves the exponent alone, which is what one expects if the
CLE is the correct diffusion limit with a slightly wrong amplitude. But it was one
step size and 3 sigma.

THE TEST, and it discriminates cleanly. Euler-Maruyama is weak-order 1, so a
DISCRETISATION excess is proportional to dt and extrapolates to zero. A PREFACTOR
excess is a property of the Langevin limit itself and is flat in dt, surviving
dt -> 0. Three step sizes a factor 4 apart, 100,000 trials per cell, against the
exact CME.

PREDICTIONS, written before running:

  P1  If the excess is real it is FLAT in dt and extrapolates to a nonzero
      constant as dt -> 0.
  P2  It should also be roughly constant across Omega and eps, since a prefactor
      multiplies p rather than shifting its exponent -- and the exponent is
      already known to match.
  P3  The exponent stays at 1.00 +- 0.02 at every dt. If the exponent moves with
      dt instead, the excess is discretisation and 21.3a's first draft was right
      after all for the wrong measurement.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.approximations import run_batch
from crnl.networks.am_reversible import am_reversible
from crnl.vectorized import compile_network
from experiments.approximation_hierarchy import _setup, p_cme

GAMMA, THETA = 0.30, 0.80


def _stop(thr):
    return lambda m: np.abs(m[:, 0] - m[:, 1]) >= thr


def cell(omega, eps_frac, dt, trials, seed):
    n0, thr, _ = _setup(GAMMA, omega, eps_frac, THETA)
    comp = compile_network(am_reversible(GAMMA), float(omega))
    out = run_batch(comp, n0, np.random.default_rng(seed), step=dt,
                    stop=_stop(thr), trials=trials, t_max=4000.0, poisson=False)
    fin = out["n"]
    done = _stop(thr)(fin)
    k = int(done.sum())
    if k == 0:
        return float("nan"), float("nan"), 0
    p = float((fin[done][:, 0] <= fin[done][:, 1]).mean())
    return p, float(np.sqrt(p * (1 - p) / k)), k


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omegas", type=int, nargs="+", default=[40, 70, 100, 140])
    ap.add_argument("--eps-fracs", type=float, nargs="+", default=[0.25, 0.40])
    ap.add_argument("--dts", type=float, nargs="+", default=[0.02, 0.005, 0.00125])
    ap.add_argument("--trials", type=int, default=100_000)
    ap.add_argument("--seed", type=int, default=20260729)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/cle_prefactor.json"))
    args = ap.parse_args()

    t0 = time.time()
    rows = []
    exact = {(om, ef): p_cme(GAMMA, om, ef, THETA)
             for om in args.omegas for ef in args.eps_fracs}
    print(f"gamma={GAMMA} theta={THETA} trials={args.trials:,}")
    for ef in args.eps_fracs:
        print(f"\n=== eps/delta* = {ef}   (ratio CLE/CME per cell)")
        print(f"{'Omega':>6} {'CME':>11} " +
              " ".join(f"{'dt=' + str(d):>17}" for d in args.dts))
        for om in args.omegas:
            pc = exact[(om, ef)]
            line = f"{om:>6} {pc:>11.4e} "
            for dt in args.dts:
                p, se, k = cell(om, ef, dt, args.trials, args.seed + om)
                rows.append({"omega": om, "eps_frac": ef, "dt": dt, "p_cme": pc,
                             "p_cle": p, "se": se, "n_done": k,
                             "ratio": p / pc, "z": (p - pc) / se})
                line += f" {p/pc:>9.4f} (z{(p-pc)/se:>+5.1f})"
            print(line)

    print(f"\n=== P1: is the excess flat in dt?  (pooled over Omega, per eps)")
    print(f"{'eps':>6} " + " ".join(f"{'dt=' + str(d):>16}" for d in args.dts))
    for ef in args.eps_fracs:
        line = f"{ef:>6.2f} "
        for dt in args.dts:
            sel = [r for r in rows if r["eps_frac"] == ef and r["dt"] == dt]
            rr = np.array([r["ratio"] for r in sel])
            # se on the mean ratio, propagated from each cell's binomial se
            ses = np.array([r["se"] / r["p_cme"] for r in sel])
            m = rr.mean()
            sem = float(np.sqrt((ses ** 2).sum()) / len(sel))
            line += f" {m:>8.4f}+-{sem:<6.4f}"
        print(line)

    print(f"\n=== P3: does the exponent move with dt?")
    print(f"{'eps':>6} {'CME':>9} " + " ".join(f"{'dt=' + str(d):>11}" for d in args.dts))
    for ef in args.eps_fracs:
        om = np.array(args.omegas, dtype=float)
        pc = np.array([exact[(o, ef)] for o in args.omegas])
        sc = np.polyfit(om, -np.log(pc), 1)[0]
        line = f"{ef:>6.2f} {sc:>9.5f} "
        for dt in args.dts:
            sel = sorted([r for r in rows if r["eps_frac"] == ef and r["dt"] == dt],
                         key=lambda r: r["omega"])
            p = np.array([r["p_cle"] for r in sel])
            g = np.isfinite(p) & (p > 0)
            s = np.polyfit(om[g], -np.log(p[g]), 1)[0] if g.sum() >= 3 else np.nan
            line += f" {s/sc:>11.4f}"
        print(line)

    print("\nA ratio that shrinks toward 1 as dt falls is discretisation "
          "(Euler-Maruyama is\nweak-order 1). One that holds is the Langevin "
          "limit's own prefactor.")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
