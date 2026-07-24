"""Predict the barrier before measuring it: c(eps) from the saddle geometry.

Every other experiment here *fits* the barrier c out of data. This one derives it
first and then checks the data lands on the predicted line -- a far sharper test of
understanding than fitting after the fact (design.md §8).

THE DERIVATION (full version in design.md §9)

Reduce AM to its decision coordinate delta = x - y. From the mass-action ODEs,

    d(delta)/dt = delta * (1 - x - y) = delta * b

so at the symmetric saddle (x=y=b=1/3) the decision direction is linearly unstable
with rate lambda = 1/3 -- exactly the saddle's positive eigenvalue.

Only the two recruitment reactions move delta, each by +-1/Omega, and at the saddle
each fires at rate Omega/9. The diffusion coefficient of delta is therefore

    D_delta = 1/2 * [ (1/Omega)^2 * Omega/9  +  (1/Omega)^2 * Omega/9 ] = 1/(9 Omega)

Near the saddle delta is a linearly *unstable* Ornstein-Uhlenbeck process,
d(delta) = lambda*delta*dt + sqrt(2 D_delta) dW, whose eventual sign is Gaussian
with variance sigma^2 = D_delta / lambda = 1/(3 Omega). The wrong-basin (splitting)
probability is then

    P(error) = Phi(-delta_0/sigma) ~ exp( -delta_0^2 / (2 sigma^2) )
             = exp( -(3/2) * delta_0^2 * Omega )

    ==>   c(eps) = (3/2) * eps^2        [KAPPA_THEORY = 1.5, quadratic in bias]

This experiment measures c at several eps -- each on its own Omega grid, sized to
that barrier so the fit lands in the error-rich band -- and fits c ~ kappa*eps^p.

Headline result (committed figure): p = 2.08 (quadratic confirmed), and c/eps^2
drifts DOWN toward the predicted 1.5 as eps -> 0, reaching 1.59 at the smallest,
cleanest point (eps=0.04, R2=0.999) -- within 6% of a first-principles prediction.
The upward drift at larger eps is the expected breakdown of the small-bias
linearization (plus bias amplification during the B-buildup transient).

    python -m experiments.quasipotential --quick
    python -m experiments.quasipotential --eps 0.04 0.06 0.08 0.10 0.14 0.20
"""

from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from experiments.restoration_wall import run_sweep, fit_exponent

#: Analytic prediction c(eps) = KAPPA_THEORY * eps^2, from the saddle geometry.
KAPPA_THEORY = 1.5


def omega_grid_for(eps, kappa_guess=1.8, factors=(0.3, 0.45, 0.65, 0.9, 1.3, 1.9)):
    """Omega grid centred where the wall is visible for this eps.

    error ~ exp(-c*Omega) with c ~ kappa*eps^2, so error ~ 0.1 near
    Omega = 2.3/c. Scaling the grid per eps is what keeps every fit inside the
    error-rich band instead of saturating (small eps) or vanishing (large eps).
    """
    om0 = 2.3 / (kappa_guess * eps ** 2)
    return sorted(set(int(round(om0 * f)) for f in factors))


def measure_c(eps, trials, seed, jobs, kappa_guess=1.8):
    grid = omega_grid_for(eps, kappa_guess)
    res = run_sweep(grid, trials, eps, seed, jobs)
    fit = fit_exponent(res)
    return {
        "eps": eps,
        "grid": grid,
        "c": fit.c if fit else float("nan"),
        "r2": fit.r2 if fit else float("nan"),
        "band": [fit.omega_lo, fit.omega_hi] if fit else None,
    }


def fit_power_law(eps, c):
    """Fit c = kappa * eps^p; returns (p, kappa) using only finite positive c."""
    eps = np.asarray(eps, dtype=float)
    c = np.asarray(c, dtype=float)
    good = np.isfinite(c) & (c > 0)
    if good.sum() < 2:
        return float("nan"), float("nan")
    p, logk = np.polyfit(np.log(eps[good]), np.log(c[good]), 1)
    return float(p), float(np.exp(logk))


def make_figure(rows, p, kappa, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    eps = np.array([r["eps"] for r in rows])
    c = np.array([r["c"] for r in rows])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.loglog(eps, c, "o", ms=8, color="#1f77b4", label="measured c(ε)")
    xx = np.geomspace(eps.min(), eps.max(), 50)
    ax1.loglog(xx, kappa * xx ** p, "-", color="#d62728",
               label=f"fit  c = {kappa:.2f}·ε^{p:.2f}")
    ax1.loglog(xx, KAPPA_THEORY * xx ** 2, "--", color="black",
               label=f"THEORY  c = {KAPPA_THEORY}·ε²")
    ax1.set_xlabel("initial bias ε")
    ax1.set_ylabel("restoration barrier c(ε)")
    ax1.set_title("Quasipotential prediction vs measurement")
    ax1.legend(fontsize=9)
    ax1.grid(True, which="both", alpha=0.3)

    # the sharper view: the prefactor itself, which must approach 3/2 as eps->0
    ax2.semilogx(eps, c / eps ** 2, "o-", ms=7, color="#1f77b4",
                 label="measured c(ε)/ε²")
    ax2.axhline(KAPPA_THEORY, color="black", ls="--",
                label=f"theory κ = {KAPPA_THEORY}")
    ax2.set_xlabel("initial bias ε")
    ax2.set_ylabel("c(ε) / ε²")
    ax2.set_title("Prefactor converges to theory as ε→0\n"
                  "(drift at large ε = linearization breakdown)")
    ax2.legend(fontsize=9)
    ax2.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote figure -> {out_path}")


def main():
    p_ = argparse.ArgumentParser(description=__doc__)
    p_.add_argument("--eps", type=float, nargs="+",
                    default=[0.04, 0.06, 0.08, 0.10, 0.14, 0.20])
    p_.add_argument("--trials", type=int, default=9000)
    p_.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p_.add_argument("--seed", type=int, default=0)
    p_.add_argument("--quick", action="store_true")
    p_.add_argument("--out", default=os.path.join(
        os.path.dirname(__file__), "quasipotential.png"))
    p_.add_argument("--data", default=os.path.join(
        os.path.dirname(__file__), os.pardir, "results",
        "quasipotential_ceps.json"))
    args = p_.parse_args()

    if args.quick:
        args.eps, args.trials = [0.10, 0.14, 0.20], 2500

    print(f"quasipotential test  (theory: c = {KAPPA_THEORY}·ε², "
          f"trials={args.trials}, jobs={args.jobs})")
    t0 = time.time()
    rows = []
    for eps in args.eps:
        r = measure_c(eps, args.trials, args.seed, args.jobs)
        rows.append(r)
        print(f"  eps={eps:.3f}  c={r['c']:.5g}  R2={r['r2']:.3f}  "
              f"c/eps²={r['c']/eps**2:.3f}  (theory {KAPPA_THEORY})", flush=True)

    p, kappa = fit_power_law([r["eps"] for r in rows], [r["c"] for r in rows])
    print(f"\nfit  c = {kappa:.3f}·ε^{p:.3f}   "
          f"(theory: p=2, kappa={KAPPA_THEORY})")

    os.makedirs(os.path.dirname(os.path.abspath(args.data)), exist_ok=True)
    with open(args.data, "w") as fh:
        json.dump({str(r["eps"]): r for r in rows}, fh, indent=2)
    print(f"wrote data -> {args.data}   (total {time.time()-t0:.0f}s)")
    make_figure(rows, p, kappa, args.out)


if __name__ == "__main__":
    main()
