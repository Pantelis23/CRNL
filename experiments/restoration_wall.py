"""The restoration wall: measure exp(-c(eps) * Omega) from data (design.md §4).

Fix a small fractional bias eps (default 51/49 between X and Y). Sweep Omega,
run many Gillespie trials to absorption at each, bin outcomes into
{X-wins, Y-wins, all-blank}, and watch the error fraction fall exponentially in
Omega. The deterministic ODE from the same 51/49 start glides to the X rail at
every Omega -- that curve is the lie; the stochastic error fraction is the truth;
c(eps) * Omega is the restoration wall as a number.

Key protocol choices, straight from §4:
  * Hold the *fraction* constant across Omega, never the absolute count
    difference -- only that keeps the starting distance from the separatrix
    Omega-independent, which is what makes the exponent clean.
  * Report TWO curves: raw error (Y / all trials) and conditional-on-decision
    (Y / (X + Y)). At low Omega the all-blank bin is non-negligible, so they
    differ materially; collapsing them hides the finite-count effect. The clean
    exponential lives in the conditional fraction; the raw fraction additionally
    carries the all-blank Omega-scaling.
  * Fit on the linear region only. The error-rich window overlaps the region
    where an algebraic prefactor still curves the log plot; fit near the top of
    the low range and treat a curved low-Omega tail as prefactor, not exponent.

Usage:
    python -m experiments.restoration_wall            # sensible default sweep
    python -m experiments.restoration_wall --trials 20000 --jobs 8
    python -m experiments.restoration_wall --quick    # fast smoke run

Reports eps alongside c(eps) (the exponent is a function of the bias, not a
universal constant) and writes a figure to experiments/restoration_wall.png.
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass, asdict

import numpy as np

from crnl import approximate_majority, gillespie, seed_for
from crnl import classify


# --------------------------------------------------------------------------- #
# One Omega worth of trials                                                    #
# --------------------------------------------------------------------------- #

@dataclass
class OmegaResult:
    omega: int
    trials: int
    x_wins: int
    y_wins: int  # the error (with an X-favoured start)
    blank: int

    @property
    def raw_error(self) -> float:
        """Y / all trials -- carries the all-blank scaling too."""
        return self.y_wins / self.trials

    @property
    def conditional_error(self) -> float:
        """Y / (X + Y) -- the clean exponential lives here."""
        decided = self.x_wins + self.y_wins
        return self.y_wins / decided if decided else float("nan")

    @property
    def blank_fraction(self) -> float:
        return self.blank / self.trials


def initial_counts(omega: int, bias: float) -> np.ndarray:
    """Integer counts holding the X:Y *fraction* fixed at (1+bias)/2 : (1-bias)/2.

    bias = 0.02 gives the 51/49 split of §4. B starts at 0; committed molecules
    fill Omega. Rounding is absorbed into B so the total is exactly Omega.
    """
    fx = (1.0 + bias) / 2.0
    nX = int(round(fx * omega))
    nY = int(round((1.0 - fx) * omega))
    nB = omega - nX - nY
    return np.array([nX, nY, nB])


def run_omega(omega: int, trials: int, bias: float, base_seed: int) -> OmegaResult:
    net = approximate_majority()
    n0 = initial_counts(omega, bias)
    x = y = b = 0
    for t in range(trials):
        res = gillespie(net, n0, omega, seed_for(omega, t, base=base_seed))
        outcome = classify.classify_am_outcome(res)
        if outcome == "X":
            x += 1
        elif outcome == "Y":
            y += 1
        elif outcome == "B":
            b += 1
        # 'undecided' cannot occur for AM (every trajectory absorbs); ignored.
    return OmegaResult(omega, trials, x, y, b)


def run_blank_probe(omega: int, trials: int, base_seed: int) -> OmegaResult:
    """Balanced-start trials to expose the all-blank bin at very low Omega.

    All-blank (0,0,Omega) is the finite-count outcome the deterministic repeller
    denies. It is only non-negligible at the low-Omega end and from a near-
    balanced start (the last committed pair annihilates via r1 before a lead
    amplifies), so we probe it separately from the biased restoration sweep.
    """
    net = approximate_majority()
    nX = omega // 2
    n0 = np.array([nX, omega - nX, 0])
    x = y = b = 0
    for t in range(trials):
        res = gillespie(net, n0, omega, seed_for(omega, t, base=base_seed + 101))
        outcome = classify.classify_am_outcome(res)
        if outcome == "X":
            x += 1
        elif outcome == "Y":
            y += 1
        elif outcome == "B":
            b += 1
    return OmegaResult(omega, trials, x, y, b)


def _worker(args):
    return run_omega(*args)


def _blank_worker(args):
    return run_blank_probe(*args)


def run_blank_sweep(omegas, trials, base_seed, jobs):
    tasks = [(int(w), trials, base_seed) for w in omegas]
    if jobs and jobs > 1:
        import multiprocessing as mp

        with mp.Pool(jobs) as pool:
            results = pool.map(_blank_worker, tasks)
    else:
        results = [_blank_worker(t) for t in tasks]
    return sorted(results, key=lambda r: r.omega)


def run_sweep(omegas, trials, bias, base_seed, jobs):
    tasks = [(int(w), trials, bias, base_seed) for w in omegas]
    if jobs and jobs > 1:
        # chunk trials across processes per Omega for better load balance
        import multiprocessing as mp

        with mp.Pool(jobs) as pool:
            results = pool.map(_worker, tasks)
    else:
        results = [_worker(t) for t in tasks]
    return sorted(results, key=lambda r: r.omega)


# --------------------------------------------------------------------------- #
# Fit the exponent on the clean linear region                                 #
# --------------------------------------------------------------------------- #

@dataclass
class Fit:
    c: float  # the barrier / noise margin: error ~ exp(-c * Omega)
    intercept: float
    omega_lo: int
    omega_hi: int
    n_points: int
    r2: float


def fit_exponent(results, min_events: int = 15, max_frac: float = 0.35) -> Fit | None:
    """Fit log(conditional_error) = intercept - c * Omega on the linear band.

    The band that is simultaneously error-rich and cleanly straight can be
    narrow (§4). We keep only Omega points where the error is (a) resolved --
    at least ``min_events`` error events, so the log is not dominated by
    counting noise -- and (b) not saturated -- conditional error below
    ``max_frac``, avoiding the flat low-Omega top. Everything dropped is logged,
    never silently truncated (§5).
    """
    xs, ys, kept = [], [], []
    dropped = []
    for r in results:
        ce = r.conditional_error
        if r.y_wins >= min_events and 0 < ce < max_frac:
            xs.append(r.omega)
            ys.append(np.log(ce))
            kept.append(r.omega)
        else:
            reason = (
                "too-few-events" if r.y_wins < min_events else "saturated/zero"
            )
            dropped.append((r.omega, reason, r.y_wins))
    if dropped:
        print("  fit drops (Omega, reason, y_events):")
        for om, reason, ye in dropped:
            print(f"    Omega={om:<4} {reason:<16} y_events={ye}")
    if len(xs) < 2:
        print("  not enough resolved points to fit an exponent.")
        return None

    xs = np.array(xs, dtype=float)
    ys = np.array(ys, dtype=float)
    slope, intercept = np.polyfit(xs, ys, 1)
    pred = slope * xs + intercept
    ss_res = np.sum((ys - pred) ** 2)
    ss_tot = np.sum((ys - ys.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return Fit(
        c=-slope,
        intercept=intercept,
        omega_lo=min(kept),
        omega_hi=max(kept),
        n_points=len(xs),
        r2=r2,
    )


# --------------------------------------------------------------------------- #
# Plot                                                                         #
# --------------------------------------------------------------------------- #

def make_figure(results, fit, bias, out_path, blank_results=None):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    omegas = np.array([r.omega for r in results])
    raw = np.array([r.raw_error for r in results])
    cond = np.array([r.conditional_error for r in results])
    blank = np.array([r.blank_fraction for r in results])
    n = results[0].trials

    # binomial standard error on the conditional fraction
    dec = np.array([r.x_wins + r.y_wins for r in results])
    cond_se = np.sqrt(np.clip(cond * (1 - cond), 0, None) / np.clip(dec, 1, None))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # -- left: the restoration wall, log scale --
    ax1.set_yscale("log")
    ax1.errorbar(
        omegas, np.clip(cond, 1e-6, None), yerr=cond_se, fmt="o", ms=5,
        color="#1f77b4", label="conditional  Y/(X+Y)", capsize=2, zorder=3,
    )
    ax1.plot(
        omegas, np.clip(raw, 1e-6, None), "s", ms=4, mfc="none",
        color="#7f7f7f", label="raw  Y/all", zorder=2,
    )
    if fit is not None:
        grid = np.linspace(omegas.min(), omegas.max(), 100)
        ax1.plot(
            grid, np.exp(fit.intercept - fit.c * grid), "-", color="#d62728",
            lw=2,
            label=f"fit  exp(-c·Ω),  c(ε)={fit.c:.4f}\n"
                  f"Ω∈[{fit.omega_lo},{fit.omega_hi}], R²={fit.r2:.3f}",
            zorder=4,
        )
    # ODE error is exactly 0 at every Omega; a log axis cannot draw 0, so mark
    # it as a floor annotation rather than a phantom line.
    ax1.plot([], [], color="#2ca02c", lw=1.5, ls="--",
             label="deterministic ODE: error = 0 at every Ω (off the log floor)")
    ax1.set_xlabel("population scale  Ω")
    ax1.set_ylabel("error fraction (log)")
    ax1.set_title(f"Restoration wall  ε={bias:.3f}  (start "
                  f"{(1 + bias) / 2:.0%}/{(1 - bias) / 2:.0%},  {n} trials/Ω)")
    ax1.legend(fontsize=8, loc="upper right")
    ax1.grid(True, which="both", alpha=0.25)

    # -- right: the all-blank bin the deterministic repeller denies --
    if blank_results:
        bo = np.array([r.omega for r in blank_results])
        bf = np.array([r.blank_fraction for r in blank_results])
        ax2.plot(bo, bf, "D-", color="#9467bd", ms=6)
        ax2.set_title("All-blank outcome, balanced start\n"
                      "(ODE calls (0,0,1) a repeller; finite Ω lands there "
                      "anyway)")
    else:
        ax2.plot(omegas, blank, "D-", color="#9467bd", ms=5)
        ax2.set_title("All-blank outcome (ODE says (0,0,1) is a repeller;\n"
                      "finite Ω lands there anyway)")
    ax2.set_xlabel("population scale  Ω")
    ax2.set_ylabel("all-blank fraction  B/all")
    ax2.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"\nwrote figure -> {out_path}")


# --------------------------------------------------------------------------- #
# CLI                                                                          #
# --------------------------------------------------------------------------- #

def deterministic_contrast(bias):
    """Integrate the ODE from the same biased start; it glides to X every time."""
    net = approximate_majority()
    fx = (1.0 + bias) / 2.0
    final = classify.find_stable_endpoint(net, [fx, 1 - fx, 0.0])
    return final


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bias", type=float, default=0.10,
                   help="fractional bias eps (0.10 = 55/45; 0.02 = the literal "
                        "51/49 of §4). Reported alongside c(eps).")
    p.add_argument("--omegas", type=int, nargs="+", default=None,
                   help="explicit Omega sweep; default is an auto low-window sweep")
    p.add_argument("--trials", type=int, default=8000,
                   help="Gillespie trials per Omega (§4 wants 1e4-1e5)")
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument("--seed", type=int, default=0, help="base seed for replay")
    p.add_argument("--quick", action="store_true",
                   help="tiny fast sweep for a smoke test")
    p.add_argument("--out", default=os.path.join(
        os.path.dirname(__file__), "restoration_wall.png"))
    args = p.parse_args()

    if args.quick:
        omegas = args.omegas or [20, 30, 40, 50, 60, 70, 80]
        trials = 2000
    else:
        # observable window: low Omega where errors actually occur (§4).
        omegas = args.omegas or [20, 30, 40, 50, 60, 70, 80, 100, 120]
        trials = args.trials

    print(f"CRNL restoration wall  (eps={args.bias}, {trials} trials/Omega, "
          f"jobs={args.jobs})")
    if args.bias <= 0.03:
        print("  note: at the literal 51/49 (eps=0.02) the quasipotential "
              "barrier c(eps) is\n"
              "  intrinsically tiny, so c*Omega stays < 1 across the whole "
              "observable window and\n"
              "  the wall does not clear the algebraic-prefactor crossover "
              "until Omega ~ thousands,\n"
              "  where errors fall to ~e^-25 and read as exactly zero (the §4 "
              "squeeze). Use the\n"
              "  default eps=0.10 to see a clean exp(-c*Omega); 51/49 is kept "
              "for reproducing the crossover.")
    print(f"start fraction: X={(1 + args.bias) / 2:.4f}  "
          f"Y={(1 - args.bias) / 2:.4f}")

    det = deterministic_contrast(args.bias)
    print(f"deterministic ODE from the same start settles at "
          f"X={det[0]:.4f} Y={det[1]:.4f} B={det[2]:.4f}  "
          f"-> the ODE 'lie': error 0 at every Omega\n")

    results = run_sweep(omegas, trials, args.bias, args.seed, args.jobs)

    print(f"{'Omega':>6} {'X':>7} {'Y(err)':>7} {'B':>6} "
          f"{'raw':>9} {'cond':>9} {'blank':>8}")
    for r in results:
        print(f"{r.omega:>6} {r.x_wins:>7} {r.y_wins:>7} {r.blank:>6} "
              f"{r.raw_error:>9.4g} {r.conditional_error:>9.4g} "
              f"{r.blank_fraction:>8.4g}")

    print("\nfitting exp(-c(eps)*Omega) on the clean linear band:")
    fit = fit_exponent(results)
    if fit is not None:
        print(f"\n  c(eps={args.bias}) = {fit.c:.4f}   "
              f"(barrier per unit Omega; the noise margin)")
        print(f"  fit window Omega in [{fit.omega_lo}, {fit.omega_hi}], "
              f"{fit.n_points} points, R^2={fit.r2:.3f}")
        print(f"  => at Omega=1000 this predicts error ~ "
              f"exp(-{fit.c * 1000:.0f}); the wall is exponential in Omega.")

    # dedicated low-Omega balanced sweep to expose the all-blank bin (§3.4)
    blank_omegas = [6, 8, 10, 12, 14, 16, 20]
    blank_trials = max(trials, 20000) if not args.quick else 6000
    print(f"\nall-blank probe (balanced start, low Omega, "
          f"{blank_trials} trials/Omega):")
    blank_results = run_blank_sweep(blank_omegas, blank_trials, args.seed, args.jobs)
    for r in blank_results:
        print(f"  Omega={r.omega:>3}  blank_fraction={r.blank_fraction:.4g}  "
              f"(B={r.blank})")

    make_figure(results, fit, args.bias, args.out, blank_results=blank_results)


if __name__ == "__main__":
    main()
