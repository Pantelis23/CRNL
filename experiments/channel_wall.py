"""Where the restoration wall gives way to the channel floor.

FINDINGS 11.1 found that cost per bit RISES with population and that the
efficient frontier's marginal cost explodes. The reason turns out to be a
crossover between two regimes that this project had measured separately and never
connected:

  * FINDINGS 1-2, the RESTORATION WALL: error falls like exp(-kappa e^2 Omega),
    so more molecules buy exponentially better reliability;
  * FINDINGS 11, the CHANNEL FLOOR: injected noise sigma_ch flips the sign
    outright with an Omega-independent probability, so more molecules buy
    nothing.

A saddle-point argument gives ONE expression covering both, with no fitted
parameter (see crnl/information.predicted_exponent):

    -ln p  ~  kappa Omega delta*^2 / (1 + 2 kappa Omega sigma^2)

with kappa(gamma) = lambda(gamma) / (2 D_0(gamma)) = (3/2)(1 - 2 gamma)/(1 + gamma).
The crossover is at Omega_x = 1/(2 kappa sigma^2).

NOTE ON WHICH KAPPA. FINDINGS 12 was produced with kappa = (3/2)(1 - 2 gamma),
which scales the restoring gain with gamma but leaves the diffusion at its
gamma = 0 value. FINDINGS 15 corrects that -- the reverse reactions are extra
jumps along the decision mode, so D_0(gamma) = (1+gamma)/9 -- and refitting these
same 216 cells with the corrected value lifts the pooled collapse from R^2 = 0.933
to 0.960. `information.wall_coefficient` now returns the corrected value, so a
fresh run of this experiment will NOT reproduce FINDINGS 12's printed slopes;
`information.wall_coefficient_gain_only` does.

This experiment measures p(Omega, sigma_ch) directly and tests that collapse.

    python -m experiments.channel_wall
    python -m experiments.channel_wall --quick
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

from crnl.information import (
    crossover_omega, flip_probability, predicted_exponent, wall_coefficient,
)
from crnl.networks.am_reversible import GAMMA_C, delta_star

#: p below this is at the resolution limit of the I(D) fit and is not data.
P_FLOOR = 1e-15


def sweep(gammas, omegas, noise_fracs, t_stage, depth, chunk):
    rows = []
    for gamma in gammas:
        for nf in noise_fracs:
            for omega in omegas:
                p = flip_probability(gamma, omega, t_stage, nf, depth, chunk)
                rows.append({
                    "gamma": gamma, "omega": omega, "noise_frac": nf,
                    "t_stage": t_stage, "p_flip": p,
                    "delta_star": delta_star(gamma),
                    "kappa": wall_coefficient(gamma),
                    "omega_crossover": crossover_omega(gamma, nf),
                    "predicted_exponent": predicted_exponent(gamma, omega, nf),
                    "usable": bool(np.isfinite(p) and p > P_FLOOR),
                })
    return rows


def collapse_fit(rows):
    """Regress -ln p on the parameter-free predicted exponent."""
    use = [r for r in rows if r["usable"]]
    if len(use) < 4:
        return None
    x = np.array([r["predicted_exponent"] for r in use])
    y = np.array([-np.log(r["p_flip"]) for r in use])
    slope, intercept = np.polyfit(x, y, 1)
    r2 = 1.0 - np.sum((y - (slope * x + intercept)) ** 2) / np.sum(
        (y - y.mean()) ** 2)
    return {"slope": float(slope), "intercept": float(intercept),
            "r2": float(r2), "n": len(use)}


def make_figure(rows, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
    gammas = sorted({r["gamma"] for r in rows})
    ref = gammas[0]
    nfs = sorted({r["noise_frac"] for r in rows})
    colors = plt.cm.plasma(np.linspace(0.05, 0.85, len(nfs)))

    for c, nf in zip(colors, nfs):
        sel = [r for r in rows if r["gamma"] == ref and r["noise_frac"] == nf
               and r["usable"]]
        if not sel:
            continue
        sel.sort(key=lambda r: r["omega"])
        ax1.semilogy([r["omega"] for r in sel], [r["p_flip"] for r in sel],
                     "o-", color=c, lw=1.8, ms=4,
                     label=f"σ/δ*={nf:g}  (Ω×={sel[0]['omega_crossover']:.0f})")
        ax1.axvline(sel[0]["omega_crossover"], color=c, ls=":", lw=1, alpha=0.6)
    ax1.set_xlabel("population Ω"); ax1.set_ylabel("per-stage flip probability")
    ax1.set_title(f"Wall → floor (γ={ref})\ndotted: predicted Ω× = 1/(2κσ²)")
    ax1.legend(fontsize=7.5); ax1.grid(alpha=0.25, which="both")

    marks = "os^Dv<>"
    for m, g in zip(marks, gammas):
        sel = [r for r in rows if r["gamma"] == g and r["usable"]]
        if not sel:
            continue
        ax2.plot([r["predicted_exponent"] for r in sel],
                 [-np.log(r["p_flip"]) for r in sel], m, ms=5, alpha=0.8,
                 label=f"γ={g:.2f}")
    lim = max([-np.log(r["p_flip"]) for r in rows if r["usable"]] + [1.0])
    ax2.plot([0, lim], [0, lim], "k--", lw=1, label="slope 1 (no free parameter)")
    ax2.set_xlabel("predicted  κΩδ*²/(1+2κΩσ²)")
    ax2.set_ylabel("measured  −ln p")
    ax2.set_title("One formula, both regimes, every γ and σ")
    ax2.legend(fontsize=8); ax2.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote figure -> {out_path}")


def main():
    here = os.path.dirname(__file__)
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gammas", type=float, nargs="+", default=[0.05, 0.15, 0.30, 0.45])
    p.add_argument("--omegas", type=int, nargs="+",
                   default=[4, 8, 12, 16, 24, 32, 48, 64, 96])
    p.add_argument("--noise-fracs", type=float, nargs="+",
                   default=[0.10, 0.15, 0.20, 0.28, 0.35, 0.45])
    p.add_argument("--t-stage", type=float, default=16.0)
    p.add_argument("--depth", type=int, default=40)
    p.add_argument("--chunk", type=int, default=64)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--out", default=os.path.join(here, "channel_wall.png"))
    p.add_argument("--data", default=os.path.join(
        here, os.pardir, "results", "channel_wall.json"))
    args = p.parse_args()

    if args.quick:
        args.gammas, args.omegas = [0.05, 0.30], [8, 16, 32]
        args.noise_fracs = [0.15, 0.35]

    gammas = [g for g in args.gammas if 0.0 < g < GAMMA_C]
    print(f"Channel floor vs restoration wall  (t_stage={args.t_stage}, "
          f"depth={args.depth}, exact)")

    rows = sweep(gammas, args.omegas, args.noise_fracs, args.t_stage,
                 args.depth, args.chunk)

    for g in gammas:
        print(f"\nγ={g:.2f}   δ*={delta_star(g):.4f}   κ={wall_coefficient(g):.3f}")
        print(f"{'σ/δ*':>6} {'Ω× pred':>8} " +
              " ".join(f"{o:>9}" for o in args.omegas))
        for nf in args.noise_fracs:
            sel = {r["omega"]: r for r in rows
                   if r["gamma"] == g and r["noise_frac"] == nf}
            cells = []
            for o in args.omegas:
                r = sel.get(o)
                cells.append("%9.2e" % r["p_flip"] if r and r["usable"]
                             else "        -")
            print(f"{nf:>6.2f} {crossover_omega(g, nf):>8.1f} " + " ".join(cells))

    fits = {"pooled": collapse_fit(rows)}
    for g in gammas:
        fits[f"gamma={g}"] = collapse_fit([r for r in rows if r["gamma"] == g])

    print("\nCOLLAPSE onto the parameter-free prediction  −ln p ≈ "
          "κΩδ*²/(1+2κΩσ²)")
    print(f"{'subset':>14} {'slope':>8} {'intercept':>10} {'R²':>8} {'n':>4}")
    for name, f in fits.items():
        if f is None:
            print(f"{name:>14}   (too few usable points)")
            continue
        print(f"{name:>14} {f['slope']:>8.3f} {f['intercept']:>10.2f} "
              f"{f['r2']:>8.4f} {f['n']:>4}")
    print("\nSlope 1 with a small intercept would mean the saddle point is exact; "
          "it is not\n(it drops the prefactor and the Gaussian-tail correction), "
          "so read the R² —\nthe claim is that ONE expression collapses both "
          "regimes, not that it is exact.")

    dropped = [r for r in rows if not r["usable"]]
    if dropped:
        print(f"\n{len(dropped)} of {len(rows)} cells unusable (p below the "
              f"{P_FLOOR:g} resolution of the I(D) fit) — reported, never fitted.")

    os.makedirs(os.path.dirname(os.path.abspath(args.data)), exist_ok=True)
    with open(args.data, "w") as fh:
        json.dump({"t_stage": args.t_stage, "depth": args.depth,
                   "fits": fits, "rows": rows}, fh, indent=2)
    print(f"wrote data -> {args.data}")
    make_figure(rows, args.out)


if __name__ == "__main__":
    main()
