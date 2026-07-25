"""Part C: what does a restoring stage cost, and where does it stop working?

FINDINGS 7 measured WHY restoration matters -- a restoring cascade carries a bit
to depth 45 where analog passthrough is a coin flip by 22 -- but priced nothing,
because irreversible AM's dissipation is formally infinite. This puts a number on
it, exactly (crnl/cascade_exact.py: no trials, no sampling error).

WHAT THIS EXPERIMENT WILL NOT CLAIM, and why the tables below look defensive.
Two earlier designs each produced a clean, plausible, WRONG headline:

  1. Stopping a stage at 0.7*delta_star and emitting +-1 made the stop predicate
     fire on the INITIAL state: 83-96% of stages ran zero reactions and the
     harness did the restoring for free. The measured "cost per stage" was a
     duty cycle.
  2. Fixing that, but comparing against a control confined to +-1 while scaling
     the noise by delta_star(gamma), produced "restoration requires a minimum
     Omega". That was the control: its dynamic range was absolute while its
     noise was in landscape units, and the mismatch grows with gamma. Under a
     matched-rail control the crossover vanishes entirely.

So this script reports BOTH control conventions on every row and refuses to
name a winner inside a tie band. A verdict that disagrees between the two
columns is a statement about the comparator, not about chemistry.

    python -m experiments.dissipation_cascade
    python -m experiments.dissipation_cascade --quick
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

from crnl.cascade_exact import (
    TIE_BAND, run_cascade, run_control, verdict,
)
from crnl.networks.am_reversible import GAMMA_C, delta_star


def one_cell(gamma: float, omega: int, t_stage: float, depth: int,
             noise_frac: float, chunk: int) -> dict:
    """One (gamma, Omega, t_stage) cell, with both controls."""
    r = run_cascade(gamma, omega, t_stage, depth, noise_frac, chunk)
    d = r["delta_star"]
    s = r["sigma_ch"]
    wide = run_control(omega, depth, s, d, rail=1.0)
    matched = run_control(omega, depth, s, d, rail=d)
    chem = r["p_correct"][-1]
    r["control_wide"] = wide["p_correct"]
    r["control_matched"] = matched["p_correct"]
    r["verdict_wide"] = verdict(chem, wide["p_correct"][-1])
    r["verdict_matched"] = verdict(chem, matched["p_correct"][-1])
    r["conventions_agree"] = r["verdict_wide"] == r["verdict_matched"]
    return r


def print_table(rows, depth):
    print(f"\n{'gamma':>6} {'Omega':>6} {'t_stg':>6} {'d*':>6} {'dS/stg':>8} "
          f"{'chem':>7} {'ctrl±1':>7} {'v':>6} {'ctrl±d*':>8} {'v':>6}  agree")
    for r in rows:
        print(f"{r['gamma']:>6.2f} {r['omega']:>6} {r['t_stage']:>6.1f} "
              f"{r['delta_star']:>6.3f} {r['ds_per_stage']:>8.2f} "
              f"{r['p_correct'][-1]:>7.4f} {r['control_wide'][-1]:>7.4f} "
              f"{r['verdict_wide']:>6} {r['control_matched'][-1]:>8.4f} "
              f"{r['verdict_matched']:>6}  "
              f"{'yes' if r['conventions_agree'] else '** NO **'}")
    disputed = [r for r in rows if not r["conventions_agree"]]
    if disputed:
        print(f"\n{len(disputed)} of {len(rows)} cells DISAGREE between the two "
              "control conventions. Those cells say nothing about chemistry:")
        for r in disputed:
            print(f"    gamma={r['gamma']:.2f} Omega={r['omega']} "
                  f"t={r['t_stage']}: {r['verdict_wide']} vs "
                  f"{r['verdict_matched']}")
    print(f"\n(verdicts use a tie band of {TIE_BAND}; both arms decay toward 0.5 "
          "with depth, so raw comparisons return verdicts on noise)")


def make_figure(rows, depth, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    depths = np.arange(1, depth + 1)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16.5, 5))

    # panel 1: survival vs depth at the reference Omega, both controls
    ref = sorted({r["omega"] for r in rows})[len(set(r["omega"] for r in rows)) // 2]
    sel = [r for r in rows if r["omega"] == ref]
    colors = plt.cm.viridis(np.linspace(0.1, 0.85, len(sel)))
    for c, r in zip(colors, sel):
        ax1.plot(depths, r["p_correct"], "-", color=c, lw=2,
                 label=f"γ={r['gamma']:.2f} (A={r['affinity']:.2f})")
    if sel:
        ax1.plot(depths, sel[0]["control_wide"], ":", color="#d62728", lw=1.8,
                 label="control, rails ±1")
        ax1.plot(depths, sel[0]["control_matched"], "--", color="#7f7f7f", lw=1.8,
                 label="control, rails ±δ*  (fair)")
    ax1.axhline(0.5, color="black", lw=0.7, ls=":")
    ax1.set_xlabel("cascade depth (stages)"); ax1.set_ylabel("P(bit correct)")
    ax1.set_title(f"Survival vs depth (Ω={ref})\ntwo controls, because they disagree")
    ax1.legend(fontsize=7.5); ax1.grid(alpha=0.25)

    # panel 2: the cost
    for c, r in zip(colors, sel):
        ax2.plot(depths, r["cum_ds"], "-", color=c, lw=2, label=f"γ={r['gamma']:.2f}")
    ax2.set_xlabel("cascade depth"); ax2.set_ylabel("cumulative ⟨ΔS⟩  [k_B]")
    ax2.set_title("The bill: free energy to carry one bit\n"
                  "(either control is the zero line)")
    ax2.legend(fontsize=8); ax2.grid(alpha=0.25)

    # panel 3: fidelity vs cost per stage -- the honest headline
    for c, r in zip(colors, sel):
        ax3.plot(r["ds_per_stage"], r["p_correct"][-1], "o", color=c, ms=11,
                 label=f"γ={r['gamma']:.2f}")
    if sel:
        ax3.axhline(sel[0]["control_matched"][-1], color="#7f7f7f", ls="--",
                    lw=1.5, label="fair control (free)")
    ax3.set_xlabel("dissipation per stage  [k_B T]")
    ax3.set_ylabel(f"P(bit correct) at depth {depth}")
    ax3.set_title("HEADLINE: restoration costs MORE where it works LESS\n"
                  "(cost rises toward γ_c as fidelity collapses)")
    ax3.legend(fontsize=8); ax3.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote figure -> {out_path}")


def main():
    here = os.path.dirname(__file__)
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gammas", type=float, nargs="+",
                   default=[0.05, 0.15, 0.30, 0.45])
    p.add_argument("--omegas", type=int, nargs="+", default=[30, 60, 120])
    p.add_argument("--t-stage", type=float, nargs="+", default=[8.0])
    p.add_argument("--depth", type=int, default=30)
    p.add_argument("--noise-frac", type=float, default=0.35,
                   help="sigma_ch as a fraction of delta_star(gamma). In "
                        "LANDSCAPE units: delta_star(0)=1 exactly, so FINDINGS "
                        "7 is the gamma->0 member of this family")
    p.add_argument("--chunk", type=int, default=64)
    p.add_argument("--allow-large", action="store_true",
                   help="permit Omega > 120 (Omega=240 costs ~3 min and ~400 MB "
                        "per cell even chunked)")
    p.add_argument("--quick", action="store_true")
    p.add_argument("--out", default=os.path.join(here, "dissipation_cascade.png"))
    p.add_argument("--data", default=os.path.join(
        here, os.pardir, "results", "dissipation_cascade.json"))
    args = p.parse_args()

    if args.quick:
        args.gammas, args.omegas, args.depth = [0.15, 0.45], [30], 12

    too_big = [o for o in args.omegas if o > 120]
    if too_big and not args.allow_large:
        p.error(f"Omega {too_big} exceeds 120; pass --allow-large (each such "
                "cell costs minutes and hundreds of MB)")

    gammas = [g for g in args.gammas if 0.0 < g < GAMMA_C]
    skipped = [g for g in args.gammas if g not in gammas]
    if skipped:
        print(f"skipped gamma = {skipped}: no bistable landscape (need "
              f"0 < gamma < {GAMMA_C})")

    print(f"Part C: the price of a restoring stage  (depth={args.depth}, "
          f"sigma_ch={args.noise_frac}·δ*(γ), exact -- no sampling)")

    rows = []
    for t_stage in args.t_stage:
        for omega in args.omegas:
            for gamma in gammas:
                rows.append(one_cell(gamma, omega, t_stage, args.depth,
                                     args.noise_frac, args.chunk))
    print_table(rows, args.depth)

    os.makedirs(os.path.dirname(os.path.abspath(args.data)), exist_ok=True)
    with open(args.data, "w") as fh:
        json.dump({"depth": args.depth, "noise_frac": args.noise_frac,
                   "tie_band": TIE_BAND, "t_stages": args.t_stage,
                   "omegas": args.omegas, "skipped_gammas": skipped,
                   "rows": rows}, fh, indent=2)
    print(f"wrote data -> {args.data}")
    make_figure(rows, args.depth, args.out)


if __name__ == "__main__":
    main()
