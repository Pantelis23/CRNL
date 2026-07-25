"""What does it cost to deliver one bit to depth D? No comparator required.

Every verdict in FINDINGS 10 needed a passive control, and that control's dynamic
range was a free parameter which had already manufactured one withdrawn claim.
This experiment removes the comparator: it measures the mutual information
between the input bit and the depth-D output, and divides the cumulative
dissipation by it. The result is k_B T per bit delivered -- no rails to choose,
no tie band, and directly comparable to k_B T ln 2.

THE DEPTH IS PART OF THE QUESTION. At depth 1 the measure is degenerate: a stage
with t_stage -> 0 does nothing, costs nothing, and still scores well because one
channel application barely damages a bit sitting on a rail (0.89 kT/bit at
t=0.05 versus 20.2 at t=16). Once the bit must survive a depth at which a passive
channel would lose it, doing nothing becomes the worst strategy rather than the
best, and an interior optimum in t_stage appears. So every number here is quoted
with its depth, and the depth profile is written to the JSON.

    python -m experiments.bit_cost
    python -m experiments.bit_cost --quick
"""

from __future__ import annotations

import argparse
import json
import math
import os

import numpy as np

from crnl.information import cost_per_bit
from crnl.networks.am_reversible import GAMMA_C, cycle_affinity, delta_star
from crnl.networks.am_reversible import am_reversible, reverse_pairing

LN2 = math.log(2.0)


def sweep(gammas, omegas, t_stages, depth, noise_frac, chunk):
    rows = []
    for omega in omegas:
        for gamma in gammas:
            net = am_reversible(gamma)
            A = cycle_affinity(net, reverse_pairing(net))
            for t in t_stages:
                prof = cost_per_bit(gamma, omega, t, depth, noise_frac, chunk)
                last = prof[-1]
                rows.append({
                    "gamma": gamma, "omega": omega, "t_stage": t,
                    "affinity": A, "delta_star": delta_star(gamma),
                    "depth": depth, "I_bits": last["I_bits"],
                    "ds_total": last["ds"], "kT_per_bit": last["kT_per_bit"],
                    "landauer_ratio": last["kT_per_bit"] / LN2,
                    "profile": [{"depth": r["depth"], "I_bits": r["I_bits"],
                                 "ds": r["ds"], "kT_per_bit": r["kT_per_bit"]}
                                for r in prof],
                })
    return rows


def make_figure(rows, depth, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16.5, 5))
    omegas = sorted({r["omega"] for r in rows})
    gammas = sorted({r["gamma"] for r in rows})
    ref_om = omegas[len(omegas) // 2]

    # 1: depth profile of cost per bit, best t_stage per gamma
    colors = plt.cm.viridis(np.linspace(0.1, 0.85, len(gammas)))
    for c, g in zip(colors, gammas):
        cand = [r for r in rows if r["omega"] == ref_om and r["gamma"] == g
                and math.isfinite(r["kT_per_bit"])]
        if not cand:
            continue
        best = min(cand, key=lambda r: r["kT_per_bit"])
        d = [p["depth"] for p in best["profile"]]
        y = [p["kT_per_bit"] for p in best["profile"]]
        ax1.semilogy(d, y, "-", color=c, lw=2,
                     label=f"γ={g:.2f} (t*={best['t_stage']:g})")
    ax1.axhline(LN2, color="black", ls=":", lw=1.2, label="k_B T ln 2")
    ax1.set_xlabel("cascade depth D"); ax1.set_ylabel("k_B T per bit delivered")
    ax1.set_title(f"Cost per bit grows superlinearly with depth (Ω={ref_om})")
    ax1.legend(fontsize=8); ax1.grid(alpha=0.25, which="both")

    # 2: information retained
    for c, g in zip(colors, gammas):
        cand = [r for r in rows if r["omega"] == ref_om and r["gamma"] == g
                and math.isfinite(r["kT_per_bit"])]
        if not cand:
            continue
        best = min(cand, key=lambda r: r["kT_per_bit"])
        ax2.plot([p["depth"] for p in best["profile"]],
                 [p["I_bits"] for p in best["profile"]], "-", color=c, lw=2,
                 label=f"γ={g:.2f}")
    ax2.set_xlabel("cascade depth D"); ax2.set_ylabel("I(bit ; output)  [bits]")
    ax2.set_title("What actually survives\n(no control needed to say this)")
    ax2.legend(fontsize=8); ax2.grid(alpha=0.25)

    # 3: the headline -- cheapest bit vs population
    for om, mark in zip(omegas, "os^Dv"):
        pts = [(r["gamma"], r["kT_per_bit"]) for r in rows
               if r["omega"] == om and math.isfinite(r["kT_per_bit"])]
        by_g = {}
        for g, v in pts:
            by_g[g] = min(v, by_g.get(g, float("inf")))
        if by_g:
            gs = sorted(by_g)
            ax3.semilogy(gs, [by_g[g] for g in gs], mark + "-", lw=1.8,
                         label=f"Ω={om}")
    ax3.axhline(LN2, color="black", ls=":", lw=1.2, label="k_B T ln 2")
    ax3.set_xlabel("γ  (weaker drive →)")
    ax3.set_ylabel(f"k_B T per bit at depth {depth}")
    ax3.set_title("HEADLINE: weak drive is not cheap — it delivers nothing\n"
                  "and a bigger population costs more per bit, not less")
    ax3.legend(fontsize=8); ax3.grid(alpha=0.25, which="both")

    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote figure -> {out_path}")


def main():
    here = os.path.dirname(__file__)
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gammas", type=float, nargs="+",
                   default=[0.05, 0.15, 0.30, 0.45])
    p.add_argument("--omegas", type=int, nargs="+", default=[30, 60, 120])
    p.add_argument("--t-stages", type=float, nargs="+",
                   default=[4.0, 16.0, 64.0])
    p.add_argument("--depth", type=int, default=30,
                   help="the depth the bit must reach. NOT optional and never "
                        "1: at depth 1 the measure rewards a stage that does "
                        "nothing (see module docstring)")
    p.add_argument("--noise-frac", type=float, default=0.35)
    p.add_argument("--chunk", type=int, default=64)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--out", default=os.path.join(here, "bit_cost.png"))
    p.add_argument("--data", default=os.path.join(
        here, os.pardir, "results", "bit_cost.json"))
    args = p.parse_args()

    if args.quick:
        args.gammas, args.omegas = [0.05, 0.30], [30]
        args.t_stages, args.depth = [4.0, 16.0], 12

    if args.depth < 5:
        p.error(f"--depth {args.depth} is too shallow to be well posed: below "
                "~5 the cheapest 'restoration' is one that does not restore")

    gammas = [g for g in args.gammas if 0.0 < g < GAMMA_C]
    print(f"Cost per bit delivered to depth {args.depth}  "
          f"(σ_ch = {args.noise_frac}·δ*(γ), exact — no sampling, no control)")

    rows = sweep(gammas, args.omegas, args.t_stages, args.depth,
                 args.noise_frac, args.chunk)

    print(f"\n{'gamma':>6} {'Omega':>6} {'t_stg':>6} {'I(bits)':>8} "
          f"{'total ΔS':>9} {'kT/bit':>10} {'/ln2':>9}")
    for r in rows:
        ratio = ("%.0fx" % r["landauer_ratio"]) if math.isfinite(
            r["landauer_ratio"]) else "inf"
        kpb = ("%.1f" % r["kT_per_bit"]) if math.isfinite(
            r["kT_per_bit"]) else "inf"
        print(f"{r['gamma']:>6.2f} {r['omega']:>6} {r['t_stage']:>6.1f} "
              f"{r['I_bits']:>8.4f} {r['ds_total']:>9.1f} {kpb:>10} {ratio:>9}")

    finite = [r for r in rows if math.isfinite(r["kT_per_bit"])]
    if finite:
        b = min(finite, key=lambda r: r["kT_per_bit"])
        print(f"\nCHEAPEST BIT at depth {args.depth}: {b['kT_per_bit']:.1f} k_BT "
              f"(γ={b['gamma']}, Ω={b['omega']}, t_stage={b['t_stage']:g}, "
              f"I={b['I_bits']:.4f} bits)")
        print(f"  = {b['landauer_ratio']:.0f}x k_B T ln 2. Landauer bounds "
              "ERASURE, not transmission, so this is a scale comparison and not "
              "a claim that the bound is approached.")

    os.makedirs(os.path.dirname(os.path.abspath(args.data)), exist_ok=True)
    with open(args.data, "w") as fh:
        json.dump({"depth": args.depth, "noise_frac": args.noise_frac,
                   "omegas": args.omegas, "t_stages": args.t_stages,
                   "rows": rows}, fh, indent=2)
    print(f"wrote data -> {args.data}")
    make_figure(rows, args.depth, args.out)


if __name__ == "__main__":
    main()
