"""Part 0: does a landscape exist at all, and at what drive?

Verification of the closed-form result, not a search. Reversible AM's symmetric
point sits at (1/3,1/3,1/3) for EVERY gamma, and the decision mode there has
eigenvalue

    lambda_antisym(gamma) = (1 - 2 gamma)/3

which is +1/3 at gamma=0 (irreversible AM's saddle, design.md 2.3) and vanishes at

    gamma_c = 1/2      A(gamma_c) = 3 ln 2 = 2.0794...

Above gamma_c the two attractors have merged into the symmetric point: a single
minimum, no threshold, nothing to restore toward -- no population size Omega can
help. So "no restoration at equilibrium" is a statement about the DETERMINISTIC
landscape, with an exact threshold, and A(gamma_c) is the minimum drive that buys
a landscape at all.

A rejected alternative, recorded so it is not retried: bisecting gamma on whether
mirror-image ODE runs separate is biased UPWARD by critical slowing down (the
estimate drifts 0.50071 -> 0.50208 as the separation tolerance tightens from 1e-3
to 1e-6, i.e. tightening the tolerance makes the answer worse).

    python -m experiments.reversible_landscape
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

from crnl.networks.am_reversible import (
    GAMMA_C,
    am_reversible,
    cycle_affinity,
    delta_star,
    fixed_points,
    lambda_antisym,
    lambda_sym,
    reverse_pairing,
)
from crnl.deterministic import jacobian

PITCHFORK_AMPLITUDE = 4 * np.sqrt(2) / 3      # 1.8856180832


def scan(gammas):
    rows = []
    for g in gammas:
        net = am_reversible(g)
        J = jacobian(net, [1 / 3, 1 / 3, 1 / 3])
        fps = fixed_points(g)
        att = [f for f in fps if f["kind"] == "attractor"]
        rows.append({
            "gamma": float(g),
            "affinity": (cycle_affinity(net, reverse_pairing(net))
                         if g > 0 else None),   # gamma=0 is infinite drive; None keeps the JSON strict
            "lambda_antisym": lambda_antisym(g),
            "lambda_antisym_numeric": float(J[0, 0] - J[0, 1]),
            "lambda_sym": lambda_sym(g),
            "delta_star": delta_star(g),
            "n_fixed_points": len(fps),
            "x_plus": att[0]["x"] if att else None,
            "y_plus": att[0]["y"] if att else None,
            "b_star": att[0]["b"] if att else None,
            "bistable": bool(att),
        })
    return rows


def make_figure(rows, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g = np.array([r["gamma"] for r in rows])
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 5))

    ax1.plot(g, [r["lambda_antisym"] for r in rows], "-", lw=2,
             color="#1f77b4", label="closed form (1−2γ)/3")
    ax1.plot(g, [r["lambda_antisym_numeric"] for r in rows], "o", ms=4,
             color="#d62728", label="numeric Jacobian")
    ax1.axhline(0, color="grey", lw=0.8)
    ax1.axvline(GAMMA_C, color="black", ls="--", lw=1, label=f"γ_c={GAMMA_C}")
    ax1.set_xlabel("γ  (reverse/forward rate ratio)")
    ax1.set_ylabel("restoring gain  λ")
    ax1.set_title("Decision-mode gain vanishes at γ_c")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.25)

    ax2.plot(g, [r["b_star"] if r["b_star"] is not None else np.nan for r in rows],
             "-", lw=2, color="#2ca02c", label="b* = γ/(1+γ)")
    ax2.plot(g, [r["x_plus"] if r["x_plus"] is not None else np.nan for r in rows],
             "-", lw=2, color="#1f77b4", label="x* (winner)")
    ax2.plot(g, [r["y_plus"] if r["y_plus"] is not None else np.nan for r in rows],
             "-", lw=2, color="#ff7f0e", label="y* (loser)")
    ax2.axvline(GAMMA_C, color="black", ls="--", lw=1)
    ax2.set_xlabel("γ")
    ax2.set_ylabel("attractor composition")
    ax2.set_title("Rails move inward; the thermal\npopulation is mainly B")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.25)

    mask = g < GAMMA_C
    ax3.plot(g[mask], [r["delta_star"] for r, m in zip(rows, mask) if m],
             "-", lw=2, color="#9467bd", label="δ*(γ)")
    gg = np.linspace(max(g.min(), 0.30), GAMMA_C, 200)
    ax3.plot(gg, PITCHFORK_AMPLITUDE * np.sqrt(GAMMA_C - gg), "--",
             color="#d62728", label=f"(4√2/3)·√(γ_c−γ)")
    ax3.axvline(GAMMA_C, color="black", ls="--", lw=1)
    ax3.set_xlabel("γ")
    ax3.set_ylabel("attractor separation δ*")
    ax3.set_title("Pitchfork: the two rails merge at γ_c")
    ax3.legend(fontsize=8)
    ax3.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote figure -> {out_path}")


def main():
    here = os.path.dirname(__file__)
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--points", type=int, default=121)
    p.add_argument("--gamma-max", type=float, default=0.8)
    p.add_argument("--out", default=os.path.join(here, "reversible_landscape.png"))
    p.add_argument("--data", default=os.path.join(
        here, os.pardir, "results", "reversible_landscape.json"))
    args = p.parse_args()

    gammas = np.linspace(0.0, args.gamma_max, args.points)
    rows = scan(gammas)

    print(f"gamma_c = {GAMMA_C}   A(gamma_c) = 3*ln2 = {3*np.log(2):.10f}")
    print(f"{'gamma':>7} {'lambda':>10} {'delta*':>9} {'b*':>8} {'A':>9} bistable")
    for r in rows[::max(1, len(rows) // 12)]:
        b = "" if r["b_star"] is None else f"{r['b_star']:.4f}"
        A = "inf" if r["affinity"] is None else f"{r['affinity']:.4f}"
        print(f"{r['gamma']:>7.3f} {r['lambda_antisym']:>10.5f} "
              f"{r['delta_star']:>9.5f} {b:>8} {A:>9} {r['bistable']}")

    worst = max(abs(r["lambda_antisym"] - r["lambda_antisym_numeric"]) for r in rows)
    print(f"\nclosed form vs numeric Jacobian: worst deviation {worst:.2e}")

    os.makedirs(os.path.dirname(os.path.abspath(args.data)), exist_ok=True)
    with open(args.data, "w") as fh:
        json.dump({"gamma_c": GAMMA_C, "affinity_at_gamma_c": 3 * np.log(2),
                   "pitchfork_amplitude": PITCHFORK_AMPLITUDE, "scan": rows}, fh,
                  indent=2)
    print(f"wrote data -> {args.data}")
    make_figure(rows, args.out)


if __name__ == "__main__":
    main()
