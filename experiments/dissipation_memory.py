"""Part B: what does REMEMBERING a decision require?

At finite gamma a decided state is only metastable: the reverse reactions
regenerate B and let the loser take over, so the memory has a finite lifetime tau.
For each gamma we solve the exact CME for the mean first-passage time from the
positive attractor to the opposite side, and the stationary entropy-production
rate sigma.

WHY THE QUESTION IS PHRASED AS "WHAT DRIVE", NOT "WHAT POWER":

  * sigma -> 0 in BOTH limits: at gamma=1 by detailed balance, and as gamma->0
    because sigma = A*J and the cycle flux J collapses faster than A = -3 ln gamma
    grows. It peaks in the middle.
  * sigma and tau move in OPPOSITE directions across the bistable range, so
    "power buys retention" is false. But note sigma is a rate and tau a time --
    that pairing is not a correlation. The dimensionless product sigma*tau (total
    dissipation per lifetime) is monotone and reads the other way.

So the honest headline is tau against the DRIVE (gamma, or the affinity A), with
sigma reported alongside as the transient price and sigma*tau as the cost-like
observable. The zero-power memory limit is the textbook ideal zero-leak ratchet
(gamma->0 is a singular limit where the states become absorbing, and gamma->0 IS
A->infinity) -- not a new result.

The CME is the primary instrument here because it is EXACT and one solve gives the
whole first-passage field -- not because it is faster. (An earlier docstring claimed
SSA would need "hundreds of hours" at Omega=120; measured, one flip at gamma=0.35
costs 5.5 minutes. SSA only becomes hopeless at gamma <= 0.30, which is where this
solve's own validity guard rejects the answer.) SSA cross-checks are Plan 2.

    python -m experiments.dissipation_memory
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

from crnl.cme import ep_rate, first_passage
from crnl.networks.am_reversible import (
    GAMMA_C,
    am_reversible,
    cycle_affinity,
    delta_star,
    fixed_points,
    reverse_pairing,
)
from crnl.stochastic import seed_for
from crnl.thermo import gillespie_instrumented
from crnl.vectorized import compile_network

STEPS_PER_SECOND = 84_000.0      # measured on this engine, reversible AM
STEPS_PER_UNIT_TIME = 0.4        # x Omega; measured 0.38-0.48 across (gamma, Omega)


def ssa_flip_rate(gamma: float, omega: int, tau_cme: float, max_seconds: float,
                  seed: int = 0, theta_frac: float = 0.7,
                  target_flips: int = 30) -> dict:
    """Measure the flip rate by SSA and compare with 1/tau from the CME.

    GATED ON THE CME, BY DESIGN: the exact tau prices the run BEFORE paying for
    it. And the budget is ENFORCED, not merely accepted -- an earlier draft took
    max_seconds as an argument and never used it, so raising it to 3600 admitted
    a 1.8 h run with max_steps as the only backstop.

    CONVENTION: FlipCounter counts ONE-WAY crossings, so flips/T -> 1/tau.
    Using 1/(2 tau) (the ROUND-TRIP rate) halves tau_ssa and was measured to
    give ratios of 0.37-0.58 with 52-82 flips -- a false alarm outside sampling
    noise. A residual ratio of ~0.8 IS expected and is not a bug: the flip clock
    starts at the ARM (0.7 delta*) while the CME's tau starts at the ATTRACTOR
    (measured 0.77-0.85 for that offset alone).
    """
    arm = theta_frac * delta_star(gamma)
    net = am_reversible(gamma)
    pairing = reverse_pairing(net)
    comp = compile_network(net, float(omega))

    att = [f for f in fixed_points(gamma) if f["kind"] == "attractor"]
    hi = max(att, key=lambda f: f["x"])
    n_x = int(round(hi["x"] * omega))
    n_b = int(round(hi["b"] * omega))
    start = np.array([n_x, omega - n_x - n_b, n_b], dtype=np.int64)

    t_needed = tau_cme * target_flips
    # hard cap from the budget, so the run cannot outlive it
    max_steps = int(max_seconds * STEPS_PER_SECOND)
    r = gillespie_instrumented(comp, start, seed_for(omega, 0, base=seed),
                               pairing, flip_arm=arm, omega=omega,
                               t_max=t_needed, max_steps=max_steps,
                               species=list(net.species))
    rate = r.flips / r.t_final if r.t_final > 0 else float("nan")
    return {
        "flips": r.flips, "t_total": r.t_final, "steps": r.steps,
        "flip_rate": rate,
        "tau_ssa": (1.0 / rate) if rate > 0 else float("nan"),
        "enough": r.flips >= 10, "budget_hit": r.steps >= max_steps,
        "t_target": t_needed,
    }


def run_point(gamma: float, omega: int, theta_frac: float = 0.7) -> dict:
    net = am_reversible(gamma)
    pairing = reverse_pairing(net)
    A = cycle_affinity(net, pairing)

    # start at the positive attractor, rounded onto the integer lattice
    att = [f for f in fixed_points(gamma) if f["kind"] == "attractor"]
    hi = max(att, key=lambda f: f["x"])
    n_x = int(round(hi["x"] * omega))
    n_b = int(round(hi["b"] * omega))
    n_y = omega - n_x - n_b
    start = np.array([n_x, n_y, n_b], dtype=np.int64)

    # a flip is crossing to the far side of the landscape
    theta = theta_frac * delta_star(gamma)

    def absorbing(n):
        return (int(n[0]) - int(n[1])) <= -theta * omega

    fp = first_passage(net, omega, float(omega), start, absorbing, pairing)
    sigma = ep_rate(net, omega, float(omega), pairing)
    tau = fp["mean_time"]
    return {
        "gamma": gamma,
        "affinity": A,
        "omega": omega,
        "start": start.tolist(),
        "theta": theta,
        "tau": tau,
        "sigma": sigma,
        "sigma_tau": sigma * tau if fp["valid"] else float("nan"),
        "flip_rate": 1.0 / (2.0 * tau) if fp["valid"] and tau > 0 else float("nan"),
        "valid": fp["valid"],
        "residual": fp["residual"],
    }


def make_figure(rows, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ok = [r for r in rows if r["valid"]]
    by_omega = {}
    for r in ok:
        by_omega.setdefault(r["omega"], []).append(r)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 5))
    colors = plt.cm.viridis(np.linspace(0, 0.8, len(by_omega)))

    for c, (om, rs) in zip(colors, sorted(by_omega.items())):
        g = [r["gamma"] for r in rs]
        ax1.semilogy(g, [r["tau"] for r in rs], "o-", color=c, label=f"Ω={om}")
        ax3.semilogy([r["affinity"] for r in rs], [r["tau"] for r in rs],
                     "o-", color=c, label=f"Ω={om}")
    ax1.axvline(GAMMA_C, color="black", ls="--", lw=1, label=f"γ_c={GAMMA_C}")
    ax1.set_xlabel("γ"); ax1.set_ylabel("memory lifetime τ")
    ax1.set_title("Retention vs drive\n(weaker drive → shorter memory)")
    ax1.legend(fontsize=8); ax1.grid(alpha=0.25, which="both")

    for c, (om, rs) in zip(colors, sorted(by_omega.items())):
        g = [r["gamma"] for r in rs]
        ax2.plot(g, [r["sigma"] for r in rs], "s-", color=c, label=f"σ, Ω={om}")
    ax2.set_xlabel("γ"); ax2.set_ylabel("stationary dissipation rate σ")
    ax2.set_title("σ → 0 at BOTH limits\n(detailed balance; and flux collapse)")
    ax2.legend(fontsize=8); ax2.grid(alpha=0.25)

    ax3.set_xlabel("affinity A = −3 ln γ  [k_B T per cycle]")
    ax3.set_ylabel("memory lifetime τ")
    ax3.set_title("HEADLINE: retention is bought by DRIVE\n(not by dissipation rate)")
    ax3.legend(fontsize=8); ax3.grid(alpha=0.25, which="both")

    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote figure -> {out_path}")


def main():
    here = os.path.dirname(__file__)
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--omegas", type=int, nargs="+", default=[30, 60, 120])
    p.add_argument("--gammas", type=float, nargs="+",
                   default=[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.49])
    p.add_argument("--theta-frac", type=float, default=0.7)
    p.add_argument("--out", default=os.path.join(here, "dissipation_memory.png"))
    p.add_argument("--data", default=os.path.join(
        here, os.pardir, "results", "dissipation_memory.json"))
    p.add_argument("--ssa", action="store_true")
    p.add_argument("--max-seconds", type=float, default=60.0,
                   help="per-point wall-clock ceiling, ENFORCED via max_steps")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    print("Part B: exact cost of remembering (CME; SSA cross-checks are Plan 2)")
    rows = []
    for omega in args.omegas:
        print(f"-- Omega={omega} --")
        print(f"{'gamma':>6} {'A':>7} {'tau':>12} {'sigma':>10} "
              f"{'sigma*tau':>12} {'flip rate':>11} valid")
        for g in args.gammas:
            if g >= GAMMA_C:
                continue
            r = run_point(g, omega, args.theta_frac)
            rows.append(r)
            print(f"{r['gamma']:>6.3f} {r['affinity']:>7.3f} {r['tau']:>12.4g} "
                  f"{r['sigma']:>10.4g} {r['sigma_tau']:>12.4g} "
                  f"{r['flip_rate']:>11.4g} {r['valid']}")

    dropped = [(r["omega"], r["gamma"]) for r in rows if not r["valid"]]
    if dropped:
        print(f"\nDROPPED (solve invalid, e.g. tau<0 or >1e14): {dropped}")

    if args.ssa:
        print(f"\nSSA cross-check of tau (budget {args.max_seconds:.0f}s/point; "
              "flip rate compared against 1/tau -- see ssa_flip_rate docstring). "
              "A ratio near 0.8 is the arm-vs-attractor offset, not an error.")
        print(f"{'Omega':>6} {'gamma':>6} {'tau CME':>11} {'tau SSA':>11} "
              f"{'flips':>6} {'ratio':>7}")
        for r in rows:
            if not r["valid"]:
                continue
            predicted = (r["tau"] * 30 * STEPS_PER_UNIT_TIME * r["omega"]
                         / STEPS_PER_SECOND)
            if predicted > args.max_seconds:
                print(f"{r['omega']:>6} {r['gamma']:>6.2f} {r['tau']:>11.4g} "
                      f"{'skipped':>11} {'-':>6} {'-':>7}  "
                      f"(~{predicted:.0f}s predicted)")
                continue
            s = ssa_flip_rate(r["gamma"], r["omega"], r["tau"], args.max_seconds,
                              args.seed, args.theta_frac)
            r["ssa"] = {k: v for k, v in s.items()}     # no wall-clock in the JSON
            ratio = s["tau_ssa"] / r["tau"] if s["enough"] else float("nan")
            print(f"{r['omega']:>6} {r['gamma']:>6.2f} {r['tau']:>11.4g} "
                  f"{s['tau_ssa']:>11.4g} {s['flips']:>6} {ratio:>7.2f}"
                  + ("" if s["enough"] else "  (too few flips: not a measurement)")
                  + ("  BUDGET HIT" if s["budget_hit"] else ""))

    os.makedirs(os.path.dirname(os.path.abspath(args.data)), exist_ok=True)
    with open(args.data, "w") as fh:
        json.dump({"theta_frac": args.theta_frac, "rows": rows}, fh, indent=2)
    print(f"wrote data -> {args.data}")
    make_figure(rows, args.out)


if __name__ == "__main__":
    main()
