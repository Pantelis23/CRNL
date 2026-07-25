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
                  target_flips: int = 30, n_seeds: int = 8) -> dict:
    """Measure the flip rate by SSA and compare with 1/tau from the CME.

    GATED ON THE CME, BY DESIGN: the exact tau prices the run BEFORE paying for
    it. And the budget is ENFORCED, not merely accepted -- an earlier draft took
    max_seconds as an argument and never used it, so raising it to 3600 admitted
    a 1.8 h run with max_steps as the only backstop.

    CONVENTION: FlipCounter counts ONE-WAY crossings, so flips/T -> 1/tau.
    Using 1/(2 tau) (the ROUND-TRIP rate) halves tau_ssa and was measured to
    give ratios of 0.37-0.58 with 52-82 flips -- a false alarm outside sampling
    An arm-vs-attractor offset was PREDICTED (the flip clock starts at the arm,
    the CME's tau at the attractor) and estimated at 0.77-0.85 from direct
    MFPT-from-arm runs. Seed-averaged measurement does NOT show it: six points
    give 0.946/0.980/0.964/1.043/1.038/0.868, mean 0.97, every one within ~1.5
    SEM of 1.0. The prediction confused two different quantities -- an MFPT from
    the arm is not the mean time between crossings of a long stationary
    trajectory, which is dominated by the full dwell near the attractor. Expect
    ~1.0; a systematic 0.8 would now be evidence of a real problem.

    AVERAGES OVER n_seeds INDEPENDENT TRAJECTORIES, and this is not optional.
    A single trajectory carrying ~30 flips has a measured spread of sd = 0.16 to
    0.32 in the ratio tau_SSA/tau_CME, so an individual point lands anywhere in
    [0.64, 1.58] with no bug present at all -- a first run of this sweep produced
    exactly two such outliers (1.50 at Omega=60 gamma=0.45, 1.58 at Omega=120
    gamma=0.49) and they were pure sampling noise: the SAME points over 8 seeds
    give 0.950 +- 0.091 and 1.047 +- 0.112. Reporting one trajectory would mean
    publishing the tail of a distribution as if it were a measurement, and would
    make any real factor-level discrepancy indistinguishable from luck.
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
    # hard cap from the budget, SHARED across the seeds so the total run cannot
    # outlive it -- a per-seed cap would silently cost n_seeds times the budget.
    max_steps = max(1, int(max_seconds * STEPS_PER_SECOND / n_seeds))

    taus, flips_all, steps_all, budget_hit = [], [], [], False
    for s in range(n_seeds):
        r = gillespie_instrumented(comp, start, seed_for(omega, s, base=seed),
                                   pairing, flip_arm=arm, omega=omega,
                                   t_max=t_needed, max_steps=max_steps,
                                   species=list(net.species))
        flips_all.append(r.flips)
        steps_all.append(r.steps)
        budget_hit |= r.steps >= max_steps
        if r.flips > 0 and r.t_final > 0:
            taus.append(r.t_final / r.flips)          # = 1 / rate

    total_flips = int(sum(flips_all))
    if not taus:
        nan = float("nan")
        return {"flips": total_flips, "n_seeds": n_seeds, "tau_ssa": nan,
                "tau_ssa_sem": nan, "flip_rate": nan, "enough": False,
                "budget_hit": budget_hit, "t_target": t_needed,
                "steps": int(sum(steps_all))}
    taus = np.array(taus)
    mean = float(taus.mean())
    sem = float(taus.std(ddof=1) / np.sqrt(len(taus))) if len(taus) > 1 else 0.0
    return {
        "flips": total_flips, "n_seeds": n_seeds, "seeds_used": len(taus),
        "steps": int(sum(steps_all)),
        "tau_ssa": mean, "tau_ssa_sem": sem,
        "flip_rate": 1.0 / mean if mean > 0 else float("nan"),
        # 10 flips was the old single-trajectory bar; with n_seeds trajectories
        # the meaningful quantity is the pooled count.
        "enough": total_flips >= 60 and len(taus) >= 3,
        "budget_hit": budget_hit, "t_target": t_needed,
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
        # ONE-WAY crossing rate, matching FlipCounter's convention and the SSA
        # cross-check above. 1/(2 tau) is the ROUND-TRIP rate; using it here made
        # the exact and sampled columns disagree by exactly 2 while both were
        # internally correct.
        "flip_rate": 1.0 / tau if fp["valid"] and tau > 0 else float("nan"),
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
    p.add_argument("--ssa-seeds", type=int, default=8,
                   help="independent trajectories per SSA point. NOT optional at "
                        "1: a single trajectory has sd ~0.26 in the ratio, so it "
                        "publishes the tail of a distribution as a measurement")
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
              "Expect ratio ~1.0; the predicted arm-vs-attractor offset does "
              "not survive seed-averaged measurement.")
        print(f"{'Omega':>6} {'gamma':>6} {'tau CME':>11} {'tau SSA':>20} "
              f"{'flips':>6} {'ratio':>16}")
        for r in rows:
            if not r["valid"]:
                continue
            # x n_seeds: the budget must price ALL the trajectories, not one
            predicted = (r["tau"] * 30 * STEPS_PER_UNIT_TIME * r["omega"]
                         / STEPS_PER_SECOND)
            if predicted > args.max_seconds:
                print(f"{r['omega']:>6} {r['gamma']:>6.2f} {r['tau']:>11.4g} "
                      f"{'skipped':>11} {'-':>6} {'-':>7}  "
                      f"(~{predicted:.0f}s predicted)")
                continue
            s = ssa_flip_rate(r["gamma"], r["omega"], r["tau"], args.max_seconds,
                              args.seed, args.theta_frac, n_seeds=args.ssa_seeds)
            r["ssa"] = {k: v for k, v in s.items()}     # no wall-clock in the JSON
            ratio = s["tau_ssa"] / r["tau"] if s["enough"] else float("nan")
            r_sem = s["tau_ssa_sem"] / r["tau"] if s["enough"] else float("nan")
            print(f"{r['omega']:>6} {r['gamma']:>6.2f} {r['tau']:>11.4g} "
                  f"{s['tau_ssa']:>13.4g}±{s['tau_ssa_sem']:<6.3g} "
                  f"{s['flips']:>6} {ratio:>9.3f}±{r_sem:<6.3f}"
                  + ("" if s["enough"] else "  (too few flips: not a measurement)")
                  + ("  BUDGET HIT" if s["budget_hit"] else ""))

    os.makedirs(os.path.dirname(os.path.abspath(args.data)), exist_ok=True)
    with open(args.data, "w") as fh:
        json.dump({"theta_frac": args.theta_frac, "rows": rows}, fh, indent=2)
    print(f"wrote data -> {args.data}")
    make_figure(rows, args.out)


if __name__ == "__main__":
    main()
