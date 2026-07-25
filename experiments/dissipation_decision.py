"""Part A: what does a decision cost, in free energy?

Sweep gamma over the bistable range and, for each, solve the exact CME for a
biased start run to |delta| >= theta_decide(gamma). Report the error probability
against the dissipation, decomposed exactly (docs/.../dissipation-design.md 2.5):

    dS_total = ln[W(n_stop)/W(n_0)]  +  (A/3) * (M_forward - M_reverse)
               ^^^ boundary term         ^^^ cycle term

Two independently interpretable Hill/Schnakenberg terms. The boundary term is
small and nearly gamma-independent; the cycle term carries the physics, and is
the x-axis of the headline plot.

TWO PROTOCOL POINTS THAT ARE NOT OPTIONAL, both learned the hard way, and both
with DIFFERENT fixes:

  * The decision THRESHOLD must scale with the landscape. A FIXED theta = 0.5
    is unreachable above gamma ~ 0.417 (delta*(0.49) = 0.187), which turns
    "deciding" into "fluctuating past the attractor" and inflates the measured
    dissipation by an order of magnitude -- a protocol artifact that looks
    exactly like physics. So theta is `theta_frac * delta_star(gamma)`.
  * The initial BIAS must not jitter freely on the integer lattice -- one
    molecule of bias is worth ~20 k_B of dissipation. So the bias is *also*
    expressed as a fraction of delta_star(gamma) (keeping the protocol
    difficulty comparable across gamma), but rounded to an exact integer
    count and the REALISED fraction is reported alongside the requested one,
    so any residual lattice jitter is visible rather than hidden.
  * The start is B(0) = 0 with the bias carried by the committed species, per
    design.md 4 and restoration_wall.py. The alternative ((1/3,1/3,1/3) + bias)
    differs by 36% in P(error).

Landauer: the cost is extensive in Omega and grows as ln(1/gamma), roughly
(A/3)*O(Omega) ~ 1.9*Omega*ln(1/gamma) k_B T. So any "orders of magnitude above
k_B T ln 2" statement is an instance at a stated (Omega, gamma), not a bound --
Landauer does not apply to this protocol; the point is the scale gap.

    python -m experiments.dissipation_decision
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
    initial_counts,
    reverse_pairing,
    theta_decide,
)
from crnl.thermo import decompose, ln_multinomial


def run_gamma(gamma: float, omega: int, bias_frac: float, theta_frac: float) -> dict | None:
    """One exact solve at fixed gamma.

    Both the bias and the threshold are expressed as fractions of
    delta_star(gamma) and then quantised to integer molecule counts -- the
    landscape shrinks as gamma rises, so a fixed count for either one either
    sits outside the landscape (threshold) or jitters uncontrollably on the
    integer lattice relative to it (bias). Returns None (skip) rather than a
    row when the quantised protocol is degenerate at this gamma: threshold
    outside the landscape, or bias/threshold collision (already-absorbed
    start).
    """
    d_star = delta_star(gamma)
    theta_counts = max(2, int(round(theta_frac * d_star * omega)))
    bias_counts = max(1, int(round(bias_frac * d_star * omega)))

    if not (bias_counts < theta_counts):
        print(f"SKIPPED gamma={gamma:.3f}: bias_counts={bias_counts} >= "
              f"theta_counts={theta_counts} (start already absorbed)")
        return None
    if not (theta_counts / omega < d_star):
        print(f"SKIPPED gamma={gamma:.3f}: theta={theta_counts / omega:.4f} "
              f">= delta_star={d_star:.4f} (threshold outside the landscape)")
        return None

    net = am_reversible(gamma)
    pairing = reverse_pairing(net)
    A = cycle_affinity(net, pairing)
    n0 = initial_counts(omega, gamma, count_diff=bias_counts)

    def absorbing(n):
        return abs(int(n[0]) - int(n[1])) >= theta_counts

    fp = first_passage(net, omega, float(omega), n0, absorbing, pairing)
    sigma = ep_rate(net, omega, float(omega), pairing)

    # Use the EXACT expected boundary term from the solve. A hand-built
    # representative stopping state gets the sign wrong, because real absorbing
    # states carry a large blank population that an n_B = 0 state does not.
    dec = decompose(n0, None, fp["net_reaction_firings"], A,
                    boundary=fp["boundary"])
    return {
        "gamma": gamma,
        "affinity": A,
        "theta_counts": theta_counts,
        "theta": theta_counts / omega,
        "delta_star": d_star,
        "theta_over_delta_star": (theta_counts / omega) / d_star,
        "bias_counts": bias_counts,
        "realised_bias": float((n0[0] - n0[1]) / omega),
        "bias_over_delta_star": float((n0[0] - n0[1]) / omega) / d_star,
        "p_error": 1.0 - fp["split"],
        "mean_time": fp["mean_time"],
        "sigma": sigma,
        "net_reaction_firings": fp["net_reaction_firings"],
        "boundary": dec["boundary"],
        "cycle": dec["cycle"],
        "total": dec["total"],
        "valid": fp["valid"],
        "residual": fp["residual"],
    }


def make_figure(rows, omega, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ok = [r for r in rows if r["valid"]]
    g = np.array([r["gamma"] for r in ok])
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 5))

    ax1.semilogy(g, [r["p_error"] for r in ok], "o-", color="#d62728")
    ax1.set_xlabel("γ"); ax1.set_ylabel("P(error)")
    ax1.set_title(f"Accuracy falls as the drive weakens (Ω={omega})")
    ax1.grid(alpha=0.25, which="both")

    ax2.plot(g, [r["cycle"] for r in ok], "o-", color="#1f77b4",
             label="cycle term (A/3)·⟨M⟩")
    ax2.plot(g, [r["boundary"] for r in ok], "s-", color="#7f7f7f",
             label="boundary term ln W ratio")
    ax2.plot(g, [r["total"] for r in ok], "^--", color="black", label="total ⟨ΔS⟩")
    ax2.set_xlabel("γ"); ax2.set_ylabel("entropy production  [k_B]")
    ax2.set_title("Exact decomposition\n(no cancellation, both terms meaningful)")
    ax2.legend(fontsize=8); ax2.grid(alpha=0.25)

    ax3.semilogy([r["cycle"] for r in ok], [r["p_error"] for r in ok],
                 "o-", color="#1f77b4")
    ax3.axvline(np.log(2), color="black", ls=":", lw=1, label="k_B T ln 2 (scale ref)")
    ax3.set_xlabel("dissipation, cycle term  [k_B T]")
    ax3.set_ylabel("P(error)")
    ax3.set_title("HEADLINE: what accuracy costs\n(free energy spent per decision)")
    ax3.legend(fontsize=8); ax3.grid(alpha=0.25, which="both")

    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote figure -> {out_path}")


def main():
    here = os.path.dirname(__file__)
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--omega", type=int, default=120)
    p.add_argument("--gammas", type=float, nargs="+",
                   default=[0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.49])
    p.add_argument("--bias-frac", type=float, default=0.2,
                   help="initial |n_X - n_Y| as a fraction of delta_star(gamma), "
                        "quantised to an integer molecule count per gamma")
    p.add_argument("--theta-frac", type=float, default=0.7,
                   help="decision threshold as a fraction of delta_star(gamma), "
                        "quantised to an integer molecule count per gamma")
    p.add_argument("--out", default=os.path.join(here, "dissipation_decision.png"))
    p.add_argument("--data", default=os.path.join(
        here, os.pardir, "results", "dissipation_decision.json"))
    args = p.parse_args()

    print(f"Part A: exact cost of a decision  (Omega={args.omega}, "
          f"bias_frac={args.bias_frac}, theta_frac={args.theta_frac}, "
          "both scaled per-gamma by delta_star(gamma))")
    results = [(g, run_gamma(g, args.omega, args.bias_frac, args.theta_frac))
               for g in args.gammas if g < GAMMA_C]
    skipped = [g for g, r in results if r is None]
    rows = [r for _, r in results if r is not None]

    print(f"{'gamma':>6} {'A':>7} {'th/d*':>7} {'bias/d*':>8} {'P(err)':>9} "
          f"{'<T>':>9} {'cycle':>10} {'bound':>8} {'total':>10} valid")
    for r in rows:
        print(f"{r['gamma']:>6.3f} {r['affinity']:>7.3f} "
              f"{r['theta_over_delta_star']:>7.3f} "
              f"{r['bias_over_delta_star']:>8.4f} {r['p_error']:>9.4g} "
              f"{r['mean_time']:>9.3g} {r['cycle']:>10.2f} {r['boundary']:>8.2f} "
              f"{r['total']:>10.2f} {r['valid']}")

    if skipped:
        print(f"\nSKIPPED (degenerate quantised protocol): gamma = {skipped}")

    dropped = [r["gamma"] for r in rows if not r["valid"]]
    if dropped:
        print(f"\nDROPPED (solve invalid): gamma = {dropped}")

    os.makedirs(os.path.dirname(os.path.abspath(args.data)), exist_ok=True)
    with open(args.data, "w") as fh:
        json.dump({"omega": args.omega, "bias_frac": args.bias_frac,
                   "theta_frac": args.theta_frac, "skipped_gammas": skipped,
                   "rows": rows}, fh, indent=2)
    print(f"wrote data -> {args.data}")
    make_figure(rows, args.omega, args.out)


if __name__ == "__main__":
    main()
