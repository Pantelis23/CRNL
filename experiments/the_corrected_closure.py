"""§107 (T-CASC-w) -- §103's closure with all four of §106's corrections applied.

§106 left the arc with a precise scope statement and three known defects, two of which it
fixed on paper without folding them back into the model. This does that, and sweeps Omega to
see whether the closure ratio flattens.

THE FOUR CORRECTIONS.
  (1) INDEXING (§106.3). escape_rate(Om, x_up) is the rate of a stage whose UPSTREAM sits at
      x_up. §102's predict() keyed each stage to its OWN operating point. Stage 1 has no
      upstream: its rate is escape_rate(Om, r3), which equals the free gap exactly.
  (2) TWO-STATE OCCUPANCY (§106.2). A free stage has no absorbing boundary, so
      P(low at t) = pi_low * (1 - exp(-lambda t)), not 1 - exp(-lambda t). pi_low runs
      0.9057 -> 0.5247 over Omega = 14-70.
  (3) THE AVERAGING (§102.1, §106's binding constraint). The escape rate must be averaged over
      the fluctuating input, and the correct fast-limit average is over the ACTION, not the
      input: with A(x) = -ln k(x)/Omega, exp(-Omega<A>) = exp(<ln k>) -- the GEOMETRIC mean of
      the rate. §102.1 used k(<x>), which is a different and smaller number.
  (4) p_transmit from §105's derived head start rather than §100's measured constant.

WHAT WAS SCOUTED, and is therefore not a prediction (rule 2). Correction (3) was checked
against §102.1's measured effective rate before this file was written: the geometric mean gives
1.0927 / 0.9372 / 0.9344 times k_eff at Omega = 14 / 30 / 55, against k(<x>)'s 0.822 / 0.757 /
0.792 and the frozen limit's 2.0 / 3.6 / 7.0. So the geometric mean is known to be the right
average to within ~9%, with no visible Omega drift. What is NOT known is what it does to the
closure.

PREDICTIONS.

  P1  WIRING. With correction (1), stage 1's rate must equal the free spectral gap to solver
      tolerance at every Omega. Stage 1 has no upstream and therefore no averaging: correction
      (3) must not touch it.

  P2  THE TEST. §106's closure ratio ran 1.1845 -> 1.8158 over Omega = 14-70, trending
      monotonically away from 1. **I predict the corrected model lands inside [0.7, 1.4] at
      every Omega AND stops trending monotonically.** The second half matters more than the
      first: a model that is uniformly wrong by a constant factor is a model with a missing
      prefactor, while one that drifts with Omega is a model with a missing mechanism. Judged
      on the trend, not on any single cell (rule 20).

  P3  THE ABLATION, so credit lands where it is due. Report the closure ratio under
      uncorrected, +indexing, +two-state, +geometric, +derived p_transmit, cumulatively. I
      predict (3) carries most of the improvement, because §106 localised the entire residual
      to `pure` -- stage 2's own escape -- and (3) is the only correction that touches it.
      (2) should be nearly invisible in the ratio, since §106 showed pi_low cancels.

  P4  DIRECTION. The geometric mean exceeds the rate at the mean, so `pure` grows and the ratio
      FALLS. §106's ratios are all above 1, so this must move them toward 1 rather than past
      it. If the corrected ratio undershoots badly at small Omega, the averaging has been
      over-applied and the fast limit is not the right one there -- which is testable, since
      §102.1's position was furthest from the fast end at Omega = 14.

WHAT THIS CANNOT SETTLE. The geometric mean is the fast-limit average. §102.1 measured the
system at position 0.09-0.22 from the fast end, i.e. NEAR but not AT that limit, and the
position drifts with Omega. So a residual of order the position is expected and is not a
defect; what would be a defect is a residual that grows with Omega.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

import experiments.chemical_cascade as cc
from experiments.chain_without_a_joint_solve import chain_operating_points
from experiments.depth_compounding import C, R2, R3
from experiments.escape_accounts_for_it import escape_rate
from experiments.jensen_shift import lna_width
from experiments.predicting_transmission import descent_rate
from experiments.the_head_start import bistability_edge, p_headstart, rate_weighted_head_start

P_TRANSMIT_MEASURED = 0.9376


def pi_low(om, x_up, cap_mult=1.25):
    """Stationary weight below the saddle for a stage whose upstream sits at x_up."""
    cap = int(np.ceil(cap_mult * R3 * om))
    lp = np.zeros(cap + 1)
    acc = 0.0
    for n in range(1, cap + 1):
        l, _ = cc.rates_stage(float(n - 1), x_up * om, om, C, R3, False, "hill")
        _, u = cc.rates_stage(float(n), x_up * om, om, C, R3, False, "hill")
        acc += np.log(l) - np.log(u)
        lp[n] = acc
    w = np.exp(lp - lp.max())
    w /= w.sum()
    return float(w[np.arange(cap + 1) < R2 * om].sum())


def input_law(om, mu, cap_mult=1.25):
    """The LNA law of a stage sitting at mu -- the predicted input to the next stage."""
    cap = int(np.ceil(cap_mult * R3 * om))
    xs = np.arange(int(np.ceil(R2 * om)) + 1, cap + 1) / om
    sd = lna_width(mu, om)
    w = np.exp(-0.5 * ((xs - mu) / sd) ** 2)
    return xs, w / w.sum()


def geometric_rate(om, mu):
    """exp(<ln k>) over the input law -- the action-averaged (fast-limit) escape rate."""
    xs, w = input_law(om, mu)
    ks = np.array([escape_rate(om, x) for x in xs])
    return float(np.exp((w * np.log(ks)).sum()))


def closure(om, t=2.0, indexing=False, two_state=False, geometric=False, derived_pt=False):
    """§103's closure with any subset of §106/§107's corrections switched on."""
    mus, _ = chain_operating_points(om, 2)
    if indexing:
        inputs = [R3, mus[0]]              # stage i is driven by stage i-1; stage 1 by the rail
    else:
        inputs = [mus[0], mus[1]]          # as coded in §102

    ks = []
    for i, x in enumerate(inputs):
        if geometric and i > 0:            # stage 1 has no upstream, so no averaging (P1)
            ks.append(geometric_rate(om, x))
        else:
            ks.append(escape_rate(om, x))

    if two_state:
        los = [pi_low(om, x) * (1 - np.exp(-k * t)) for k, x in zip(ks, inputs)]
    else:
        los = [1 - np.exp(-k * t) for k in ks]

    if derived_pt:
        kd, _ = descent_rate(om)
        p_t = p_headstart(kd, t, rate_weighted_head_start(om, bistability_edge()))
    else:
        p_t = P_TRANSMIT_MEASURED

    contam = los[0] * p_t
    pure = (1 - los[0]) * los[1]
    return contam / pure, ks, p_t


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/the_corrected_closure.json"))
    args = ap.parse_args()

    src = json.loads(pathlib.Path("results/where_the_expansion_frays.json").read_text())
    meas = {r["omega"]: r["meas_ratio"] for r in src["cells"]}

    from experiments.what_reflection_costs import spectral_gap
    print("P1 -- with the indexing fixed, stage 1's rate is the free gap")
    for om in sorted(meas):
        _, ks, _ = closure(om, indexing=True, geometric=True)
        g, _ = spectral_gap(om, False)
        assert abs(ks[0] / g - 1) < 1e-9, (om, ks[0], g)
    print("  exact at every Omega, and the geometric average does not touch it")

    stages = [
        ("uncorrected (§103)", {}),
        ("+ indexing", {"indexing": True}),
        ("+ two-state", {"indexing": True, "two_state": True}),
        ("+ geometric", {"indexing": True, "two_state": True, "geometric": True}),
        ("+ derived p_t", {"indexing": True, "two_state": True, "geometric": True,
                           "derived_pt": True}),
    ]

    print("\nP2/P3 -- closure ratio pred/meas, corrections applied cumulatively")
    header = f"{'Omega':>7}" + "".join(f"{n:>21}" for n, _ in stages)
    print(header)
    rows = []
    for om in sorted(meas):
        vals = [closure(om, **kw)[0] / meas[om] for _, kw in stages]
        rows.append({"omega": om, "ratios": vals})
        print(f"{om:>7}" + "".join(f"{v:>21.4f}" for v in vals), flush=True)
        args.out.write_text(json.dumps(rows, indent=2, default=float))

    print("\n  trend of each column across Omega (max/min, and monotone?)")
    for j, (name, _) in enumerate(stages):
        col = [r["ratios"][j] for r in rows]
        mono = col == sorted(col) or col == sorted(col, reverse=True)
        print(f"  {name:>21}: span {max(col)/min(col):.3f}x   "
              f"{'MONOTONE' if mono else 'not monotone'}   "
              f"in [0.7,1.4] at {sum(1 for v in col if 0.7 <= v <= 1.4)}/{len(col)} cells")

    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
