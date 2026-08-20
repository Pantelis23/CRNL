"""§103 (T-CASC-r) -- the depth ceiling from one element's rate functions, no joint solve.

§102 showed that escape rates AT the operating points account for §101's contaminated/pure
split to tens of percent. But it read those operating points off §101's joint CME, so it did
not show the chain is computable without one. This closes that gap, or fails to.

Every ingredient here is single-element and 1-D:
  * stage 1's law            -- upstream_qsd(om), the QSD of ONE free stage above its saddle
  * stage k+1's mean         -- <F(x_k)>, the exact static-transfer average (§95/§96)
  * stage k+1's width        -- lna_width(mean, om), LNA about the shifted point (§96)
  * every escape rate        -- escape_rate(om, x_up), ONE stage with its upstream pinned (§102)
Nothing below touches a joint generator. The joint solves appear only on the MEASURED side,
as §101's stored numbers.

A MODELLING CHOICE, stated because it is not forced. F(x_up) is NaN where the downstream has
no high rail at all -- below that input the downstream is monostable-low and "the operating
point of a stage that has not failed" is not defined there. Those inputs are dropped and the
law renormalised. That is the same conditioning §101 applies on the measured side (counts
above the saddle), so the two are comparable, but it is an approximation on both sides.

PREDICTIONS, WRITTEN BEFORE RUNNING.

  P1  STAGE 1. The predicted mean is the QSD mean of one free stage; the measured one is the
      free chain's stage-1 mean at t0 = 2.0 conditioned on it not having failed. These are
      the same object computed two ways and differ only because the chain has not fully
      relaxed to the QSD. Predict agreement WITHIN 1%: 3.0280 (Om=30), 2.8042 (Om=14).

  P2  STAGES 2 AND 3. §96 predicted the reflected chain's operating point to 0.12%, but the
      conditioning here differs (a free stage's law is truncated by escape, a reflected one's
      is not), so I predict a looser 2%: measured 2.9759 (Om=30 D=2), 2.6698 and 2.6113
      (Om=14 D=3).

  P3  THE SPLIT, against §102's own pre-registered gate so the two are comparable. Feed those
      predicted operating points through §102's rate curve and compute contaminated/pure.
      **Predict it survives the factor-of-two gate, with a residual LARGER than §102's
      1.04-1.36.** The reason is P2 plus §102's P2: escape_rate is steeply convex -- ~9.5x
      across the operating range at Om=30 -- so a 2% error in an operating point is amplified
      into a much bigger error in a rate. I do NOT predict the sign of the residual; the
      conditioning drops inputs at the low end, which raises the predicted mean and so lowers
      the predicted rate, while ignoring the failure channel's back-action pushes the other
      way. Both directions are reported.

  P4  THE POINT OF THE WHOLE THING. If P3 holds, then for this cascade the contaminated
      channel -- which §101 showed is the MAJORITY of the error by D = 3 -- is computable from
      a single element's rate functions, and the depth ceiling with it. If P3 fails, the
      operating points genuinely need the joint law and §102's agreement was borrowing
      information from the joint solve it read them off.

WHAT THIS CANNOT SETTLE. p_transmit is still §100's measured 0.9376 rather than a predicted
quantity, so one number in the model is empirical. It is a probability bounded by 1 and
measured between 0.73 and 0.98 (§100.2), so it cannot absorb an order of magnitude -- but the
model is not parameter-free and does not claim to be.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from experiments.depth_compounding import R2, R3
from experiments.escape_accounts_for_it import P_TRANSMIT, escape_rate
from experiments.jensen_shift import F, lna_width
from experiments.margin_law import upstream_qsd

# §101's measured free-chain operating points at t0 = 2.0 -- the target, not an input.
MEASURED = {
    (30, 2): {"mus": [3.0280, 2.9759], "contam": 7.2892e-03, "pure": 1.0875e-02},
    (14, 2): {"mus": [2.8042, 2.6698], "contam": 1.3825e-01, "pure": 1.5383e-01},
    (14, 3): {"mus": [2.8042, 2.6698, 2.6113], "contam": 2.7544e-01, "pure": 1.2520e-01},
    # §101 stored the split for this cell but not its stage means, so it tests P3 only.
    # It is the D = 3 cell at the barrier where the Omega-expansion is healthy.
    (30, 3): {"mus": None, "contam": 1.6518e-02, "pure": 1.0409e-02},
}


def lattice(om, cap_mult=1.25):
    """The concentrations a stage can occupy above its saddle."""
    cap = int(np.ceil(cap_mult * R3 * om))
    n = np.arange(int(np.ceil(R2 * om)) + 1, cap + 1)
    return n / om


def transfer(xs, px):
    """<F(x)> over a law, dropping inputs where the downstream has no high rail."""
    fs = np.array([F(x) for x in xs])
    ok = np.isfinite(fs)
    if not ok.any():
        return np.nan, 0.0
    w = px[ok] / px[ok].sum()
    return float((w * fs[ok]).sum()), float(px[ok].sum())


def chain_operating_points(om, D):
    """Every stage's operating point from single-element quantities only.

    THE INTRINSIC TERM, and §103.1 is why it is here. The map is neutral at the rail --
    F(r3) = r3 exactly, which is the condition §91 built the coupling to satisfy -- so r3 is
    a FIXED POINT of F and iterating <F(x)> alone can only converge toward the rail. It can
    never degrade, and the measured chain degrades. What supplies the degradation is each
    stage's own depression below its deterministic rail at finite Omega, `d_intr`, which §96
    adds to <F> and the first version of this file omitted.

    §96 validated d_intr at D = 2 only. Applying it once per stage is the natural
    generalisation and is NOT separately validated here; §103's P2 is the test of it.
    """
    xs, px = upstream_qsd(om)
    mus = [float((px * xs).sum())]
    d_intr = mus[0] - R3                     # stage 1's own depression below its rail
    cur_xs, cur_px = xs, px
    kept = [1.0]
    for _ in range(1, D):
        mu, frac = transfer(cur_xs, cur_px)
        mu += d_intr
        kept.append(frac)
        sd = lna_width(mu, om)
        cur_xs = lattice(om)
        w = np.exp(-0.5 * ((cur_xs - mu) / sd) ** 2)
        cur_px = w / w.sum()
        mus.append(mu)
    return mus, kept


def split_from(om, mus, t, p_transmit=P_TRANSMIT):
    ks = [escape_rate(om, x) for x in mus]
    surv = float(np.prod([np.exp(-k * t) for k in ks[:-1]]))
    return ks, (1.0 - surv) * p_transmit, surv * (1.0 - np.exp(-ks[-1] * t))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/chain_without_a_joint_solve.json"))
    args = ap.parse_args()
    t = 2.0
    out = []

    print("P1/P2 -- operating points, predicted from single-element quantities")
    print(f"{'Om':>4}{'D':>3}{'stage':>7}{'predicted':>12}{'measured':>11}{'rel err':>10}")
    for (om, D), ref in sorted(MEASURED.items()):
        mus, kept = chain_operating_points(om, D)
        if ref["mus"] is None:
            print(f"{om:>4}{D:>3}{'--':>7}{'':>12}{'(means not stored; P3 only)':>38}")
        else:
            for i, (p, m) in enumerate(zip(mus, ref["mus"])):
                print(f"{om:>4}{D:>3}{i+1:>7}{p:>12.4f}{m:>11.4f}{(p-m)/m:>9.2%}")
        out.append({"omega": om, "D": D, "pred_mus": mus, "meas_mus": ref["mus"],
                    "kept": kept})

    print("\nP3 -- the split, against §102's pre-registered factor-of-two gate")
    print(f"{'Om':>4}{'D':>3}{'contam/pure pred':>18}{'measured':>11}{'pred/meas':>11}"
          f"{'§102 was':>10}")
    was = {(30, 2): 1.359, (14, 2): 1.159, (14, 3): 1.040, (30, 3): float("nan")}
    for r in out:
        om, D = r["omega"], r["D"]
        ref = MEASURED[(om, D)]
        ks, cp, pp = split_from(om, r["pred_mus"], t)
        pred, meas = cp / pp, ref["contam"] / ref["pure"]
        r.update({"ks": ks, "pred_ratio": pred, "meas_ratio": meas,
                  "pred_over_meas": pred / meas})
        gate = "" if 0.5 < pred / meas < 2.0 else "   OUTSIDE THE GATE"
        print(f"{om:>4}{D:>3}{pred:>18.4f}{meas:>11.4f}{pred/meas:>11.4f}"
              f"{was[(om, D)]:>10.3f}{gate}")

    inside = [0.5 < r["pred_over_meas"] < 2.0 for r in out]
    worse = [abs(np.log(r["pred_over_meas"])) > abs(np.log(was[(r["omega"], r["D"])]))
             for r in out if np.isfinite(was[(r["omega"], r["D"])])]
    print(f"\n  inside the gate: {sum(inside)}/{len(out)}"
          f"   residual larger than §102's: {sum(worse)}/{len(worse)}")
    print("  (P3 predicted both: survives the gate, with a larger residual than §102)")

    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
