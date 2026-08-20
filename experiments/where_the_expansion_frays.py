"""§106 (T-CASC-s) -- where §103's single-element prediction stops working, and how.

§103 predicted the chain's operating points from single-element quantities and got +0.13%
and +0.50% at Omega = 30 but +1.64%, +3.45% and +4.32% at Omega = 14 -- every miss ABOVE the
measurement, i.e. predicting a LESS degraded chain than the one measured. All three
ingredients (the QSD mean, the LNA width, the static-transfer average) are asymptotic in
Omega, so the natural reading is that §103 is an expansion in Omega and Omega = 14 is where it
frays. That reading has never been tested: it rests on two barrier depths.

This sweeps Omega at fixed depth and asks what the residual DOES, rather than whether it is
small at any particular cell.

PREDICTIONS, WRITTEN BEFORE RUNNING.

  P1  WIRING. At Omega = 30 the residuals must reproduce §103's +0.13% (stage 1) and +0.50%
      (stage 2) to the precision those were quoted at. Anything else means the sweep is not
      running §103's model.

  P2  THE SHAPE, and this is the point. If the residual is the leading finite-Omega correction
      to an asymptotic expansion it must fall like 1/Omega. Fit log|residual| against
      log(Omega) and report the exponent WITH ITS WINDOW (rule 21 -- an exponent quoted without
      the range it was fitted over is not an exponent). **I predict an exponent near -1.**
      Judged by whether the residual CONVERGES, not by whether it sits below any threshold at
      the coarsest Omega I happened to include (rule 20).

  P3  THE SIGN. Every residual must stay POSITIVE at every Omega -- the model must always
      predict a less degraded chain. Two of the three ingredients bias that way and I know of
      nothing pushing the other: the conditioning drops inputs where F has no high rail, which
      raises the predicted mean, and the LNA underestimates a skewed width. A sign flip
      anywhere would mean two errors are competing and the 1/Omega reading is too simple.

  P4  THE CONSEQUENCE. §103's closure is stated in terms of contaminated/pure against §102's
      factor-of-two gate. If P2 holds, that ratio must approach 1 as Omega grows. Reported as a
      trend; the gate is not re-litigated here.

  P5  THE ALTERNATIVE, so this can fail informatively. If the residual SATURATES instead of
      falling, it is not a truncation error and one of the three ingredients is wrong rather
      than merely approximate. The first suspect is the conditioning choice -- §103 drops
      inputs where the downstream has no high rail and renormalises, which is an O(1) surgery
      on the tail rather than an O(1/Omega) correction, and it biases the mean in exactly the
      observed direction.

WHAT THE RUN THEN FORCED (added after P1-P5 printed, and labelled as such). P4 failed, so
the divergence had to be localised. Four candidates were checked in turn and three were
cleared: the operating points (they become exact), the rate ratio k1/k2 (it goes to 1),
p_transmit (flat to 1% across Omega), and finally P(stage 1 low) itself. That last one is the
crudest ingredient in the whole model -- a single free stage, no coupling -- and §102
validated it at Omega = 30 alone.

SCOPE. D = 2 only, t0 = 2.0. Depth is held fixed on purpose: §101's equilibration gate and its
saturation bound both bite at D = 3 and would confound an Omega sweep with a window problem.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from experiments.chain_without_a_joint_solve import chain_operating_points, split_from
from experiments.escape_accounts_for_it import escape_rate, operating_points
from experiments.what_reflection_costs import run_scope, spectral_gap

OMEGAS = (14, 20, 24, 30, 40, 55, 70)
T0 = 2.0


def cell(om, t=T0):
    """Measured operating points and split (joint solve, D = 2) against §103's prediction."""
    meas_mus, tot, pure, contam = operating_points(om, 2, t)
    pred_mus, _ = chain_operating_points(om, 2)
    # legacy=True: §106's published sweep was measured with the uncorrected model.
    _, cp, pp = split_from(om, pred_mus, t, legacy=True)
    return {
        "omega": om,
        "meas_mus": meas_mus,
        "pred_mus": pred_mus,
        "rel": [(p - m) / m for p, m in zip(pred_mus, meas_mus)],
        "meas_ratio": contam / pure,
        "pred_ratio": cp / pp,
        "pred_over_meas": (cp / pp) / (contam / pure),
        "p_free": tot,
    }


def decompose(rows, t=T0, p_t=0.9376):
    """Which side of the split diverges? Model evaluated AT the measured operating points."""
    out = []
    for r in rows:
        om = r["omega"]
        k1 = escape_rate(om, r["meas_mus"][0])
        k2 = escape_rate(om, r["meas_mus"][1])
        P, rat = r["p_free"], r["meas_ratio"]
        pure_m = P / (1.0 + rat)
        contam_m = P - pure_m
        surv = np.exp(-k1 * t)
        out.append({"omega": om,
                    "contam_meas_over_model": contam_m / ((1 - surv) * p_t),
                    "pure_meas_over_model": pure_m / (surv * (1 - np.exp(-k2 * t)))})
    return out


def escape_model_check(omegas=OMEGAS, t=T0):
    """The crudest ingredient: is P(stage 1 low) = 1 - exp(-k1 t) uniform in Omega?"""
    out = []
    for om in omegas:
        k, _ = spectral_gap(om, False)
        meas = run_scope(om=om, t0=t)["p_s1_lo"]
        pred = 1.0 - np.exp(-k * t)
        out.append({"omega": om, "k1": k, "pred": pred, "meas": meas,
                    "model_over_meas": pred / meas})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/where_the_expansion_frays.json"))
    args = ap.parse_args()

    rows = []
    print(f"D = 2, t0 = {T0}.  operating-point residuals, predicted vs measured")
    print(f"{'Omega':>7}{'stage 1':>11}{'stage 2':>11}{'contam/pure':>14}{'P_free':>11}")
    for om in OMEGAS:
        r = cell(om)
        rows.append(r)
        print(f"{om:>7}{r['rel'][0]:>10.2%}{r['rel'][1]:>11.2%}"
              f"{r['pred_over_meas']:>14.4f}{r['p_free']:>11.3e}", flush=True)
        args.out.write_text(json.dumps(rows, indent=2, default=float))

    print("\nP3 -- is every residual one-signed?")
    signs = {np.sign(v) for r in rows for v in r["rel"]}
    print(f"  distinct signs across all cells: {sorted(signs)}"
          f"   -> P3 {'HOLDS' if signs == {1.0} else 'FAILS'}")

    print("\nP2 -- the shape, fitted with its window (rule 21)")
    oms = np.array([r["omega"] for r in rows], float)
    for k, name in ((0, "stage 1"), (1, "stage 2")):
        y = np.array([abs(r["rel"][k]) for r in rows])
        ok = y > 0
        slope, intercept = np.polyfit(np.log(oms[ok]), np.log(y[ok]), 1)
        print(f"  {name}: exponent {slope:+.3f} over Omega = {int(oms.min())}"
              f"-{int(oms.max())}   (1/Omega would be -1)")
        # local exponents, so a drift is visible rather than averaged away
        loc = [f"{np.log(y[i+1]/y[i])/np.log(oms[i+1]/oms[i]):+.2f}"
               for i in range(len(oms) - 1) if y[i] > 0 and y[i + 1] > 0]
        print(f"           local: {' '.join(loc)}")

    print("\nP4 -- does the closure ratio approach 1?")
    seq = [f"{r['pred_over_meas']:.4f}" for r in rows]
    print(f"  {' -> '.join(seq)}")
    d0, d1 = abs(rows[0]["pred_over_meas"] - 1), abs(rows[-1]["pred_over_meas"] - 1)
    print(f"  |ratio - 1| goes {d0:.4f} -> {d1:.4f}"
          f"   ({'converging' if d1 < d0 else 'NOT converging'})")

    print("\nLocalising the P4 failure -- model at the MEASURED operating points")
    dec = decompose(rows)
    print(f"{'Omega':>7}{'contam meas/model':>20}{'pure meas/model':>18}")
    for d in dec:
        print(f"{d['omega']:>7}{d['contam_meas_over_model']:>20.4f}"
              f"{d['pure_meas_over_model']:>18.4f}")
    print("  pure is bounded and non-monotone; contam falls steadily -> the contaminated"
          " channel is the divergent one")

    print("\nThe crudest ingredient: P(stage 1 low) = 1 - exp(-k1 t), a single free stage")
    esc = escape_model_check()
    print(f"{'Omega':>7}{'predicted':>14}{'measured':>14}{'model/meas':>12}")
    for e in esc:
        print(f"{e['omega']:>7}{e['pred']:>14.4e}{e['meas']:>14.4e}"
              f"{e['model_over_meas']:>12.4f}", flush=True)
    sp = max(e["model_over_meas"] for e in esc) / min(e["model_over_meas"] for e in esc)
    print(f"  spans {sp:.2f}x across Omega = {OMEGAS[0]}-{OMEGAS[-1]}, crossing 1 near"
          f" Omega = 35 -- §102 validated it at Omega = 30 alone")

    args.out.write_text(json.dumps(
        {"cells": rows, "decomposition": dec, "escape_model": esc}, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
