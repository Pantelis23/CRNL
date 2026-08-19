"""T-CASC-h: the missing term is CURVATURE, not gain -- Jensen through the transfer function

§94 found three readings of the width sequence and refuted all three. The sharpest anomaly was
not about depth at all: **the deterministic gain at the measured operating points is g^2 = 0.0109,
thirteen times too small to explain sigma_2/sigma_1 = 1.0855.** Something is missing at the FIRST
step, before any compounding.

**The candidate is curvature.** The stage-to-stage map x_out = F(x_up) -- the high root of the
downstream landscape given a pinned upstream -- is CONCAVE near the rail, because the Hill coupling
is close to saturation there. So by Jensen

    <F(x_up)>  <  F(<x_up>)      and the shortfall is  (1/2) F''(mu) sigma_up^2

which lowers the downstream's operating point. A lower operating point sits on a flatter part of
the landscape, so the downstream is WIDER -- and the effect is second order in sigma, which is
exactly why a first-order gain misses it.

**This is a mechanism (rule 17), so it gets an absolute test, not a fit.** F and its derivatives
come from the deterministic landscape; sigma_1 is measured; the prediction for mu_2 and sigma_2 is
then a number with nothing free in it, checked against §94's exact values 2.95635 and 0.54181.

PREDICTIONS, written before running.

  P1  GATE. F(r3) = r3 exactly (the coupling is neutral at the rail by construction, §91 P1(a)),
      and F', F'' must be converged in the finite-difference step -- reported across h, not at one
      h (rule 13: an approximation's own numerical parameter is a second axis).
  P2  **THE CONCAVITY.** F'' < 0 at the operating point. If F is convex there, Jensen pushes the
      other way and the mechanism is dead on arrival.
  P3  **THE MEAN, ABSOLUTE.** mu_2 = F(mu_1) + (1/2) F''(mu_1) sigma_1^2 + delta_intr, where
      delta_intr = mu_1 - r3 is the element's own stochastic shift, measured on stage 1 whose input
      is a noiseless chemostat. **Predicted: this lands near the measured 2.95635.** Nothing is
      fitted -- F comes from the roots, sigma_1 from §94's exact solve.
  P4  **THE WIDTH.** With the operating point shifted, the LNA width at that point predicts
      sigma_2. **Predicted: it accounts for most of the gap the linear gain missed**, i.e. lands
      far closer to 0.54181 than the gain-only 0.5019 does.
  P5  **HOW MUCH OF THE FACTOR OF 13.** Report the shortfall explicitly: gain-only, gain plus
      curvature, and measured. **If curvature closes only part of it, say what fraction and stop**
      -- a partial mechanism named as partial is worth more than a fitted one.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

import experiments.chemical_cascade as cc
from experiments.cascade_schlogl import schlogl_consts
from experiments.margin_law import R1, R2, R3

C = schlogl_consts(R1, R2, R3)
# §94's exact values, quoted as stored numbers (rule 16).
# **CORRECTED in §96.** The originals were MU1, SD1 = 3.02222, 0.49922 and MU2, SD2 = 2.95635,
# 0.54181, computed with a high-side filter that a REFLECTED stage does not need -- it merely
# dropped the boundary lattice site. Unconditioned, which is correct for a stage that cannot
# escape, the values are below. The originals stand in FINDINGS §94; these supersede them.
MU1, SD1 = 3.02117, 0.50120
MU2, SD2 = 2.95165, 0.54965


def F(x_up):
    """Downstream high rail given a pinned upstream -- the stage-to-stage map."""
    r = cc.downstream_roots(x_up, C, R3, "hill")
    return float(r[-1]) if len(r) >= 3 else np.nan


def dF(x, h):
    return (F(x + h) - F(x - h)) / (2 * h)


def d2F(x, h):
    return (F(x + h) - 2 * F(x) + F(x - h)) / h ** 2


def lna_width(x_rail, om):
    """LNA sd of a stage sitting at x_rail: sqrt(V/Omega), V = (lam+mu)/(2|f'|)."""
    k1a, k1r, k2b, k2r = C
    lam = k1a * x_rail ** 2 + k2b
    mu = k1r * x_rail ** 3 + k2r * x_rail
    fp = 2 * k1a * x_rail - 3 * k1r * x_rail ** 2 - k2r
    return float(np.sqrt((lam + mu) / (2 * abs(fp)) / om))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=30)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/jensen_shift.json"))
    args = ap.parse_args()
    om = args.omega
    cc.HILL_N, cc.HILL_K = 4.0, 1.0
    out = {}

    print("=== P1 GATE: neutrality at the rail, and are F', F'' converged in h?")
    print(f"  F(r3) = {F(R3):.10f}   r3 = {R3}   diff {abs(F(R3)-R3):.2e}")
    print(f"{'h':>10}{'F1(mu1)':>12}{'F2(mu1)':>12}")
    d1s, d2s = [], []
    for h in (1e-2, 5e-3, 2e-3, 1e-3, 5e-4):
        a, b = dF(MU1, h), d2F(MU1, h)
        d1s.append(a); d2s.append(b)
        print(f"{h:>10.0e}{a:>12.6f}{b:>12.6f}")
    conv = abs(d2s[-1] - d2s[-2]) < 0.02 * abs(d2s[-1])
    print(f"  -> P1 {'HOLDS: neutral at the rail and the second derivative is converged' if (abs(F(R3)-R3) < 1e-9 and conv) else 'FAILS'}")
    Fp, Fpp = d1s[-1], d2s[-1]

    print(f"\n=== P2: is the map CONCAVE at the operating point?")
    print(f"  F'(mu1) = {Fp:.6f}   F''(mu1) = {Fpp:.6f}")
    print(f"  -> P2 {'HOLDS: concave, so Jensen lowers the downstream operating point' if Fpp < 0 else 'FAILS: convex -- Jensen pushes the other way and the mechanism is dead'}")

    print("\n=== P3: the mean, absolute. Nothing fitted.")
    d_intr = MU1 - R3
    jensen = 0.5 * Fpp * SD1 ** 2
    pred_no_j = F(MU1) + d_intr
    pred = pred_no_j + jensen
    print(f"  element's own shift, from stage 1:  delta_intr = {d_intr:+.5f}")
    print(f"  F(mu1)                             = {F(MU1):.5f}")
    print(f"  Jensen term (1/2) F'' sigma_1^2    = {jensen:+.5f}")
    print(f"  mu_2 predicted WITHOUT Jensen      = {pred_no_j:.5f}")
    print(f"  mu_2 predicted WITH Jensen         = {pred:.5f}")
    print(f"  mu_2 measured                      = {MU2:.5f}")
    gap_no_j, gap = MU2 - pred_no_j, MU2 - pred
    frac = 1 - abs(gap) / abs(gap_no_j) if gap_no_j else np.nan
    print(f"  residual without Jensen {gap_no_j:+.5f};  with Jensen {gap:+.5f}"
          f"   -> closes {100*frac:.0f}% of it")
    out.update({"F_mu1": F(MU1), "Fp": Fp, "Fpp": Fpp, "d_intr": d_intr,
                "jensen": jensen, "pred": pred, "pred_no_j": pred_no_j,
                "gap": gap, "gap_no_j": gap_no_j})
    print(f"  -> P3 {'HOLDS: the curvature term moves the mean the right way and closes most of the gap' if frac > 0.5 else ('the curvature term has the right SIGN but closes only part of the gap -- reported as partial' if frac > 0 else 'FAILS: it moves the mean the WRONG way')}")

    print("\n=== P4/P5: the width, and how much of the factor of 13 this accounts for")
    s_gain = lna_width(F(MU1) + d_intr, om)
    s_shift = lna_width(pred, om)
    s_meas_pt = lna_width(MU2, om)
    print(f"{'operating point':>34}{'x':>10}{'LNA sd':>10}")
    print(f"{'gain only (no curvature)':>34}{pred_no_j:>10.5f}{s_gain:>10.5f}")
    print(f"{'with the Jensen shift':>34}{pred:>10.5f}{s_shift:>10.5f}")
    print(f"{'at the MEASURED mean':>34}{MU2:>10.5f}{s_meas_pt:>10.5f}")
    print(f"{'measured sigma_2':>34}{'':>10}{SD2:>10.5f}")
    print(f"  sigma_1 measured = {SD1:.5f};  LNA at stage 1's own mean ="
          f" {lna_width(MU1, om):.5f}")
    need = SD2 - s_gain
    got = s_shift - s_gain
    print(f"  gap to close from the gain-only width: {need:+.5f};"
          f"  curvature supplies {got:+.5f}  ({100*got/need:.0f}%)")
    out.update({"s_gain": s_gain, "s_shift": s_shift, "s_meas_pt": s_meas_pt,
                "frac_width": got / need if need else np.nan})
    print(f"  -> P4/P5 {'HOLD: moving the operating point by the curvature term accounts for most of the width the linear gain missed' if 0.5 < got/need < 1.6 else ('the curvature moves the width the right way but explains only part -- named as partial, not fitted' if got/need > 0 else 'FAILS: curvature moves the width the wrong way')}")
    print("  **Whatever fraction this is, it is reported and not tuned. A partial mechanism")
    print("    named as partial is worth more than a fitted one (rule 17).**")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
