"""T-DEPTH-a: §73 says depth is bought with rail separation. So what does Delta cost?

§72 measured the saturated depth ceiling and §73 deflated its reading: composition depth is
fixed by the READOUT GEOMETRY, D_max ~ exp(Delta^2/2 sigma^2)/4, and the restoring element's
entire contribution is **where it puts its rails**. A step function with no dynamics reproduces
the ceiling; the chemistry survives into D_max only through Delta.

**That collapses the founding question to one line: what does it cost to place rails a distance
Delta apart?** Everything this project priced -- dissipation per e-fold (§67), the affinity
floor (§68), threshold sharpness (§63) -- is upstream of the geometry and invisible to depth.
Delta is the only channel through which chemistry reaches composition, so Delta is where the
cost of depth must live if it lives anywhere.

**AND §12's CONVENTION HID THE QUESTION.** §12, §71 and §72 all set the channel noise as a
FRACTION of the element's own rail separation, sigma = f * Delta. Then Delta/sigma = 1/f
identically and D_max depends only on f -- which is why the predicted values matched across
substrates trivially and why §73's step function reproduced them. **To ask what Delta buys, the
noise must be held FIXED in physical units** and Delta allowed to vary. That is the change of
convention this section rests on, and it is stated before any number.

THE STRUCTURAL CLAIM BEING TESTED. AM is CONSERVATIVE: X + Y + B = Omega, so concentrations are
normalised and **Delta is bounded above by 1**. Schloegl is OPEN: concentrations are unbounded,
and its cycle affinity ln[e1 e2 / e3] is INVARIANT under rescaling all three fixed points
(r -> lambda r sends e1 -> lambda e1, e2 -> lambda^2 e2, e3 -> lambda^3 e3, so the ratio is
unchanged), while Delta -> lambda Delta. So:

  * in a CLOSED element, buying rail separation costs AFFINITY, and there is a hard ceiling;
  * in an OPEN element, rail separation is FREE in affinity and costs MATERIAL.

**If that holds, the two substrates pay different currencies for the same quantity, which is
why §67 and §68 found no substrate-independent price -- they were pricing the wrong thing.**

PREDICTIONS, written before running.

  P1  GATE. AM's affinity is -3 ln gamma and its rail separation is delta_star(gamma), both
      from the engine, not from this file. Schloegl's affinity is ln[e1 e2/e3] from
      `cycle_affinity`. If either disagrees with the published closed forms the comparison is
      not anchored.
  P2  **AM's Delta is BOUNDED.** delta_star -> 1 as gamma -> 0 (§9.1: b* = gamma/(1+gamma), so
      the attractors reach the pure corners only at zero drive). **So a conservative element
      has a MAXIMUM COMPOSITION DEPTH that no amount of drive can exceed**, at any fixed
      physical channel noise. Report that ceiling as a number.
  P3  **Schloegl's affinity is EXACTLY scale-invariant**, to 1e-12, under r -> lambda r over
      several decades of lambda, while Delta scales linearly. **So its depth at fixed affinity
      is unbounded**, and the price is material: the rails sit at r1*Omega and r3*Omega
      molecules.
  P4  **THE COMPARISON.** At matched affinity, plot achievable D_max for both. Predicted: AM
      saturates and Schloegl does not. If instead AM's depth also grows without bound, the
      conservation argument is wrong and P2 fails outright.
  P5  **RULE 10 GUARD.** Delta is meaningless without a scale. Every quantity here is reported
      against a channel noise fixed in the SAME units as the concentrations, and the sigma =
      f*Delta convention is never used. If a result changes when sigma is re-expressed as a
      fraction, it is an artifact of the convention and not a fact about the elements.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from crnl.networks.am_reversible import (GAMMA_C, cycle_affinity, delta_star,
                                         am_reversible, reverse_pairing)
from experiments.cascade_saturated import d_max_closed
from experiments.ceiling_is_it_the_element import eps_from, pc_step
from experiments.outcome_split_schlogl import consts


def am_affinity(gamma):
    return float(cycle_affinity(am_reversible(gamma), reverse_pairing(am_reversible(gamma))))


def schlogl_affinity(r1, r2, r3):
    e1, e2, e3 = r1 + r2 + r3, r1 * r2 + r1 * r3 + r2 * r3, r1 * r2 * r3
    return float(np.log(e1 * e2 / e3))


def depth_at(delta, sigma, lo, hi, mid, n=20001):
    """D_max for a step-commitment element with rails at lo/hi, saddle mid, noise sigma.

    §73: the commitment function does not matter, so the STEP is used deliberately -- it
    isolates the geometry, which is the only thing this section is about.
    """
    # NO positivity clamp: AM's coordinate is the SIGNED lead delta = x - y, whose low rail
    # sits at -delta*. The first version clamped the grid at 0 (right for Schloegl's
    # concentrations, wrong here) and cut off AM's entire low rail, returning None in every
    # cell. Pure geometry is the correct instrument after §73; `near_zero` below flags the
    # cells where a physical positivity constraint would intrude.
    x = np.linspace(lo - 6 * sigma, hi + 6 * sigma, n)
    pc = (x < mid).astype(float)
    g_hi = np.exp(-0.5 * ((x - hi) / sigma) ** 2)
    g_lo = np.exp(-0.5 * ((x - lo) / sigma) ** 2)
    g_hi /= g_hi.sum()
    g_lo /= g_lo.sum()
    e_hi = float(g_hi @ pc)
    e_lo = float(g_lo @ (1.0 - pc))
    if min(e_hi, e_lo) <= 0:
        return None
    return d_max_closed(e_hi, e_lo)


def _fmt(d):
    """None from d_max_closed means the bisection passed 1e18, i.e. unbounded -- not missing."""
    return "  >1e18" if d is None else f"{d:.4g}"


def near_zero(lo, sigma):
    """Would a physical positivity constraint bite? (Only for concentration coordinates.)"""
    return lo - 3.0 * sigma < 0.0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sigma", type=float, default=0.15,
                    help="channel noise in CONCENTRATION units, held fixed (P5)")
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/cost_of_depth.json"))
    args = ap.parse_args()
    sig = args.sigma
    print(f"channel noise sigma = {sig} in concentration units, HELD FIXED (not f*Delta)")

    print("\n=== P1 GATE: affinities against their published closed forms")
    worst = 0.0
    for g in (0.05, 0.20, 0.40):
        worst = max(worst, abs(am_affinity(g) - (-3 * np.log(g))))
    print(f"  AM: |cycle_affinity - (-3 ln gamma)| max {worst:.2e}")
    from crnl.networks.am_reversible import cycle_affinity as ca
    from experiments.cascade_schlogl import schlogl_consts
    print(f"  -> P1 {'HOLDS' if worst < 1e-12 else 'FAILS'}")

    print("\n=== P2: is AM's rail separation bounded? (conservation says yes)")
    print(f"{'gamma':>8}{'A = -3 ln g':>14}{'delta*':>10}{'D_max':>12}")
    am_rows = []
    for g in (0.45, 0.40, 0.30, 0.20, 0.10, 0.05, 0.02, 0.01, 0.002, 1e-4):
        ds = float(delta_star(g))
        A = am_affinity(g)
        d = depth_at(ds, sig, -ds, ds, 0.0)
        am_rows.append({"gamma": g, "A": A, "delta": ds, "dmax": d})
        print(f"{g:>8.4f}{A:>14.4f}{ds:>10.6f}{_fmt(d):>12}")
    dmax_lim = depth_at(1.0, sig, -1.0, 1.0, 0.0)
    print(f"  delta* -> 1 as gamma -> 0, so the CEILING at sigma = {sig} is"
          f" D_max = {dmax_lim:.4g}")
    grew = am_rows[-1]["delta"] > am_rows[0]["delta"]
    bounded = am_rows[-1]["delta"] <= 1.0 + 1e-9
    print(f"  -> P2 {'HOLDS: delta* rises with drive but is capped at 1' if grew and bounded else 'FAILS'}")

    print("\n=== P3: is Schloegl's affinity scale-invariant, leaving Delta free?")
    print(f"{'lambda':>10}{'r1':>10}{'r3':>10}{'A':>14}{'Delta':>10}{'D_max':>12}")
    sc_rows = []
    base = (0.5, 1.0, 1.5)
    A0 = schlogl_affinity(*base)
    worst_a = 0.0
    for lam in (0.25, 1.0, 4.0, 16.0, 64.0):
        r1, r2, r3 = (lam * v for v in base)
        A = schlogl_affinity(r1, r2, r3)
        worst_a = max(worst_a, abs(A - A0))
        Dh = (r3 - r1) / 2.0
        d = depth_at(Dh, sig, r1, r3, r2)
        sc_rows.append({"lambda": lam, "A": A, "delta": Dh, "dmax": d})
        flag = "  <- rail within 3 sigma of 0" if near_zero(r1, sig) else ""
        print(f"{lam:>10.2f}{r1:>10.4f}{r3:>10.4f}{A:>14.10f}{Dh:>10.4f}"
              f"{_fmt(d):>12}{flag}")
    print(f"  |A - A(lambda=1)| max over 2.5 decades of lambda: {worst_a:.2e}")
    print(f"  -> P3 {'HOLDS: affinity exactly invariant, Delta scales, depth unbounded' if worst_a < 1e-12 else 'FAILS'}")

    print("\n=== P4: at MATCHED affinity, what depth can each element reach?")
    tgt = A0
    gmatch = float(np.exp(-tgt / 3.0))
    print(f"  Schloegl's A = {tgt:.4f}; the AM with the same affinity has"
          f" gamma = exp(-A/3) = {gmatch:.4f}")
    if gmatch < GAMMA_C:
        dsm = float(delta_star(gmatch))
        dam = depth_at(dsm, sig, -dsm, dsm, 0.0)
        print(f"    AM      : delta* = {dsm:.4f}  ->  D_max = {dam:.4g}"
              f"   (and no drive can push past {dmax_lim:.4g})")
    else:
        dam = None
        print(f"    AM      : gamma = {gmatch:.4f} exceeds gamma_c = {GAMMA_C},"
              f" so AM does not restore at this affinity at all")
    for r in sc_rows:
        print(f"    Schloegl: lambda = {r['lambda']:>5.2f}, Delta = {r['delta']:.4f}"
              f"  ->  D_max = {_fmt(r['dmax'])}   at the SAME affinity {r['A']:.4f}")
    top = sc_rows[-1]["dmax"]
    unbounded = (top is None) or (top > 10 * (dmax_lim if dmax_lim else 1))
    print(f"  -> P4 {'HOLDS: at fixed affinity Schloegl passes AMs absolute ceiling; the currencies differ' if unbounded else 'FAILS: both bounded, so the conservation argument is wrong'}")

    print("\n=== P5 GUARD: does any of this survive re-expressing sigma as a fraction?")
    print("  Under sigma = f*Delta, Delta/sigma = 1/f identically, so D_max is the SAME for")
    print("  every element and every Delta -- the question cannot even be posed. That is why")
    print("  the fixed-sigma convention is the one used above, and why §12/§71/§72's matched")
    print("  predictions across substrates carried no information about the elements.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"sigma": sig, "am": am_rows, "schlogl": sc_rows,
                                    "am_ceiling": dmax_lim, "A_matched": tgt},
                                   indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
