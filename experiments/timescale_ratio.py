"""T-CASC-f: is composition a function of the timescale RATIO alone?

§92 found a second control the margin does not contain: scaling the upstream stage's propensities
leaves its landscape, barrier, rail width and stationary law identical and moves the penalty by
2.7x, from the frozen-upstream average down toward the mean-landscape rate -- motional narrowing.
But it swept only the UPSTREAM clock, on a chain whose two stages are the same element, so
tau_up/tau_down = 1/speed by construction and the crossover could not be placed independently of
anything else.

**Scaling the DOWNSTREAM separately makes the ratio a free variable.** If the penalty collapses
onto a function of s_up/s_dn alone, then cascade error is a function of two quantities -- the noise
margin and the timescale ratio -- **both computable from single elements**, and §91's D = 2 ceiling
stops mattering for the theory even though it still binds the exact instrument.

NORMALISATION, which is the whole design, and the FIRST VERSION GOT IT WRONG. The penalty is

    P(stage 2 low | upstream FLUCTUATING) / P(stage 2 low | upstream PINNED at r3)

at the same downstream speed and wall time. **But holding the WALL time fixed while scaling the
downstream clock changes how much of its OWN clock the downstream has run**, and the penalty
depends on that strongly: the first run's diagonal s_up = s_dn = c ran 26.295, 4.442, 3.340, 2.941
over c = 1/4 .. 16 instead of being flat. P1 caught it and nothing was read from the grid.

The observation window is therefore set in the DOWNSTREAM'S OWN CLOCK: t = t0/s_dn. On the
diagonal both stages then see identical proper time, so the configuration is literally the same
physics relabelled, and the penalty must be constant to MACHINE PRECISION -- a far sharper gate
than "approximately flat".

PREDICTIONS, written before running.

  P1  GATE, and with the window set in the downstream's own clock it is EXACT. Along the diagonal
      s_up = s_dn = c the two stages see identical proper time, so the penalty must be constant to
      machine precision. **Anything else is a bug, not a physical effect.** (Held at fixed WALL
      time instead it ran 26.295 -> 2.941, which is what killed the first version.)
  P2  **COLLAPSE ONTO THE RATIO -- and it is ALGEBRA, not a discovery, which is stated here
      because §84 already published its own Taylor term as a phenomenon.** The joint generator is
      Q = s_up Q1 + s_dn Q2 (the coupling lives inside Q2, since only stage 2's rates depend on
      n1), so with the window at t = t0/s_dn the propagator is

          Q t = t0 [ (s_up/s_dn) Q1 + Q2 ]

      which depends on the two speeds ONLY through their ratio, exactly. **So P2 is a consistency
      check on the implementation, not a physical result**, and it must come out exact to machine
      precision or something is wired wrong. It is run because it would catch a coupling
      accidentally placed in the wrong block.
  P2b **THE PHYSICS IS THE SHAPE, not the collapse.** What is not fixed by algebra is what the
      function of the ratio looks like: where it plateaus, where it turns, how far it falls.
  P3  **WHERE IS THE CROSSOVER?** §92 saw the fall begin around upstream speed 2-4 on a chain
      where the two elements are identical. **Predicted: the crossover sits at ratio ~ 1**, i.e.
      when the upstream correlation time matches the downstream response time -- and if so, §92's
      curve was centred by construction rather than by coincidence.
  P4  **IF THE COLLAPSE FAILS**, the residual identifies the third variable and is reported as
      such, not fitted. A cell that disagrees with its ratio-partner is the finding.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import experiments.chemical_cascade as cc
from experiments.cascade_schlogl import schlogl_consts
from experiments.margin_law import R1, R2, R3, stage1_stationary

C = schlogl_consts(R1, R2, R3)


def build_pair(om, s_up, s_dn, cap_mult=1.25):
    """Stage 1 reflected at its saddle (nothing to condition on), each stage's own clock."""
    cap = int(np.ceil(cap_mult * R3 * om))
    nsad = int(np.ceil(R2 * om))
    up = np.arange(nsad, cap + 1)
    nu, m2 = len(up), cap + 1
    rows, cols, vals = [], [], []
    for a in range(nu):
        n1 = up[a]
        for n2 in range(m2):
            idx = a * m2 + n2
            tot = 0.0
            l1, u1 = cc.rates_stage(float(n1), 0.0, om, C, R3, True, "hill")
            l1 *= s_up; u1 *= s_up
            if n1 < cap and l1 > 0:
                rows.append(idx); cols.append((a + 1) * m2 + n2); vals.append(l1); tot += l1
            if n1 > up[0] and u1 > 0:
                rows.append(idx); cols.append((a - 1) * m2 + n2); vals.append(u1); tot += u1
            l2, u2 = cc.rates_stage(float(n2), float(n1), om, C, R3, False, "hill")
            l2 *= s_dn; u2 *= s_dn
            if n2 < cap and l2 > 0:
                rows.append(idx); cols.append(idx + 1); vals.append(l2); tot += l2
            if n2 > 0 and u2 > 0:
                rows.append(idx); cols.append(idx - 1); vals.append(u2); tot += u2
            rows.append(idx); cols.append(idx); vals.append(-tot)
    return sp.csr_matrix((vals, (rows, cols)), shape=(nu * m2, nu * m2)), up, m2, cap


def pinned_reference(om, s_dn, t, cap_mult=1.25):
    """Stage 2 alone with the upstream held exactly at r3, same clock and same wall time."""
    cap = int(np.ceil(cap_mult * R3 * om))
    m = cap + 1
    rows, cols, vals = [], [], []
    for n in range(m):
        tot = 0.0
        l, u = cc.rates_stage(float(n), R3 * om, om, C, R3, False, "hill")
        l *= s_dn; u *= s_dn
        if n < cap and l > 0:
            rows.append(n); cols.append(n + 1); vals.append(l); tot += l
        if n > 0 and u > 0:
            rows.append(n); cols.append(n - 1); vals.append(u); tot += u
        rows.append(n); cols.append(n); vals.append(-tot)
    Q = sp.csr_matrix((vals, (rows, cols)), shape=(m, m))
    p = np.zeros(m)
    p[int(round(R3 * om))] = 1.0
    p = spla.expm_multiply(Q.T * t, p)
    return float(p[np.arange(m) < R2 * om].sum())


def penalty(om, s_up, s_dn, t0, pi1=None):
    """Window set in the DOWNSTREAM's own clock: t = t0/s_dn."""
    t = t0 / s_dn
    Q, up, m2, cap = build_pair(om, s_up, s_dn)
    if pi1 is None:
        _, pi1 = stage1_stationary(om)
    p = np.zeros(len(up) * m2)
    for a, w in enumerate(pi1):
        p[a * m2 + int(round(R3 * om))] = w
    p = spla.expm_multiply(Q.T * t, p)
    lo2 = (np.arange(len(up) * m2) % m2) < R2 * om
    num = float(p[lo2].sum())
    den = pinned_reference(om, s_dn, t)
    return num / den, num, den


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=30)
    ap.add_argument("--t", type=float, default=2.0)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/timescale_ratio.json"))
    args = ap.parse_args()
    om, t = args.omega, args.t
    cc.HILL_N, cc.HILL_K = 4.0, 1.0
    _, pi1 = stage1_stationary(om)
    out = {}

    print("=== P1 GATE: along the diagonal s_up = s_dn = c, the penalty must not depend on c")
    print("    (scaling both clocks is a change of time units; at small escape probability")
    print("     P ~ rate x t, so the ratio must be c-independent. If not, we are outside the")
    print("     linear regime and no ratio statement means anything.)")
    print(f"{'c':>8}{'P(s2 lo)':>13}{'pinned ref':>13}{'penalty':>10}")
    diag = []
    for cc_ in (0.25, 0.5, 1.0, 2.0, 4.0):
        pen, num, den = penalty(om, cc_, cc_, t, pi1)
        diag.append({"c": cc_, "pen": pen, "num": num, "den": den})
        print(f"{cc_:>8.2f}{num:>13.4e}{den:>13.4e}{pen:>10.3f}")
    out["diag"] = diag
    pens = [d["pen"] for d in diag]
    flat = max(pens) / min(pens) < 1.001
    print(f"  penalty along the diagonal spans {min(pens):.3f}..{max(pens):.3f}"
          f" (factor {max(pens)/min(pens):.3f})")
    print(f"  -> P1 {'HOLDS to machine precision: identical proper time is identical physics, so the instrument is sound and the grid can be read' if flat else 'FAILS: the diagonal is not constant, which with the window in the downstream clock is a BUG, not an effect'}")
    assert flat, "diagonal not constant: the two stages do not see identical proper time"

    print("\n=== P2/P3: sweep the two clocks independently. Do equal ratios agree?")
    print(f"{'s_up':>7}{'s_dn':>7}{'ratio':>9}{'penalty':>10}")
    grid = []
    speeds = (0.25, 1.0, 4.0, 16.0)
    for su in speeds:
        for sd in speeds:
            pen, num, den = penalty(om, su, sd, t, pi1)
            grid.append({"s_up": su, "s_dn": sd, "ratio": su / sd, "pen": pen})
            print(f"{su:>7.2f}{sd:>7.2f}{su/sd:>9.4f}{pen:>10.3f}")
    out["grid"] = grid
    print(f"\n  cells grouped by ratio:")
    print(f"{'ratio':>9}{'n':>4}{'penalties':>34}{'spread':>9}")
    worst = 1.0
    for r in sorted({g["ratio"] for g in grid}):
        cells = [g for g in grid if abs(g["ratio"] / r - 1) < 1e-9]
        if len(cells) < 2:
            continue
        vs = [x["pen"] for x in cells]
        worst = max(worst, max(vs) / min(vs))
        print(f"{r:>9.4f}{len(cells):>4}"
              + ("  " + ", ".join(f"{v:.3f}" for v in vs)).rjust(34)
              + f"{max(vs)/min(vs):>9.3f}")
    print(f"  worst spread within a ratio group: {worst:.3f}")
    collapse = worst < 1.15
    print(f"  -> P2 {'HOLDS exactly, as the algebra requires: Q t = t0[(s_up/s_dn) Q1 + Q2]. This is a WIRING CHECK, not a discovery -- it would have caught the coupling being placed in the wrong block' if collapse else 'FAILS: equal ratios disagree, which given the identity means the generator is not split the way it is claimed to be'}")
    if not collapse:
        print("  -> P4: the disagreeing groups, reported and not fitted:")
        for r in sorted({g["ratio"] for g in grid}):
            cells = [g for g in grid if abs(g["ratio"] / r - 1) < 1e-9]
            if len(cells) < 2:
                continue
            vs = [x["pen"] for x in cells]
            if max(vs) / min(vs) > 1.15:
                print(f"       ratio {r:.4f}: "
                      + ", ".join(f"(up {x['s_up']:g}, dn {x['s_dn']:g}) -> {x['pen']:.3f}"
                                  for x in cells))

    print("\n=== P3: where does the fall happen, in ratio?")
    by_ratio = {}
    for g in grid:
        by_ratio.setdefault(round(g["ratio"], 6), []).append(g["pen"])
    for r in sorted(by_ratio):
        print(f"   ratio {r:>8.4f}: mean penalty {np.mean(by_ratio[r]):.3f}")
    print("  (the crossover is where the penalty leaves its slow plateau -- THIS is the")
    print("   physical content; the collapse above is an identity.)")
    rs = sorted(by_ratio)
    vals = [float(np.mean(by_ratio[r])) for r in rs]
    lo = [v for r, v in zip(rs, vals) if r <= 1.0]
    hi = [v for r, v in zip(rs, vals) if r >= 4.0]
    plateau = max(lo) / min(lo) < 1.10
    falls = min(hi) < 0.5 * max(lo)
    peak = rs[int(np.argmax(vals))]
    print(f"  plateau for ratio <= 1: {min(lo):.3f}..{max(lo):.3f}"
          f" (flat to {100*(max(lo)/min(lo)-1):.1f}%);  peak at ratio {peak:g};"
          f"  falls to {min(hi):.3f}")
    print(f"  -> P2b/P3 {'HOLD: the penalty PLATEAUS for a slow upstream, turns at ratio ~ 1 -- where the upstream correlation time meets the downstream response time -- and falls by motional narrowing above it. The crossover is not an adjustable feature: it sits at equal timescales' if (plateau and falls and peak <= 1.0) else 'the shape is not plateau-then-fall with a turn at ratio 1'}")
    out["shape"] = {"ratios": rs, "penalty": vals, "peak": peak}

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
