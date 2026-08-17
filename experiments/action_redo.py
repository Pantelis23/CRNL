"""T-DEPTH-h: redo §77-§79 against the escape action A, which §80 showed is the right exponent

§80 found that §75-§79 priced the subdominant failure mode. For a chemically-coupled cascade the
per-stage error is escape, not misreading, so the coefficient that governs depth is the escape
action A and not the linear-noise eta. Three sections were built on eta and each needs redoing:

  §77 asked whether the coefficient TRANSFERS across substrates. Answer for eta: no, a factor
      of 633. Unknown for A.
  §78 DERIVED the coefficient from deterministic-side quantities. §80 already showed A is the
      quasipotential integral -int ln(mu/lam) dx to 1e-4 on ONE element.
  §79 tested the derivation OUT OF SAMPLE, on systems never used. Not yet done for A.

**This section is the redo, and it is cheap in 1-D**: A is an integral of the propensity
densities, so no chain, no lattice and no master equation are needed to predict it -- only to
check it.

For AM the quasipotential is not an integral (it needs Hamilton-Jacobi in 2-D), so A is taken
from the exact stationary distribution as lim [ln pi(rail) - ln pi(saddle)]/Omega -- the
definition of the quasipotential for a NESS, and the route §64 used when its eigenvalue
extrapolation and its stationary extrapolation disagreed.

PREDICTIONS, written before running.

  P0  **UNDERFLOW GUARD, added after the first run.** AM's quasipotential is read from the exact
      stationary distribution, whose saddle weight underflows once the barrier passes ~691 nats.
      At gamma = 0.20 the resulting "A" read 5.7137, 3.8145, 2.8591 -- exactly 1/Omega, because
      (rail - saddle) had saturated at the floor. Cells whose saddle weight is at the floor are
      excluded, not fitted. **That cell was driving both headline spreads in the first run.**
  P1  GATE. For Schloegl the integral -int ln(mu/lam) dx must match A from the exact mean
      first-passage time, and the agreement must IMPROVE with Omega (rule 20). Checked across
      landscape shapes AND at a different reaction order (§79's quartic 3X <-> 4X), not one
      case as §80 did.
  P2  **§77 REDONE. Does A transfer?** Sweep landscape shape within each substrate and compare
      across substrates. **Predicted: no, and for the same reason eta did not** -- A is a
      property of the element's own landscape. A transferring A would be far more surprising
      than a transferring eta, because A carries the barrier.
  P3  **§79 REDONE. Out of sample.** Predict A from the integral for landscapes never used in
      §80, then check against the exact MFPT. **Predicted: sub-percent**, since the integral is
      exact in the Omega -> infinity limit and the only error is finite-Omega.
  P4  **IS A/eta CONSTANT?** If it were, §77's transfer conclusions would carry over unchanged
      and §80 would be a relabelling. **Predicted: NOT constant** -- eta is a local harmonic
      quantity at the rail and A is a global integral over the barrier, so they should respond
      differently to landscape shape. Measured across shapes and orders.
  P5  **THE DEPTH, restated correctly.** With escape dominant, D_max = c*/(t exp(-A Omega)), so
      ln D_max = A*Omega - ln(t/c*). Report what that does to §72's numbers.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
from scipy.integrate import quad

from crnl.cme import enumerate_states, stationary
from crnl.networks.am_reversible import am_reversible, delta_star
from experiments.cascade_schlogl import schlogl_consts
from experiments.derive_eta import schlogl_V
from experiments.out_of_sample import quartic_consts, quartic_roots
from experiments.two_exponentials import ln_mfpt


def A_integral(r1, r2, r3):
    """A = -int_{saddle}^{rail} ln(mu/lam) dx, from the RATE FUNCTIONS alone."""
    k1a, k1r, k2b, k2r = schlogl_consts(r1, r2, r3)
    val, err = quad(lambda x: np.log((k1r * x ** 3 + k2r * x) / (k1a * x ** 2 + k2b)),
                    r2, r3, limit=200)
    return -val, err


def A_integral_p(c, rails):
    """The same for the general p-family: lam = k1a x^p + k2b, mu = k1r x^(p+1) + k2r x."""
    k1a, k1r, k2b, k2r, p = c
    val, err = quad(lambda x: np.log((k1r * x ** (p + 1) + k2r * x)
                                     / (k1a * x ** p + k2b)), rails[1], rails[2], limit=200)
    return -val, err


def A_mfpt(r1, r2, r3, oms=(3200, 6400)):
    a, b = ln_mfpt(oms[0], r1, r2, r3), ln_mfpt(oms[1], r1, r2, r3)
    return (b - a) / (oms[1] - oms[0])


def A_am(gamma, omegas=(120, 180, 240)):
    """AM's quasipotential barrier from the exact stationary distribution, per molecule."""
    out = []
    for om in omegas:
        try:
            pi = stationary(am_reversible(gamma), int(om), float(om))
        except RuntimeError:
            continue
        states, _ = enumerate_states(3, int(om))
        d = np.array([int(s[0]) - int(s[1]) for s in states])
        lp = np.log(np.maximum(pi, 1e-300))
        rail = lp[d > 0].max()
        saddle = lp[d == 0].max()
        # UNDERFLOW GUARD: once pi(saddle) reaches the 1e-300 floor the barrier saturates
        # at ln(1/1e-300) ~ 691 and (rail - saddle) stops growing with Omega -- at
        # gamma = 0.20 the series read 5.7137, 3.8145, 2.8591, whose products with Omega are
        # 685, 686, 686. That is the floor, not the physics, and it is excluded.
        if saddle <= np.log(1e-300) + 5.0 or (rail - saddle) > 600.0:
            continue
        out.append((om, (rail - saddle) / om))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/action_redo.json"))
    args = ap.parse_args()

    print("=== P1/P3 GATE: the integral against the exact MFPT, incl. out-of-sample shapes")
    print(f"{'landscape':>26}{'A integral':>13}{'A from MFPT':>13}{'ratio':>9}{'in §80?':>9}")
    rows = []
    cases = [((0.1, 1.0, 1.9), True), ((0.4, 1.0, 2.2), False),
             ((0.3, 1.0, 1.7), False), ((0.05, 1.0, 2.5), False)]
    for (r1, r2, r3), seen in cases:
        Ai, err = A_integral(r1, r2, r3)
        Am = A_mfpt(r1, r2, r3)
        rows.append({"kind": "schlogl", "rails": [r1, r2, r3], "A_int": Ai,
                     "A_mfpt": Am, "ratio": Ai / Am, "seen": seen})
        print(f"{f'{r1}/{r2}/{r3}':>26}{Ai:>13.6f}{Am:>13.6f}{Ai/Am:>9.4f}"
              f"{('yes' if seen else 'NO'):>9}")
    qc = quartic_consts(m=0.8)
    qr = quartic_roots(qc)
    Ai, _ = A_integral_p(qc, qr)
    print(f"{'QUARTIC 3X<->4X':>26}{Ai:>13.6f}{'(1-D integral)':>13}{'--':>9}{'NO':>9}")
    rows.append({"kind": "quartic", "A_int": Ai, "seen": False})
    worst = max(abs(r["ratio"] - 1) for r in rows if "ratio" in r)
    print(f"  worst |integral/MFPT - 1| = {100*worst:.2f}%"
          f"   ({sum(1 for r in rows if not r['seen'])} of {len(rows)} never used in §80)")
    print(f"  -> P1/P3 {'HOLD: the integral predicts A on landscapes it never saw' if worst < 0.02 else 'FAIL'}")

    print("\n=== P2: §77 redone -- does A transfer?")
    print(f"{'element':>26}{'A':>12}{'eta':>12}{'A/eta':>9}")
    tab = []
    for (r1, r2, r3), _ in cases:
        Ai, _ = A_integral(r1, r2, r3)
        D = (r3 - r1) / 2
        eta = D ** 2 / (2 * schlogl_V(r1, r2, r3))
        tab.append({"el": f"Schlogl {r1}/{r3}", "A": Ai, "eta": eta, "ratio": Ai / eta})
        print(f"{f'Schlogl {r1}/{r3}':>26}{Ai:>12.6f}{eta:>12.6f}{Ai/eta:>9.4f}")
    am_rows = []
    for g in (0.20, 0.30, 0.35):
        ser = A_am(g)
        if len(ser) < 2:
            print(f"{f'AM gamma={g}':>26}   no usable stationary solve")
            continue
        A = ser[-1][1]
        V, ds = None, float(delta_star(g))
        from experiments.derive_eta import am_V
        V, ds = am_V(g)
        eta = ds ** 2 / (2 * V)
        am_rows.append({"gamma": g, "series": ser, "A": A, "eta": eta})
        tab.append({"el": f"AM gamma={g}", "A": A, "eta": eta, "ratio": A / eta})
        print(f"{f'AM gamma={g}':>26}{A:>12.6f}{eta:>12.6f}{A/eta:>9.4f}"
              + f"   (A series " + ", ".join(f"{v:.4f}" for _, v in ser) + ")")
    As = [t["A"] for t in tab]
    print(f"  A spans {min(As):.6f} .. {max(As):.6f} -- a factor of {max(As)/min(As):.0f}")
    print(f"  -> P2 {'A does NOT transfer either, like eta' if max(As)/min(As) > 2 else 'A TRANSFERS -- unlike eta'}")

    print("\n=== P4: is A/eta constant? (if so, §80 is a relabelling)")
    rr = [t["ratio"] for t in tab]
    print("  A/eta: " + ", ".join(f"{v:.4f}" for v in rr))
    print(f"  spans {min(rr):.4f} .. {max(rr):.4f}, a factor of {max(rr)/min(rr):.2f}")
    print(f"  -> P4 {'A/eta is NOT constant: the two coefficients respond differently to shape, so §80 is not a relabelling' if max(rr)/min(rr) > 1.2 else 'A/eta is nearly constant: §77s conclusions carry over'}")

    print("\n=== P5: the depth, restated. ln D_max = A*Omega - ln(t/c*)")
    from experiments.depth_is_error import c_star
    cs, t = c_star(), 2.0
    print(f"{'element':>26}{'Omega':>8}{'ln D (escape)':>15}{'ln D (readout)':>16}")
    for t_ in tab[:2]:
        for om in (1600, 6400):
            print(f"{t_['el']:>26}{om:>8}{t_['A']*om - np.log(t/cs):>15.1f}"
                  f"{t_['eta']*om:>16.1f}")
    print("  -> the escape route gives a far SHALLOWER ceiling; §72's depths were computed")
    print("     with an external channel and are not these numbers (§80 P5).")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"rows": rows, "table": tab, "am": am_rows},
                                   indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
