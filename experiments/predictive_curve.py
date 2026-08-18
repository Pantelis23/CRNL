"""T15-n.3: predict the escape curve from the tilted generator, with no stationary solve

§86 showed the 1-D reductions integrate the right formula along the WRONG CURVE: the
quasipotential is realised on the ridge argmax_b pi(u,b), not on the deterministic slow manifold,
and integrating along the ridge closes the deficit (0.9027 -> 1.0088, 0.9348 -> 1.0107).
**But the ridge was read off the exact stationary distribution, so §86 is a diagnosis and
predicts nothing.** T15-n.3 asks whether the correct curve is computable in advance.

THE ROUTE, which the §85.3 literature check named. Write the WKB Hamiltonian of the full 2-D
jump process in the lead/total coordinates u = x - y, s = x + y (b = 1 - s):

    H = xy(e^-2ps - 1) + bx(e^(pu+ps) - 1) + by(e^(-pu+ps) - 1)
        + g b^2(e^2ps - 1) + g x^2(e^(-pu-ps) - 1) + g y^2(e^(pu-ps) - 1)

s is the fast variable. The slow-fast reduction freezes (u, pu) and sends the FAST HAMILTONIAN
FLOW to its fixed point -- not to the deterministic steady state:

    dH/dps = 0      (s-dot = 0)          and      dH/ds = 0      (ps-dot = 0)

Then H = 0 fixes pu(u), and the action is A = int pu du. **The deterministic slow manifold is
exactly the pu = ps = 0 branch of this**, which is why §84/§85 got it: they solved the right
equations at zero momentum. The escape path carries pu != 0, and s shifts with it. **If that
shift reproduces §86's measured displacement, the reduction becomes predictive and the ridge stops
being an oracle.**

Three equations, three unknowns (s, ps, pu) at each u. No master equation, no stationary solve,
no lattice.

PREDICTIONS, written before running.

  P1  GATE, and it is exact. At pu = ps = 0 the Hamiltonian must vanish identically for every
      (u, s) -- that is the statement that the zero-momentum sheet is the deterministic dynamics.
      And dH/dps = 0 at zero momentum must reproduce slow_manifold(u) to machine precision. If
      the p = 0 branch does not recover §85's curve, the Hamiltonian is written wrong and nothing
      below counts.
  P2  **THE PREDICTION. The solved curve must be DISPLACED from the deterministic manifold in the
      same direction and by comparable magnitude as §86's measured ridge** -- upward (more blank),
      by ~0.003-0.005 at gamma = 0.40 and ~0.001-0.002 at gamma = 0.44. §86's ridge values are
      quoted as stored numbers (rule 16) and not recomputed.
  P3  **THE ACTION, in absolute terms.** A = int pu du against the exact first-passage action.
      **Predicted: it lands near 1.00**, like §86's ridge integration did, and decisively better
      than the deterministic manifold's 0.9027 / 0.9348. If it reproduces the deficit instead of
      curing it, the fast-variable momentum is not what the elimination was losing.
  P4  **THE POINT: this uses no stationary solve.** Report that the whole computation is root
      finding on a 3-equation system. A prediction that needs the exact answer is a diagnosis.
  P5  **RULE 9.** Both gammas, and the sign must not be fitted -- it is whatever the equations
      give. If the predicted displacement has the WRONG SIGN, §86's agreement was a coincidence
      of the marginal and the ridge is not the instanton's curve.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, fsolve

from experiments.where_the_deficit_lives import MEASURED, slow_manifold, u_star

# §86's measured ridge displacement, quoted as stored numbers (rule 16).
RIDGE_DISP = {0.40: 0.00254, 0.44: 0.00104}
DET_RATIO = {0.40: 0.9027, 0.44: 0.9348}


def H(u, s, pu, ps, g):
    """WKB Hamiltonian of the 2-D jump process in (u, s) with b = 1 - s."""
    x, y, b = (s + u) / 2.0, (s - u) / 2.0, 1.0 - s
    return (x * y * (np.exp(-2 * ps) - 1.0)
            + b * x * (np.exp(pu + ps) - 1.0)
            + b * y * (np.exp(-pu + ps) - 1.0)
            + g * b * b * (np.exp(2 * ps) - 1.0)
            + g * x * x * (np.exp(-pu - ps) - 1.0)
            + g * y * y * (np.exp(pu - ps) - 1.0))


def _grad(u, s, pu, ps, g, h=1e-6):
    dps = (H(u, s, pu, ps + h, g) - H(u, s, pu, ps - h, g)) / (2 * h)
    dsv = (H(u, s + h, pu, ps, g) - H(u, s - h, pu, ps, g)) / (2 * h)
    return dps, dsv


def fast_fixed(u, pu, g, guess=None):
    """At FIXED pu, send the fast pair (s, ps) to the fixed point of the fast Hamiltonian flow.

    Two equations, two unknowns. pu is a parameter here, deliberately: see `pu_at`.
    """
    x0 = guess if guess is not None else np.array([1.0 - slow_manifold(u, g), 0.0])

    def eqs(v):
        return list(_grad(u, v[0], pu, v[1], g))

    sol, _, ier, _ = fsolve(eqs, x0, full_output=True)
    if ier != 1 or not (0.0 < sol[0] < 1.0):
        return None
    return float(sol[0]), float(sol[1])


def pu_at(u, g, guess=None, pu_max=1.5):
    """The NON-TRIVIAL zero-energy momentum at u.

    **The three-equation form of this collapsed, and the reason is in P1.** H(u, s, 0, 0) = 0
    IDENTICALLY IN s, so the whole p = 0 sheet is a two-parameter family of exact roots of
    (dH/dps = 0, dH/ds = 0, H = 0), and fsolve lands on it from any starting point: the first
    run returned pu = 0 and ps = 0 at every u, an action of exactly 0.000000, and a
    "displacement" that was just wherever the solver stopped drifting. Every equation was
    satisfied and the answer meant nothing.

    So pu is not solved for jointly. The fast pair is solved at FIXED pu, and the energy
    condition H = 0 is then a scalar equation in pu whose trivial root at 0 is excluded by the
    bracket.
    """
    def F(pv):
        r = fast_fixed(u, pv, g, guess)
        if r is None:
            return np.nan
        return H(u, r[0], pv, r[1], g)

    # SIGN. The escape momentum is NEGATIVE in this parametrisation -- the first scan swept
    # pu > 0 only and found nothing at all. The 1-D reduced chain says |pu| ~ |ln(mu/lam)|
    # (0.083 at gamma = 0.40, u = u*/2), and the non-trivial root sits near -0.07 there.
    # The action is therefore A = -int pu du, positive.
    grid = -np.geomspace(1e-4, pu_max, 80)
    prev_p, prev_f = None, None
    for pv in grid:
        fv = F(pv)
        if not np.isfinite(fv):
            continue
        if prev_f is not None and np.sign(fv) != np.sign(prev_f):
            try:
                root = brentq(F, prev_p, pv, xtol=1e-13)
            except (ValueError, RuntimeError):
                return None
            r = fast_fixed(u, root, g, guess)
            if r is None or abs(root) < 1e-4:
                return None
            return {"u": u, "s": r[0], "b": 1.0 - r[0], "ps": r[1], "pu": float(root),
                    "b_det": slow_manifold(u, g)}
        prev_p, prev_f = pv, fv
    return None


def solve_at(u, g, guess=None):
    return pu_at(u, g, guess)


def curve(g, n=60):
    """The predicted escape curve and momentum, swept in u by continuation."""
    us = u_star(g)
    out, guess = [], None
    for u in np.linspace(us * 0.02, us * 0.98, n):
        r = solve_at(u, g, guess)
        if r is None:
            continue
        guess = np.array([r["s"], r["ps"]])
        out.append(r)
    return out


def action(g, n=60):
    """A = int pu du over the predicted curve."""
    rows = curve(g, n)
    if len(rows) < 5:
        return None, rows
    uu = np.array([r["u"] for r in rows])
    pp = np.array([r["pu"] for r in rows])
    us = u_star(g)
    uu = np.concatenate(([0.0], uu, [us]))
    pp = np.concatenate(([0.0], pp, [0.0]))
    return float(-np.trapezoid(pp, uu)), rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/predictive_curve.json"))
    args = ap.parse_args()
    out = {}

    print("=== P1 GATE: is the zero-momentum sheet the deterministic dynamics?")
    worst_h, worst_m = 0.0, 0.0
    for g in (0.40, 0.44):
        us = u_star(g)
        for u in np.linspace(0.05 * us, 0.95 * us, 9):
            for s in (0.5, 0.62, 0.7):
                worst_h = max(worst_h, abs(H(u, s, 0.0, 0.0, g)))
            sdet = 1.0 - slow_manifold(u, g)
            worst_m = max(worst_m, abs(_grad(u, sdet, 0.0, 0.0, g)[0]))
    print(f"  max |H(u, s, 0, 0)| over 54 points          = {worst_h:.3e}")
    print(f"  max |dH/dps at p=0 on the slow manifold|    = {worst_m:.3e}")
    ok = worst_h < 1e-12 and worst_m < 1e-8
    print(f"  -> P1 {'HOLDS: the p = 0 branch IS the deterministic slow manifold, so §84/§85 solved the right equations at zero momentum' if ok else 'FAILS -- the Hamiltonian is written wrong'}")
    assert ok, "the zero-momentum sheet is not the deterministic dynamics"

    print("\n=== P2/P5: the PREDICTED curve against §86's measured ridge displacement")
    print(f"{'gamma':>7}{'u/u*':>8}{'b predicted':>13}{'b det':>10}{'displacement':>14}"
          f"{'pu':>9}{'ps':>10}")
    p2 = []
    for g in (0.40, 0.44):
        rows = curve(g, n=40)
        if not rows:
            print(f"{g:>7}   no solution found")
            continue
        us = u_star(g)
        sel = [rows[int(f * (len(rows) - 1))] for f in (0.15, 0.35, 0.55, 0.75, 0.9)]
        for r in sel:
            print(f"{g:>7}{r['u']/us:>8.2f}{r['b']:>13.5f}{r['b_det']:>10.5f}"
                  f"{r['b']-r['b_det']:>14.5f}{r['pu']:>9.4f}{r['ps']:>10.5f}")
        disp = float(np.mean([r["b"] - r["b_det"] for r in rows]))
        p2.append({"gamma": g, "disp_pred": disp, "disp_meas": RIDGE_DISP[g]})
        print(f"   mean predicted displacement {disp:+.5f}"
              f"   vs §86 MEASURED {RIDGE_DISP[g]:+.5f}"
              f"   ratio {disp/RIDGE_DISP[g]:.3f}")
    out["p2"] = p2
    if len(p2) == 2:
        signok = all(r["disp_pred"] > 0 for r in p2)
        magok = all(0.3 < r["disp_pred"] / r["disp_meas"] < 3.0 for r in p2)
        shrinks = p2[1]["disp_pred"] < p2[0]["disp_pred"]
        print(f"  -> P2/P5 {'HOLD: predicted displacement has the RIGHT SIGN, the right order, and shrinks toward gamma_c as the measured one does' if (signok and magok and shrinks) else ('FAILS on SIGN -- the ridge is not the instantons curve and §86s agreement was a coincidence of the marginal' if not signok else 'partial: sign right, magnitude or trend off')}")

    print("\n=== P3: the action, absolute, against the exact first passage")
    print(f"{'gamma':>7}{'measured':>11}{'on b_det':>11}{'ratio':>8}{'PREDICTED':>12}{'ratio':>8}")
    p3 = []
    for g in (0.40, 0.44):
        A, rows = action(g, n=60)
        m = MEASURED[g]
        if A is None:
            print(f"{g:>7}   no curve")
            continue
        p3.append({"gamma": g, "meas": m, "A_pred": A, "r_pred": A / m,
                   "r_det": DET_RATIO[g]})
        print(f"{g:>7}{m:>11.6f}{DET_RATIO[g]*m:>11.6f}{DET_RATIO[g]:>8.4f}"
              f"{A:>12.6f}{A/m:>8.4f}")
    out["p3"] = p3
    if len(p3) == 2:
        better = all(abs(1 - r["r_pred"]) < abs(1 - r["r_det"]) for r in p3)
        close = all(abs(1 - r["r_pred"]) < 0.03 for r in p3)
        print(f"  -> P3 {'HOLDS: the tilted-generator action beats the deterministic manifold and lands within 3%. The reduction is PREDICTIVE' if (better and close) else ('improves on the manifold but does not reach 3% -- the momentum is part of what was lost, not all of it' if better else 'FAILS: no better than the deterministic manifold, so the fast momentum is not what the elimination lost')}")

    print("\n=== P4: what this cost")
    print("  three equations (dH/dps = 0, dH/ds = 0, H = 0) solved by continuation in u.")
    print("  No master equation, no stationary solve, no lattice, no Omega. §86 needed the")
    print("  exact answer to find the curve; this needs only the rate functions.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
