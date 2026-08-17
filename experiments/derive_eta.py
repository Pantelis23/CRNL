"""T-DEPTH-e: derive eta. And T-DEPTH-e's own kill test named the wrong theory.

§77 measured eta = d ln(1/eps)/d Omega, nats of reliability per molecule, and left its
landscape dependence undertermined: 0.6813, 1.8346, 9.8813 at gamma = 0.30, 0.20, 0.05, plainly
not linear. **It proposed deriving eta from WKB as the barrier action per molecule. That is the
wrong theory for this eps.**

WKB gives the ESCAPE probability -- the chance the element spontaneously crosses its own saddle
-- and its exponent is the quasipotential barrier. But §75/§77's eps is not an escape: it is the
GAUSSIAN READOUT of the rail's own fluctuation, eps = Phi(-Delta/sigma), which is what a
chemically-coupled cascade actually applies (§75). The relevant object is therefore not the
barrier but the LINEAR-NOISE VARIANCE at the rail:

    sigma_x^2 = V / Omega    (van Kampen; V from the Lyapunov equation at the fixed point)
    ln(1/eps) ~ Delta^2 / (2 sigma^2) = Delta^2 Omega / (2V)
    **eta = Delta^2 / (2V)**

**If that is right, eta is computable from the DETERMINISTIC field and its linearisation
alone** -- no master equation, no stationary solve, no simulation, no entropy. The whole
founding question would then reduce to a quantity available from the ODE.

For a 1-D birth-death chain the Lyapunov equation is scalar: with drift f(x) and diffusion
b(x) = (lambda + mu)/Omega, the stationary variance at a stable fixed point is
V = b(x*) / (2 |f'(x*)|). For AM it is the 2 x 2 problem on the conserved simplex,
J C + C J^T + D = 0 with D = S diag(a) S^T, and V is the variance of the LEAD coordinate
x - y, i.e. u^T C u with u = (1, -1, 0).

PREDICTIONS, written before running. **P2 is absolute: no fit, no free parameter.**

  P1  GATE. The LNA prediction sigma^2 = V/Omega must reproduce the rail widths §75 measured
      from the exact stationary distribution, converging in Omega (rule 20 -- convergence, not
      a tolerance). Two independent routes to sigma; if they disagree the LNA is not the right
      description and P2 is meaningless.
  P2  **THE TEST. eta = Delta^2/(2V) against §77's measured 0.6813, 1.8346, 9.8813 at
      gamma = 0.30, 0.20, 0.05, and against Schloegl's 0.059574 and 0.015617.** Nothing is
      fitted: Delta is delta_star(gamma) or the root spacing, V comes from the Lyapunov
      equation, and the comparison is a ratio to a stored number (rule 16).
  P3  **THE CONSEQUENCE, if P2 holds.** eta -- and therefore reliability, and therefore depth,
      and therefore everything §73-§77 reduced the founding question to -- **is a property of
      the deterministic field and its Jacobian.** The master equation is needed to VALIDATE
      that, not to compute it.
  P4  **WHAT WOULD REFUTE IT.** A systematic gap growing toward gamma_c, where the rail
      flattens and the LNA's harmonic assumption fails. Predicted: the LNA is worst at
      gamma = 0.30 (nearest gamma_c = 0.5, shallowest rail) and best at gamma = 0.05. If the
      error instead grows toward SMALL gamma the diagnosis is wrong.
  P5  **RULE 9.** Check across substrates and across landscape shapes, not one axis: AM over
      gamma, Schloegl over root spacing AND over lambda, where §75's collapse demands
      eta ~ lambda exactly.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
from scipy.linalg import solve_lyapunov

from crnl.networks.am_reversible import am_reversible, delta_star
from experiments.cascade_schlogl import rates, schlogl_consts
from experiments.chemical_channel_noise import am_rail_width, rail_width

# §77's measured values, quoted as stored numbers (rule 16)
AM_ETA = {0.30: 0.6813, 0.20: 1.8346, 0.05: 9.8813}
SCH_ETA = {0.9: 0.059574, 0.6: 0.015617}


def schlogl_V(r1, r2, r3, omega=100000):
    """Scalar Lyapunov: V = b(x*)/(2|f'(x*)|) at the HIGH rail, in concentration units."""
    c = schlogl_consts(r1, r2, r3)
    k1a, k1r, k2b, k2r = c

    def f(x):
        return k1a * x ** 2 - k1r * x ** 3 + k2b - k2r * x

    def fp(x):
        return 2 * k1a * x - 3 * k1r * x ** 2 - k2r

    b = (k1a * r3 ** 2 + k1r * r3 ** 3 + k2b + k2r * r3)      # (lambda + mu)/Omega at x*
    assert abs(f(r3)) < 1e-10, f(r3)
    return b / (2 * abs(fp(r3)))


def am_V(gamma):
    """Lead-coordinate variance coefficient at AM's rail, from J C + C J^T + D = 0."""
    net = am_reversible(gamma)
    ds = float(delta_star(gamma))
    b = gamma / (1.0 + gamma)                    # §9.1: b* = gamma/(1+gamma), exact
    x = np.array([(1 - b + ds) / 2, (1 - b - ds) / 2, b])
    S = net.stoichiometry_matrix().astype(float)
    v = net.fluxes(x)
    assert np.abs(S @ v).max() < 1e-9, np.abs(S @ v).max()

    h = 1e-6
    J = np.zeros((3, 3))
    for j in range(3):
        xp, xm = x.copy(), x.copy()
        xp[j] += h; xm[j] -= h
        J[:, j] = (S @ net.fluxes(xp) - S @ net.fluxes(xm)) / (2 * h)
    D = S @ np.diag(v) @ S.T

    # the conserved direction makes J singular; work in the 2-D subspace orthogonal to (1,1,1)
    Q, _ = np.linalg.qr(np.array([[1.0, 1, 1]]).T, mode="complete")
    P = Q[:, 1:]                                   # basis of the conserved-total subspace
    Jr, Dr = P.T @ J @ P, P.T @ D @ P
    Cr = solve_lyapunov(Jr, -Dr)
    C = P @ Cr @ P.T
    u = np.array([1.0, -1.0, 0.0])
    return float(u @ C @ u), ds


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/derive_eta.json"))
    args = ap.parse_args()

    print("=== P1 GATE: does the LNA reproduce the exact rail widths, converging in Omega?")
    print(f"{'element':>18}{'Omega':>8}{'sigma exact':>13}{'sqrt(V/Om)':>12}{'ratio':>9}")
    conv = {}
    for g in (0.20, 0.05):
        V, ds = am_V(g)
        ser = []
        for om in (60, 120, 240, 400):
            try:
                sd, _ = am_rail_width(g, om)
            except RuntimeError:
                continue
            pred = np.sqrt(V / om)
            ser.append(abs(sd / pred - 1))
            print(f"{f'AM gamma={g}':>18}{om:>8}{sd:>13.6f}{pred:>12.6f}{sd/pred:>9.4f}")
        conv[f"AM {g}"] = ser
    for spread in (0.9,):
        r1, r3 = 1.0 - spread, 1.0 + spread
        V = schlogl_V(r1, 1.0, r3)
        ser = []
        for om in (1600, 6400, 25600):
            r = rail_width(om, r1, 1.0, r3)
            pred = np.sqrt(V / om)
            ser.append(abs(r["sd_exact"] / pred - 1))
            print(f"{f'Schlogl s={spread}':>18}{om:>8}{r['sd_exact']:>13.6f}{pred:>12.6f}"
                  f"{r['sd_exact']/pred:>9.4f}")
        conv[f"Sch {spread}"] = ser
    ok = all(len(v) >= 2 and v[-1] < v[0] for v in conv.values())
    for k, v in conv.items():
        print(f"  {k}: |ratio-1| " + " -> ".join(f"{x:.4f}" for x in v))
    print(f"  -> P1 {'HOLDS: the LNA discrepancy shrinks with Omega on every element' if ok else 'FAILS'}")

    print("\n=== P2: eta = Delta^2/(2V), against §77's measured values. No fit.")
    print(f"{'element':>18}{'Delta':>9}{'V':>12}{'eta predicted':>15}"
          f"{'eta measured':>14}{'ratio':>9}")
    rows = []
    for g, meas in sorted(AM_ETA.items()):
        V, ds = am_V(g)
        pred = ds ** 2 / (2 * V)
        rows.append({"el": f"AM gamma={g}", "pred": pred, "meas": meas,
                     "ratio": pred / meas})
        print(f"{f'AM gamma={g}':>18}{ds:>9.5f}{V:>12.6f}{pred:>15.4f}{meas:>14.4f}"
              f"{pred/meas:>9.4f}")
    for spread, meas in sorted(SCH_ETA.items()):
        r1, r3 = 1.0 - spread, 1.0 + spread
        V = schlogl_V(r1, 1.0, r3)
        D = spread
        pred = D ** 2 / (2 * V)
        rows.append({"el": f"Schlogl s={spread}", "pred": pred, "meas": meas,
                     "ratio": pred / meas})
        print(f"{f'Schlogl s={spread}':>18}{D:>9.5f}{V:>12.6f}{pred:>15.6f}{meas:>14.6f}"
              f"{pred/meas:>9.4f}")
    worst = max(abs(r["ratio"] - 1) for r in rows)
    print(f"  worst |predicted/measured - 1| = {100*worst:.2f}%")
    print(f"  -> P2 {'HOLDS: eta is derived, not fitted' if worst < 0.10 else 'FAILS'}")

    print("\n=== P4: is the error worst where the LNA should be worst (rails shallow)?")
    am = [r for r in rows if r["el"].startswith("AM")]
    print("  |ratio-1| by gamma: " + ", ".join(
        f"{r['el'].split('=')[1]}: {100*abs(r['ratio']-1):.2f}%" for r in am))
    trend = abs(am[-1]["ratio"] - 1) < abs(am[0]["ratio"] - 1)
    print(f"  -> {'as predicted: worst at the gamma nearest gamma_c, where the rail is shallowest' if not trend else 'the error grows toward SMALL gamma -- the diagnosis in P4 is wrong'}")

    print("\n=== P5 (rule 9): does eta ~ lambda exactly, as §75's collapse demands?")
    for lam in (1.0, 2.0, 4.0):
        V = schlogl_V(lam * 0.1, lam * 1.0, lam * 1.9)
        D = lam * 0.9
        print(f"  lambda={lam:>4}: V = {V:.6f}, eta = {D**2/(2*V):.6f}")
    v1 = 0.81 / (2 * schlogl_V(0.1, 1.0, 1.9))
    v4 = (4 * 0.9) ** 2 / (2 * schlogl_V(0.4, 4.0, 7.6))
    print(f"  eta(lambda=4)/eta(lambda=1) = {v4/v1:.4f}   (§75's collapse requires 4)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"rows": rows, "conv": conv}, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
