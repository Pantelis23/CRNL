"""Is §15's closed form right for the wrong reason? — the Gaussian truncation, tested

§35 established the largest open quantitative fact here: §15's parameter-free prediction
is **7.5-15.5% too steep** against the asymptotic collapse rate, at every gamma. §28.3
named the 1-D slaved reduction as the suspect and could never test it absolutely.

**THERE IS A SECOND APPROXIMATION IN THAT FORMULA AND NOBODY HAS NAMED IT.** `V_exact`
computes `V = int mu/D`, which is the quasipotential of a DIFFUSION -- the CLE / van
Kampen / Gaussian approximation to the master equation. But a master equation's true
large-deviation function comes from the full WKB Hamiltonian

    H(x, p) = sum_r a_r(x) * (exp(p . S_r) - 1)

and truncating that at second order in `p` is exactly what the diffusion approximation
does. It is a long-known result that this gets rare-event exponents wrong by an O(1)
factor, and this project has been quoting the truncated version as a first-principles
prediction since §15.

**FOR delta THE REDUCTION IS EXACT AND THE COMPARISON IS ONE LINE.** With U = (1,-1,0),
only four of the six reactions move delta, by exactly +-1:

    lambda (up)   = f2 + r3        mu_down (down) = f3 + r2

so `mu = lambda - mu_down` and `D = lambda + mu_down`, and

    Gaussian:   rate = -2 * int (lambda - mu_down)/(lambda + mu_down)  d(delta)
    exact WKB:  rate = -1 * int ln( lambda / mu_down )                 d(delta)

These agree to first order about the saddle, where lambda ~ mu_down, and diverge as the
ratio leaves 1. **Since ln(r) > 2(r-1)/(r+1) for every r > 1, the WKB barrier is
STRICTLY LARGER than the Gaussian one** -- so the exact-Hamiltonian prediction should be
even STEEPER than §15's, which is already too steep.

**IF THAT HOLDS, §15's ACCURACY IS PARTLY ACCIDENTAL.** Two approximations sit in that
formula with OPPOSITE signs: the Gaussian truncation makes the barrier too small
relative to true WKB, and the 1-D slaving makes it too large relative to the true 2-D
minimum-action path (a minimum over all paths cannot exceed the value along the slaved
one). §15 keeps both and lands 7.5-15.5% high. Remove only the Gaussian one and it
should get WORSE -- which would mean the 1-D reduction error is substantially larger
than §28.3's estimate, and that the two have been cancelling.

PREDICTIONS, written before running:

  P0  GATE. My lambda/mu_down decomposition must reproduce `mu_D`'s own `mu` and `D` to
      machine precision, and my Gaussian quadrature must reproduce `V_exact` to the same.
      If either fails, the decomposition is wrong and nothing below is admissible.
  P1  |rate_WKB| > |rate_Gauss| > |c_measured| at every gamma -- the exact Hamiltonian
      moves the prediction AWAY from the truth. This is the sign test and it follows
      from ln(r) > 2(r-1)/(r+1) alone, so failing it means the decomposition or the
      conventions are wrong, not the physics.
  P2  THE CONSEQUENCE, and the reason this is worth running. If P1 holds, the 1-D
      slaving error is not 7.5-15.5% but |rate_WKB|/|c| - 1, i.e. LARGER, and §15's
      closed form owes part of its agreement to a cancellation between two errors of
      opposite sign. **That reframes the project's central closed form from "accurate"
      to "accurate by cancellation", which is a different claim with a different
      lifetime.**
  P3  The WKB-Gaussian gap grows as the barrier grows (smaller gamma, larger x0),
      because the ratio lambda/mu_down travels further from 1. So the cancellation is
      gamma-dependent -- which would EXPLAIN why the residual is gamma-dependent at all,
      something §28.3 could only describe with a straight line and could not account for.
  P4  If instead |rate_WKB| lands CLOSER to `c` than the Gaussian does, the Gaussian
      truncation was the dominant error all along, the 1-D reduction is nearly exact,
      and §28.3's attribution to the slaved manifold is wrong. That is the opposite
      outcome and it would be a cleaner result than P2 -- so both are worth the run.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from crnl.networks.am_reversible import am_reversible, delta_star
from experiments.collapse_slope_absolute import U, V_exact, mu_D, slaved_state

# measured asymptotic rates from §35 (prefactor ansatz), and §28.3's effective slopes
MEASURED = {0.20: -0.049064, 0.25: -0.035156, 0.30: -0.023643, 0.35: -0.014195}


def up_down(net, delta: float):
    """(lambda, mu_down): the +1 and -1 jump rates of delta on the slaved manifold."""
    st = slaved_state(net, delta)
    if st is None:
        return None
    f = net.fluxes(st)
    S = net.stoichiometry_matrix()
    up = dn = 0.0
    for r in range(net.n_reactions):
        j = float(U @ S[:, r])
        if j > 0.5:
            up += float(f[r])
        elif j < -0.5:
            dn += float(f[r])
    return up, dn


def rates(net, x0: float, n: int = 4001):
    xs = np.linspace(1e-6, x0, n)
    g, w = [], []
    for x in xs:
        ud = up_down(net, float(x))
        if ud is None or min(ud) <= 0:
            return None
        up, dn = ud
        g.append((up - dn) / (up + dn))
        w.append(np.log(up / dn))
    return (-2.0 * float(np.trapezoid(g, xs)), -1.0 * float(np.trapezoid(w, xs)))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.20, 0.25, 0.30, 0.35])
    ap.add_argument("--eps-frac", type=float, default=0.35)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/wkb_vs_gaussian.json"))
    args = ap.parse_args()

    print("=== P0 gate: does the up/down split reproduce mu_D and V_exact?")
    worst_md = worst_v = 0.0
    for g in args.gammas:
        net = am_reversible(g)
        for d in (0.05, 0.15, 0.30, 0.50):
            md, ud = mu_D(net, d), up_down(net, d)
            if md is None or ud is None:
                continue
            worst_md = max(worst_md, abs(md[0] - (ud[0] - ud[1])) / abs(md[0]),
                           abs(md[1] - (ud[0] + ud[1])) / abs(md[1]))
        x0 = args.eps_frac * delta_star(g)
        mine = rates(net, x0)
        worst_v = max(worst_v, abs(mine[0] - (-2.0 * V_exact(net, x0)))
                      / abs(2.0 * V_exact(net, x0)))
    print(f"  worst relative mismatch, mu and D:  {worst_md:.3e}")
    print(f"  worst relative mismatch, Gaussian rate vs -2*V_exact: {worst_v:.3e}")
    ok = worst_md < 1e-12 and worst_v < 1e-6
    print(f"  -> P0 {'HOLDS' if ok else 'FAILS'}")
    if not ok:
        raise SystemExit("gate failed; nothing below is admissible")

    print(f"\n=== P1/P2: the Gaussian truncation, removed")
    print(f"{'gamma':>7}{'measured c':>13}{'Gaussian':>12}{'exact WKB':>12}"
          f"{'G/meas':>9}{'WKB/meas':>10}{'WKB/G':>9}")
    rows = []
    for g in args.gammas:
        net = am_reversible(g)
        x0 = args.eps_frac * delta_star(g)
        rg, rw = rates(net, x0)
        c = MEASURED[g]
        rows.append({"gamma": g, "measured": c, "gauss": rg, "wkb": rw,
                     "g_over_meas": rg / c, "wkb_over_meas": rw / c,
                     "wkb_over_g": rw / rg})
        print(f"{g:>7.2f}{c:>13.6f}{rg:>12.6f}{rw:>12.6f}"
              f"{rg/c:>9.4f}{rw/c:>10.4f}{rw/rg:>9.4f}")

    print(f"\n  P1: is |WKB| > |Gaussian| > |measured| everywhere?")
    p1 = all(abs(r["wkb"]) > abs(r["gauss"]) > abs(r["measured"]) for r in rows)
    print(f"    {'YES' if p1 else 'NO'}")
    closer = [r for r in rows if abs(r["wkb"] / r["measured"] - 1)
              < abs(r["gauss"] / r["measured"] - 1)]
    print(f"  P4: cells where the exact WKB is CLOSER to the truth than the Gaussian: "
          f"{len(closer)}/{len(rows)}")

    print(f"\n=== P2: what the 1-D slaving error must then be")
    print(f"{'gamma':>7}{'§28.3 said':>12}{'§35 asympt':>12}"
          f"{'implied 1-D error':>20}{'Gaussian defect':>17}")
    for r in rows:
        implied = abs(r["wkb"] / r["measured"]) - 1.0
        gauss_def = 1.0 - abs(r["gauss"] / r["wkb"])
        print(f"{r['gamma']:>7.2f}{'-':>12}{abs(r['g_over_meas'])-1:>12.4f}"
              f"{implied:>20.4f}{gauss_def:>17.4f}")

    print(f"\n=== P3: does the WKB-Gaussian gap grow as the barrier grows?")
    gaps = np.array([r["wkb_over_g"] for r in rows])
    gs = np.array([r["gamma"] for r in rows])
    print(f"  WKB/Gaussian by gamma: "
          + "  ".join(f"{g:.2f}:{x:.4f}" for g, x in zip(gs, gaps)))
    print(f"  monotone decreasing in gamma? "
          f"{'yes' if np.all(np.diff(gaps) < 0) else 'NO'}"
          f"   (deeper barrier -> larger gap)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
