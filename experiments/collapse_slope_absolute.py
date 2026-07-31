"""T14-c: predict §27's collapse slope from closed forms — absolutely, not by fitting

§27 measured `ln P(error) = -0.024904 * Omega + ...` over 6.53 decades, exactly. That
is a FITTED slope. Rule 16: a law that is only ever fitted is never tested. §15 gives
`kappa(gamma) = (3/2)(1-2g)/(1+g)` and `delta_star(gamma)` in closed form, and
`breaking_diffusion` gives `D_0`, so the slope is computable with **no free
parameter** and the ratio predicted/measured is the test.

THE PREDICTION. For a 1-D diffusion `d(delta) = mu dt + sqrt(D/Omega) dW` with
absorbing barriers, the splitting probability is governed by the scale function
`S'(x) = exp(-2 Omega V(x))` with

    V(x) = integral from 0 to x of mu(d) / D(d) dd

Laplace: the numerator is dominated by its endpoint at the start `x0` and the
denominator by the saddle at 0, so

    d(ln P) / d(Omega)  =  -2 V(x0)                     [EXACT form]

and near the saddle, where `mu ~ lambda*delta` and `D ~ D(0)`, that collapses to

    d(ln P) / d(Omega)  =  -kappa * x0^2                [QUADRATIC form]

**The normalisation is the trap here and is worth stating.** `breaking_mode` is a
UNIT vector, so `D_0 = sum_r (v.S_r)^2 f_r` is computed with `v = (1,-1,0)/sqrt(2)`,
whereas the coordinate `delta = x - y` corresponds to the UNNORMALISED `u = (1,-1,0)`.
So `D_delta = 2 D_0`, and `kappa = lambda/(2 D_0) = lambda/D_delta` means the
near-saddle exponent is `V ~ lambda x^2/(2 D_delta) = kappa x^2 / 2`, giving
`-2V = -kappa x0^2` and NOT `-2 kappa x0^2`. Getting that factor wrong turns a 12%
agreement into a 2.3x failure, and it would have been reported as physics. The
quadratic limit of the exact integral is checked numerically below as a guard.

The 1-D reduction slaves the pool: for each `delta`, `s = x + y` is solved so that
`ds/dt = 0`, and `mu`, `D` are evaluated there. That is the same reduction §15 used
and the one §24.1a showed carries a subleading error that shrinks with the timescale
separation.

PREDICTIONS, written before running:

  P1  The EXACT scale-function integral beats the quadratic, because the quadratic
      linearises a drift evaluated 35% of the way to the attractor.
  P2  The quadratic OVER-estimates the slope magnitude (predicts too-steep decay).
      §22.4 established `kappa*delta^2` is stiffer than the exact barrier away from
      the saddle, so its exponent should be too large. If it comes out too SHALLOW,
      the §22.4 reading does not transfer to this observable and that is worth more
      than the agreement.
  P3  The exact form lands within ~10% of measured at gamma = 0.30. A hand estimate
      put the quadratic near 12%, so this is not a blind prediction for that arm --
      stated so the quadratic's agreement is not presented as a clean forecast.
  P4  WHAT WOULD FALSIFY THE WHOLE PICTURE: a ratio that DRIFTS systematically with
      gamma. A constant offset means a missing prefactor; a drift means the
      `kappa delta*^2` reading of the barrier is wrong in a gamma-dependent way, and
      §15's closed form would not survive it.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
from scipy.optimize import brentq

from crnl.information import wall_coefficient
from crnl.networks.am_reversible import am_reversible, delta_star
from crnl.networks.n_winner_reversible import breaking_diffusion, lambda_breaking
from experiments.approximation_hierarchy import _setup, p_cme

U = np.array([1.0, -1.0, 0.0])          # the delta = x - y direction, UNNORMALISED


def slaved_state(net, delta: float, lo: float = 1e-9, hi: float = 1.0 - 1e-9):
    """(x, y, b) at separation `delta` with the pool on its own nullcline."""
    S = net.stoichiometry_matrix()

    def ds_dt(s):
        x, y = 0.5 * (s + delta), 0.5 * (s - delta)
        if x < 0 or y < 0 or 1.0 - s < 0:
            return np.nan
        d = S @ net.fluxes(np.array([x, y, 1.0 - s]))
        return d[0] + d[1]

    a, b = max(abs(delta) + 1e-9, lo), hi
    fa, fb = ds_dt(a), ds_dt(b)
    if not (np.isfinite(fa) and np.isfinite(fb)) or fa * fb > 0:
        return None
    s = brentq(ds_dt, a, b, xtol=1e-14)
    return np.array([0.5 * (s + delta), 0.5 * (s - delta), 1.0 - s])


def mu_D(net, delta: float):
    """Drift and diffusion of `delta` on the slaved manifold (unnormalised mode)."""
    st = slaved_state(net, delta)
    if st is None:
        return None
    S = net.stoichiometry_matrix()
    f = net.fluxes(st)
    mu = float(U @ (S @ f))
    D = float(sum((U @ S[:, r]) ** 2 * f[r] for r in range(net.n_reactions)))
    return mu, D


def V_exact(net, x0: float, n: int = 4001) -> float:
    """V(x0) = int_0^x0 mu/D, by direct quadrature on the slaved manifold."""
    xs = np.linspace(1e-6, x0, n)
    vals = []
    for x in xs:
        md = mu_D(net, float(x))
        if md is None or md[1] <= 0:
            return float("nan")
        vals.append(md[0] / md[1])
    return float(np.trapezoid(np.array(vals), xs))


def measured_slope(gamma: float, eps: float, theta: float, omegas) -> dict:
    """eps-controlled fit, as §27: the integer lattice makes realised eps wobble."""
    ds = delta_star(gamma)
    om, lp, er = [], [], []
    for o in omegas:
        try:
            p = p_cme(gamma, o, eps, theta)
        except Exception:
            continue
        if not np.isfinite(p) or p <= 0:
            continue
        n0, _thr, _ = _setup(gamma, o, eps, theta)
        om.append(float(o))
        lp.append(np.log(p))
        er.append(int(n0[0] - n0[1]) / (ds * o))
    om, lp, er = np.array(om), np.array(lp), np.array(er)
    if len(om) < 4:
        return {"slope": float("nan"), "n": len(om)}
    A = np.vstack([om, om * (er - er.mean()), np.ones_like(om)]).T
    c, *_ = np.linalg.lstsq(A, lp, rcond=None)
    r = lp - A @ c
    return {"slope": float(c[0]), "r2": float(1 - r.var() / lp.var()),
            "n": len(om), "eps_mean": float(er.mean()),
            "decades": float(np.log10(np.exp(lp.max() - lp.min())))}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.15, 0.30, 0.45])
    ap.add_argument("--eps-frac", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--omegas", type=int, nargs="+",
                    default=[40, 80, 120, 160, 200, 260, 340, 420, 500])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/collapse_slope_absolute.json"))
    args = ap.parse_args()

    print(f"eps/delta* = {args.eps_frac}, theta = {args.theta}; "
          f"exact CME sweep per gamma\n")
    print(f"{'gamma':>6} {'measured':>11} {'quadratic':>11} {'ratio':>7} "
          f"{'exact int':>11} {'ratio':>7} {'decades':>8} {'R^2':>9}")
    rows = []
    for g in args.gammas:
        net = am_reversible(g)
        ds = delta_star(g)
        x0 = args.eps_frac * ds
        kap = wall_coefficient(g)

        m = measured_slope(g, args.eps_frac, args.theta, args.omegas)
        quad = -kap * x0 ** 2
        vex = V_exact(net, x0)
        exact = -2.0 * vex

        # guard: the quadratic limit of the exact integral must reproduce -kappa x^2
        xs = 1e-3 * ds
        guard = (-2.0 * V_exact(net, xs, 401)) / (-kap * xs ** 2)

        rq = quad / m["slope"] if m["slope"] else float("nan")
        rx = exact / m["slope"] if m["slope"] else float("nan")
        print(f"{g:>6.2f} {m['slope']:>11.6f} {quad:>11.6f} {rq:>7.3f} "
              f"{exact:>11.6f} {rx:>7.3f} {m.get('decades', float('nan')):>8.2f} "
              f"{m.get('r2', float('nan')):>9.6f}")
        rows.append({"gamma": g, "kappa": kap, "delta_star": ds, "x0": x0,
                     "measured": m["slope"], "quadratic": quad, "exact": exact,
                     "ratio_quad": rq, "ratio_exact": rx,
                     "near_saddle_guard": guard, **m})

    print(f"\n  near-saddle guard (exact integral / quadratic as x->0, must be ~1): "
          + "  ".join(f"{r['near_saddle_guard']:.4f}" for r in rows))
    for key, nm in (("ratio_quad", "quadratic"), ("ratio_exact", "exact integral")):
        v = np.array([r[key] for r in rows])
        if np.isfinite(v).all():
            print(f"  {nm:>15} predicted/measured: mean {v.mean():.3f}, "
                  f"spread {v.min():.3f}-{v.max():.3f}  "
                  f"(P4: a DRIFT with gamma falsifies; a constant offset is a prefactor)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
