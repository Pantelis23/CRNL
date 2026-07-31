"""T14-c-iii + the eps axis: how far does §15's closed form actually reach?

§28.2 tested the parameter-free prediction of the collapse slope at five gammas on
MATCHED grids and found every ratio above 1, drifting monotonically toward 1 as the
timescale separation rises: `excess = 5.987 * sep^(-2.031)`, R^2 = 0.887. Two things
were left open and this run does both.

  (i) T14-c-iii -- SMALL GAMMA. The power law is carried by five points over a factor
      of 3 in `sep`, and the first two gammas differ by 0.006 in excess where sep^-2
      predicts ~37% between them. At gamma = 0.05 and 0.10 the separation falls to
      3.67 and 4.50, so `sep^-2` demands a LARGE excess -- ratios near 1.44 and 1.30.
      If the excess saturates near §28.2's 0.13-0.14 instead, the reduction error has
      a floor and the inverse-square reading is wrong.

  (ii) THE EPS AXIS, which the closed form has never faced. Every test so far used
      eps/delta* = 0.35. The 1-D reduction is a property of the slaved MANIFOLD, not
      of where on it the trajectory starts, so if the residual really is the reduction
      then the ratio at fixed gamma should be roughly EPS-INDEPENDENT. If instead it
      moves strongly with eps, the residual is not (only) the reduction, and §28.2's
      attribution is incomplete.

SELF-CALIBRATING GRIDS. §28.2 took its Omega endpoints from previously measured
slopes, which does not extend to cells never measured. Here each cell BISECTS on
Omega to find where P ~ 1e-2 and P ~ 1e-6, then places 12 cells equally spaced
between them. Same rule everywhere, no hand-picked lists, and the achieved P window
is reported so matching can be checked rather than assumed -- §28.1's whole error was
an unmatched grid nobody had checked.

Cells whose required Omega exceeds OMEGA_CAP are SKIPPED AND REPORTED, not quietly
dropped: at gamma = 0.35, eps = 0.20 the collapse is so shallow that four decades
would need Omega ~ 2000 (~2M states).

PREDICTIONS, written before running:

  P1  At gamma = 0.05 and 0.10 the excess GROWS as sep^-2 demands (ratios ~1.44,
      ~1.30). This is the prediction most likely to fail, because §28.2's own first
      two points already sit flatter than the law.
  P2  At fixed gamma the ratio is roughly eps-independent -- within a few percent
      between eps = 0.35 and 0.50 -- because the reduction is a property of the
      manifold. Strong eps-dependence would mean §28.2's attribution is incomplete.
  P3  Every ratio stays above 1, continuing §28.2 and §22.4. A ratio below 1 anywhere
      that is NOT a grid artifact would mean the prediction is too shallow somewhere,
      which no account here allows.
  P4  If P1 fails and P2 holds, the honest reading is that the residual IS the
      reduction but does not follow a clean power law in `sep` -- and §28.2's exponent
      should be withdrawn as a description while its attribution survives.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from crnl.information import wall_coefficient
from crnl.networks.am_reversible import am_reversible, delta_star
from experiments.approximation_hierarchy import _setup, p_cme
from experiments.collapse_slope_absolute import V_exact

P_HI, P_LO, N_CELLS, OMEGA_CAP = 1e-2, 1e-6, 12, 900


def find_omega(gamma, eps, theta, target, lo=20, hi=OMEGA_CAP):
    """Smallest Omega with P(Omega) <= target. None if beyond the cap."""
    def P(o):
        try:
            v = p_cme(gamma, int(o), eps, theta)
        except Exception:
            return np.nan
        return v if (np.isfinite(v) and v > 0) else np.nan
    phi = P(hi)
    if not np.isfinite(phi) or phi > target:
        return None
    while hi - lo > 4:
        mid = (lo + hi) // 2
        pm = P(mid)
        if not np.isfinite(pm):
            lo = mid
        elif pm > target:
            lo = mid
        else:
            hi = mid
    return int(hi)


def fit_cell(gamma, eps, theta):
    o_hi = find_omega(gamma, eps, theta, P_HI)
    o_lo = find_omega(gamma, eps, theta, P_LO)
    if o_hi is None or o_lo is None or o_lo <= o_hi:
        return None
    grid = sorted({int(round(o)) for o in np.linspace(o_hi, o_lo, N_CELLS)})
    ds = delta_star(gamma)
    om, lp, er, ps = [], [], [], []
    for o in grid:
        try:
            p = p_cme(gamma, o, eps, theta)
        except Exception:
            continue
        if not np.isfinite(p) or p <= 0 or p >= 1.0:
            continue
        n0, _t, _ = _setup(gamma, o, eps, theta)
        om.append(float(o)); lp.append(np.log(p)); ps.append(p)
        er.append(int(n0[0] - n0[1]) / (ds * o))
    om, lp, er = np.array(om), np.array(lp), np.array(er)
    if len(om) < 6:
        return None
    A = np.vstack([om, om * (er - er.mean()), np.ones_like(om)]).T
    c, *_ = np.linalg.lstsq(A, lp, rcond=None)
    r = lp - A @ c
    dof = max(len(om) - 3, 1)
    se = float(np.sqrt((np.linalg.pinv(A.T @ A) * (r @ r) / dof)[0, 0]))
    return {"slope": float(c[0]), "se": se, "r2": float(1 - r.var() / lp.var()),
            "n_cells": len(om), "o_lo": int(om.min()), "o_hi": int(om.max()),
            "decades": float(np.log10(max(ps) / min(ps))),
            "p_hi": float(max(ps)), "p_lo": float(min(ps))}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35])
    ap.add_argument("--eps-fracs", type=float, nargs="+", default=[0.35, 0.50])
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/collapse_slope_grid.json"))
    args = ap.parse_args()

    print(f"self-calibrating matched grids: every cell spans P = {P_HI:g} -> {P_LO:g}"
          f", {N_CELLS} points, Omega cap {OMEGA_CAP}\n")
    print(f"{'gamma':>6}{'eps':>6}{'sep':>7}{'Omega':>13}{'dec':>6}"
          f"{'measured':>12}{'predicted':>11}{'ratio':>8}{'excess':>9}{'R^2':>9}")
    rows, skipped = [], []
    for eps in args.eps_fracs:
        for g in args.gammas:
            f = fit_cell(g, eps, args.theta)
            if f is None:
                skipped.append((g, eps))
                print(f"{g:>6.2f}{eps:>6.2f}   SKIPPED (needs Omega > {OMEGA_CAP})")
                continue
            net = am_reversible(g)
            pred = -2.0 * V_exact(net, eps * delta_star(g))
            ratio = pred / f["slope"]
            sep = 3 * (1 + 2 * g) / (1 - 2 * g)
            rows.append({"gamma": g, "eps": eps, "sep": sep, "predicted": pred,
                         "ratio": ratio, "excess": ratio - 1.0,
                         "kappa": wall_coefficient(g), **f})
            print(f"{g:>6.2f}{eps:>6.2f}{sep:>7.2f}{f['o_lo']:>7}-{f['o_hi']:<5}"
                  f"{f['decades']:>6.2f}{f['slope']:>12.6f}{pred:>11.6f}"
                  f"{ratio:>8.3f}{ratio-1:>9.4f}{f['r2']:>9.6f}")

    if skipped:
        print(f"\n  skipped {len(skipped)} cells (reported, not dropped): {skipped}")

    print(f"\n=== P1: does the excess keep growing as sep^-2 at small gamma?")
    print(f"  §28.2 law: excess = 5.987 * sep^-2.031")
    for r in sorted(rows, key=lambda x: (x["eps"], x["gamma"])):
        pl = 5.987 * r["sep"] ** -2.031
        print(f"   gamma={r['gamma']:.2f} eps={r['eps']:.2f} sep={r['sep']:>5.2f}  "
              f"excess={r['excess']:>7.4f}   law predicts {pl:>7.4f}   "
              f"ratio {r['excess']/pl if pl else float('nan'):>6.2f}")

    print(f"\n=== P2: is the ratio eps-independent at fixed gamma?")
    for g in args.gammas:
        cell = {r["eps"]: r["ratio"] for r in rows if r["gamma"] == g}
        if len(cell) >= 2:
            v = list(cell.values())
            print(f"   gamma={g:.2f}: " + "  ".join(f"eps={k}:{x:.3f}" for k, x in cell.items())
                  + f"   spread {100*(max(v)-min(v))/np.mean(v):>5.1f}%")

    below = [r for r in rows if r["ratio"] < 1]
    print(f"\n=== P3: ratios below 1: {len(below)}" +
          (f"  -> {[(r['gamma'], r['eps'], round(r['ratio'],3)) for r in below]}" if below else "  (none)"))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
