"""T14-c-ii: is §28.1's +-10% scatter physics, or unmatched Omega grids?

§28.1 predicted the collapse slope from §15's closed forms at five gammas and got
ratios 1.155, 1.118, 0.952, 1.061, 1.044 -- within +-15% everywhere, but NOT monotone
in gamma, with the sign flipping at gamma = 0.25. That killed §28's 1-D-reduction
story (a reduction error tracking timescale separation cannot change sign) and left
the residual as unexplained scatter.

**A confound was introduced there and is the first thing to eliminate:** each gamma
used a DIFFERENT Omega list, hand-picked to reach a usable depth. So the cells are not
matched across gamma, each carries its own realised-eps rounding pattern, and that
wobble alone moved §27's raw local slopes by 40%. +-10% of scatter from unmatched
grids is entirely plausible without any physics at all.

THE MATCHING RULE, one rule for every gamma:

    every gamma spans the SAME probability range, 1e-2 down to 1e-6, with ln P
    equally spaced -- which, because the collapse is exponential in Omega, means
    equally spaced Omega between the two endpoints.

Endpoints come from §28.1's own measured slope and intercept per gamma, so the grid
is chosen by the data rather than by hand. Same dynamic range, same cell count, same
P window: the only difference between two gammas is gamma.

The precision floor is not needed here -- 1e-6 sits nine orders above the
cancellation limit that contaminated §28's gamma = 0.15 sweep -- but the eps control
in the fit is kept, since realised eps still wobbles cell to cell within a gamma.

PREDICTIONS, written before running:

  P1  The scatter SHRINKS materially -- to roughly +-5% or better -- and the sign
      flip at gamma = 0.25 disappears. That would make the residual an artifact of
      grid choice, and §15's closed form exact to the precision of the test.
  P2  If instead the scatter survives at +-10% with gamma = 0.25 still below 1, it is
      physics: `kappa delta*^2` carries a residual gamma-dependence that §15's closed
      form does not capture. That is the outcome that costs something, and it is the
      reason the grid rule is fixed before the run.
  P3  The mean ratio stays above 1 (prediction slightly too steep), continuing §22.4's
      finding that `kappa delta^2` is stiffer than the exact barrier. A mean that
      moves to 1.00 under matched grids would say the whole offset was grid choice
      too, which I do not expect.
  P4  Whatever happens, the shallower window (4 decades against §28.1's 5.5-9.9) makes
      each individual slope LESS precisely determined, not more. If the scatter grows,
      that is the likeliest reason and it is not evidence of physics -- so the
      per-gamma R^2 and slope standard error are reported alongside.
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

# gamma -> (slope, ln P at Omega = 40) from §28.1's measured sweeps
ANCHOR = {
    0.15: (-0.066953, np.log(1.52e-2)),
    0.20: (-0.050710, np.log(2.38e-2)),
    0.25: (-0.041823, np.log(5.26e-2)),
    0.30: (-0.024630, np.log(9.83e-2)),
    0.35: (-0.014605, np.log(1.58e-1)),
}
P_HI, P_LO, N_CELLS = 1e-2, 1e-6, 12


def omega_grid(gamma: float) -> list[int]:
    slope, ln_p40 = ANCHOR[gamma]
    o_hi = 40 + (np.log(P_HI) - ln_p40) / slope
    o_lo = 40 + (np.log(P_LO) - ln_p40) / slope
    return sorted({int(round(o)) for o in np.linspace(o_hi, o_lo, N_CELLS)})


def fit(gamma, eps, theta, omegas):
    ds = delta_star(gamma)
    om, lp, er, ps = [], [], [], []
    for o in omegas:
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
    if len(om) < 5:
        return None
    A = np.vstack([om, om * (er - er.mean()), np.ones_like(om)]).T
    c, res_, rank, _ = np.linalg.lstsq(A, lp, rcond=None)
    r = lp - A @ c
    dof = max(len(om) - 3, 1)
    cov = np.linalg.pinv(A.T @ A) * (r @ r) / dof
    return {"slope": float(c[0]), "se": float(np.sqrt(cov[0, 0])),
            "r2": float(1 - r.var() / lp.var()), "n_cells": len(om),
            "decades": float(np.log10(max(ps) / min(ps))),
            "p_hi": float(max(ps)), "p_lo": float(min(ps)),
            "omegas": [int(x) for x in om]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eps-frac", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/collapse_slope_matched.json"))
    args = ap.parse_args()

    print(f"matched grids: every gamma spans P = {P_HI:g} -> {P_LO:g}, "
          f"{N_CELLS} cells, ln P equally spaced\n")
    print(f"{'gamma':>6}{'Omega range':>14}{'cells':>6}{'decades':>8}"
          f"{'measured':>12}{'+-se':>9}{'predicted':>11}{'ratio':>8}{'R^2':>9}")
    rows = []
    for g in sorted(ANCHOR):
        grid = omega_grid(g)
        f = fit(g, args.eps_frac, args.theta, grid)
        if f is None:
            print(f"{g:>6.2f}   no usable cells"); continue
        net = am_reversible(g); x0 = args.eps_frac * delta_star(g)
        pred = -2.0 * V_exact(net, x0)
        ratio = pred / f["slope"]
        rows.append({"gamma": g, "predicted": pred, "ratio": ratio,
                     "sep": 3 * (1 + 2 * g) / (1 - 2 * g),
                     "kappa": wall_coefficient(g), **f})
        print(f"{g:>6.2f}{f['omegas'][0]:>7}-{f['omegas'][-1]:<6}{f['n_cells']:>6}"
              f"{f['decades']:>8.2f}{f['slope']:>12.6f}{f['se']:>9.6f}"
              f"{pred:>11.6f}{ratio:>8.3f}{f['r2']:>9.6f}")

    r = np.array([x["ratio"] for x in rows])
    print(f"\n  ratios: " + "  ".join(f"{x:.3f}" for x in r))
    print(f"  mean {r.mean():.4f}   sd {r.std(ddof=1):.4f}   "
          f"spread {r.min():.3f}-{r.max():.3f}")
    print(f"  §28.1 unmatched grids gave: mean 1.0660  sd 0.0770  "
          f"spread 0.952-1.155")
    print(f"  any ratio below 1? {'YES at gamma=' + str([x['gamma'] for x in rows if x['ratio'] < 1]) if (r < 1).any() else 'no'}")
    print(f"  monotone in gamma? {'yes' if np.all(np.diff(r) < 0) or np.all(np.diff(r) > 0) else 'NO'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
