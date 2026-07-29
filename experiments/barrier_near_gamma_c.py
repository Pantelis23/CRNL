"""T4's kill test: how does the barrier vanish as the landscape dies?

T4 guessed the restoration barrier falls like `(gamma_c - gamma)`, one vanishing
factor from `delta*^2`. FINDINGS 12 implied TWO -- `kappa * delta*^2` with both
`kappa` and `delta*^2` linear in `(gamma_c - gamma)` -- hence `(gamma_c - gamma)^2`,
and the population needed for fixed reliability diverging like `1/(gamma_c-gamma)^2`.
That has only ever been INFERRED from a collapse, never measured by sweeping gamma.

This measures it directly with the exact quasipotential (FINDINGS 15), whose usable
window -- large gamma, small barrier -- happens to be exactly this region.

WHY IT IS A 2-FOR-1. FINDINGS 15's correction rests on `kappa = lambda/(2 D_0)` with
`lambda = (1-2 gamma)/3` vanishing LINEARLY at gamma_c. If the barrier exponent is
not 2, either delta*'s square-root pitchfork or lambda's linear vanishing is wrong,
and 15 breaks with T4. THEORIES says they stand or fall together; this is that test.

PREDICTION, written before running: **exponent 2.** Near a pitchfork the
quasipotential is quartic, `W = -a x^2/2 + b x^4/4` with `a` proportional to
`(gamma_c - gamma)`, so the barrier is `a^2/(4b)` -- quadratic. A measured 1
restores T4's original guess and falsifies the kappa correction; anything else
means the pitchfork normal form does not describe this landscape.

THE INSTRUMENT'S WINDOW IS WHAT SETS Omega, and it closes from both sides
(FINDINGS 15.1): the probability floor caps the resolvable barrier at ~30/Omega, so
Omega must be SMALL, while the lattice must resolve a landscape of width delta*, so
Omega must be LARGE. As gamma -> gamma_c the barrier shrinks fast (predicted
delta*^4) and delta* shrinks slowly (delta*^2), so the window OPENS -- which is why
this region is the instrument's best and the rest of the project's worst.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.networks.am_reversible import GAMMA_C, am_reversible, delta_star
from crnl.quasipotential import P_FLOOR, barriers

SITES_MIN = 60          # lattice steps across the landscape width
FLOOR_SAFETY = 25.0     # keep Omega*dW below this (the hard cap is ~30)


def omega_window(gamma: float, barrier_guess: float) -> tuple[int, int]:
    ds = delta_star(gamma)
    lo = int(np.ceil(SITES_MIN / ds))
    hi = int(FLOOR_SAFETY / max(barrier_guess, 1e-12))
    return lo, hi


def barrier_at(gamma: float, omega: int) -> float:
    b = barriers(am_reversible(gamma), int(omega), float(omega))
    return 0.5 * (b["dW_x"] + b["dW_y"])      # equal by symmetry; averaged anyway


def converged_barrier(gamma: float, omegas) -> tuple[float, list]:
    """dW extrapolated to 1/Omega -> 0, the way FINDINGS 15 extrapolated kappa."""
    got = []
    for om in omegas:
        try:
            got.append((om, barrier_at(gamma, om)))
        except (ValueError, RuntimeError):
            continue
    if len(got) < 2:
        return (got[0][1] if got else float("nan")), got
    x = np.array([1.0 / o for o, _ in got])
    y = np.array([v for _, v in got])
    return float(np.polyfit(x, y, 1)[1]), got


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.40, 0.42, 0.44, 0.455, 0.47, 0.48, 0.487])
    ap.add_argument("--n-omega", type=int, default=3)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/barrier_near_gamma_c.json"))
    args = ap.parse_args()

    t0 = time.time()
    print(f"{'gamma':>7} {'gc-g':>8} {'delta*':>8} {'Omegas used':>22} "
          f"{'dW(Om->inf)':>12} {'Om*dW':>7}")
    rows = []
    for g in args.gammas:
        ds = delta_star(g)
        # first guess from the predicted form, then open the window around it
        guess = 1.5 * (1 - 2 * g) / (1 + g) * ds ** 2
        lo, hi = omega_window(g, guess)
        if lo >= hi:
            print(f"{g:>7.4f}   window closed (lattice {lo} >= floor {hi})")
            continue
        oms = np.unique(np.geomspace(lo, min(hi, 6 * lo), args.n_omega).astype(int))
        dw, got = converged_barrier(g, oms)
        if not np.isfinite(dw) or dw <= 0 or not got:
            print(f"{g:>7.4f}   no usable measurement")
            continue
        print(f"{g:>7.4f} {GAMMA_C - g:>8.4f} {ds:>8.4f} "
              f"{str([o for o, _ in got]):>22} {dw:>12.3e} "
              f"{got[-1][0] * got[-1][1]:>7.2f}")
        rows.append({"gamma": g, "gap": GAMMA_C - g, "delta_star": ds,
                     "dW": dw, "points": got})

    if len(rows) >= 4:
        x = np.log(np.array([r["gap"] for r in rows]))
        y = np.log(np.array([r["dW"] for r in rows]))
        p, c = np.polyfit(x, y, 1)
        res = y - (p * x + c)
        print(f"\n  dW = {np.exp(c):.4f} * (gamma_c - gamma)^{p:.4f}"
              f"     R^2 = {1 - res.var()/y.var():.6f}")
        print(f"  residuals (ln): " + " ".join(f"{v:+.4f}" for v in res))
        print(f"\n  PREDICTED exponent 2 (kappa and delta*^2 each vanish linearly)")
        print(f"  T4's original guess was 1.  Measured: {p:.4f}")
        # local slopes -- an exponent that drifts is not an exponent
        print(f"\n  local slopes between consecutive gamma:")
        for i in range(len(rows) - 1):
            s = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
            print(f"    {rows[i]['gamma']:.3f} -> {rows[i+1]['gamma']:.3f}: {s:.4f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
