"""T14-d, second attempt: extract the prefactor from the LOCAL SLOPE, not from ln P.

§35.2 failed to resolve `b` and reported the null. Two flaws in that design, both mine:

  **(1) EVERY WINDOW SHARED THE SAME LOWER EDGE.** Sweep B varied only the upper edge
  of [150, X], so every fit contained the most-contaminated low-Omega data and none of
  them was a convergence test. A quantity that drifts as the upper edge moves out is
  not thereby converging -- it is being dragged by a fixed anchor. Windows must slide,
  not stretch.

  **(2) FITTING ln P IS BADLY CONDITIONED FOR THIS QUESTION.** Over any window,
  `ln Omega` is nearly linear in `Omega`, so the c-term and the b-term are strongly
  collinear; adding `1/Omega` made it worse, which is exactly why §35.2's four-term fit
  scattered b from +0.06 to -0.74 while the three-term fit spanned only -0.37..-0.45.
  That was read as "the missing-term reading is unsafe". It is better read as: the
  basis was wrong.

**THE RIGHT INSTRUMENT. If `ln P = -c*Omega + b*ln(Omega) + a`, then**

    s(Omega) = d(ln P)/d(Omega) = -c + b/Omega

**so the LOCAL SLOPE is LINEAR IN 1/Omega, with intercept -c and slope b.** Two
parameters, orthogonal in a way `ln Omega` and `Omega` never are, the constant `a`
differentiated away entirely -- and, decisively, **the linearity is itself a test of
the functional form**. A straight line means the log ansatz is right and b is read off
directly. Curvature means there is more structure, and its SIGN says what kind.

THE eps LATTICE HAS TO BE HANDLED FIRST, because local slopes amplify it. Realised eps
wobbles ~10% cell to cell and §27 measured that alone bouncing raw local slopes by 40%.
So every ln P is corrected to a common eps using the sensitivity from a global
eps-controlled fit BEFORE any local slope is taken -- otherwise this instrument is
strictly worse than the one it replaces.

PREDICTIONS, written before running:

  P1  The eps-corrected local slope is LINEAR in 1/Omega with high R^2. That would
      confirm the three-term form outright and make b a directly-read quantity rather
      than a fitted-and-confounded one.
  P2  b read this way is consistent across gamma -- i.e. §35.2's apparent
      gamma-dependence was the shared-lower-edge artifact. **This is what I expect**,
      and it is the specific claim the user's hunch is aimed at.
  P3  The intercept -c must reproduce §35's c to the ~0.2% that §35.2 already measured
      across windows. If it does not, this instrument disagrees with the established
      one and neither is usable until that is settled.
  P4  IF s vs 1/Omega is CURVED, the log ansatz is incomplete and the sign of the
      curvature discriminates: upward curvature means a positive 1/Omega^2 term (a
      further algebraic correction), while a systematic sag at LARGE 1/Omega (small
      Omega) means a second exponential contribution dying out. **A second exponential
      is physically expected** -- trajectories that cross without relaxing versus those
      that relax to the attractor first have different rates -- but the attractor route
      is ~8x steeper here, so it should die long before the window opens, and seeing it
      would be a surprise worth chasing.
  P5  Sliding windows of FIXED length at increasing position give a b that converges as
      the window moves out. §35.2's stretching windows could not show convergence even
      in principle, and reporting both makes the difference visible rather than argued.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from experiments.deep_tail import solve_direct


def curve(gamma, eps, theta, omegas):
    rows = []
    for om in omegas:
        r = solve_direct(gamma, int(om), eps, theta)
        if np.isfinite(r["p"]) and 0.0 < r["p"] < 1.0:
            rows.append(r)
    return rows


def eps_corrected(rows):
    """Remove the lattice wobble before differentiating (§27)."""
    om = np.array([r["omega"] for r in rows], float)
    lp = np.log([r["p"] for r in rows])
    er = np.array([r["eps_realised"] for r in rows])
    ec = er - er.mean()
    A = np.vstack([om, np.log(om), om * ec, np.ones_like(om)]).T
    c, *_ = np.linalg.lstsq(A, lp, rcond=None)
    return om, lp - c[2] * (om * ec), float(c[2])


def local_slopes(om, lp):
    """Centred differences; endpoints one-sided."""
    s = np.gradient(lp, om)
    return s


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.25, 0.30, 0.35])
    ap.add_argument("--omegas", type=int, nargs="+",
                    default=list(range(150, 2000, 100)))
    ap.add_argument("--eps-frac", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--win-len", type=int, default=8,
                    help="points per sliding window (fixed length, P5)")
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/local_slope_law.json"))
    args = ap.parse_args()

    t0 = time.time()
    print("s(Omega) = d(lnP)/dOmega = -c + b/Omega  ->  b is the SLOPE of s against 1/Omega")
    out = {}
    for g in args.gammas:
        rows = curve(g, args.eps_frac, args.theta, args.omegas)
        om, lp, k_eps = eps_corrected(rows)
        s = local_slopes(om, lp)
        x = 1.0 / om
        print(f"\n=== gamma = {g}   {len(om)} cells, Omega {int(om.min())}..{int(om.max())}"
              f", eps-sensitivity {k_eps:.4f}   ({time.time()-t0:.0f}s)")
        print(f"  {'Omega':>7}{'1/Omega':>10}{'ln P':>12}{'local slope':>14}")
        for o, xx, l, ss in zip(om, x, lp, s):
            print(f"  {int(o):>7}{xx:>10.5f}{l:>12.3f}{ss:>14.6f}")

        # P1/P3: the whole-range line
        p1 = np.polyfit(x, s, 1)
        r1 = s - np.polyval(p1, x)
        R1 = 1 - r1.var() / s.var()
        # P4: is there curvature?
        p2 = np.polyfit(x, s, 2)
        r2 = s - np.polyval(p2, x)
        R2 = 1 - r2.var() / s.var()
        print(f"  linear   s = {p1[1]:.6f} + {p1[0]:.4f}/Omega     R^2 = {R1:.6f}"
              f"   rms {np.sqrt((r1**2).mean()):.3e}")
        print(f"  quadratic  extra 1/Omega^2 coeff {p2[0]:>12.2f}   R^2 = {R2:.6f}"
              f"   rms {np.sqrt((r2**2).mean()):.3e}")
        print(f"  -> b = {p1[0]:.4f}   c = {-p1[1]:.6f}"
              f"   curvature {'MATTERS' if R2 - R1 > 0.02 else 'negligible'}")

        # P5: sliding fixed-length windows, position varying
        L = args.win_len
        print(f"  sliding windows of {L} points (fixed length, position varies):")
        print(f"  {'Omega window':>16}{'b':>10}{'c':>12}{'R^2':>10}")
        slid = []
        for i in range(0, len(om) - L + 1):
            xa, sa = x[i:i + L], s[i:i + L]
            pp = np.polyfit(xa, sa, 1)
            rr = sa - np.polyval(pp, xa)
            win = f"{int(om[i])}-{int(om[i+L-1])}"
            print(f"  {win:>16}{pp[0]:>10.4f}{-pp[1]:>12.6f}"
                  f"{1 - rr.var()/sa.var():>10.6f}")
            slid.append({"lo": int(om[i]), "hi": int(om[i + L - 1]),
                         "b": float(pp[0]), "c": float(-pp[1])})
        bs = np.array([w["b"] for w in slid])
        print(f"    b across sliding windows: {bs.min():.4f}..{bs.max():.4f}"
              f"   spread {100*(bs.max()-bs.min())/abs(bs.mean()):.2f}%")

        out[str(g)] = {"omega": om.tolist(), "lnP": lp.tolist(),
                       "slope": s.tolist(), "b_line": float(p1[0]),
                       "c_line": float(-p1[1]), "r2_line": float(R1),
                       "r2_quad": float(R2), "quad_coeff": float(p2[0]),
                       "sliding": slid, "eps_sensitivity": k_eps}

    print(f"\n=== P2: is b consistent across gamma once read this way?")
    print(f"{'gamma':>7}{'b (whole range)':>18}{'c':>12}{'R^2 linear':>12}")
    for g, v in out.items():
        print(f"{float(g):>7.2f}{v['b_line']:>18.4f}{v['c_line']:>12.6f}"
              f"{v['r2_line']:>12.6f}")
    bb = np.array([v["b_line"] for v in out.values()])
    print(f"  spread across gamma: {100*(bb.max()-bb.min())/abs(bb.mean()):.2f}%"
          f"   (§35.2's ln-P fits gave 11.6% unmatched, 8.0% matched)")
    print(f"  mean b = {bb.mean():.4f}   against WKB's -0.5 -> "
          f"{abs(bb.mean()/-0.5):.4f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
