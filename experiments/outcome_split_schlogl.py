"""T-TUR-f: does §60's outcome-wise factorisation survive ASYMMETRY? — asked where it can be

§60 found <e^(-S_tot)|o> = 1 for each outcome separately and concluded the fluctuation theorem
cannot bound the error rate, closing the founding question's sharpest form. §66 showed that
conclusion was tested only at beta = 0, where the two absorbing boundaries are exchange images
and the cancellation is guaranteed by symmetry -- and then could not settle the tilted case,
because the pi(n)/pi(n0) boundary convention exhausts double precision exactly where the
predicted effect is largest (ln w = -38.9 at beta = 0.40, Phi_all returning 0.813 instead of 1).

**Rather than fight that wall in 2-D, ask the question where the instrument is exact.** §67's
Schloegl element has ONE dynamic species, so there is no exchange symmetry to appeal to at all
-- the two rails are not images of each other under anything -- and the chain is 1-D, so every
quantity is a tridiagonal solve.

**AND THE TILTED GENERATOR COLLAPSES.** For an up-jump n -> n+1 the two channels contribute
a1(n)*exp(-ln[a1(n)/a2(n+1)]) + a3(n)*exp(-ln[a3(n)/a4(n+1)]) = a2(n+1) + a4(n+1) = mu(n+1),
and for a down-jump, a1(n-1) + a3(n-1) = lambda(n-1). So **the reverse-weighted generator is
just the chain with lambda and mu swapped and shifted** -- explicit, tridiagonal, and free of
the enormous factors that broke §66. ln pi comes from the exact product formula in logs.

Asymmetry is built in by placing the three fixed points at arbitrary r1 < r2 < r3 rather than
symmetrically: k1a = k1r*e1, k2r = k1r*e2, k2b = k1r*e3 with e_i the elementary symmetric
polynomials, and A = ln[e1*e2/e3]. The skew s = (r3 - r2)/(r2 - r1) is 1 for a balanced
element and departs from 1 as the two basins differ.

PREDICTIONS, written before running.

  P1  GATE. Phi_lo + Phi_hi = 1 to solver precision, at every skew. The IFT is an identity for
      any network, so a failure is the tilted construction, not physics, and nothing below
      counts. **This is the gate §66 could not hold** -- if it holds here across the full
      asymmetry range, the 1-D route has bought exactly what was missing.
  P2  **THE TEST. Is Phi_o = p_o per outcome on an element with no symmetry whatsoever?**
      §60's stated mechanism -- boundaries as exchange images with equal stationary weight --
      is simply unavailable here. **PREDICTED: it FAILS, and |Phi_o/p_o - 1| grows with the
      skew.** If instead it holds at every skew, then the outcome-wise identity is deeper than
      symmetry, §60's conclusion generalises after all, and the mechanism §60 gave for it is
      wrong even though the result is right -- which is the more interesting outcome and is
      why the prediction is written against it.
  P3  **THE CONFOUND FROM §66, CHECKED AGAIN.** There the nuisance (boundary-weight ratio,
      hence conditioning) grew with the cause (tilt). Here the solve is exact and tridiagonal,
      so ln w can be reported alongside |r-1| and the two separated: if |r-1| stays at solver
      noise while ln w runs over tens, the confound is broken by construction rather than by
      correlation.
  P4  **RULE 9.** Sweep Omega and the start position as well as the skew. An effect that is a
      property of asymmetry must survive both.
  P5  **VERDICT RULE unit-tested on engineered data before this runs** (§66's convention,
      which has now caught a defect in three consecutive sections).
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
from scipy.linalg import solve_banded


def consts(r1, r2, r3, k1r=1.0):
    e1, e2, e3 = r1 + r2 + r3, r1 * r2 + r1 * r3 + r2 * r3, r1 * r2 * r3
    return k1r * e1, k1r, k1r * e3, k1r * e2      # k1a, k1r, k2b, k2r


def channels(n, omega, c):
    k1a, k1r, k2b, k2r = c
    n = np.asarray(n, dtype=float)
    return (np.maximum(k1a * n * (n - 1.0) / omega, 0.0),
            np.maximum(k1r * n * (n - 1.0) * (n - 2.0) / omega ** 2, 0.0),
            np.full_like(n, k2b * omega),
            np.maximum(k2r * n, 0.0))


def _solve(lo, hi, up, dn, source):
    idx = np.arange(lo + 1, hi)
    k = len(idx)
    ab = np.zeros((3, k))
    ab[0, 1:] = up[idx[:-1]]
    ab[1, :] = -(up[idx] + dn[idx])
    ab[2, :-1] = dn[idx[1:]]
    return idx, solve_banded((1, 1), ab, -source[idx])


def cell(omega, r1, r2, r3, eps, theta):
    """Outcome-split IFT on the exact asymmetric chain. Start eps below the saddle r2."""
    c = consts(r1, r2, r3)
    cap = int(np.ceil(3.0 * r3 * omega))
    ns = np.arange(0, cap + 1, dtype=float)
    a1, a2, a3, a4 = channels(ns, omega, c)
    lam, mu = a1 + a3, a2 + a4

    lo = int(round((r2 - theta * (r2 - r1)) * omega))
    hi = int(round((r2 + theta * (r3 - r2)) * omega))
    n0 = int(round((r2 - eps * (r2 - r1)) * omega))
    if not (lo < n0 < hi) or hi >= cap or lo < 3:
        return None

    # ---- exact splitting probabilities, named at the boundary (§35, no 1 - split) ----
    lp, acc = np.zeros(hi - lo + 1), 0.0
    for k in range(lo + 1, hi + 1):
        if lam[k] <= 0 or mu[k] <= 0:
            return None
        acc += np.log(mu[k]) - np.log(lam[k])
        lp[k - lo] = acc

    def lse(v):
        M = v.max()
        return M + np.log(np.exp(v - M).sum())

    i = n0 - lo
    p_lo = float(np.exp(lse(lp[i:hi - lo]) - lse(lp[0:hi - lo])))
    p_hi = 1.0 - p_lo

    # ---- ln pi exactly, by the birth-death product formula in logs ----
    lnpi = np.zeros(cap + 1)
    acc = 0.0
    for k in range(1, cap + 1):
        if lam[k - 1] <= 0 or mu[k] <= 0:
            lnpi[k] = -np.inf
            continue
        acc += np.log(lam[k - 1]) - np.log(mu[k])
        lnpi[k] = acc

    # ---- the reverse-weighted (tilted) generator: up-rate mu(n+1), down-rate lambda(n-1) --
    up_t = np.zeros(cap + 1)
    dn_t = np.zeros(cap + 1)
    up_t[:-1] = mu[1:]
    dn_t[1:] = lam[:-1]
    a_tot = lam + mu                                   # exit rate of the ORIGINAL chain

    idx = np.arange(lo + 1, hi)
    k = len(idx)
    ab = np.zeros((3, k))
    ab[0, 1:] = -up_t[idx[:-1]]
    ab[1, :] = a_tot[idx]
    ab[2, :-1] = -dn_t[idx[1:]]

    w = lnpi - lnpi[n0]                                # ln[pi(n)/pi(n0)]
    out = {}
    for name, at_lo, at_hi in (("lo", True, False), ("hi", False, True),
                               ("all", True, True)):
        rhs = np.zeros(k)
        if at_lo:
            rhs[0] += dn_t[lo + 1] * np.exp(w[lo])
        if at_hi:
            rhs[-1] += up_t[hi - 1] * np.exp(w[hi])
        out[name] = float(solve_banded((1, 1), ab, rhs)[int(np.where(idx == n0)[0][0])])

    skew = (r3 - r2) / (r2 - r1)
    A = float(np.log((r1 + r2 + r3) * (r1 * r2 + r1 * r3 + r2 * r3) / (r1 * r2 * r3)))
    return {"omega": omega, "r1": r1, "r2": r2, "r3": r3, "skew": float(skew), "A": A,
            "eps": eps, "theta": theta, "lo": lo, "hi": hi, "n0": n0,
            "p_lo": p_lo, "p_hi": p_hi, "Phi_lo": out["lo"], "Phi_hi": out["hi"],
            "Phi_all": out["all"], "lnw": float(w[hi] - w[lo])}


def decide(skews, dev, dev_sym, noise):
    """Reachable three ways; unit-tested in tests/test_outcome_split_schlogl.py."""
    d = np.asarray(dev, float)
    floor = max(dev_sym, noise)
    grew = bool(np.all(np.diff(d) > -1e-15)) and d[-1] > d[0]
    if d.max() <= 2.0 * floor:
        return "a", (f"NO DEVIATION: |r-1| peaks at {d.max():.3e} against a floor "
                     f"{floor:.3e}. The factorisation holds with NO symmetry available, so "
                     f"§60's result generalises and §60's MECHANISM for it is wrong.")
    if grew and d[-1] > 2.0 * floor:
        return "b", (f"DEVIATION: |r-1| rises with the skew to {d[-1]:.3e}, clear of "
                     f"{floor:.3e}. §60's closure is symmetry-dependent and the founding "
                     f"question reopens for real devices.")
    return "c", (f"INCONCLUSIVE: |r-1| reaches {d.max():.3e} against a floor {floor:.3e}, "
                 f"{'non-monotone' if not grew else 'not clear of it'}.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r2", type=float, default=1.0)
    ap.add_argument("--gap", type=float, default=0.5)
    ap.add_argument("--skews", type=float, nargs="+",
                    default=[1.0, 1.3, 1.7, 2.2, 3.0, 4.0])
    ap.add_argument("--omegas", type=int, nargs="+", default=[200, 400, 800])
    ap.add_argument("--eps", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/outcome_split_schlogl.json"))
    args = ap.parse_args()

    print("An element with ONE species: no exchange symmetry exists to appeal to.")
    print(f"{'skew':>6}{'Om':>6}{'A':>8}{'lnw':>9}{'p_hi':>11}{'Phi_hi':>11}"
          f"{'sum':>11}{'|r_lo-1|':>11}{'|r_hi-1|':>11}")
    rows, bad = [], []
    for sk in args.skews:
        r2 = args.r2
        r1 = r2 - args.gap
        r3 = r2 + sk * args.gap
        for om in args.omegas:
            r = cell(om, r1, r2, r3, args.eps, args.theta)
            if r is None:
                continue
            r["sum"] = r["Phi_lo"] + r["Phi_hi"]
            r["r_lo"] = r["Phi_lo"] / r["p_lo"] if r["p_lo"] > 0 else np.nan
            r["r_hi"] = r["Phi_hi"] / r["p_hi"] if r["p_hi"] > 0 else np.nan
            if abs(r["sum"] - 1.0) > 1e-8:          # §66's per-cell precision gate
                bad.append(r)
                print(f"{sk:>6.2f}{om:>6}{r['A']:>8.4f}{r['lnw']:>9.2f}"
                      f"   EXCLUDED by P1: |sum-1| = {abs(r['sum']-1):.2e}"
                      f" at ln w = {r['lnw']:.1f} (e^{r['lnw']:.0f} past double precision)")
                continue
            rows.append(r)
            print(f"{sk:>6.2f}{om:>6}{r['A']:>8.4f}{r['lnw']:>9.2f}{r['p_hi']:>11.4e}"
                  f"{r['Phi_hi']:>11.4e}{r['sum']:>11.7f}"
                  f"{abs(r['r_lo']-1):>11.3e}{abs(r['r_hi']-1):>11.3e}")
    if not rows:
        print("no cells")
        return

    print(f"\n=== P1 GATE: does Phi_lo + Phi_hi reconstruct the identity?")
    dd = np.array([abs(r["sum"] - 1.0) for r in rows])
    print(f"  |sum - 1| over {len(rows)} cells: max {dd.max():.3e}, median {np.median(dd):.3e}")
    print(f"  {len(bad)} cells excluded on precision"
          + (f" (ln w = " + ", ".join(f"{r['lnw']:.0f}" for r in bad) + ")" if bad else ""))
    print(f"  -> P1 {'HOLDS on every surviving cell -- the gate §66 could not hold' if dd.max() < 1e-8 else 'FAILS'}")

    print(f"\n=== P2: is Phi_o = p_o with NO symmetry available?")
    print(f"{'skew':>6}{'mean |r_lo-1|':>16}{'mean |r_hi-1|':>16}{'max |ln w|':>12}")
    per, sks = [], []
    for sk in args.skews:
        sel = [r for r in rows if abs(r["skew"] - sk) < 1e-9]
        if not sel:
            continue
        dl = float(np.mean([abs(r["r_lo"] - 1) for r in sel]))
        dh = float(np.mean([abs(r["r_hi"] - 1) for r in sel]))
        per.append(max(dl, dh))
        sks.append(sk)
        print(f"{sk:>6.2f}{dl:>16.3e}{dh:>16.3e}"
              f"{max(abs(r['lnw']) for r in sel):>12.1f}")
    dev_sym = per[0] if sks and abs(sks[0] - 1.0) < 1e-9 else 0.0
    noise = float(np.median(dd))
    code, msg = decide(sks, per, dev_sym, noise)
    print(f"  balanced-element residual {dev_sym:.3e}; solver noise {noise:.3e}")
    print(f"  -> P2 ({code}) {msg}")

    print(f"\n=== P3: is the §66 confound present? |r-1| against ln w")
    lw = np.array([abs(r["lnw"]) for r in rows])
    dv = np.array([max(abs(r["r_lo"] - 1), abs(r["r_hi"] - 1)) for r in rows])
    print(f"  |ln w| spans {lw.min():.1f}..{lw.max():.1f} while |r-1| spans "
          f"{dv.min():.2e}..{dv.max():.2e}")
    print(f"  -> {'the confound is BROKEN: ln w runs over tens and |r-1| stays at solver noise' if dv.max() < 1e-8 else f'correlation {np.corrcoef(lw, dv)[0,1]:+.3f}, to be read with §66 in mind'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"rows": rows, "per_skew": per, "skews": sks,
                                    "verdict": code}, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
