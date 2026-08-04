"""T14 UNLOCKED: the collapse to 1e-30, and the asymptotic rate against §15's closed form.

THEORIES T14 has said since §21 that the sampling floor is the binding constraint:
"the founding claim concerns a switch that errs at 1e-15 while every measured number
sits between 1e-1 and 1e-2 ... Exact CME reaches small probabilities but its state space
grows as ~Omega^2/2; sampling handles any Omega but floors at 1/N. Large Omega AND small
probability is reachable by neither."

**The probability half of that was an implementation artifact, not a limit.** `p_cme`
returns `1 - split` -- the error probability as a difference of two numbers near 1 -- so
it dies to catastrophic cancellation around 1e-12, and §28 lost its gamma = 0.15 cells
to exactly that. But `splitting_probability` takes the favoured-set predicate, so naming
the WRONG outcome as favoured solves for the small number DIRECTLY, with no subtraction
anywhere. Measured: identical to the subtractive route to 7-8 digits in the overlap, and
it keeps going -- **P = 6.35e-33 at Omega = 2000 in 115 s**, 25 orders below anything
this project has measured.

**IT IS VALIDATED COMPONENTWISE, not by a norm.** A norm residual is dominated by the
large components and would not notice a garbage small one. Each row of the transient
generator has at most ~7 nonzeros, so the true residual is computable by exact summation
per row; one refinement step then gives the componentwise relative correction. At
Omega = 2000, h(start) = 6.35e-33 with |delta/h| = 1.0e-13 at the start state and
1.2e-13 as the MAXIMUM over every positive component. That is the M-matrix property --
the LU solve carries no subtractive cancellation, so relative accuracy survives to
arbitrarily small values. The check is re-run here rather than cited.

WHAT THE UNLOCK IS FOR. Every collapse slope this project has published -- §12, §15,
§27's 6.53 decades, §28's whole absolute-test arc -- is a **finite-Omega effective
slope** fitted over a window of a few decades. The local slope is not constant: over
Omega = 400..2000 at gamma = 0.25 it drifts from -0.036816 to -0.036096, monotonically.
So `P ~ A(Omega) exp(-c*Omega)` with an algebraic prefactor, and a straight-line fit
returns `c` contaminated by `A`. With 30+ decades the two separate for the first time,
and §15's closed form can be tested against the ASYMPTOTIC rate rather than against an
effective slope that happens to be measured where the contamination is largest.

PREDICTIONS, written before running:

  P1  The componentwise refinement stays below 1e-10 at every cell, including the
      deepest. If it does not, nothing else in the section is admissible.
  P2  The local slope drifts monotonically with Omega at every gamma, so no single
      number is "the" collapse slope. §27's "collapse holds over 6.53 decades" was
      measured over a window too short to see it, and remains true as stated -- a
      straight line does fit a short window -- while being the wrong asymptotic object.
  P3  THE TEST, AND THE UNCOMFORTABLE OUTCOME IS THE LIKELY ONE. Fitting with a
      prefactor term and extrapolating, the asymptotic `c` should differ from §28.3's
      effective slope, and because the drift moves the slope AWAY from the prediction
      (|slope| falls, prediction is larger in magnitude), the disagreement with
      `-2*V_exact` should be **LARGER** than §28.3's 6.3% at gamma = 0.25, not smaller.
      **That would qualify §28's headline**: its agreement was flattered by finite-Omega
      contamination. If instead `c` lands closer to the closed form, the residual §28
      attributed to the 1-D reduction was largely finite-Omega and the closed form is
      better than §28 thought -- also a real result, in the opposite direction.
  P4  RULE 15. Three ansaetze are fitted and ALL are reported, never only the flattering
      one: pure exponential, plus algebraic prefactor `a*ln(Omega)`, plus `d/Omega`
      correction. If they disagree on `c` by more than a few percent the asymptotic rate
      is UNRESOLVED and the section says so rather than picking one.
  P5  eps-controlled fits throughout (§27: the integer lattice moves realised eps ~10%
      and that alone bounced raw local slopes 40%). The Omega = 1400 point in the first
      probe sat visibly off the drift curve for exactly this reason.
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import time

import numpy as np
import scipy.sparse.linalg as spla

from crnl.cme import enumerate_states, generator
from crnl.networks.am_reversible import am_reversible, delta_star
from experiments.approximation_hierarchy import _setup
from experiments.collapse_slope_absolute import V_exact


def solve_direct(gamma: float, omega: int, eps: float, theta: float,
                 refine: bool = False) -> dict:
    """P(error) solved for DIRECTLY -- the wrong outcome is the favoured set."""
    net = am_reversible(gamma)
    n0, thr, _ = _setup(gamma, omega, eps, theta)
    states, index = enumerate_states(net.n_species, int(omega))
    absorb = np.array([abs(int(s[0]) - int(s[1])) >= thr for s in states])
    fav = np.array([int(s[1]) > int(s[0]) for s in states])[absorb].astype(float)
    Q = generator(net, int(omega), float(omega))
    tr = np.where(~absorb)[0]
    tmap = {int(i): r for r, i in enumerate(tr)}
    A = Q[tr][:, tr].tocsr()
    b = -(Q[tr][:, np.where(absorb)[0]].tocsr() @ fav)
    h = spla.spsolve(A, b)
    si = tmap[index[tuple(np.asarray(n0, dtype=np.int64))]]
    out = {"omega": omega, "p": float(h[si]),
           "eps_realised": float(int(n0[0] - n0[1]) / (delta_star(gamma) * omega)),
           "states": int(A.shape[0])}
    if refine:
        r = np.empty_like(h)
        ip, ix, da = A.indptr, A.indices, A.data
        for i in range(A.shape[0]):
            r[i] = math.fsum([-da[k] * h[ix[k]] for k in range(ip[i], ip[i + 1])]
                             + [b[i]])
        d = spla.spsolve(A, r)
        pos = h > 0
        out["refine_start"] = float(abs(d[si] / h[si])) if h[si] else float("nan")
        out["refine_max"] = float(np.abs(d[pos] / h[pos]).max())
    return out


def fit_all(om: np.ndarray, lp: np.ndarray, er: np.ndarray) -> dict:
    """Three ansaetze, all eps-controlled, all reported (rule 15)."""
    ec = er - er.mean()
    out = {}
    designs = {
        "pure":      [om, om * ec, np.ones_like(om)],
        "prefactor": [om, np.log(om), om * ec, np.ones_like(om)],
        "inverse":   [om, 1.0 / om, om * ec, np.ones_like(om)],
    }
    for name, cols in designs.items():
        A = np.vstack(cols).T
        c, *_ = np.linalg.lstsq(A, lp, rcond=None)
        res = lp - A @ c
        out[name] = {"c": float(c[0]),
                     "extra": float(c[1]) if name != "pure" else None,
                     "r2": float(1 - res.var() / lp.var()),
                     "rms": float(np.sqrt((res ** 2).mean()))}
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.20, 0.25, 0.30, 0.35])
    ap.add_argument("--omegas", type=int, nargs="+",
                    default=[150, 300, 450, 600, 750, 900, 1050, 1200, 1350, 1500])
    ap.add_argument("--eps-frac", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--refine-at", type=int, nargs="+", default=[600, 1500])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/deep_tail.json"))
    args = ap.parse_args()

    t0 = time.time()
    print("P(error) solved DIRECTLY (wrong outcome as favoured set) -- no subtraction")
    allr = {}
    for g in args.gammas:
        print(f"\n=== gamma = {g}   eps = {args.eps_frac}")
        print(f"{'Omega':>6}{'states':>9}{'P(error)':>15}{'log10 P':>10}"
              f"{'local slope':>13}{'eps real':>10}{'refine':>10}{'sec':>7}")
        rows, prev = [], None
        for om in args.omegas:
            t1 = time.time()
            r = solve_direct(g, om, args.eps_frac, args.theta,
                             refine=om in args.refine_at)
            dt = time.time() - t1
            p = r["p"]
            sl = ((np.log(p) - np.log(prev[1])) / (om - prev[0])
                  if prev and p > 0 and prev[1] > 0 else float("nan"))
            rows.append(r)
            rf = f"{r['refine_max']:.1e}" if "refine_max" in r else ""
            print(f"{om:>6}{r['states']:>9}{p:>15.5e}"
                  f"{np.log10(p) if p > 0 else float('nan'):>10.2f}{sl:>13.6f}"
                  f"{r['eps_realised']:>10.4f}{rf:>10}{dt:>7.1f}")
            prev = (om, p)

        good = [r for r in rows if r["p"] > 0]
        om = np.array([r["omega"] for r in good], float)
        lp = np.log([r["p"] for r in good])
        er = np.array([r["eps_realised"] for r in good])
        fits = fit_all(om, lp, er)
        net = am_reversible(g)
        pred = -2.0 * V_exact(net, args.eps_frac * delta_star(g))
        decades = (max(lp) - min(lp)) / np.log(10)
        print(f"  spans {decades:.2f} decades   "
              f"(§27's best was 6.53, §28.3's grids 4)")
        print(f"  {'ansatz':>10}{'c':>12}{'extra':>11}{'R^2':>10}{'rms':>9}"
              f"{'pred/c':>9}")
        for k, v in fits.items():
            ex = f"{v['extra']:.4f}" if v["extra"] is not None else "-"
            print(f"  {k:>10}{v['c']:>12.6f}{ex:>11}{v['r2']:>10.6f}"
                  f"{v['rms']:>9.4f}{pred/v['c']:>9.4f}")
        cs = np.array([v["c"] for v in fits.values()])
        spread = (cs.max() - cs.min()) / abs(cs.mean())
        print(f"  closed form -2*V_exact = {pred:.6f};  spread across ansaetze "
              f"{100*spread:.2f}%"
              f"   -> asymptotic rate {'RESOLVED' if spread < 0.03 else 'UNRESOLVED'}")
        allr[str(g)] = {"rows": rows, "fits": fits, "predicted": float(pred),
                        "decades": float(decades), "spread": float(spread)}

    print(f"\n=== P3: the closed form against the ASYMPTOTIC rate, per gamma")
    print(f"{'gamma':>7}{'pure c':>11}{'prefac c':>11}{'closed form':>13}"
          f"{'pred/pure':>11}{'pred/prefac':>13}{'§28.3 ratio':>13}")
    prior = {0.20: 1.110, 0.25: 1.063, 0.30: 1.034, 0.35: 1.004}
    for g, v in allr.items():
        pu, pf = v["fits"]["pure"]["c"], v["fits"]["prefactor"]["c"]
        p28 = prior.get(float(g))
        print(f"{float(g):>7.2f}{pu:>11.6f}{pf:>11.6f}{v['predicted']:>13.6f}"
              f"{v['predicted']/pu:>11.4f}{v['predicted']/pf:>13.4f}"
              f"{(f'{p28:.3f}' if p28 else '-'):>13}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(allr, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
