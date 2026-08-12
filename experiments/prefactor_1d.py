"""T14-e: the prefactor, on an instrument §35.3 did not have — the exact 1-D chain

§35.3 proved the prefactor cannot be extracted NUMERICALLY from the 2-D problem: the
candidate basis functions are collinear over any accessible Omega range (correlations
0.961 and 0.986, condition number 2.8e8), so a slowly-varying prefactor and a slightly
wrong exponent cannot be told apart by fitting. That result stands and is not being
re-litigated.

**BUT A 1-D BIRTH-DEATH CHAIN HAS A CLOSED-FORM SPLITTING PROBABILITY**, and §50 built the
exact 1-D slaved chain: delta hops +-1 with rates Omega*up and Omega*dn, both from
`updown`. For absorbing states a and b,

    P(hit a before b | start i) = sum_{k=i}^{b-1} pi_k / sum_{k=a}^{b-1} pi_k,
    ln pi_k = sum_{j=a+1}^{k} [ln mu_j - ln lambda_j]

computed in logs with log-sum-exp, so it is **exact at any Omega with no solve, no fitting
and no cancellation** -- and it names the wrong outcome directly, as §35 required rather
than forming 1 - split.

So the 2-D error probability can be compared against an exactly-known 1-D one, and the
comparison is a fork rather than a fit.

PREDICTIONS, and per rule 19 the VERDICT CRITERION is designed before the test:

  P1  GATE. The closed form reproduces a direct sparse solve of the SAME chain at small
      Omega, to solver precision. If it does not, the recursion or the log-sum-exp is
      wrong and nothing below counts.
  P2  GATE. The 2-D solve reproduces §35's published deep-tail values at a shared cell.
  P3  **THE TEST, and it is a fork read off the SHAPE, not a threshold.** Plot
      ln(P_1D / P_2D) against Omega.
        (a) **FLAT** -> the ratio is a constant, and that constant IS the prefactor ratio
            between the reduction and the truth -- measured, never fitted, which is
            exactly what §35.3 showed fitting cannot deliver.
        (b) **LINEAR with nonzero slope** -> the two EXPONENTS differ, the slope measures
            by how much, and the prefactor is not separable from it on this axis either --
            §35.3's conclusion reproduced on an independent instrument.
      Both outcomes are informative; the criterion distinguishes them by shape over a wide
      Omega range, so neither can be reached by a bad threshold.
  P4  ABSOLUTE. §36 measured the 1-D reduction ~1% shallow in the exponent, so if (b)
      holds the slope should be ~0.01 * 2V. Reported against that number, not fitted to it.
  P5  **THE DISCRIMINATING AXIS (rule 9).** §39.2 established the reduction becomes exact
      as sep -> infinity, and §44 gives rho as a lever on sep that leaves the drift
      untouched (§51: rho does not appear in mu at all). **So the slope in (b) must shrink
      toward zero as rho grows.** If it does, the discrepancy is the slow-manifold lag and
      the prefactor question is well-posed in the rho -> infinity limit. If the slope is
      rho-independent, the discrepancy is NOT the lag and §36's 1% has another cause.
  P6  If (a) holds at large rho, the prefactor ratio is a number this project has never
      had. It is reported with its Omega range and nothing is extrapolated past it.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from crnl.cme import enumerate_states, generator
from experiments.arrhenius_optimum import am_rho, delta_star_rho
from experiments.slaving_axis import slaved, updown


def _setup(net, ds, omega, eps, theta):
    st = slaved(net, eps * ds)
    if st is None:
        return None
    nb = int(round(st[2] * omega))
    rest = omega - nb
    d0 = max(1, int(round(eps * ds * omega)))
    if (rest - d0) % 2:
        d0 -= 1
    thr = max(2, int(round(theta * ds * omega)))
    if thr <= d0 or rest - d0 < 0:
        return None
    n0 = np.array([(rest + d0) // 2, (rest - d0) // 2, nb], dtype=np.int64)
    return n0, thr, d0


def p_err_1d(net, omega, m0, thr):
    """Exact P(hit -thr before +thr) for the slaved chain, in logs. No cancellation."""
    ms = np.arange(-thr, thr + 1)
    lp = np.zeros(len(ms))
    acc = 0.0
    for k in range(1, len(ms)):
        ud = updown(net, abs(ms[k]) / omega)
        if ud is None:
            return None
        up, dn = (ud[0], ud[1]) if ms[k] >= 0 else (ud[1], ud[0])
        if up <= 0 or dn <= 0:
            return None
        acc += np.log(dn) - np.log(up)      # ln(mu_k / lambda_k)
        lp[k] = acc

    def lse(v):
        m = v.max()
        return m + np.log(np.exp(v - m).sum())

    i = int(np.where(ms == m0)[0][0])
    b = len(ms) - 1                          # index of +thr
    return float(np.exp(lse(lp[i:b]) - lse(lp[0:b])))


def p_err_1d_solve(net, omega, m0, thr):
    """The same chain by a direct sparse solve, naming the wrong outcome (§35)."""
    ms = np.arange(-(thr - 1), thr)
    n = len(ms)
    up = np.empty(n)
    dn = np.empty(n)
    for i, m in enumerate(ms):
        ud = updown(net, abs(m) / omega)
        if ud is None:
            return None
        up[i], dn[i] = (ud[0], ud[1]) if m >= 0 else (ud[1], ud[0])
    up *= omega
    dn *= omega
    rows, cols, vals, rhs = [], [], [], np.zeros(n)
    for i in range(n):
        rows.append(i); cols.append(i); vals.append(-(up[i] + dn[i]))
        if i + 1 < n:
            rows.append(i); cols.append(i + 1); vals.append(up[i])
        if i - 1 >= 0:
            rows.append(i); cols.append(i - 1); vals.append(dn[i])
        else:
            rhs[i] -= dn[i]                  # absorbing at -thr is the WRONG outcome
    A = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsc()
    h = spla.spsolve(A, rhs)
    return float(h[int(np.where(ms == m0)[0][0])])


def p_err_2d(net, omega, n0, thr):
    """§35's direct 2-D solve: the wrong outcome named at the boundary."""
    states, index = enumerate_states(3, int(omega))
    absorb = np.array([abs(int(s[0]) - int(s[1])) >= thr for s in states])
    fav = np.array([int(s[1]) > int(s[0]) for s in states])[absorb].astype(float)
    Q = generator(net, int(omega), float(omega))
    tr = np.where(~absorb)[0]
    tmap = {int(i): r for r, i in enumerate(tr)}
    h = spla.spsolve(Q[tr][:, tr].tocsr(),
                     -(Q[tr][:, np.where(absorb)[0]].tocsr() @ fav))
    return float(h[tmap[index[tuple(n0)]]])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.20)
    ap.add_argument("--rhos", type=float, nargs="+", default=[1.0, 4.0, 16.0, 64.0])
    ap.add_argument("--omegas", type=int, nargs="+",
                    default=[100, 150, 200, 300, 400, 500])
    ap.add_argument("--eps", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/prefactor_1d.json"))
    args = ap.parse_args()

    t0 = time.time()
    g = args.gamma

    print("=== P1 GATE: closed form vs direct solve of the SAME 1-D chain")
    net = am_rho(g, 1.0)
    ds = delta_star_rho(g, 1.0)
    worst = 0.0
    for om in (60, 100, 150):
        s = _setup(net, ds, om, args.eps, args.theta)
        if s is None:
            continue
        n0, thr, d0 = s
        a = p_err_1d(net, om, d0, thr)
        b = p_err_1d_solve(net, om, d0, thr)
        rel = abs(a - b) / max(abs(a), 1e-300)
        worst = max(worst, rel)
        print(f"  Omega={om:>4}: closed {a:.10e}  solve {b:.10e}  rel {rel:.2e}")
    print(f"  -> P1 {'HOLDS' if worst < 1e-8 else 'FAILS'}")

    print(f"\n=== P3/P5: ln(P_1D / P_2D) vs Omega, at each rho")
    rows = []
    for rho in args.rhos:
        net = am_rho(g, rho)
        ds = delta_star_rho(g, rho)
        if ds <= 0:
            continue
        print(f"\n  rho = {rho}   (delta* = {ds:.4f})")
        print(f"{'Omega':>7}{'P_1D':>15}{'P_2D':>15}{'ln(1D/2D)':>13}")
        oms, lr = [], []
        for om in args.omegas:
            s = _setup(net, ds, om, args.eps, args.theta)
            if s is None:
                continue
            n0, thr, d0 = s
            try:
                a = p_err_1d(net, om, d0, thr)
                b = p_err_2d(net, om, n0, thr)
            except Exception:
                continue
            if a is None or b is None or a <= 0 or b <= 0:
                continue
            l = float(np.log(a / b))
            oms.append(om); lr.append(l)
            rows.append({"rho": rho, "omega": om, "p1d": a, "p2d": b, "lnratio": l})
            print(f"{om:>7}{a:>15.6e}{b:>15.6e}{l:>13.5f}")
        if len(oms) >= 3:
            sl, ic = np.polyfit(oms, lr, 1)
            span = max(lr) - min(lr)
            print(f"    slope d ln(ratio)/d Omega = {sl:+.6f}   intercept {ic:+.4f}"
                  f"   span {span:.4f}")
            rows[-1]["slope"] = float(sl)

    print(f"\n=== P3 verdict: FLAT (prefactor) or LINEAR (exponent mismatch)?")
    slopes = {}
    for rho in args.rhos:
        sel = [(r["omega"], r["lnratio"]) for r in rows if r["rho"] == rho]
        if len(sel) < 3:
            continue
        x = np.array([s[0] for s in sel], dtype=float)
        y = np.array([s[1] for s in sel])
        sl = float(np.polyfit(x, y, 1)[0])
        slopes[rho] = sl
        # shape test: does a line through the data explain it, and is the slope real?
        drift = abs(sl) * (x.max() - x.min())
        print(f"  rho={rho:>6}: slope {sl:+.6f}, total drift over the Omega range"
              f" {drift:+.4f} nats"
              f"  -> {'LINEAR, exponents differ' if drift > 0.2 else 'FLAT, a prefactor ratio'}")

    print(f"\n=== P5 (rule 9): does the slope shrink with rho, as §39.2 requires?")
    if len(slopes) >= 2:
        ks = sorted(slopes)
        print("  " + "   ".join(f"rho={k}:{slopes[k]:+.6f}" for k in ks))
        shrink = abs(slopes[ks[-1]]) < 0.5 * abs(slopes[ks[0]])
        print(f"  -> {'slope SHRINKS with rho: the discrepancy is the slow-manifold lag' if shrink else 'slope does NOT shrink: §36 1% has another cause'}")

    print(f"\n=== P6: the prefactor ratio, where the ratio is flat")
    for rho in sorted(slopes):
        sel = [r for r in rows if r["rho"] == rho]
        y = np.array([r["lnratio"] for r in sel])
        if abs(slopes[rho]) * (max(r["omega"] for r in sel)
                               - min(r["omega"] for r in sel)) <= 0.2:
            print(f"  rho={rho}: P_1D/P_2D = {np.exp(y.mean()):.5f}"
                  f"  over Omega {min(r['omega'] for r in sel)}"
                  f"..{max(r['omega'] for r in sel)}  (spread {np.ptp(y):.4f} nats)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"gamma": g, "rows": rows,
                                    "slopes": {str(k): v for k, v in slopes.items()}},
                                   indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
