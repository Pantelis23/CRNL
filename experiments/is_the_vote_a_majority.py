"""§111 -- does §34's crossover law transfer off AM, and is its k-independence a symmetry accident?

§34 derived the concatenation crossover in closed form and tested it absolutely across three
gamma on AM:

    Omega_x = [ (m-1)(a - ln tau) - ln C(k,m) ] / [ c (k - m) ],   m = ceil((k+1)/2)

and reported a structurally surprising confirmation: **the leading term contains no k**,
because m-1 = k-m = (k-1)/2 for EVERY odd k. Margins differ by 2x between k = 3 and 5; the
crossings nearly coincide.

**That identity is a statement about MAJORITY.** m = ceil((k+1)/2) is the number of failures
that flips a majority vote, and the pool merge is a majority vote only if the two rails sit
symmetrically about the saddle -- which is an AM fact, not a restoration fact. AM is exchange
symmetric by construction. **Schlogl is not.** Its rails are r1 = 0.15 and r3 = 3.1827 about a
saddle at r2 = 1.0, so the count-weighted mean of a merged pool is pulled far to the high side.

Pooling k tanks of volume Omega gives volume k*Omega holding the summed counts, so the merged
concentration is the mean of the k committed concentrations:

    x_merged(j) = [ (k-j) r3 + j r1 ] / k        for j tanks committed low

which falls below the saddle only when j > k (r3 - r2)/(r3 - r1) = 0.71972 k. That is not a
majority. Nothing in §34's derivation breaks -- m is simply a different integer -- but the
identity that made the answer k-independent does.

PREDICTIONS, WRITTEN BEFORE RUNNING.

  P1  WIRING. The Schlogl hold lifetime must be exponential in volume: ln T(N) = c*N + a with
      a straight-line fit good to a few percent over the solvable range. If it is not
      exponential, §34's derivation does not apply here for a reason that has nothing to do
      with symmetry and the rest of this section is void.

  P2  THE MERGE IS NOT A MAJORITY, and this is measured, not asserted. Computing the merged
      tank's own commitment probability p_merge(j) -- an actual Schlogl tank of volume k*Omega
      started at x_merged(j), no free comparison anywhere -- it must be ~0 for j below the
      threshold and ~1 above, with the step at

          m(k) = floor(0.71972 k) + 1   =   3, 4, 6   for k = 3, 5, 7

      against AM's majority m = 2, 3, 4. **The k = 3 case is the sharp one: m = k, so a merged
      Schlogl trio fails only if ALL THREE tanks fail.** If instead the step lands at the
      majority value the asymmetry does not survive the merge and P3 is moot.

  P3  THE STRUCTURAL CLAIM. §34's k-independence is an accident of symmetry. The ratio
      (m-1)/(k-m) is 1 for every odd k on AM; on Schlogl it is UNDEFINED at k = 3 (m = k, so
      the denominator vanishes and re-merging has no crossover of that form at all), 3 at
      k = 5, and 5 at k = 7. So predicted crossovers must be strongly k-DEPENDENT here --
      the opposite of §34's headline -- and the k = 3 protocol should have no crossover.

  P4  THE ABSOLUTE TEST where the formula is defined (k = 5, 7): predicted Omega_x against the
      continuously-interpolated measured crossover, with c and a taken from the HOLD protocol
      alone and no crossover measurement entering the prediction (§34's standard, §30.1's
      reminder). Reported with its residual; no tolerance gate, since §34's own P2 failed over
      its full grid for a diagnosed reason (the tau << T linearisation) and the same bias
      applies here.

WHAT THIS CANNOT SETTLE. Two substrates is not a law. If the transfer holds it says §34's
DERIVATION is substrate-independent while its k-independence is not; if it fails it says the
derivation needs something AM supplies and Schlogl does not. Either way m is computed from the
rails, so this is a test of the derivation, not of the number 0.71972.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import experiments.chemical_cascade as cc
from experiments.depth_compounding import C, R1, R2, R3

VOLUMES = (10, 14, 18, 22, 26, 30, 36, 42)
KS = (3, 5, 7)
# The exact MFPT solve returns NEGATIVE lifetimes above N ~ 150 -- the same instrument
# ceiling §34 hit on AM at N <= 72. Every volume used for a crossover must satisfy
# k*N <= MFPT_CEILING, and hold_lifetime refuses to return a non-positive number.
MFPT_CEILING = 150


def merge_fraction():
    """The fraction of failed tanks at which a merged pool drops below the saddle."""
    return (R3 - R2) / (R3 - R1)


def m_of_k(k):
    """Failures needed to flip the MERGE, from the rails -- majority only if symmetric."""
    return int(math.floor(merge_fraction() * k)) + 1


def x_merged(k, j):
    return ((k - j) * R3 + j * R1) / k


def _generator(N, cap_mult=1.25):
    cap = int(np.ceil(cap_mult * R3 * N))
    m = cap + 1
    rows, cols, vals = [], [], []
    for s in range(m):
        tot = 0.0
        lam, mu = cc.rates_stage(float(s), 0.0, N, C, R3, True, "hill")
        if s < cap and lam > 0:
            rows.append(s); cols.append(s + 1); vals.append(lam); tot += lam
        if s > 0 and mu > 0:
            rows.append(s); cols.append(s - 1); vals.append(mu); tot += mu
        rows.append(s); cols.append(s); vals.append(-tot)
    return sp.csr_matrix((vals, (rows, cols)), shape=(m, m)), m, cap


def hold_lifetime(N):
    """MFPT from the high rail to below the saddle -- one tank of volume N holding a bit."""
    Q, m, cap = _generator(N)
    a = int(np.ceil(R2 * N)) - 1
    idx = np.arange(a, cap + 1)
    A = np.zeros((len(idx), len(idx)))
    rhs = -np.ones(len(idx))
    Qd = Q.toarray()
    for i, s in enumerate(idx):
        if s == a:
            A[i, i] = 1.0
            rhs[i] = 0.0
            continue
        A[i, i] = Qd[s, s]
        if s + 1 <= cap:
            A[i, i + 1] = Qd[s, s + 1]
        if s - 1 >= a:
            A[i, i - 1] = Qd[s, s - 1]
    T = np.linalg.solve(A, rhs)
    val = float(T[list(idx).index(int(round(R3 * N)))])
    if not (val > 0.0):
        raise FloatingPointError(
            f"MFPT solve returned {val:.4e} at N={N}: the dense solve has lost conditioning. "
            f"§111 caps volumes at k*N <= {MFPT_CEILING} for exactly this reason.")
    return val


def merge_commit(k, j, N, t_settle=None):
    """P(merged tank of volume k*N ends HIGH), started at the pooled concentration.

    An actual Schlogl tank: no sign(), no comparison, the chemistry decides.
    """
    kN = k * N
    Q, m, cap = _generator(kN)
    n0 = int(round(x_merged(k, j) * kN))
    n0 = min(max(n0, 0), cap)
    p = np.zeros(m)
    p[n0] = 1.0
    if t_settle is None:
        t_settle = 20.0 / abs(-3 * C[1] * R3 ** 2 + 2 * C[0] * R3 - C[3])
    p = spla.expm_multiply(Q.T * t_settle, p)
    return float(p[np.arange(m) > R2 * kN].sum())


def remerge_lifetime(T, k, tau):
    m = m_of_k(k)
    return T ** m / (math.comb(k, m) * tau ** (m - 1))


def predicted_crossover(c, a, k, tau):
    m = m_of_k(k)
    if k == m:
        return float("inf")          # no crossover of this form: the denominator vanishes
    return ((m - 1) * (a - math.log(tau)) - math.log(math.comb(k, m))) / (c * (k - m))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/is_the_vote_a_majority.json"))
    args = ap.parse_args()

    print(f"merge fraction (r3-r2)/(r3-r1) = {merge_fraction():.5f}"
          f"   (AM's symmetric value would be 0.5)")
    print(f"{'k':>4}{'m (Schlogl)':>14}{'m (majority)':>14}{'(m-1)/(k-m)':>14}")
    for k in KS:
        m = m_of_k(k)
        maj = (k + 1) // 2
        r = "undefined" if k == m else f"{(m - 1) / (k - m):.4f}"
        print(f"{k:>4}{m:>14}{maj:>14}{r:>14}")

    print("\nP1 -- is the hold lifetime exponential in volume?")
    Ts = []
    for N in VOLUMES:
        T = hold_lifetime(N)
        Ts.append(T)
        print(f"  N = {N:>3}   T = {T:.6e}   ln T = {math.log(T):.4f}", flush=True)
    Ns = np.array(VOLUMES, float)
    lnT = np.log(np.array(Ts))
    c, a = np.polyfit(Ns, lnT, 1)
    resid = lnT - (c * Ns + a)
    r2 = 1 - (resid ** 2).sum() / ((lnT - lnT.mean()) ** 2).sum()
    print(f"  ln T = {c:.6f} N + {a:.4f}   R2 = {r2:.6f}   max|resid| = {abs(resid).max():.4f}")

    print("\nP2 -- where does the merge actually flip? (a real tank, no free comparison)")
    print("  reported as a SOFT threshold: the deterministic prediction is where x_merged")
    print("  crosses the saddle, but a pool landing within a fluctuation of it is a coin flip,")
    print("  so the margin x_merged - r2 is quoted with every cell (rule 20).")
    rows = []
    for k in KS:
        N = 14
        ps = [merge_commit(k, j, N) for j in range(k + 1)]
        margins = [x_merged(k, j) - R2 for j in range(k + 1)]
        step = next((j for j, p in enumerate(ps) if p < 0.5), None)
        det = m_of_k(k)
        soft = [j for j, p in enumerate(ps) if 0.05 < p < 0.95]
        rows.append({"k": k, "N": N, "p_high": ps, "margins": margins,
                     "hard_step": step, "predicted_m": det,
                     "majority_m": (k + 1) // 2, "soft_cells": soft})
        print(f"  k={k}: predicted m = {det}   (majority would be {(k + 1) // 2})")
        for j, (p, mg) in enumerate(zip(ps, margins)):
            flag = "  <-- SOFT" if 0.05 < p < 0.95 else ""
            print(f"        j={j}: P(high) = {p:.4f}   x_merged - r2 = {mg:+.4f}{flag}")
        print(f"        first j with P < 0.5: {step}", flush=True)

    print("\nP3/P4 -- crossovers")
    tau = 1.0 / abs(-3 * C[1] * R3 ** 2 + 2 * C[0] * R3 - C[3])
    print(f"  tau = t_relax = {tau:.4f}")
    xo = []
    for k in KS:
        pred = predicted_crossover(c, a, k, tau)
        vols = [N for N in VOLUMES if k * N <= MFPT_CEILING]
        lo = []
        for N in vols:
            lo.append(math.log(hold_lifetime(k * N))
                      - math.log(remerge_lifetime(hold_lifetime(N), k, tau)))
        sign_change = [i for i in range(len(lo) - 1) if lo[i] * lo[i + 1] < 0]
        meas = None
        if sign_change:
            i = sign_change[0]
            f = -lo[i] / (lo[i + 1] - lo[i])
            meas = vols[i] + f * (vols[i + 1] - vols[i])
        xo.append({"k": k, "predicted": pred, "measured": meas, "volumes": vols,
                   "ln_ratio": lo,
                   "ratio": (pred / meas) if (meas and np.isfinite(pred)) else None})
        pr = "inf (no crossover)" if not np.isfinite(pred) else f"{pred:.3f}"
        ms = "none in range" if meas is None else f"{meas:.3f}"
        print(f"  k={k}: predicted Omega_x = {pr:>18}   measured = {ms}"
              f"   (volumes {vols[0]}-{vols[-1]}, capped by k*N <= {MFPT_CEILING})", flush=True)
        print(f"        ln(L_hold/L_remerge) = "
              + ", ".join(f"{v:+.2f}" for v in lo))

    args.out.write_text(json.dumps(
        {"merge_fraction": merge_fraction(), "c": c, "a": a, "r2": r2,
         "hold": [{"N": int(n), "T": T} for n, T in zip(VOLUMES, Ts)],
         "merge": rows, "crossovers": xo}, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
