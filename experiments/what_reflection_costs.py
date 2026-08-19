"""§100 -- what the reflecting boundary buys and what it hides.

§92.1(b) fixed a real conditioning bias by REFLECTING stage 1 at its saddle, and
§93-§98 all inherit that construction. A reflecting wall is exactly the family rule 10
warns about -- the harness doing something the chemistry cannot -- so it is used here
as a declared control rather than smuggled in, and its cost has never been measured.
Two questions, both answered on the same run and paired within it (rule 18).

QUESTION 1 -- SCOPE. With stage 1 reflected, stage 1 CANNOT FAIL. So the penalty that
§92-§98 measure is the transfer of stage 1's *rail fluctuations* into stage 2, not the
accumulation of stage 1's *errors*. §91 (unreflected) measured the latter. The depth
ceiling D_max = c*/(penalty x eps) needs both, and it has been built from the reflected
penalty alone.

QUESTION 2 -- T-CASC-l. §99(b) opened a suspect: §93 placed the composition maximum at
"matched relaxation times", but §93's knob is a GLOBAL speed scaling on stage 1, which
multiplies its relaxation time AND its escape time by the same factor. That axis cannot
say which timescale was matched (rule 9). Omega separates them: an escape time grows like
exp(A*Omega), a relaxation time is a macroscopic quantity and is Omega-INDEPENDENT.

PREDICTIONS, WRITTEN BEFORE RUNNING.

  P1  WIRING. The reflected two-stage chain at s_up = s_dn = 1, Omega = 30, t0 = 2.0
      reproduces §93's penalty of 4.4419. Same code path, so this must agree to solver
      tolerance; a mismatch means the harness moved under §93 and everything below is void.

  P2  DIRECTION, and it is certain in sign. Freeing stage 1 can only ADD stage-2 error,
      never remove it: reflection deletes trajectories in which stage 1 fell, and every
      one of those raises P(stage 2 low) or leaves it alone. So P_free >= P_refl. The
      measurement is the SIZE, and the size is the scope statement -- not a pass/fail.
      Reported as a ratio, with no tolerance attached (rule 20): if the excess is a
      fraction of a percent, reflection is a harmless control and §92-§98's scope is
      untouched; if it is comparable to P_refl itself, the arc's ceiling is CONDITIONAL
      on the upstream surviving and must say so.

  P3  DECOMPOSITION. The excess should factor as P(stage 1 low) x P(stage 2 low | stage 1
      low). I do NOT predict the second factor is near 1. The hill coupling saturates --
      that is its entire point, and §91 measured it holding a 3.39 sigma margin -- so a
      dip in stage 1 need not move stage 2. Both factors are reported separately, because
      a large product from a large P(s1 low) and a large product from perfect transmission
      are different physical statements.

  P4  T-CASC-l, THE DISCRIMINATOR. The spectral gap of stage 1's generator, at Omega =
      14, 30, 55:
        - REFLECTED: stage 1 has no escape channel at all, so its only correlation time is
          the intra-basin relaxation time. The gap must approach the macroscopic
          |f'(r3)| = 6.6195 and must NOT vary exponentially with Omega. Judged by whether
          the sequence CONVERGES toward that value (rule 20), not by a tolerance at any
          one Omega.
        - FREE: the gap is the escape rate and must collapse roughly like exp(-A*Omega),
          A = 0.190241, i.e. by ~2 orders of magnitude over this range.
      If both hold, §93's crossover matched a RELAXATION time because the construction left
      stage 1 no other timescale to match -- which settles T-CASC-l without appealing to
      the type-I/type-II classification at all, and simultaneously narrows what §93 claims.
      If the reflected gap instead tracks exp(-A*Omega), the reflection is not doing what
      the docstrings say and §92-§98 need re-reading.

WHAT THIS CANNOT SETTLE. It says nothing about whether the type-I/type-II distinction
would matter for a cascade whose upstream IS allowed to fail. That is a different chain
and it is not built here.
"""

from __future__ import annotations

import json
import pathlib

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import experiments.chemical_cascade as cc
from experiments.cascade_schlogl import schlogl_consts
from experiments.margin_law import R1, R2, R3, stage1_stationary
from experiments.timescale_ratio import build_pair, pinned_reference

C = schlogl_consts(R1, R2, R3)
OMEGA, T0 = 30, 2.0
A_UP = 0.190241                       # §94's escape action for this element
GAP_OMEGAS = (14, 30, 55)


def relax_rate():
    """Macroscopic relaxation rate at the high rail -- no Omega in it anywhere."""
    k1a, k1r, k2b, k2r = C
    return abs(-3.0 * k1r * R3**2 + 2.0 * k1a * R3 - k2r)


def build_free(om, s_up=1.0, s_dn=1.0, cap_mult=1.25):
    """Both stages FREE over 0..cap. Same rates and same clock as build_pair, no wall."""
    cap = int(np.ceil(cap_mult * R3 * om))
    m = cap + 1
    rows, cols, vals = [], [], []
    for n1 in range(m):
        for n2 in range(m):
            idx = n1 * m + n2
            tot = 0.0
            l1, u1 = cc.rates_stage(float(n1), 0.0, om, C, R3, True, "hill")
            l1 *= s_up; u1 *= s_up
            if n1 < cap and l1 > 0:
                rows.append(idx); cols.append(idx + m); vals.append(l1); tot += l1
            if n1 > 0 and u1 > 0:
                rows.append(idx); cols.append(idx - m); vals.append(u1); tot += u1
            l2, u2 = cc.rates_stage(float(n2), float(n1), om, C, R3, False, "hill")
            l2 *= s_dn; u2 *= s_dn
            if n2 < cap and l2 > 0:
                rows.append(idx); cols.append(idx + 1); vals.append(l2); tot += l2
            if n2 > 0 and u2 > 0:
                rows.append(idx); cols.append(idx - 1); vals.append(u2); tot += u2
            rows.append(idx); cols.append(idx); vals.append(-tot)
    return sp.csr_matrix((vals, (rows, cols)), shape=(m * m, m * m)), m, cap


def stage1_generator(om, reflected, cap_mult=1.25):
    """Stage 1 alone, with or without the wall at its saddle."""
    cap = int(np.ceil(cap_mult * R3 * om))
    lo = int(np.ceil(R2 * om)) if reflected else 0
    ns = np.arange(lo, cap + 1)
    k = len(ns)
    rows, cols, vals = [], [], []
    for i, n in enumerate(ns):
        tot = 0.0
        l, u = cc.rates_stage(float(n), 0.0, om, C, R3, True, "hill")
        if i < k - 1 and l > 0:
            rows.append(i); cols.append(i + 1); vals.append(l); tot += l
        if i > 0 and u > 0:
            rows.append(i); cols.append(i - 1); vals.append(u); tot += u
        rows.append(i); cols.append(i); vals.append(-tot)
    return sp.csr_matrix((vals, (rows, cols)), shape=(k, k)), ns


def spectral_gap(om, reflected):
    """Slowest non-zero relaxation rate of stage 1: -Re(lambda_1) of the generator."""
    Q, ns = stage1_generator(om, reflected)
    k = Q.shape[0]
    ev = np.linalg.eigvals(Q.toarray()) if k < 400 else \
        spla.eigs(Q.T.tocsc(), k=3, which="SM", return_eigenvectors=False)
    ev = np.sort(-np.real(ev))            # 0 first, then the gap
    return float(ev[1]), k


def run_scope(om=OMEGA, t0=T0):
    """Paired within one run: reflected vs free, identical rates, clock and seed law."""
    t = t0
    _, pi1 = stage1_stationary(om)

    Qr, up, m2, cap = build_pair(om, 1.0, 1.0)
    pr = np.zeros(len(up) * m2)
    for a, w in enumerate(pi1):
        pr[a * m2 + int(round(R3 * om))] = w
    pr = spla.expm_multiply(Qr.T * t, pr)
    lo_r = (np.arange(len(up) * m2) % m2) < R2 * om
    p_refl = float(pr[lo_r].sum())

    Qf, m, _ = build_free(om)
    pf = np.zeros(m * m)
    nsad = int(np.ceil(R2 * om))
    for a, w in enumerate(pi1):                       # same seed law, embedded at n1 = up[a]
        pf[(nsad + a) * m + int(round(R3 * om))] = w
    pf = spla.expm_multiply(Qf.T * t, pf)
    grid = pf.reshape(m, m)
    n1lo = np.arange(m) < R2 * om
    p_free = float(grid[:, np.arange(m) < R2 * om].sum())
    p_s1_lo = float(grid[n1lo, :].sum())
    p_both = float(grid[np.ix_(n1lo, np.arange(m) < R2 * om)].sum())
    cond = p_both / p_s1_lo if p_s1_lo > 0 else float("nan")

    den = pinned_reference(om, 1.0, t)
    return {"omega": om, "t": t0, "p_refl": p_refl, "p_free": p_free,
            "pen_refl": p_refl / den, "pen_free": p_free / den, "den": den,
            "excess": p_free - p_refl, "ratio": p_free / p_refl,
            "p_s1_lo": p_s1_lo, "p_s2lo_given_s1lo": cond, "product": p_s1_lo * cond}


def main():
    print(f"macroscopic relaxation rate |f'(r3)| = {relax_rate():.4f}  (no Omega in it)\n")

    s = run_scope()
    print(f"P1  reflected penalty at s_up = s_dn = 1, Omega = {s['omega']}, t0 = {s['t']}: "
          f"{s['pen_refl']:.4f}   (§93 measured 4.4419)")

    print(f"\nP2  stage-2 error, paired in one run:")
    print(f"      reflected upstream  P = {s['p_refl']:.6e}   penalty {s['pen_refl']:.4f}")
    print(f"      free upstream       P = {s['p_free']:.6e}   penalty {s['pen_free']:.4f}")
    print(f"      excess {s['excess']:.6e}   ratio free/reflected = {s['ratio']:.4f}")

    print(f"\nP3  decomposition of the excess:")
    print(f"      P(stage 1 low)               = {s['p_s1_lo']:.6e}")
    print(f"      P(stage 2 low | stage 1 low) = {s['p_s2lo_given_s1lo']:.4f}")
    print(f"      product                      = {s['product']:.6e}")

    print(f"\nP4  stage-1 spectral gap -- the T-CASC-l discriminator")
    print(f"{'Omega':>7}{'reflected':>13}{'states':>8}{'free':>14}{'states':>8}"
          f"{'exp(-A*Om)':>13}")
    gaps = []
    for om in GAP_OMEGAS:
        gr, kr = spectral_gap(om, True)
        gf, kf = spectral_gap(om, False)
        gaps.append({"omega": om, "reflected": gr, "free": gf})
        print(f"{om:>7}{gr:>13.4f}{kr:>8}{gf:>14.4e}{kf:>8}"
              f"{np.exp(-A_UP * om):>13.4e}")
    r = [g["reflected"] for g in gaps]
    f = [g["free"] for g in gaps]
    print(f"\n      reflected gap varies by {max(r)/min(r):.3f}x over Omega {GAP_OMEGAS[0]}"
          f"-{GAP_OMEGAS[-1]}; target |f'(r3)| = {relax_rate():.4f}")
    print(f"      free gap varies by {max(f)/min(f):.3e}x; exp(-A*Omega) varies by "
          f"{np.exp(-A_UP*GAP_OMEGAS[0])/np.exp(-A_UP*GAP_OMEGAS[-1]):.3e}x")

    out = pathlib.Path("results/what_reflection_costs.json")
    out.write_text(json.dumps({"scope": s, "gaps": gaps,
                               "relax_rate": relax_rate()}, indent=2))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
