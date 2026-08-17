"""T-CASC-a: the depth ceiling in the SATURATED regime, with no matrix exponential at all

§71 could not reach the Omega-independent regime §12's ceiling formula describes, and §71.2
showed the wall is structural rather than a matter of machine size: propensities are
extensive, so ||Qt|| ~ Omega, and every Krylov exponential costs O(||Qt||) while dense
scaling-and-squaring costs O(n^3). Both wall out near 1e4 states. §71.2 proposed the second
eigenvalue of the stage map as the reachable instrument -- **but that still needs exp(Qt)v,
so it inherits the same wall.**

**The right reduction is cheaper and is exactly the regime in question.** §12's ceiling is a
statement about the limit where each stage FULLY RESTORES: the per-stage error then saturates
at an Omega-independent floor set by the channel. In that limit the chemistry contributes one
number per rail -- the probability of committing to the wrong rail -- and

    per-stage error  eps = sum_n  q(n) * p_cross(n)

where q is the channel's output distribution from a rail and **p_cross is the exact
birth-death splitting probability**, which §61 and §69 already compute in closed form from the
scale function, in logs, with no solve. **Cost O(cap) per cell and no exponential anywhere**,
so Omega is limited only by memory.

The cascade is then a two-state Markov chain. The element is ASYMMETRIC, so the two rails have
DIFFERENT error probabilities and this is a binary ASYMMETRIC channel -- eps_hi != eps_lo --
which the two-state chain handles exactly and a naive (1-2eps)^D would not.

**What this measures and what it does not.** It is the t_stage -> infinity limit. That is the
regime §12's formula describes and the one §71 failed to reach, but it is an idealisation: a
real cascade with finite stage time restores less, so this is an UPPER bound on D_max. Stated
here rather than discovered later.

PREDICTIONS, written before running.

  P1  GATE. As t_stage grows, §71's finite-time cascade must approach this limit. Compare
      D_max from `run` at increasing t against the saturated value at the same Omega. If they
      do not converge, the two are not measuring the same thing and nothing below counts.
  P2  GATE. p_cross must be a probability in [0,1], equal 1/2 at the saddle by construction of
      the splitting problem, and go to 0 at the high rail and 1 at the low rail.
  P3  **THE TEST. Does D_max saturate in Omega?** §12 says it must: the per-stage error has an
      Omega-independent floor. §71 saw it still climbing at Omega = 1300 with 84.8% spread.
      Here Omega can run to 1e5, so the question is answerable rather than merely open.
  P4  **THE COMPARISON, absolute (rule 16).** Ratio of measured D_max to exp(Delta^2/2 sigma^2)/4,
      against AM's PUBLISHED ratios (3.00, 3.38, 3.33 at sigma/Delta = 0.45, 0.35, 0.28).
      **Gated on P3** -- if D_max has not saturated the ratio is withheld, as §71 withheld it.
  P3' **SATURATION TEST, second version.** The first demanded a spread under 5% across the
      last three Omega. D_max is an INTEGER depth, so at D_max ~ 7 one unit of quantisation is
      14% and the test could never pass there however saturated the answer was -- the f = 0.45
      column read 6, 7, 7, 7, 7, 7, 8 and was declared "still moving". A criterion must be
      satisfiable in the regime it is applied to (rule 19). It is now
      spread <= max(5%, 1.5 units / mean), so quantisation cannot masquerade as drift.
  P5  **RULE 9.** Sweep the rail geometry as well as Omega: a ceiling that is a property of
      sigma/Delta must not move when Delta is held and the rails are placed asymmetrically.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from experiments.cascade_schlogl import rates, schlogl_consts

AM_RATIO = {0.45: 9.0 / 3.0, 0.35: 50.0 / 14.8, 0.28: 489.0 / 147.0}


def p_cross(omega, r1, r2, r3, cap_mult=1.6):
    """P(hit the LOW rail before the HIGH rail | start n), exactly, in logs.

    The scale function of a birth-death chain: pi_k = prod (mu/lambda), and
    P(low first | n) = sum_{k=n}^{hi-1} pi_k / sum_{k=lo}^{hi-1} pi_k.
    """
    c = schlogl_consts(r1, r2, r3)
    cap = int(np.ceil(cap_mult * r3 * omega))
    lam, mu = rates(omega, c, cap)
    lo, hi = int(round(r1 * omega)), int(round(r3 * omega))
    lp = np.full(cap + 1, -np.inf)
    acc = 0.0
    lp[lo] = 0.0
    for k in range(lo + 1, hi + 1):
        if lam[k] <= 0 or mu[k] <= 0:
            return None, None, None
        acc += np.log(mu[k]) - np.log(lam[k])
        lp[k] = acc

    def lse(v):
        v = v[np.isfinite(v)]
        if v.size == 0:
            return -np.inf
        M = v.max()
        return M + np.log(np.exp(v - M).sum())

    # Suffix log-sum-exp in O(cap): the per-n loop was O(cap^2) and could not reach the
    # Omega this experiment exists to test.
    seg = lp[lo:hi].copy()
    M = seg[np.isfinite(seg)].max()
    e = np.where(np.isfinite(seg), np.exp(seg - M), 0.0)
    suff = np.cumsum(e[::-1])[::-1]          # suff[i] = sum_{k>=i} e[k]
    out = np.zeros(cap + 1)
    out[:lo + 1] = 1.0
    out[hi:] = 0.0
    out[lo:hi] = suff / suff[0]
    return out, lo, hi


def eps_pair(omega, r1, r2, r3, sigma, cap_mult=1.6):
    """Per-stage error from each rail: channel smear, then exact commitment."""
    pc, lo, hi = p_cross(omega, r1, r2, r3, cap_mult)
    if pc is None:
        return None
    cap = len(pc) - 1
    x = np.arange(cap + 1) / omega
    g_hi = np.exp(-0.5 * ((x - r3) / sigma) ** 2)
    g_lo = np.exp(-0.5 * ((x - r1) / sigma) ** 2)
    g_hi /= g_hi.sum()
    g_lo /= g_lo.sum()
    e_hi = float(g_hi @ pc)                 # started high, committed low = error
    e_lo = float(g_lo @ (1.0 - pc))         # started low, committed high = error
    return e_hi, e_lo


def mutual_info_depth(e_hi, e_lo, depth):
    """Exact I(input ; output) after `depth` stages of the binary ASYMMETRIC channel."""
    T = np.array([[1.0 - e_lo, e_lo], [e_hi, 1.0 - e_hi]])   # rows: true low, true high
    P = np.eye(2)
    out = []

    def h(v):
        v = np.asarray(v, float)
        v = v[v > 0]
        return float(-(v * np.log2(v)).sum())

    for _ in range(depth):
        P = P @ T
        marg = 0.5 * (P[0] + P[1])
        out.append(h(marg) - 0.5 * (h(P[0]) + h(P[1])))
    return np.array(out)


def mutual_info_at(e_hi, e_lo, depth):
    """I at a single depth in CLOSED FORM, so no O(depth) iteration.

    T has eigenvalues 1 and lam = 1 - e_hi - e_lo, so T^D = pi (+) lam^D (I - pi), i.e.
    T^D[i,j] = pi_j + lam^D (delta_ij - pi_j). Added because the iterative version costs
    O(depth) per evaluation and §74 needs depths of 1e10, where it simply does not return.
    Gated against `mutual_info_depth` at moderate depth.
    """
    tot = e_hi + e_lo
    if tot <= 0:
        return 1.0
    lam = 1.0 - tot
    pi = np.array([e_hi / tot, e_lo / tot])
    ld = lam ** depth if abs(lam) < 1 else 1.0
    P = np.array([pi + ld * (np.array([1.0, 0.0]) - pi),
                  pi + ld * (np.array([0.0, 1.0]) - pi)])

    def h(v):
        v = np.asarray(v, float)
        v = v[v > 0]
        return float(-(v * np.log2(v)).sum())

    return h(0.5 * (P[0] + P[1])) - 0.5 * (h(P[0]) + h(P[1]))


def d_max_closed(e_hi, e_lo, level=0.5):
    """Depth where I falls through `level`, by bisection on the closed form. O(log)."""
    if max(e_hi, e_lo) <= 0 or (e_hi + e_lo) >= 1.0:
        return None
    lo, hi = 1.0, 2.0
    while mutual_info_at(e_hi, e_lo, hi) >= level:
        hi *= 2.0
        if hi > 1e18:
            return None
    while hi - lo > max(1e-9, 1e-9 * hi):
        mid = 0.5 * (lo + hi)
        if mutual_info_at(e_hi, e_lo, mid) >= level:
            lo = mid
        else:
            hi = mid
    return float(hi)


def d_max_saturated(e_hi, e_lo, level=0.5, cap_depth=10 ** 7):
    """Depth where I falls through `level`. Bisection on the closed-form chain."""
    if max(e_hi, e_lo) <= 0:
        return None
    lo, hi = 1, 2
    while hi < cap_depth:
        if mutual_info_depth(e_hi, e_lo, hi)[-1] < level:
            break
        lo, hi = hi, hi * 2
    if hi >= cap_depth:
        return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mutual_info_depth(e_hi, e_lo, mid)[-1] >= level:
            lo = mid
        else:
            hi = mid
    return float(hi)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r1", type=float, default=0.1)
    ap.add_argument("--r2", type=float, default=1.0)
    ap.add_argument("--r3", type=float, default=1.9)
    ap.add_argument("--omegas", type=int, nargs="+",
                    default=[400, 900, 1800, 3600, 7200, 14400, 28800])
    ap.add_argument("--fracs", type=float, nargs="+", default=[0.45, 0.35, 0.28])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/cascade_saturated.json"))
    args = ap.parse_args()
    D = (args.r3 - args.r1) / 2.0
    print(f"rails {args.r1}/{args.r3}, saddle {args.r2}, Delta = {D}")

    print("\n=== P2 GATE: is p_cross a proper splitting probability?")
    pc, lo, hi = p_cross(600, args.r1, args.r2, args.r3)
    sad = int(round(args.r2 * 600))
    print(f"  p_cross at low rail {pc[lo]:.6f}, saddle {pc[sad]:.6f}, high rail {pc[hi]:.6f}")
    print(f"  range [{pc.min():.3e}, {pc.max():.6f}], monotone decreasing:"
          f" {bool(np.all(np.diff(pc[lo:hi+1]) <= 1e-12))}")
    ok2 = (abs(pc[lo] - 1) < 1e-12 and abs(pc[hi]) < 1e-12
           and abs(pc[sad] - 0.5) < 0.05 and pc.min() >= -1e-12 and pc.max() <= 1 + 1e-12)
    print(f"  -> P2 {'HOLDS' if ok2 else 'FAILS'}")

    print("\n=== P1 GATE: does §71's finite-t cascade approach this limit as t grows?")
    from experiments.cascade_schlogl import run, d_max as dmax_ft
    om, f = 400, 0.45
    sat = d_max_saturated(*eps_pair(om, args.r1, args.r2, args.r3, f * D))
    print(f"  saturated D_max at Omega={om}, f={f}: {sat}")
    for t in (2.0, 6.0, 18.0):
        Is, _ = run(om, args.r1, args.r2, args.r3, t, 60, f * D)
        print(f"    t_stage={t:>5}: finite-time D_max = {dmax_ft(Is)}")
    print("  -> P1: the finite-time ceiling must rise toward the saturated value as t grows")

    print(f"\n=== P3/P4: does D_max saturate in Omega?")
    print(f"{'f':>6}" + "".join(f"{f'W={o}':>12}" for o in args.omegas))
    rows, sat_ok = [], {}
    for fr in args.fracs:
        sigma = fr * D
        ds = []
        for om in args.omegas:
            ep = eps_pair(om, args.r1, args.r2, args.r3, sigma)
            if ep is None:
                ds.append(None)
                continue
            d = d_max_saturated(*ep)
            ds.append(d)
            rows.append({"frac": fr, "omega": om, "e_hi": ep[0], "e_lo": ep[1], "dmax": d})
        print(f"{fr:>6.2f}" + "".join(f"{d:>12.1f}" if d else f"{'--':>12}" for d in ds))
        v = [d for d in ds[-3:] if d]
        if v:
            spread = (max(v) - min(v)) / np.mean(v)
            tol = max(0.05, 1.5 / np.mean(v))     # D_max is an integer: see P3'
            sat_ok[fr] = len(v) >= 3 and spread <= tol
            print(f"       last three spread {100*spread:>6.2f}% against tolerance"
                  f" {100*tol:>5.2f}% (quantisation-aware)"
                  f"   -> {'SATURATED' if sat_ok[fr] else 'still moving'}")
        else:
            sat_ok[fr] = False

    print(f"\n=== P4: ratio to exp(Delta^2/2 sigma^2)/4, gated on P3")
    for fr in args.fracs:
        pred = float(np.exp(1.0 / (2 * fr ** 2)) / 4.0)
        v = [r["dmax"] for r in rows if r["frac"] == fr and r["dmax"]]
        if not v:
            continue
        if not sat_ok.get(fr):
            print(f"  f={fr}: WITHHELD -- not saturated")
            continue
        ratio = v[-1] / pred
        am = AM_RATIO.get(fr)
        rel = abs(ratio - am) / am
        print(f"  f={fr}: measured {v[-1]:.1f}, predicted {pred:.2f}, ratio {ratio:.2f}"
              f"  against AM's {am:.2f} -- {100*rel:.0f}%"
              f"  -> {'TRANSFERS' if rel <= 0.5 else 'does NOT transfer'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"rows": rows, "Delta": D,
                                    "saturated": {str(k): bool(v) for k, v in sat_ok.items()}},
                                   indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
