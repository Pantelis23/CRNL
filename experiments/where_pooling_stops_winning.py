"""§112 (T16-c, T16-d) -- the merge threshold is rail geometry, and it decides which protocol wins.

§111 found that a merged Schlogl pool crosses the saddle only when j > 0.71972 k, not at a
majority, and that at k = 3 this means UNANIMITY: re-merging never loses, at any volume. It
tested one element and left the number 0.71972 looking like a Schlogl fact.

It is not. For any bistable one-species element, pooling k tanks of volume Omega conserves
counts, so the merged concentration is the count-weighted mean

    x_merged(k, j) = [ (k-j) r3 + j r1 ] / k ,

which falls below the saddle exactly when j > f k, with

    f = (r3 - r2) / (r3 - r1)          -- three roots, no dynamics, no fitting

so the number of failures that flips the pool is m(k) = floor(f k) + 1. Substituting into §34's
crossover law Omega_x = [(m-1)(a - ln tau) - ln C(k,m)] / [c (k - m)] and expanding gives three
regimes, and the boundaries between them are pure rail geometry:

    f >= 1 - 1/k   ->  m = k.  Denominator vanishes. ln(L_hold/L_remerge) = (k-1)(ln tau - a),
                       CONSTANT in Omega: re-merging wins at every volume, forever.
    f <  1/k       ->  m = 1.  ln(L_hold/L_remerge) = ln T(k*Omega) - ln T(Omega) + ln k > 0 and
                       rising: pooling wins at every volume, forever.
    otherwise      ->  a finite crossover, at the volume §34 wrote down.

**So §32's finite crossover between pooling and re-merging is not generic.** It exists only in a
window of rail geometry, and which protocol wins outside that window is fixed by one number
computed from three roots. AM sits at f = 1/2 by construction, in the middle of the window for
every k >= 3 -- which is why §32 found a crossover at all.

The elements here are constructed: r1 = 0.15 and r3 = 3.0 are held fixed and the saddle r2 is
moved, which sweeps f from 0.21 to 0.84 while `schlogl_consts` returns the exact rate constants
for those three roots (verified in P1). They are real chemistry, not a reparametrisation of one
element: c and a differ by 9x across the set.

TWO INSTRUMENT CHANGES, both removing something the harness was doing (rule 10):

  (i)  §111's dense MFPT solve returned NEGATIVE lifetimes above N ~ 150 and was capped there.
       A 1-D birth-death MFPT has an exact all-positive downward recursion,
       d_{n-1} = (lam_n d_n + 1)/mu_n with d_{cap-1} = 1/mu_cap and T = sum d_n, which in log
       domain cannot lose conditioning at any volume. The cap goes away. P1 checks it.
  (ii) §111 measured commitment by evolving expm_multiply for a settle time of 20 tau. That is
       a free knob, and it is wrong for a shallow well, where the tank escapes during the
       settle. In 1-D the commitment probability is the SPLITTING probability between the two
       rails, which is exact and has no time parameter at all.
  (iii) ADDED AFTER THE PREDICTIONS WERE WRITTEN, because rule 13 says to sweep an
       approximation's own numerical parameter: the reflecting box at cap_mult = 1.25, inherited
       from the cascade code via §111, is NOT converged -- ln T is off by up to 0.249 there and
       rises monotonically with the box, since the wall truncates excursions above the rail and
       can only shorten the escape. 3.0 and 4.0 agree exactly on every element, so §112 runs at
       3.0. P0 measures it. This LOOSENS §111's headline from 0.31% to 1.77% and makes §112's own
       crossover residuals worse (median 3.3% -> 8.4%); both are reported. The splitting
       probability is exactly box-independent, so P2 and P2b are untouched.

PREDICTIONS, WRITTEN BEFORE RUNNING.

  P1  THE INSTRUMENT. On §111's own rails the log-domain recursion must reproduce its dense
      solve to better than 1e-9 in ln T for N <= 60, and stay finite, positive and asymptotically
      linear past N = 150 where the dense solve fails. §111's published headline is recomputed on
      the new instrument and printed BESIDE the original; the original stands either way (rule
      7). If the two disagree at small N the new instrument is wrong and nothing below is
      readable.

  P2  m IS RAIL GEOMETRY, ACROSS SIX ELEMENTS. For each element and k in {3, 5, 7}, the measured
      step in P(commit high | j) -- a real merged tank of volume k*Omega, splitting probability
      between its own two rails, no free comparison -- must land at m = floor(f k) + 1, NOT at
      the majority ceil((k+1)/2). Reported per cell with the margin x_merged - r2. If the step
      tracks the majority instead of f, m is not rail geometry and P3-P5 are void.

  P2b T16-c, CLOSED EITHER WAY. §111 left a straddling cell (k = 7, j = 5) whose pool sits 1.6%
      ABOVE the saddle yet gave P(high) = 0.4930. The deterministic premise says a positive
      margin must drive P(high) -> 1 as Omega grows, and a negative one -> 0. Every soft cell
      (0.05 < P < 0.95) found in P2 is swept in Omega. Criterion is the LIMIT, not a tolerance
      (rule 20): P must move monotonically toward the sign of its own margin. **If a soft cell
      with a positive margin instead goes to 0, the merged tank's commitment is not decided by
      which basin its mean lies in, and every m in this section and in §111 is wrong for a
      reason that has nothing to do with symmetry.**

  P3  THE PHASE BOUNDARY, MEASURED WITHIN ONE ELEMENT. The unanimity condition f >= 1 - 1/k is
      k-dependent, so a single element can sit on both sides of it. E-084 (f = 0.8421) has
      1 - 1/3 = 0.667 and 1 - 1/5 = 0.8 below it but 1 - 1/7 = 0.857 above it. So the SAME
      element, same rails, same c and a, must show ln(L_hold/L_remerge) FLAT in Omega at k = 3
      and k = 5, and RISING THROUGH ZERO at k = 7. Nothing cross-run enters (rule 18). The flat
      value must equal (k-1)(ln tau - a) with c and a from the hold protocol alone.

  P4  THE OTHER EDGE. E-021 (f = 0.2105 < 1/3) has m = 1 at k = 3: one failed tank flips the
      pool. Then ln(L_hold/L_remerge) = ln T(3 Omega) - ln T(Omega) + ln 3, positive and RISING
      at every Omega -- pooling wins everywhere, no crossover on that side either. The
      prediction is parameter-free; what is measured is whether the merged tank agrees that one
      failure suffices, which is P2's E-021 row at k = 3.

  P5  WHAT THE STRADDLE COSTS. E-062 (f = 0.62, m = 2) and E-070 (f = 0.70, m = 3) sit either
      side of 2/3 at k = 3 and differ in r2 by 0.228. They must behave QUALITATIVELY
      differently -- finite crossover vs none. Both curves reported in full.

  P5b THE FIT WINDOW, stated because it is the one soft number here. ln T = c N + a has
      curvature: the local slope on §111's rails drifts 0.18983 -> 0.19020 between N = 150 and
      N = 1000. c and a therefore depend on the window and so does every absolute Omega_x. Two
      disjoint windows are fitted per element and the spread in Omega_x quoted (rule 15, rule
      21). **If the window spread exceeds the residual being claimed, the absolute prediction is
      unresolved and says so.**

  P6  ARITHMETIC, NOT A MEASUREMENT, and labelled so because rule 21's second half was bought
      with a section that dressed the expansion of its own formula as a discovery. §34's
      k-independence needs (m-1)/(k-m) = 1, i.e. m = (k+1)/2, i.e. floor(f k) = (k-1)/2, i.e.

          1/2 - 1/(2k)  <=  f  <  1/2 + 1/(2k)

      -- a band around the symmetric point that shrinks like 1/k, so the headline holds for ALL
      k only at f = 1/2 exactly. Tabulated to bound how much asymmetry §34 tolerates. This is
      expansion of a formula already published in §34 plus §111's m. It is not a finding.

A CRITERION CHANGED AFTER SEEING OUTPUT, DISCLOSED HERE BECAUSE RULE 19 IS ABOUT EXACTLY THIS.
P2 as first coded read the step at a single Omega and scored `first j with P < 0.5`. Run once, it
printed MISMATCH for E-021 at k = 5, where the cell in question has a margin of +0.03 and gave
P = 0.454 -- a hard gate on a quantity a fluctuation wide, the rule-20 error this section was
written to avoid. The scoring was replaced with the limit test P2b had ALREADY specified in
writing: every unsaturated cell is swept in Omega and must move toward the sign of its own
margin, and the step is read at the largest volume, where cells have saturated. The replacement
is strictly more conservative -- it can no longer pass on a coin-flip cell -- and it was made
before P3-P6 were read. The original P2 wording is left above unedited (rule 3).

WHAT THIS CANNOT SETTLE. Every element is Schlogl-form: one species, cubic drift, rails moved by
relocating the saddle. f is a single number summarising three roots; two elements with the same
f but different well SHAPES share m but not c or a, so their crossover volumes differ. This
tests that m -- and therefore which protocol wins in the limits -- is fixed by f, not that f is
the only thing that matters. It says nothing about multi-species elements, where the merged
pool's state need not be the count-weighted mean of the committed ones.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib

import numpy as np

import experiments.chemical_cascade as cc
from experiments.cascade_schlogl import schlogl_consts

R1 = 0.15
R3 = 3.0
KS = (3, 5, 7)
# The reflecting box is an approximation and its width is a second axis (rule 13). §111 and the
# cascade code inherited 1.25, which is NOT converged: ln T still moves by up to 0.47 between
# cap_mult 1.15 and 2.0, always upward, because the wall truncates excursions above the rail and
# so shortens the escape. Successive differences shrink by 3-4 orders per step and 3.0 vs 4.0 is
# exactly zero on every element, so §112 runs at 3.0 and the box contributes nothing. P0 measures it.
CAP_MULT = 3.0

# r1 and r3 fixed, the saddle moved: f = (r3 - r2)/(r3 - r1) runs 0.21 -> 0.84.
ELEMENTS = (
    ("E-021", 2.4000),
    ("E-045", 1.7175),
    ("E-050", 1.5750),
    ("E-062", 1.2330),
    ("E-070", 1.0050),
    ("E-084", 0.6060),
)


class Element:
    """A bistable Schlogl element specified by its three rails."""

    def __init__(self, name, r1, r2, r3, cap_mult=CAP_MULT):
        self.name, self.r1, self.r2, self.r3 = name, r1, r2, r3
        self.cap_mult = cap_mult
        self.c = schlogl_consts(r1, r2, r3)
        if not all(v > 0 for v in self.c):
            raise ValueError(f"{name}: rails {(r1, r2, r3)} give a non-positive rate constant")
        k1a, k1r, k2b, k2r = self.c
        self.tau = 1.0 / abs(-3 * k1r * r3 ** 2 + 2 * k1a * r3 - k2r)

    @property
    def f(self):
        return (self.r3 - self.r2) / (self.r3 - self.r1)

    def m_of_k(self, k):
        """Failures needed to flip a merged pool of k tanks. Majority only if f = 1/2."""
        return int(math.floor(self.f * k)) + 1

    def x_merged(self, k, j):
        return ((k - j) * self.r3 + j * self.r1) / k

    def roots(self):
        k1a, k1r, k2b, k2r = self.c
        z = np.roots([-k1r, k1a, -k2r, k2b])
        return np.sort([v.real for v in z if abs(v.imag) < 1e-9])

    def _rates(self, N):
        """Vectorised copy of `cc.rates_stage(..., first=True)`, asserted equal to it in tests.

        Written out only because the converged box needs ~40k sites at the largest volumes and
        the scalar loop is too slow there. The chemistry is not restated -- it is the same two
        expressions, and `test_vectorised_rates_match_the_shared_helper` pins them elementwise.
        """
        cap = int(np.ceil(self.cap_mult * self.r3 * N))
        k1a, k1r, k2b, k2r = self.c
        n = np.arange(cap + 1, dtype=float)
        mu = k1r * n * (n - 1.0) * (n - 2.0) / N ** 2 + k2r * n
        lam = k1a * n * (n - 1.0) / N + k2b * N
        np.clip(lam, 0.0, None, out=lam)
        np.clip(mu, 0.0, None, out=mu)
        lam[cap] = 0.0                      # reflecting at the top of the box
        return lam, mu, cap

    def log_hold_lifetime(self, N):
        """ln MFPT, high rail -> below the saddle, by an exact all-positive recursion.

        d_n = T_{n+1} - T_n obeys d_{n-1} = (lam_n d_n + 1)/mu_n downward from
        d_{cap-1} = 1/mu_cap, and T_{n0} = sum_{n=a}^{n0-1} d_n. Every term is a sum of
        positive quantities, so there is no cancellation at any volume -- which is what
        §111's dense solve lost above N ~ 150 (it returned NEGATIVE lifetimes).
        """
        lam, mu, cap = self._rates(N)
        a = int(np.ceil(self.r2 * N)) - 1
        n0 = int(round(self.r3 * N))
        if not (a < n0 <= cap):
            raise ValueError(f"{self.name}: N={N} too small -- saddle site {a} vs rail {n0}")
        ld = np.full(cap, -np.inf)
        ld[cap - 1] = -math.log(mu[cap])
        for n in range(cap - 1, a, -1):
            ld[n - 1] = np.logaddexp(math.log(lam[n]) + ld[n], 0.0) - math.log(mu[n])
        seg = ld[a:n0]
        mx = seg.max()
        return float(mx + math.log(np.exp(seg - mx).sum()))

    def commit_high(self, k, j, N):
        """P(a merged pool of k tanks of volume N settles HIGH), exactly.

        The splitting probability between the element's own two rails, started at the pooled
        concentration. §111 evolved expm_multiply for 20 tau; that settle time is a free knob
        and it is wrong for a shallow well, where the tank escapes while it is settling. In
        1-D the answer is exact and has no time parameter:

            P_n(hit b before a) = (sum_{i=a}^{n-1} rho_i) / (sum_{i=a}^{b-1} rho_i),
            rho_i = prod_{s=a+1}^{i} mu_s / lam_s.
        """
        kN = k * N
        lam, mu, cap = self._rates(kN)
        a = max(int(round(self.r1 * kN)), 0)
        b = min(int(round(self.r3 * kN)), cap)
        n = int(round(self.x_merged(k, j) * kN))
        n = min(max(n, a), b)
        lr = np.zeros(b - a + 1)                     # ln rho_i, i = a..b
        for i in range(a + 1, b + 1):
            lr[i - a] = lr[i - a - 1] + math.log(mu[i]) - math.log(lam[i])
        mx = lr.max()
        w = np.exp(lr - mx)
        num = w[: n - a].sum()
        den = w[: b - a].sum()
        return float(num / den)

    def ln_remerge(self, lnT, k, m):
        """ln of §32's re-merge protocol lifetime: T^m / (C(k,m) tau^(m-1))."""
        return m * lnT - math.log(math.comb(k, m)) - (m - 1) * math.log(self.tau)

    def ln_ratio(self, N, k):
        """ln(L_hold / L_remerge): positive means POOLING wins, negative means re-merging."""
        m = self.m_of_k(k)
        return self.log_hold_lifetime(k * N) - self.ln_remerge(self.log_hold_lifetime(N), k, m)

    def fit(self, Ns):
        lnT = np.array([self.log_hold_lifetime(int(N)) for N in Ns])
        Ns = np.asarray(Ns, float)
        c, a = np.polyfit(Ns, lnT, 1)
        resid = lnT - (c * Ns + a)
        return float(c), float(a), float(abs(resid).max())

    def predicted_crossover(self, c, a, k):
        m = self.m_of_k(k)
        if k == m:
            return float("inf")
        return ((m - 1) * (a - math.log(self.tau))
                - math.log(math.comb(k, m))) / (c * (k - m))

    def grid(self, lo=0.6, hi=12.0, n=14, c=None):
        """Volumes spaced so that ln T runs over a fixed span, normalised across elements."""
        if c is None:
            c, _, _ = self.fit([max(10, int(round(x / 0.2))) for x in (2, 4, 6)])
        out = []
        for x in np.linspace(lo, hi, n):
            N = int(round(x / c))
            N = max(N, 8)
            if not out or N > out[-1]:
                out.append(N)
        return out


def _p0(els, report):
    """RULE 13: the reflecting box is an approximation and its width is a second axis.

    Not pre-registered -- it was run because rule 13 says to, after the predictions were written
    and while the suite was running. It changed the section: §111 and the cascade code inherited
    cap_mult = 1.25, which is not converged.
    """
    print("P0 -- box-width convergence (rule 13), NOT pre-registered")
    caps = (1.25, 1.5, 2.0, 3.0, 4.0)
    print(f"{'element':>9}{'N':>5}" + "".join(f"{'cap=' + str(c):>13}" for c in caps)
          + f"{'1.25 vs 3.0':>14}")
    rows, worst = [], 0.0
    for el in els:
        for N in (20, 60):
            vs = [Element(el.name, el.r1, el.r2, el.r3, cap_mult=c).log_hold_lifetime(N)
                  for c in caps]
            d = abs(vs[0] - vs[-1])
            worst = max(worst, d)
            rows.append({"element": el.name, "N": N, "caps": list(caps), "lnT": vs,
                         "err_at_1_25": d})
            print(f"{el.name:>9}{N:>5}" + "".join(f"{v:>13.7f}" for v in vs) + f"{d:>14.2e}")
    conv = max(abs(r["lnT"][-2] - r["lnT"][-1]) for r in rows)
    # >= not >: the tail of the sequence is exactly zero, and 0 > 0 is False
    shrink = all(all(abs(r["lnT"][i + 1] - r["lnT"][i]) >= abs(r["lnT"][i + 2] - r["lnT"][i + 1])
                     for i in range(len(caps) - 2)) for r in rows)
    print(f"   worst error at the INHERITED cap_mult = 1.25: {worst:.3f} in ln T")
    print(f"   successive differences shrink on every element: {shrink}")
    print(f"   worst difference between cap_mult 3.0 and 4.0: {conv:.2e}  -> converged")
    print(f"   §112 therefore runs at cap_mult = {CAP_MULT}")

    # the splitting probability does NOT depend on the box -- both rails sit inside any of them
    d = max(abs(Element(el.name, el.r1, el.r2, el.r3, cap_mult=1.25).commit_high(k, j, 30)
                - Element(el.name, el.r1, el.r2, el.r3, cap_mult=3.0).commit_high(k, j, 30))
            for el in els for k in KS for j in range(k + 1))
    print(f"   commit_high is box-INDEPENDENT: worst difference over 90 cells = {d:.1e}")
    print("   -> P2 and P2b below are unaffected by any of this; only the MFPT moves")
    report["p0"] = {"rows": rows, "worst_at_1_25": worst, "converged_2_vs_3": conv,
                    "commit_high_box_dependence": d, "cap_mult_used": CAP_MULT}


def _p1(report):
    """The instrument, against §111's dense solve and past its ceiling."""
    from experiments.depth_compounding import R2 as S_R2, R3 as S_R3
    from experiments.is_the_vote_a_majority import (
        MFPT_CEILING, VOLUMES, hold_lifetime as dense, m_of_k as s_m, merge_fraction,
    )
    # cap_mult=1.25 here on purpose: this compares two SOLVERS, so the box must be held at
    # §111's value or the gate would fail for a reason that has nothing to do with the solver.
    sch = Element("Schlogl-111", 0.15, S_R2, S_R3, cap_mult=1.25)
    print("P1 -- the instrument, on §111's own rails (box held at §111's 1.25 to compare solvers)")
    worst = 0.0
    rows = []
    for N in (10, 14, 20, 30, 42, 60):
        lg = sch.log_hold_lifetime(N)
        dn = math.log(dense(N))
        worst = max(worst, abs(lg - dn))
        rows.append({"N": N, "log_domain": lg, "dense": dn})
        print(f"   N = {N:>4}   log-domain {lg:>11.7f}   dense {dn:>11.7f}   "
              f"diff {lg - dn:+.2e}")
    print(f"   worst |diff| for N <= 60: {worst:.2e}   "
          f"-> {'HOLDS' if worst < 1e-9 else 'FAILS'} (P1's gate)")
    print(f"   past §111's cap of k*N <= {MFPT_CEILING}, where the dense solve went negative:")
    past = []
    for N in (200, 400, 700, 1000):
        lg = sch.log_hold_lifetime(N)
        past.append({"N": N, "lnT": lg})
        print(f"   N = {N:>4}   ln T = {lg:>12.5f}   finite and positive")

    # §111's headline, recomputed. The original stands beside it (rule 7).
    k = 3
    vols = [N for N in VOLUMES if k * N <= MFPT_CEILING]
    print(f"\n   §111's headline recomputed, under its OWN box and under a converged one:")
    print(f"     merge fraction  {merge_fraction():.5f}   m(3) = {s_m(3)} = k (unanimity)")
    head = []
    for cm in (1.25, 2.0):
        e = Element("S", 0.15, S_R2, S_R3, cap_mult=cm)
        c, a, _ = e.fit(VOLUMES)
        pred = 2.0 * (math.log(e.tau) - a)
        meas = [e.log_hold_lifetime(k * N)
                - e.ln_remerge(e.log_hold_lifetime(N), k, s_m(k)) for N in vols]
        head.append({"cap_mult": cm, "pred": pred, "meas": list(map(float, meas)),
                     "first": float(meas[0]), "mean": float(np.mean(meas)),
                     "resid_first_pct": abs(pred - meas[0]) / abs(pred) * 100})
        tag = "  <- §111's own box" if cm == 1.25 else "  <- converged"
        print(f"     cap={cm:<4} predicted {pred:+.4f}   first cell {meas[0]:+.4f}   "
              f"mean {np.mean(meas):+.4f}   resid(first) "
              f"{abs(pred - meas[0]) / abs(pred) * 100:.2f}%{tag}")
    print(f"     §111 PUBLISHED -3.4131 predicted / -3.4236 measured, 0.31% -- reproduced")
    print(f"     exactly under its own box. Under a converged box the same comparison reads")
    print(f"     {head[1]['pred']:+.4f} / {head[1]['first']:+.4f}, {head[1]['resid_first_pct']:.2f}%."
          f" §111's numbers stand (rule 7); its PRECISION was")
    print(f"     partly the unconverged wall.")
    pred, meas = head[0]["pred"], head[0]["meas"]
    report["p1"] = {"headline_both_boxes": head,
                    "agreement": rows, "worst_diff": worst, "past_ceiling": past,
                    "headline_pred": pred, "headline_meas": list(map(float, meas)),
                    "headline_mean": float(np.mean(meas)),
                    "published_pred": -3.4131, "published_meas": -3.4236}
    return worst < 1e-9


def _p2(els, report):
    """P2 and P2b together: the step is never scored at a single Omega (rule 20).

    The claim under test is that a merged tank commits to whichever basin its MEAN lies in.
    That is a statement about the limit, so every cell is swept in Omega and scored by whether
    P(high) moves toward the sign of its own margin. m = floor(f k) + 1 then follows by
    arithmetic, and the step is read off at the LARGEST volume, where the cells have saturated.
    """
    print("\nP2 -- where does the merge actually flip? (splitting probability, no settle time)")
    print("   m_f = floor(f k) + 1 from the rails;  m_maj = ceil((k+1)/2) is what AM would give")
    print("   scored at the largest Omega of each element's own sweep, never at a single cell")
    rows, cells = [], []
    for el in els:
        c, a, _ = el.fit(el.grid(2.0, 8.0, 4))
        N0 = max(10, min(60, int(round(6.0 / c))))
        oms = [N0, 2 * N0, 4 * N0, 8 * N0, 16 * N0]
        print(f"\n   {el.name}: r2 = {el.r2:.4f}  f = {el.f:.5f}  "
              f"c = {c:.5f}  a = {a:+.4f}  tau = {el.tau:.4f}  (Omega {oms[0]}-{oms[-1]})")
        for k in KS:
            ps = {om: [el.commit_high(k, j, om) for j in range(k + 1)] for om in oms}
            mg = [el.x_merged(k, j) - el.r2 for j in range(k + 1)]
            top = ps[oms[-1]]
            step = next((j for j, q in enumerate(top) if q < 0.5), None)
            m_f, m_maj = el.m_of_k(k), (k + 1) // 2
            ok = step == m_f
            # per-cell limit check -- every unsaturated cell must move toward its margin's sign
            for j in range(k + 1):
                series = [ps[om][j] for om in oms]
                if min(series) > 0.98 or max(series) < 0.02:
                    continue                      # already saturated, nothing to converge
                want_hi = mg[j] > 0
                moved = (series[-1] > series[0]) if want_hi else (series[-1] < series[0])
                mono = (series == sorted(series)) if want_hi else \
                       (series == sorted(series, reverse=True))
                cells.append({"element": el.name, "k": k, "j": j, "margin": mg[j],
                              "omegas": oms, "series": series, "want_high": want_hi,
                              "moved": bool(moved), "monotone": bool(mono)})
            rows.append({"element": el.name, "f": el.f, "k": k, "omegas": oms, "c": c, "a": a,
                         "p_high": {str(o): ps[o] for o in oms}, "margins": mg, "step": step,
                         "m_from_f": m_f, "m_majority": m_maj, "agrees": ok})
            print(f"     k={k}  m_f={m_f} m_maj={m_maj}  step={step}  "
                  f"{'OK' if ok else 'MISMATCH':>8}   P(high|j) at Omega={oms[-1]}: "
                  + " ".join(f"{q:.3f}" for q in top))
    agree = sum(r["agrees"] for r in rows)
    maj = sum(r["step"] == r["m_majority"] for r in rows)
    diff = [r for r in rows if r["m_from_f"] != r["m_majority"]]
    d_agree = sum(r["agrees"] for r in diff)
    print(f"\n   step = m from f      in {agree}/{len(rows)} cells")
    print(f"   step = majority      in {maj}/{len(rows)} cells")
    print(f"   of the {len(diff)} rows where the two DISAGREE, m from f is right in {d_agree}")
    print(f"   -> P2 {'HOLDS' if agree == len(rows) else 'FAILS'}")

    print("\nP2b -- T16-c: does every unsaturated cell converge toward its own margin's sign?")
    bad = [c for c in cells if not c["moved"]]
    for cl in cells:
        tag = "" if cl["moved"] else "   <-- AWAY FROM ITS MARGIN"
        print(f"   {cl['element']} k={cl['k']} j={cl['j']}  margin {cl['margin']:+.5f} "
              f"-> {'1' if cl['want_high'] else '0'}:  "
              + " ".join(f"{q:.4f}" for q in cl["series"]) + tag)
    print(f"   -> P2b {'HOLDS' if not bad else 'FAILS'}: "
          f"{len(cells) - len(bad)}/{len(cells)} cells move toward their margin's sign; "
          f"{sum(c['monotone'] for c in cells)}/{len(cells)} monotonically")
    report["p2"] = rows
    report["p2b"] = cells


def _p345(els, report):
    print("\nP3/P4/P5 -- the three regimes, each measured within one element")
    rows = []
    for el in els:
        gr = el.grid(0.6, 10.0, 14)
        lo_w, hi_w = gr[: len(gr) // 2], gr[len(gr) // 2:]
        # ln T = cN + a is the ASYMPTOTIC (WKB) form, so the primary fit is the large-volume
        # window. The whole-grid and small-volume fits are reported beside it as the envelope
        # (rule 15 -- every candidate extrapolation, not only the flattering one), and the
        # residual of each fit is printed so the reader can see which window is even linear.
        c, a, res = el.fit(hi_w)
        c_all, a_all, res_all = el.fit(gr)
        c_lo, a_lo, res_lo = el.fit(lo_w)
        print(f"\n   {el.name}: f = {el.f:.5f}   volumes {gr[0]}-{gr[-1]}")
        print(f"      PRIMARY fit, window {hi_w[0]}-{hi_w[-1]} (asymptotic): "
              f"ln T = {c:.6f} N {a:+.4f}  max|resid| = {res:.4f}")
        print(f"      alternates: whole grid {gr[0]}-{gr[-1]} c = {c_all:.6f} a = {a_all:+.4f} "
              f"resid {res_all:.4f}   |   small {lo_w[0]}-{lo_w[-1]} c = {c_lo:.6f} "
              f"a = {a_lo:+.4f} resid {res_lo:.4f}")
        for k in KS:
            m = el.m_of_k(k)
            bnd_un, bnd_pool = 1.0 - 1.0 / k, 1.0 / k
            if m == k:
                regime = f"UNANIMITY (f >= 1-1/k = {bnd_un:.4f})"
            elif m == 1:
                regime = f"ONE FLIPS IT (f < 1/k = {bnd_pool:.4f})"
            else:
                regime = "finite crossover"
            lr = [el.ln_ratio(N, k) for N in gr]
            sl = float(np.polyfit(np.array(gr, float), np.array(lr), 1)[0])
            pred_slope = c * (k - m)
            xo = el.predicted_crossover(c, a, k)
            xo_lo = el.predicted_crossover(c_lo, a_lo, k)
            xo_hi = el.predicted_crossover(c_all, a_all, k)
            sc = [i for i in range(len(gr) - 1) if lr[i] * lr[i + 1] < 0]
            meas = None
            if sc:
                # Bisect on INTEGER volumes down to adjacent sites, then interpolate across
                # that single step. Interpolating across the coarse grid instead would span
                # 15+ volumes and is the error rule 19 was written for.
                lo, hi = gr[sc[0]], gr[sc[0] + 1]
                f_lo = lr[sc[0]]
                while hi - lo > 1:
                    mid = (lo + hi) // 2
                    v = el.ln_ratio(mid, k)
                    if v * f_lo < 0:
                        hi = mid
                    else:
                        lo, f_lo = mid, v
                f_hi = el.ln_ratio(hi, k)
                meas = lo + (-f_lo / (f_hi - f_lo)) * (hi - lo)
            flat_pred = (k - 1) * (math.log(el.tau) - a) if m == k else None
            rows.append({"element": el.name, "f": el.f, "k": k, "m": m, "regime": regime,
                         "volumes": gr, "ln_ratio": lr, "slope": sl,
                         "pred_slope": pred_slope, "c": c, "a": a,
                         "resid": res, "c_lo": c_lo, "a_lo": a_lo,
                         "c_all": c_all, "a_all": a_all,
                         "pred_xover": xo, "pred_xover_lo": xo_lo, "pred_xover_hi": xo_hi,
                         "meas_xover": meas, "flat_pred": flat_pred})
            print(f"      k={k} m={m}  {regime}")
            print(f"         ln(L_hold/L_remerge) = "
                  + ", ".join(f"{v:+.2f}" for v in lr))
            print(f"         slope {sl:+.6f} vs predicted c(k-m) = {pred_slope:+.6f}", end="")
            if m == k:
                da = (np.mean(lr) - flat_pred) / (k - 1)
                rows[-1]["implied_delta_a"] = float(da)
                print(f"   |   FLAT value {np.mean(lr):+.4f} vs (k-1)(ln tau - a) = "
                      f"{flat_pred:+.4f}   resid {abs(np.mean(lr) - flat_pred):.4f}"
                      f"   -> implied error in `a` alone: {da:+.4f}")
            elif m == 1:
                print(f"   |   all positive: {all(v > 0 for v in lr)}, "
                      f"rising: {lr == sorted(lr)}  -> pooling wins everywhere")
            else:
                ms = "none in range" if meas is None else f"{meas:.2f}"
                sp = (f"{min(xo_lo, xo_hi):.2f}-{max(xo_lo, xo_hi):.2f}"
                      if np.isfinite(xo_lo) and np.isfinite(xo_hi) else "n/a")
                print(f"\n         Omega_x predicted {xo:.2f} [window spread {sp}]   "
                      f"measured {ms}", end="")
                if meas:
                    depth = math.exp(el.log_hold_lifetime(int(round(meas)))) / el.tau
                    rows[-1]["xover_depth_T_over_tau"] = depth
                    resolved = not (min(xo_lo, xo_hi) <= meas <= max(xo_lo, xo_hi))
                    rows[-1]["resolved_vs_window"] = bool(resolved)
                    print(f"   ratio {xo / meas:.4f}"
                          f"   [{'resolved' if resolved else 'UNRESOLVED: inside window spread'}]"
                          f"   T/tau there = {depth:.1f}")
                else:
                    print()
    report["p345"] = rows


def _p7(els, report):
    """POST-HOC, NOT PRE-REGISTERED (rule 2), and it changes what P5b's residuals mean.

    §34's law needs T at BOTH Omega and k*Omega, but c and a are fitted over Omega alone. For
    E-084 at k = 7 the asymptotic window is 19-29 while the pooled tank spans 56-203, entirely
    outside it, and ln T curves upward -- so ln T(k*Omega) is under-predicted and Omega_x comes
    out high. That is a suspect (rule 17), and its kill test is to refit over a window that
    COVERS Omega..k*Omega and see whether the outlier closes.

    A first account -- that k - m = 1 makes the denominator sensitive -- was killed on sight:
    E-062 at k = 5 also has k - m = 1 and agrees to 0.9%.
    """
    print("\nP7 -- POST-HOC (not pre-registered): is the residual set by the fit window?")
    print("   refit c, a over a window spanning Omega..k*Omega and re-predict")
    print(f"{'cell':>12}{'measured':>10}{'asymptotic':>12}{'covering':>10}"
          f"{'r_asym':>9}{'r_cover':>9}{'shift':>9}")
    rows = []
    for r in report["p345"]:
        if r.get("meas_xover") is None:
            continue
        el = next(e for e in els if e.name == r["element"])
        k, meas = r["k"], r["meas_xover"]
        lo, hi = int(round(meas)), int(round(k * meas))
        cov = sorted({int(round(lo * (hi / lo) ** (i / 7))) for i in range(8)})
        c_c, a_c, res_c = el.fit(cov)
        p_c = el.predicted_crossover(c_c, a_c, k)
        ra, rc = (r["pred_xover"] / meas - 1) * 100, (p_c / meas - 1) * 100
        rows.append({"element": el.name, "k": k, "measured": meas,
                     "asym": r["pred_xover"], "cover": p_c, "cover_window": cov,
                     "r_asym": ra, "r_cover": rc, "shift": rc - ra})
        print(f"{el.name + ' k=' + str(k):>12}{meas:>10.2f}{r['pred_xover']:>12.2f}"
              f"{p_c:>10.2f}{ra:>8.1f}%{rc:>8.1f}%{rc - ra:>8.1f}")
    sp = max(abs(x["shift"]) for x in rows)
    m_a = float(np.median([abs(x["r_asym"]) for x in rows]))
    m_c = float(np.median([abs(x["r_cover"]) for x in rows]))
    best = min(rows, key=lambda x: abs(x["r_asym"]))
    print(f"   window shifts the prediction by up to {sp:.1f} percentage points")
    print(f"   median |residual|: asymptotic fit {m_a:.1f}%, covering fit {m_c:.1f}%")
    print(f"   the best asymptotic cell ({best['element']} k={best['k']}) reads "
          f"{best['r_asym']:+.1f}% there and {best['r_cover']:+.1f}% under the covering window")
    print(f"   -> the absolute agreement is WINDOW-LIMITED at the ~{max(m_a, m_c):.0f}% level;"
          f" no cell is resolved below it,\n      and quoting any single cell's residual as the"
          f" result would be over-claiming.")
    report["p7_summary"] = {"max_shift_pp": sp, "median_asym": m_a, "median_cover": m_c}
    report["p7"] = rows


def _p6(report):
    print("\nP6 -- ARITHMETIC, NOT A MEASUREMENT (rule 21): how symmetric must you be for")
    print("   §34's k-independence?  (m-1)/(k-m) = 1 iff 1/2 - 1/(2k) <= f < 1/2 + 1/(2k)")
    fs = [0.40, 0.45, 0.48, 0.50, 0.52, 0.55, 0.60]
    print("      f  " + "".join(f"{'k=%d' % k:>10}" for k in (3, 5, 7, 9, 11)))
    rows = []
    for f in fs:
        cells, r = [], {"f": f}
        for k in (3, 5, 7, 9, 11):
            m = int(math.floor(f * k)) + 1
            ratio = "inf" if m == k else f"{(m - 1) / (k - m):.3f}"
            cells.append(f"{ratio:>10}")
            r[f"k{k}"] = ratio
        rows.append(r)
        band = f"  <- inside for k <= {int(math.floor(1 / (2 * abs(f - 0.5)))) if f != 0.5 else 'all k'}"
        print(f"   {f:.2f}  " + "".join(cells) + band)
    print("   (a ratio of 1.000 is §34's k-independent case; the band shrinks like 1/k, so the")
    print("    headline holds for every k only at f = 1/2 exactly -- the symmetric point)")
    report["p6"] = rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/where_pooling_stops_winning.json"))
    args = ap.parse_args()
    report = {}

    els = [Element(n, R1, r2, R3) for n, r2 in ELEMENTS]
    print("the constructed elements: r1 and r3 fixed, the saddle moved")
    print(f"{'name':>8}{'r2':>9}{'f':>10}{'m(3)':>7}{'m(5)':>7}{'m(7)':>7}"
          f"{'maj(3)':>8}{'maj(5)':>8}{'maj(7)':>8}")
    for el in els:
        rts = el.roots()
        assert len(rts) == 3 and abs(rts[1] - el.r2) < 1e-9, f"{el.name}: rails not recovered"
        print(f"{el.name:>8}{el.r2:>9.4f}{el.f:>10.5f}"
              + "".join(f"{el.m_of_k(k):>7}" for k in KS)
              + "".join(f"{(k + 1) // 2:>8}" for k in KS))

    _p0(els, report)
    ok = _p1(report)
    if not ok:
        print("\nP1 FAILED -- the instrument disagrees with §111 where §111 is sound. Stopping.")
        args.out.write_text(json.dumps(report, indent=2, default=float))
        return
    _p2(els, report)
    _p345(els, report)
    _p7(els, report)
    _p6(report)

    args.out.write_text(json.dumps(report, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
