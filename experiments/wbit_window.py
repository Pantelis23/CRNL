"""§131 -- T-WBIT-a‴: is N's slow approach to Gamma the cost of hovering in a narrow window above the saddle?

§127-§129: in N the cheapest never-dipped fall parks cell 1 on its saddle while cell 2 escapes at drive b2, so
Gamma_N = -S2(b2). The exact chain's slope approaches it slowly: Gamma - slope ~ alpha' sigma^2 with an effective
alpha' of 4.6-5.7 (§129's intervals).

FIRST CANDIDATE, KILLED BEFORE ANY RUN (rule 16). "Cell 1 hovers at a repelling saddle at a cost ln(lambda1) per
step, traded against how long cell 2's escape takes." That gives alpha' = ln(lambda1) / (2 ln(1/mu2)), with lambda1
= 1.25 and mu2 = 0.826 (the slower of cell 2's well and saddle eigenvalues at b2), so alpha' = 0.58. The measured
value is ~5, so this model is 9x too small.

SUSPECT UNDER TEST: A NARROW WINDOW. Cell 2's escape cost rises with cell 1's height u above its saddle,
dS2 = S2'(b2) 2c u (S2' ~ 1.7). For that to cost O(1) in the exponent, cell 1 must stay within w ~ sigma^2 / (2 c S2')
of the saddle while cell 2 crosses. Cell 1's own one-step noise there is 0.354 sigma >> w. So each step of hovering
costs about ln(0.354 sigma 2c S2' / sigma^2), and k steps cost k ln(c) + (terms in sigma). With alpha' ~ k/2, k ~ 10.

THE DISCRIMINATING TEST. Hold b2 at N's value (-0.10980), so Gamma_N = -S2(b2) is identical, and vary only c,
adjusting the healthy drive h = b2 + c (2 HI1 - 1) to compensate. c = 0.20, 0.23, 0.26, 0.30, 0.34 (h 0.062-0.182;
edge 0.26-0.36, all below 1/2; cell 2 bistable at h; cell 2's own route costs S2(h) = 0.487-0.860, all above
W* = 0.4355, so parking stays the cheapest route. The margin at c = 0.20 is thin: 0.05). Exact R by §129's
positive solver.

PREDICTIONS, WRITTEN BEFORE RUNNING.
  N1  THE WINDOW. At sigma = 0.15 (M = 320), d ln R / d ln c, fitted over c = 0.23-0.34 (c = 0.20 excluded for its
      thin margin), lies in [-20, -5]. Plain large deviations, Gamma alone, predicts no leading-order dependence.
      The verdicts, chosen so that each has its own data:
        slope in [-20, -5]        -> the window suspect SURVIVES;
        |slope| < 1.5             -> the window suspect is KILLED: the prefactor does not care about c;
        anything else             -> neither: the dependence exists but not at the window's size, and the
                                     suspect is killed in magnitude (rule 19).
  N2  It is a power law. The local slopes between neighbouring c agree to within a factor 1.5 of each other.
  N3  It does not fade with noise. At sigma = 0.12 the fitted |d ln R / d ln c| is at least that at 0.15. In the
      window picture each step costs ln(c/sigma) and the number of steps does not shrink as sigma falls.
  N4  GRID. R at sigma = 0.12 changes by < 1% from M = 320 to 480 at c = 0.23 and 0.34.

ADDED AFTER THE FIRST RUN, BEFORE THE EXTENSION RAN. The first run's fit (-10.4 at sigma = 0.15) is CONFOUNDED.
Holding b2 fixed ties h to c, and at low c cell 2's own route is not negligible. Its frozen-input estimate (own
escape rate at drive h / dip rate) is 1.6 / 1.6 / 0.68 / 0.058 / 0.003 of R at c = 0.20 / 0.23 / 0.26 / 0.30 / 0.34
(sigma = 0.15). The action margin used to vouch for those cells (S2(h) - W* = 0.12 at c = 0.23) did not protect
them at this sigma. N1's criterion could therefore be met by the own route rather than the window (rule 19), and
its verdict is not read. The sigma = 0.12 own-route estimates came out negative: the 1-D direct solve lost them to
subtraction, as §128's GMRES did.
EXTENSION. Cells c = 0.30, 0.32, 0.34, 0.355, 0.365 (h up to 0.2034 < DC: cell 2 still bistable at the healthy
input). A cell is CLEAN if its own-route estimate, now from a positive-arithmetic 1-D quasi-stationary rate, is
below 0.1 R. Only clean cells are read; the rule is applied before any slope is looked at.
  N1′ On the clean cells, d ln R / d ln c lies in [-20, -5] at both sigma = 0.15 and 0.12.
  N2′ Its local slopes agree within a factor 1.5.
  N3′ |slope(0.12)| >= |slope(0.15)|.

OUTCOMES (appended after running; predictions above left as written).
  N1 NOT READ (confounded; see the addendum). As printed: -10.4 at 0.15, -8.8 at 0.12.
  N2 FAILED as printed (local -15.8 / -9.0 / -7.5, spread 2.1x). The steep end is the own route.
  N3 FAILED as printed (8.8 < 10.4). It is confounded the same way and is not read.
  N4 HELD: M 320 -> 480 moves R by 0.8% (c = 0.23) and 0.76% (c = 0.34).
  The addendum's confound, re-checked (rule 14). The first own-route estimate used `flip_rate`, which counts an exit
  at p < 1/2 rather than at cell 2's saddle, and so overstated it for h > 0. The positive rate with the right exit
  (validated against the direct solve to 5e-4 and 1e-6 where both are readable) gives own/R = 0.74 / 0.45 / 0.10 /
  0.003 at c = 0.20 / 0.23 / 0.26 / 0.30 (sigma = 0.15). The confound stands, smaller than first stated.
  N1′ HELD: on c = 0.30-0.365 (all clean, own/R <= 0.003), the slope is -7.06 (sigma = 0.15) and -7.69 (0.12).
  N2′ HELD: local slopes -7.60 / -7.30 / -6.62 / -5.62 (1.35x) and -7.82 / -7.98 / -7.51 / -6.57 (1.21x).
  N3′ HELD: 7.69 >= 7.06. FINDINGS §131.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from experiments.wbit_action import escape_action
from experiments.wbit_cascade import DC, HI1, PAIRS
from experiments.wbit_qsd import exits_pos

B2 = PAIRS["N"][1]
CS = (0.20, 0.23, 0.26, 0.30, 0.34)


def h_of(c):
    return B2 + c * (2 * HI1 - 1)


def fit_slope(cs, R):
    x, y = np.log(cs), np.log(R)
    return float(np.polyfit(x, y, 1)[0]), [float((y[i + 1] - y[i]) / (x[i + 1] - x[i])) for i in range(len(x) - 1)]


def rate_pos(M, sigma, drive, tol=1e-13):
    """One cell's per-step escape rate from its high well to p < its saddle, by positive arithmetic: the
    quasi-stationary exit probability per step (power iteration of the alive kernel, no subtraction)."""
    from experiments.wbit_cell import fixed_points, kernel
    P, x = kernel(M, 1.5, 0.5, 0.5, sigma, drive)
    sad = fixed_points(1.5, 0.5, 0.5, drive)[1]
    alive = x >= sad
    Pa = P[alive][:, alive].tocsr()
    exitp = np.asarray(P[alive][:, ~alive].sum(axis=1)).ravel()
    mu = np.zeros(alive.sum())
    mu[int(np.argmin(abs(x[alive] - fixed_points(1.5, 0.5, 0.5, drive)[-1])))] = 1.0
    prev = None
    for t in range(1, 200000):
        mu = Pa.T @ mu
        mu /= mu.sum()
        if t % 50 == 0:
            r = float(mu @ exitp)
            if prev is not None and abs(r / prev - 1) < tol:
                return r
            prev = r
    return r


def own_ratio(c, sigma, R, M=3200):
    return rate_pos(M, sigma, h_of(c)) / rate_pos(M, sigma, 0.0) / R


def extension():
    out = Path("results/wbit_window.json")
    rep = json.loads(out.read_text())
    ext = (0.30, 0.32, 0.34, 0.355, 0.365)
    for s in (0.15, 0.12):
        for c in tuple(c for c in CS if c not in ext) + ext:      # the first run's cells re-judged too
            key = f"{c}:{s}:320"
            if key not in rep:
                r = exits_pos(320, s, c, B2, h=h_of(c))
                r.update({"c": c, "h": h_of(c), "edge": HI1 - (DC + h_of(c)) / (2 * c)})
                rep[key] = r
            rep[key]["own_over_R"] = own_ratio(c, s, rep[key]["R"])
            print(key, {k: rep[key][k] for k in ("R", "h", "own_over_R")}, flush=True)
            out.write_text(json.dumps(rep, indent=2))
        clean = [c for c in ext if rep[f"{c}:{s}:320"]["own_over_R"] < 0.1]
        if len(clean) >= 2:
            slope, local = fit_slope(clean, [rep[f"{c}:{s}:320"]["R"] for c in clean])
            rep[f"ext_fit:{s}"] = {"clean": clean, "slope": slope, "local": local}
            print("ext fit", s, rep[f"ext_fit:{s}"], flush=True)
        out.write_text(json.dumps(rep, indent=2))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "ext":
        return extension()
    out = Path("results/wbit_window.json")
    rep = json.loads(out.read_text()) if out.exists() else {}
    cells = [(c, s, 320) for s in (0.15, 0.12) for c in CS] + [(c, 0.12, 480) for c in (0.23, 0.34)]
    if len(sys.argv) > 1:
        cells = [x for x in cells if x[2] == int(sys.argv[1])]
    for c, s, M in cells:
        key = f"{c}:{s}:{M}"
        if key in rep:
            continue
        r = exits_pos(M, s, c, B2, h=h_of(c))
        r.update({"c": c, "h": h_of(c), "edge": HI1 - (DC + h_of(c)) / (2 * c), "S2_h": escape_action(h_of(c))})
        rep[key] = r
        print(key, {k: r[k] for k in ("R", "h", "edge", "mass_check")}, flush=True)
        out.write_text(json.dumps(rep, indent=2))
    for s in (0.15, 0.12):
        cs = [c for c in CS if c >= 0.23 and f"{c}:{s}:320" in rep]
        if len(cs) == 4:
            slope, local = fit_slope(cs, [rep[f"{c}:{s}:320"]["R"] for c in cs])
            rep[f"fit:{s}"] = {"slope": slope, "local": local}
            print("fit", s, rep[f"fit:{s}"], flush=True)
    out.write_text(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main()
