"""§122 -- the near-miss floor from single elements, attempted as a formula; the limit is universal.

THE GOAL, set before any of this ran: predict §119's floor C = P(stage 2 wrong | stage 1 right)/(k1 t)
in every cell to within 1.25x from a formula built only from each element's own rates -- no
simulation, no joint chain -- with every prediction written before its comparison. Judged near-
impossible: §120's best single-element SIMULATION only reached 0.66-1.06, and nothing explained
the Omega-dependence of p_d.

C = phi * p_d + p0/(k1 t). The predictions below were first written as docstrings of scratch scripts,
each before that script ran, and are copied here verbatim in the order they were made.

  F1 "phi = w_d/(k1 t) equals the birth-death splitting probability of returning to the saddle
     from one step below it before reaching the low rail, to within a few %."
     -> FAILED: that probability is ~0.9; phi is ~0.4. A recrossed stage 1 sits at its saddle and
        re-dips; what counts is where it ends.
  F2 "phi equals the committor from sad-1 to the HIGH rail before the LOW rail, to within a few %."
     -> HELD: 0.987-1.003 in four cells; 0.954 in Z:40, the one cell whose window saturates
        (k1 t = 0.43).
  F3 "Asymptotic formula p_asym(Om) = erf(w e^{-ls T*} sqrt(ls Om / (2 lam_s))), T* = time for
     deterministic stage 2, from its high fixed point at input mu1, with input frozen at s1, to
     reach the middle root (separatrix) of its cubic at input mu1. Predict: p_asym and the
     reduction's p_d_ode converge as Om grows (ratio -> 1), for X and b at least; for Z the formula
     is far above at small Om (it omits the deep-excursion term, which is the falling part)."
     -> FAILED both halves: X drifted AWAY (0.95 -> 0.60), and Z's formula was far BELOW, not above.
        Freezing the input at s1 and linearising to the zone's far edge is wrong for wide zones.
  F4 "p_d = erf(U* sqrt(ls Om/(2 lam(s1)))), U* = the largest offset above the saddle from which
     the deterministic two-element orbit (stage 1 its own rate equation from s1+U, stage 2 its own
     driven by it, from its high fixed point at mu1) ends with stage 2 below its separatrix.
     Predict: reduction/formula -> 1 with Om; within 15% at Om=400 for X and b."
     -> HELD: X 1.134, b 1.138 at Om = 400; ratios fall monotonically to 1 by Om = 1600-3200 in
        all three pairs (deterministic stage 2), and with stage 2's noise too.

THE FORMULA. Stage 1 near its saddle is an unstable linear diffusion, u' = ls u + noise with variance
rate 2 lam(s1)/Om, ls = k(s1 - r1)(r3 - s1) (Lean: Schlogl.drift'_rails) and lam(s1) its birth rate
there. Its escape is the deterministic orbit started at s1 + U with U ~ N(0, lam(s1)/(ls Om)). Stage 2
falls iff U < U*. Hence p_d -> erf(U* sqrt(ls Om / (2 lam(s1)))), and U* > 0 iff the zone exists.
The asymptotic regime needs stage 1's wander sigma = sqrt(lam(s1)/(ls Om)) below the zone width w:
Omega_c = lam(s1)/(ls w^2).

THE CONSEQUENCE, stated as a suspect until something independent of these reductions confirms it
(rule 17): phi -> 1/2 and p_d -> 1 for every positive zone, so phi * p_d -> 1/2 -- a universal
floor, independent of zone width and of the downstream barrier, approached logarithmically.

PREDICTIONS FOR THIS SCRIPT'S OWN RUN (it regenerates everything above; written before it runs):
  G1  U* is stable to the integration horizon (80 -> 200) to 1e-6 relative.
  G2  The regenerated tables reproduce the scratch runs within Monte Carlo error (same seeds:
      identical to the printed digits).
  G3  (the goal) C_formula = committor * erf(...) is within 1.25x of §119's exact C in every cell
      with a zone. Expected to FAIL: the formula is 3-30x low at Omega = 30-55.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib

import numpy as np
from scipy.integrate import solve_ivp

import experiments.frozen_downstream as fd
from experiments.weak_upstream import ARMS, element, edge

PAIRS = {"X": ("E-070", "E-084"), "b": ("E-070", "F-040"), "Z": ("E-062", "E-084")}
BIG = (30, 55, 100, 200, 400, 800, 1600, 3200)


def committor(A, om, start):
    """P(reach the high rail before the low rail | start), birth-death closed form, log domain."""
    lo, hi = int(round(A.r1 * om)), int(round(A.r3 * om))
    n = np.arange(0, hi + 1, dtype=float)
    lam, mu = A.rates(n, np.zeros(len(n)), om, True)
    lr = np.zeros(hi + 1)
    for j in range(lo + 1, hi + 1):
        lr[j] = lr[j - 1] + np.log(mu[j]) - np.log(lam[j])
    t = lr[lo:hi]
    w = np.exp(t - t.max())
    return float(w[: start - lo].sum() / w.sum())


def _f1(A, x):
    k1a, k1r, k2b, k2r = A.c
    return k1a * x * x + k2b - k1r * x ** 3 - k2r * x


def _mid_root(B, x_in):
    k1a, k1r, k2b, k2r = B.c
    r = sorted(z.real for z in np.roots([-k1r, B.hill(x_in) * k1a, -k2r, k2b]) if abs(z.imag) < 1e-9)
    return r[1]


def u_star(A, B, horizon=80.0):
    mu1 = A.r3
    x2h, sep = fd.high_fixed_point(B, mu1), _mid_root(B, mu1)

    def falls(U):
        rhs = lambda t, y: [_f1(A, y[0]), fd.ode_rate(B, y[1], y[0])]
        s = solve_ivp(rhs, (0, horizon), [A.r2 + U, x2h], rtol=1e-10, atol=1e-12, method="LSODA")
        return s.y[1, -1] < sep

    lo, hi = 1e-14, 0.5
    if not falls(lo):
        return 0.0
    for _ in range(80):
        m = math.sqrt(lo * hi)
        lo, hi = (m, hi) if falls(m) else (lo, m)
    return lo


def scales(A):
    k1a, k1r, k2b, k2r = A.c
    s1 = A.r2
    return k1r * (s1 - A.r1) * (A.r3 - s1), k1a * s1 ** 2 + k2b


def p_formula(A, Us, om):
    ls, lam_s = scales(A)
    return math.erf(Us * math.sqrt(ls * om / (2 * lam_s)))


EXACT100 = """
G4, written before running: exact joint CME, arm X (E-070 -> E-084), Omega = 100, t = 16.
"Predict p_d in [0.33, 0.45] (rising from 0.289/0.289/0.308 at 30/40/55), and C = P/(k1 t) above
0.136 (its Omega = 55 value). Flat or falling p_d -> the asymptotic picture does not reach the exact
chain at accessible Omega."
  -> HELD: p_d 0.385, C 0.183; and phi 0.4311 against the committor's 0.4299 (a fresh F2 cell).
"""


def exact100(out):
    from experiments.weak_upstream import run
    A, B = element("E-070"), element("E-084")
    r = run(A, B, 100)
    s = r["series"][-1]
    rep = json.loads(out.read_text()) if out.exists() else {}
    rep["exact100"] = {"p_d": s["p_d"], "w_d": s["w_d"], "p0": s["p0"], "k1": r["k1"],
                       "C": s["p_free"] / (r["k1"] * s["t"]), "phi": s["w_d"] / (r["k1"] * s["t"]),
                       "committor": committor(A, 100, int(np.ceil(A.r2 * 100)) - 1),
                       "formula_p_d": p_formula(A, u_star(A, B), 100)}
    out.write_text(json.dumps(rep, indent=2, default=float))
    print(rep["exact100"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("results/asymptotic_floor.json"))
    ap.add_argument("--exact100", action="store_true",
                    help="only the exact Omega = 100 check on arm X (~25 min); merged into --out")
    args = ap.parse_args()
    if args.exact100:
        exact100(args.out)
        return
    exact = json.loads(pathlib.Path("results/weak_upstream.json").read_text())
    rep = {"phi": {}, "u_star": {}, "conv": {}, "goal": {}}

    for key in ("X:30", "X:40", "X:55", "Z:40", "Z:55"):
        A = element(ARMS[key[0]][0]); om = int(key[2:])
        s = exact[key]["series"][-1]
        phi = s["w_d"] / (exact[key]["k1"] * s["t"])
        c = committor(A, om, int(np.ceil(A.r2 * om)) - 1)
        rep["phi"][key] = {"phi": phi, "committor": c, "ratio": phi / c}
        print(f"F2 {key}: phi {phi:.4f} committor {c:.4f} ratio {phi / c:.4f}", flush=True)
    for name in ("E-070", "E-062"):
        A = element(name)
        rep["phi"][f"{name}/limit"] = {str(om): committor(A, om, int(np.ceil(A.r2 * om)) - 1)
                                       for om in (30, 100, 400, 1600, 6400)}

    for p, (a, b) in PAIRS.items():
        A, B = element(a), element(b)
        u80, u200 = u_star(A, B, 80.0), u_star(A, B, 200.0)
        ls, lam_s = scales(A)
        w = edge(B) - A.r2
        rep["u_star"][p] = {"u_star": u80, "u_star_h200": u200, "w": w, "lambda_s": ls,
                            "lam_s1": lam_s, "omega_c": lam_s / (ls * w * w)}
        print(f"G1 {p}: U* {u80:.6e} (h=200: {u200:.6e}, rel {abs(u200 / u80 - 1):.1e}); "
              f"Omega_c {lam_s / (ls * w * w):.1f}", flush=True)
        fd.PATHS = 20000
        rows = []
        for om in BIG:
            if om >= 800:
                fd.PATHS = 40000
            r = fd.memory_reduction(A, B, om)
            pf = p_formula(A, u80, om)
            se = math.sqrt(r["p_d_ode"] * (1 - r["p_d_ode"]) / r["back"])
            rows.append({"omega": om, "p_d_ode": r["p_d_ode"], "se": se, "formula": pf,
                         "ratio": r["p_d_ode"] / pf})
            print(f"   {p} Om {om}: reduction {r['p_d_ode']:.4f}±{se:.4f} formula {pf:.4f} "
                  f"ratio {r['p_d_ode'] / pf:.3f}", flush=True)
        rep["conv"][p] = rows
        args.out.write_text(json.dumps(rep, indent=2, default=float))

    for key, v in exact.items():
        if ":" not in key or key.split(":")[0] not in "SXbZacd" or key.startswith("conv"):
            continue
        arm, om = key.split(":"); om = int(om)
        A, B = element(ARMS[arm][0]), element(ARMS[arm][1])
        s = v["series"][-1]
        C = s["p_free"] / (v["k1"] * s["t"])
        Cf = committor(A, om, int(np.ceil(A.r2 * om)) - 1) * p_formula(A, u_star(A, B), om)
        rep["goal"][key] = {"C_exact": C, "C_formula": Cf, "ratio": C / Cf if Cf > 0 else float("inf")}
        print(f"G3 {key}: C exact {C:.4g} formula {Cf:.4g} ratio {C / Cf if Cf > 0 else float('inf'):.3g}",
              flush=True)
    args.out.write_text(json.dumps(rep, indent=2, default=float))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
