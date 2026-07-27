"""Does the design rule's coefficient go to 1? -- Q4a-i.

FINDINGS 17 measured

    ln(e-/e+) at beta*  =  r(Omega) * ln(p/(1-p))       R^2 = 0.9999

with r rising 0.644 -> 0.879 across Omega = 100 -> 400. The derivation says
r -> 1. Four times in Omega could not distinguish the candidate approaches:
extrapolating under 1/Omega, Omega^-0.75, 1/sqrt(Omega) and 1/ln(Omega) gives
limits 0.947 / 1.004 / 1.120 / 1.685, and the two BEST-fitting forms overshoot.

The obstacle was cost, not principle. `biased_source` finds beta* by maximising
the exact mutual information, which spends ~40 CME solves per (Omega, p) cell, so
a prior sweep at one Omega is already 200 solves.

THE INVERSION THAT MAKES THE RANGE AFFORDABLE. beta* is defined by dI/dbeta = 0,
and that condition is a scalar equation in p:

    -p [H'(q) + H'(e+)] * e+ * dln(e+)/dbeta
  + (1-p)[H'(q) - H'(e-)] * e- * dln(e-)/dbeta   =  0

with q = p(1-e+) + (1-p)e- and H'(x) = log2((1-x)/x). So instead of asking "given
p, which beta is optimal" (an optimisation, many solves), ask "given beta, which p
makes it optimal" (a root find in p, no solves at all once the error curve is
known). ONE sweep of e±(beta) -- about ten solves -- yields the whole beta <-> p
map, hence the entire P2 line and its slope r. That is a 20x saving, and it is
what buys the Omega range.

Derivatives are taken on ln(e±), not e±, because the errors vary exponentially in
beta and differencing the raw values loses most of the precision.

Prediction, written before running: r keeps rising and the extrapolation
stabilises once the range is long enough to separate the ansaetze. If instead r
flattens below 1, the derivation is wrong and the coefficient is a real number to
be explained -- which would be the more interesting outcome.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
from scipy.optimize import brentq

from crnl.networks.am_asymmetric import beta_critical
from experiments.biased_source import errors

# Full (0,1), not (1/2,1). A tilt favouring X always matches a prior above 1/2,
# so the narrow domain worked for every real sweep -- but it put the root ON the
# bracket edge for the symmetric case and made beta < 0 unreachable, which is a
# domain bug hiding behind a sign convention rather than a real restriction.
P_LO, P_HI = 1e-9, 1.0 - 1e-9


def _hprime(x: float) -> float:
    """d/dx of the binary entropy, in bits."""
    x = min(max(x, 1e-300), 1.0 - 1e-16)
    return float(np.log2((1.0 - x) / x))


def matched_prior(e_plus: float, e_minus: float,
                  dln_plus: float, dln_minus: float) -> float | None:
    """The prior p for which this beta is the optimum. See the module docstring."""
    def f(p: float) -> float:
        q = p * (1.0 - e_plus) + (1.0 - p) * e_minus
        hq = _hprime(q)
        return (-p * (hq + _hprime(e_plus)) * e_plus * dln_plus
                + (1.0 - p) * (hq - _hprime(e_minus)) * e_minus * dln_minus)

    lo, hi = f(P_LO), f(P_HI)
    if not np.isfinite(lo) or not np.isfinite(hi) or lo * hi > 0:
        return None
    return float(brentq(f, P_LO, P_HI, xtol=1e-12))


def rule_line(gamma: float, omega: int, eps_frac: float,
              betas: np.ndarray) -> dict | None:
    """The whole P2 line at one Omega, from a single sweep of e+-(beta)."""
    t0 = time.time()
    curve = []
    for b in betas:
        r = errors(gamma, float(b), omega, eps_frac)
        if r is None or min(r) <= 0.0:
            return None
        curve.append(r)
    ep = np.array([c[0] for c in curve])
    em = np.array([c[1] for c in curve])
    lp, lm = np.log(ep), np.log(em)

    xs, ys, ps = [], [], []
    for i in range(1, len(betas) - 1):                 # central differences
        h = betas[i + 1] - betas[i - 1]
        dlp = (lp[i + 1] - lp[i - 1]) / h
        dlm = (lm[i + 1] - lm[i - 1]) / h
        p = matched_prior(ep[i], em[i], dlp, dlm)
        # Match FINDINGS 17's prior range. Beyond p ~ 0.96 the relation bends and
        # including those points drags the slope; the first run of this script
        # reached p = 0.987 and read 0.738 where the direct optimisation, over
        # p <= 0.95, read 0.763.
        if p is None or not (0.55 < p < 0.96):
            continue
        xs.append(np.log(p / (1 - p)))
        ys.append(float(lm[i] - lp[i]))
        ps.append(p)
    if len(xs) < 4:
        return None
    xs, ys = np.array(xs), np.array(ys)
    slope, intercept = np.polyfit(xs, ys, 1)
    resid = ys - (slope * xs + intercept)
    return {
        "omega": int(omega), "slope": float(slope), "intercept": float(intercept),
        "r2": float(1 - resid.var() / ys.var()),
        "n": len(xs), "p_range": [float(min(ps)), float(max(ps))],
        "e_plus_range": [float(ep.min()), float(ep.max())],
        "seconds": time.time() - t0,
    }


def extrapolate(oms: np.ndarray, rs: np.ndarray) -> list[dict]:
    """Every candidate approach, reported together -- never just the flattering one."""
    forms = {
        "1/Omega": 1.0 / oms,
        "Omega^-0.75": oms ** -0.75,
        "1/sqrt(Omega)": 1.0 / np.sqrt(oms),
        "1/ln(Omega)": 1.0 / np.log(oms),
    }
    out = []
    for name, x in forms.items():
        sl, ic = np.polyfit(x, rs, 1)
        res = rs - (sl * x + ic)
        out.append({"form": name, "limit": float(ic),
                    "r2": float(1 - res.var() / rs.var()),
                    "max_resid": float(np.abs(res).max())})
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.35)
    ap.add_argument("--eps-frac", type=float, default=0.25)
    ap.add_argument("--omegas", type=int, nargs="+",
                    default=[100, 150, 200, 300, 400, 600, 800])
    ap.add_argument("--n-beta", type=int, default=15)
    ap.add_argument("--beta-max-frac", type=float, default=0.13,
                    help="as a fraction of beta_c, AT THE REFERENCE Omega below")
    ap.add_argument("--beta-ref-omega", type=float, default=200.0)
    ap.add_argument("--beta-omega-exponent", type=float, default=0.70,
                    help="beta* shrinks with Omega (measured exponent 0.51-0.75 "
                         "in FINDINGS 17.2), so the grid must shrink with it or "
                         "the sampled priors drift to 1 as Omega grows")
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/tilt_rule_limit.json"))
    args = ap.parse_args()

    bc = beta_critical(args.gamma)
    print(f"gamma={args.gamma}  eps={args.eps_frac}*delta*  beta_c={bc:.4f}")
    print(f"{'Omega':>7} {'slope r':>9} {'intercept':>10} {'R^2':>9} {'n':>3} "
          f"{'p range':>15} {'e+ range':>21} {'s':>6}")
    rows = []
    for om in args.omegas:
        scale = (args.beta_ref_omega / om) ** args.beta_omega_exponent
        betas = np.linspace(0.0, args.beta_max_frac * bc * scale, args.n_beta)
        r = rule_line(args.gamma, int(om), args.eps_frac, betas)
        if r is None:
            print(f"{om:>7}   no usable window")
            continue
        rows.append(r)
        print(f"{om:>7} {r['slope']:>9.4f} {r['intercept']:>10.4f} {r['r2']:>9.6f} "
              f"{r['n']:>3} {r['p_range'][0]:>6.3f}-{r['p_range'][1]:<8.3f} "
              f"{r['e_plus_range'][0]:>9.2e}-{r['e_plus_range'][1]:<10.2e} "
              f"{r['seconds']:>6.0f}")

    if len(rows) >= 4:
        oms = np.array([r["omega"] for r in rows], dtype=float)
        rs = np.array([r["slope"] for r in rows])
        print(f"\nr(Omega) over a {oms.max() / oms.min():.0f}x range: "
              + ", ".join(f"{v:.4f}" for v in rs))
        print(f"{'assumed correction':>20} {'limit':>8} {'R^2':>9} {'max resid':>10}")
        for e in extrapolate(oms, rs):
            print(f"{e['form']:>20} {e['limit']:>8.4f} {e['r2']:>9.5f} "
                  f"{e['max_resid']:>10.5f}")
        print("\nAll four are printed on purpose. Reporting only the form that "
              "lands near 1\nwould be choosing the answer, which is the failure "
              "FINDINGS 17.2 flags.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
