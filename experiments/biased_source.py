"""Is symmetric restoration still optimal when the SOURCE is biased? -- Q4a.

FINDINGS 16 showed mutual information is maximised at beta = 0 for a 50/50
source, which is the easy case: symmetry alone forces beta = 0 to be stationary,
and the only content was the sign of the second derivative. A source with prior
p != 1/2 has no such symmetry, and if the optimal tilt is nonzero it is the first
DESIGN RULE this project would have -- a statement about how to build the
chemistry, not about how it behaves.

THE PREDICTIONS, derived before running and recorded here so the fit cannot be
retrofitted. Write e+ = P(err | X sent), e- = P(err | Y sent). For small errors
the information deficit is L ~ p*h(e+) + (1-p)*h(e-) with h(e) = -e log2 e, and
e± = A exp(-Omega c±) with c+' = -c-' at beta = 0. Then dL/dbeta = 0 gives

  (P1)  p * h(e+)  =  (1-p) * h(e-)
        -- at the optimum the two symbols contribute EQUALLY to the deficit.

  (P2)  ln(e-/e+)  =  ln(p/(1-p))
        -- the leading-order form of P1, dropping the log factors. Completely
        parameter-free: no Omega, no rate constants, no fitted coefficient.

  (P3)  beta* ~ ln(p/(1-p)) / (Omega * K),  K = d(c+ - c-)/dbeta at beta = 0
        -- so the optimal tilt SHRINKS LIKE 1/Omega. The design rule vanishes in
        the very limit where restoration works best, which is the interesting
        part if it holds.

  (P4)  the deficit falls by 2*sqrt(p(1-p)) -- 13% at p = 0.75, 40% at p = 0.9.

Keeping the exact objective rather than the asymptotic one: `information` below
computes I = H(q) - [p H(e+) + (1-p) H(e-)] exactly and beta* is found by
maximising THAT. P1-P4 are then checked at the measured beta*, not assumed.
Retaining the output-entropy term matters -- it contributes at first order in e,
and dropping it is what turns P1 into the cruder P2.

Note the two objectives differ: minimising total ERROR gives p e+ = (1-p) e-,
minimising the information deficit gives P1. `--objective` runs either.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
from scipy.optimize import minimize_scalar

from crnl.cme import first_passage
from crnl.networks.am_asymmetric import am_asymmetric, beta_critical
from crnl.networks.am_reversible import delta_star

THETA = 0.80


def _h2(p: float) -> float:
    p = min(max(p, 1e-300), 1.0 - 1e-16)
    return float(-p * np.log2(p) - (1 - p) * np.log2(1 - p))


def _hdef(e: float) -> float:
    """h(e) = -e log2 e, one symbol's contribution to the deficit."""
    e = max(e, 1e-300)
    return float(-e * np.log2(e))


def errors(gamma: float, beta: float, omega: int, eps_frac: float):
    """(e+, e-) exactly. Inputs carry EXACTLY +-d0 counts -- see FINDINGS 16."""
    ds = delta_star(gamma)
    d0 = max(1, int(round(eps_frac * ds * omega)))
    nb = int(round(omega * gamma / (1.0 + gamma)))
    rest = omega - nb
    if (rest - d0) % 2:
        d0 -= 1
    thr = max(2, int(round(THETA * ds * omega)))
    net = am_asymmetric(gamma, beta)

    def absorbing(s, thr=thr):
        return abs(int(s[0]) - int(s[1])) >= thr

    out = []
    for s in (+1, -1):
        nx, ny = (rest + s * d0) // 2, (rest - s * d0) // 2
        assert nx - ny == s * d0
        fp = first_passage(net, int(omega), float(omega),
                           np.array([nx, ny, nb], dtype=np.int64), absorbing, None)
        if not fp["valid"]:
            return None
        out.append(fp["split"])
    a, b = out
    return max(1.0 - a, 0.0), max(b, 0.0)


def information(p: float, e_plus: float, e_minus: float) -> float:
    """Exact I(source ; absorbed state) in bits, prior p on X."""
    q = p * (1.0 - e_plus) + (1.0 - p) * e_minus
    return _h2(q) - (p * _h2(e_plus) + (1.0 - p) * _h2(e_minus))


def optimal_beta(gamma: float, omega: int, eps_frac: float, p: float,
                 objective: str = "info") -> dict | None:
    bc = beta_critical(gamma)
    cache: dict[float, tuple] = {}

    def ee(b: float):
        key = round(b, 12)
        if key not in cache:
            cache[key] = errors(gamma, key, omega, eps_frac)
        return cache[key]

    def loss(b: float) -> float:
        r = ee(float(b))
        if r is None:
            return 1e9
        ep, em = r
        if objective == "error":
            return p * ep + (1 - p) * em
        return -information(p, ep, em)

    res = minimize_scalar(loss, bounds=(0.0, 0.60 * bc), method="bounded",
                          options={"xatol": 1e-7})
    if not res.success:
        return None
    b = float(res.x)
    r0, rs = ee(0.0), ee(b)
    if r0 is None or rs is None:
        return None
    ep, em = rs
    L0 = _h2(p) - information(p, *r0)
    Ls = _h2(p) - information(p, ep, em)
    return {
        "gamma": gamma, "omega": omega, "eps_frac": eps_frac, "p": p,
        "beta_star": b, "beta_c": bc, "beta_star_over_bc": b / bc,
        "e_plus": ep, "e_minus": em,
        "e_at_zero": r0[0],
        # P1: the two symbols' deficit contributions
        "P1_lhs": p * _hdef(ep), "P1_rhs": (1 - p) * _hdef(em),
        # P2: log-ratio of errors against the prior log-odds
        "P2_measured": float(np.log(em / ep)) if ep > 0 else float("nan"),
        "P2_predicted": float(np.log(p / (1 - p))),
        # P3: Omega * beta*
        "omega_beta": omega * b,
        # P4: deficit reduction
        "deficit_ratio": Ls / L0 if L0 > 0 else float("nan"),
        "deficit_predicted": float(2 * np.sqrt(p * (1 - p))),
    }


def report(rows: list[dict], title: str) -> None:
    print(f"\n=== {title}")
    print(f"{'Om':>5} {'p':>5} {'beta*':>9} {'b*/bc':>7} {'Om*b*':>7} "
          f"{'e+':>9} {'e-':>9} | {'P2 meas':>8} {'P2 pred':>8} | "
          f"{'P1 l/r':>7} | {'defic':>7} {'pred':>6}")
    for r in rows:
        print(f"{r['omega']:>5} {r['p']:>5.2f} {r['beta_star']:>9.6f} "
              f"{r['beta_star_over_bc']:>7.4f} {r['omega_beta']:>7.3f} "
              f"{r['e_plus']:>9.2e} {r['e_minus']:>9.2e} | "
              f"{r['P2_measured']:>8.4f} {r['P2_predicted']:>8.4f} | "
              f"{r['P1_lhs'] / r['P1_rhs']:>7.4f} | "
              f"{r['deficit_ratio']:>7.4f} {r['deficit_predicted']:>6.4f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.35)
    ap.add_argument("--eps-frac", type=float, default=0.25)
    ap.add_argument("--objective", choices=["info", "error"], default="info")
    ap.add_argument("--omegas", type=int, nargs="+",
                    default=[100, 150, 200, 300, 400])
    ap.add_argument("--priors", type=float, nargs="+",
                    default=[0.60, 0.70, 0.80, 0.90, 0.95])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/biased_source.json"))
    args = ap.parse_args()

    g, ef = args.gamma, args.eps_frac
    print(f"gamma={g}  eps={ef}*delta*  beta_c={beta_critical(g):.4f}  "
          f"objective={args.objective}")

    # P3: does Omega * beta* hold still as Omega moves, at one prior?
    om_rows = [r for om in args.omegas
               if (r := optimal_beta(g, om, ef, 0.80, args.objective))]
    report(om_rows, "P3: sweep Omega at p = 0.80  (Om*b* should be flat)")

    # P2/P1/P4: does beta* track the prior log-odds, at one Omega?
    p_rows = [r for p in args.priors
              if (r := optimal_beta(g, 200, ef, p, args.objective))]
    report(p_rows, "P1/P2/P4: sweep the prior at Omega = 200")

    if len(p_rows) >= 3:
        x = np.array([r["P2_predicted"] for r in p_rows])
        y = np.array([r["omega_beta"] for r in p_rows])
        sl, ic = np.polyfit(x, y, 1)
        res = y - (sl * x + ic)
        print(f"\nP3 shape: Omega*beta* = {sl:.4f} * ln(p/(1-p)) + {ic:.4f}   "
              f"R^2 = {1 - res.var() / y.var():.5f}")
        print("  A zero intercept with a constant slope IS the design rule.")
        print("  The slope is NOT 1/K from the attractor-to-saddle barriers: that")
        print("  gives K = 0.874 at gamma = 0.35, hence 1.143, and the measured")
        print("  slope is 1.9x larger. Those are different barriers -- the errors")
        print("  here start from a biased INPUT, not from the attractor -- and")
        print("  conflating them is what made the original Q4a prediction wrong.")
        print("  The parameter-free statement is P2 on ln(e-/e+), not this slope.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"omega_sweep": om_rows,
                                    "prior_sweep": p_rows}, indent=2))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
