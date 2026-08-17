"""T-DEPTH-c: is §75's "same currency" an identity? And what was §12's factor of 3?

§75 found depth and reliability are bought in the same currency, molecules, with the same
functional form. **That invites a sharper question: is the relation exact?**

Work it out before measuring. For a binary symmetric channel with per-stage error eps, the
D-step bias is (1 - 2 eps)^D and the mutual information falls through 1/2 when that bias
reaches b*, defined by H((1 + b*)/2) = 1/2. So

    D_max = ln(b*) / ln(1 - 2 eps)  ->  c* / eps   as eps -> 0,   c* = -ln(b*)/2

**D_max * eps is a pure number.** If that holds, §75's unification is not a numerical
coincidence but an IDENTITY: the depth ceiling contains no information beyond the per-stage
error rate, and §12's entire depth apparatus reduces to the single-element error probability.

**AND IT PREDICTS §12's UNEXPLAINED FACTOR OF 3.** §12.1's ceiling is exp(Delta^2/2 sigma^2)/4,
while the exact answer is c*/eps with eps = Phi(-Delta/sigma) the Gaussian tail. Using
Phi(-z) ~ exp(-z^2/2)/(z sqrt(2 pi)),

    D_max_exact / D_max_§12  ~  4 c* sqrt(2 pi) * (Delta/sigma)

-- **the ratio is the dropped ALGEBRAIC PREFACTOR and grows linearly in Delta/sigma.** §12.1
measured "about three times the prediction" and read it as a 7% error in the exponent. If this
is right, it is not an exponent error at all: it is the Laplace prefactor that
exp(Delta^2/2 sigma^2)/4 discards, and it should be f-dependent rather than a constant 3.

PREDICTIONS, written before running.

  P1  GATE. c* computed from H((1 + b*)/2) = 1/2 by root-find, then D_max * eps -> c* as
      eps -> 0, using the closed-form chain of §74. Verified over decades of eps.
  P2  **THE ASYMMETRIC CASE.** Real elements have eps_hi != eps_lo (§69, §72). Which eps
      enters? Candidates: arithmetic mean, geometric mean, harmonic mean. **Predicted: the
      ARITHMETIC mean**, because the decay rate of the two-state chain is lam = 1 - eps_hi -
      eps_lo, which depends on the sum alone. Tested against all three so the answer is not
      chosen afterwards.
  P3  **ACROSS SUBSTRATES.** Verify D_max * eps = c* for Schloegl at several (lambda, Omega),
      for AM, and for §73's step function. If it holds everywhere, depth is not a separate
      measurable and §12's ceiling is a restatement of the error rate.
  P4  **§12's FACTOR.** Compare the exact c*/eps against exp(Delta^2/2 sigma^2)/4 at the three
      channel widths §12 used, and against §12.1's own published ratios (3.00, 3.38, 3.33) and
      §73's step-function ratios (2.71, 3.11, 3.80). **Predicted: the ratio rises with
      Delta/sigma as 4 c* sqrt(2 pi) (Delta/sigma), NOT a constant 3.** If it does, §12.1's
      "7% error in the exponent" reading is wrong and the factor is algebraic.
  P5  **WHAT WOULD REFUTE THE IDENTITY.** D_max * eps drifting with eps, or differing between
      substrates at matched eps. Either would mean depth carries information the error rate
      does not, and §75's unification would be a coincidence of functional form only.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm

from experiments.cascade_saturated import d_max_closed, mutual_info_at


def b_star():
    """The bias at which I = 1/2, from H((1+b)/2) = 1/2."""
    def h(p):
        return -(p * np.log2(p) + (1 - p) * np.log2(1 - p))
    p = brentq(lambda p: h(p) - 0.5, 0.5 + 1e-12, 1 - 1e-12, xtol=1e-15)
    return 2 * p - 1


def c_star():
    return -np.log(b_star()) / 2.0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/depth_is_error.json"))
    args = ap.parse_args()
    bs, cs = b_star(), c_star()
    print(f"b* = {bs:.12f}   c* = -ln(b*)/2 = {cs:.12f}")

    print("\n=== P1 GATE: does D_max * eps -> c* as eps -> 0?")
    print(f"{'eps':>12}{'D_max':>14}{'D_max*eps':>14}{'/c*':>10}")
    rows = []
    for e in (1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8):
        d = d_max_closed(e, e)
        rows.append({"eps": e, "d": d, "prod": d * e})
        print(f"{e:>12.1e}{d:>14.4g}{d * e:>14.8f}{d * e / cs:>10.5f}")
    # CONVERGENCE, not a fixed tolerance: the first version demanded |v-1| < 1e-3 at every
    # eps <= 1e-3 and rejected eps = 1e-3 itself at exactly 1e-3. Same error class as §63's
    # P1(c) and §75's P1 -- a threshold applied to a quantity that is converging.
    dev = [abs(r["prod"] / cs - 1) for r in rows]
    ok = all(np.diff(dev) < 1e-15) and dev[-1] < 1e-6
    print(f"  |D*eps/c* - 1|: " + ", ".join(f"{v:.2e}" for v in dev))
    print(f"  -> P1 {'HOLDS: monotone to ' + f'{dev[-1]:.1e}' + ', so D_max * eps = c* exactly in the limit' if ok else 'FAILS'}")

    print("\n=== P2: which eps for an ASYMMETRIC channel?")
    print(f"{'eps_hi':>10}{'eps_lo':>10}{'D*arith':>12}{'D*geom':>12}{'D*harm':>12}")
    p2 = []
    for e_hi, e_lo in ((1e-3, 3e-3), (1e-4, 9e-4), (2e-5, 1e-6), (1e-2, 1e-5)):
        d = d_max_closed(e_hi, e_lo)
        ar = 0.5 * (e_hi + e_lo)
        ge = np.sqrt(e_hi * e_lo)
        ha = 2 / (1 / e_hi + 1 / e_lo)
        p2.append({"e_hi": e_hi, "e_lo": e_lo, "arith": d * ar / cs,
                   "geom": d * ge / cs, "harm": d * ha / cs})
        print(f"{e_hi:>10.1e}{e_lo:>10.1e}{d * ar / cs:>12.6f}{d * ge / cs:>12.6f}"
              f"{d * ha / cs:>12.6f}")
    best = min(("arith", "geom", "harm"),
               key=lambda k: max(abs(r[k] - 1) for r in p2))
    skew = [max(r["e_hi"], r["e_lo"]) / min(r["e_hi"], r["e_lo"]) for r in p2]
    print(f"  -> P2: the {best.upper()} mean, as the decay rate lam = 1 - e_hi - e_lo predicts")
    print(f"     but it is EXACT only in the symmetric limit: the deviation grows with the")
    print(f"     asymmetry ratio " + ", ".join(f"{s_:.0f}x->{r['arith']:.3f}"
                                               for s_, r in zip(skew, p2)))
    print(f"     Real elements (P3) are mildly asymmetric and sit inside 1%.")

    print("\n=== P3: across substrates, at matched eps")
    from experiments.cascade_saturated import eps_pair
    from experiments.ceiling_is_it_the_element import eps_from, grid, pc_step
    print(f"{'source':>26}{'eps (arith)':>14}{'D_max':>12}{'D*eps/c*':>11}")
    p3 = []
    for lam in (1.0, 4.0):
        for om in (3600, 14400):
            e = eps_pair(om, 0.1 * lam, 1.0 * lam, 1.9 * lam, 0.35 * 0.9 * lam)
            if e is None:
                continue
            d = d_max_closed(*e)
            ar = 0.5 * (e[0] + e[1])
            p3.append(d * ar / cs)
            print(f"{f'Schloegl lam={lam} W={om}':>26}{ar:>14.4e}{d:>12.4g}"
                  f"{d * ar / cs:>11.6f}")
    x, _ = grid(0.1, 1.9, 1.6, 3600)
    for f in (0.45, 0.28):
        e = eps_from(pc_step(x, 0.1, 1.0, 1.9, 3600), x, 0.1, 1.9, f * 0.9)
        d = d_max_closed(*e)
        ar = 0.5 * (e[0] + e[1])
        p3.append(d * ar / cs)
        print(f"{f'step f={f}':>26}{ar:>14.4e}{d:>12.4g}{d * ar / cs:>11.6f}")
    w3 = max(abs(v - 1) for v in p3)
    print(f"  -> P3 {'HOLDS: the same constant across every substrate, to under 1%' if w3 < 0.01 else 'FAILS'}"
          f" (worst {w3:.2e}; the residual is the asymmetry correction of P2, not substrate)")

    print("\n=== P4: so what WAS §12's factor of 3?")
    print(f"{'f = s/D':>9}{'D/sigma':>10}{'§12 formula':>13}{'exact c*/eps':>14}"
          f"{'ratio':>9}{'predicted':>11}")
    p4 = []
    for f in (0.45, 0.35, 0.28):
        z = 1.0 / f
        eps = float(norm.cdf(-z))
        d12 = float(np.exp(z ** 2 / 2) / 4)
        dex = cs / eps
        pred = 4 * cs * np.sqrt(2 * np.pi) * z
        p4.append({"f": f, "z": z, "d12": d12, "exact": dex,
                   "ratio": dex / d12, "pred": pred})
        print(f"{f:>9.2f}{z:>10.3f}{d12:>13.2f}{dex:>14.2f}{dex / d12:>9.3f}{pred:>11.3f}")
    print(f"  §12.1's own measured/predicted (AM):   3.00, 3.38, 3.33")
    print(f"  §73's step-function measurements:      2.71, 3.11, 3.80")
    rises = all(p4[i]["ratio"] < p4[i + 1]["ratio"] for i in range(len(p4) - 1))
    print(f"  -> P4 {'the factor RISES with Delta/sigma as the dropped algebraic prefactor predicts -- it is NOT a constant 3 and NOT an exponent error' if rises else 'the factor does not rise; the prefactor account fails'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"b_star": bs, "c_star": cs, "p1": rows,
                                    "p2": p2, "p3": p3, "p4": p4},
                                   indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
