"""The cost of restoration, predicted in CLOSED FORM and tested absolutely — T-COST-b

§38 priced restoration per e-fold of gain, `G = Sigma / (Omega ln(theta/eps))`, and left
a residual 7-27% preparation-dependence with negative (unphysical) fitted intercepts.
The proposed fix was a four-term fit. **Physics rules that out and gives something
better.** As eps -> theta the start IS the threshold, absorption is immediate and
Sigma -> 0, so any model must vanish with the gain: the constant and the bare-Omega term
are forbidden, and `Sigma = ln(gain) * (G*Omega + B)` is the most a fit may have.

**But ln(theta/eps) was only ever an approximation to the traversal.** It is the time to
cross under `d(delta)/dt = lambda*delta`, i.e. pure exponential growth. The real drift
saturates near the attractor, so the true traversal time is `int d(delta)/mu(delta)` and
the entropy spent is the ep rate integrated along that path:

    **Sigma_pred = Omega * int_{x0}^{x_thr}  sigma(delta) / mu(delta)  d(delta)**

with `sigma(delta) = sum_r f_r ln(f_r / f_rev)` the local entropy rate per molecule and
`mu(delta)` the drift, both evaluated on the slaved manifold from the network's own
fluxes. **This has no fitted parameter of any kind.** It is the same move §15/§28 made
for the collapse slope, applied to the cost -- rule 16's absolute test rather than a
better fit.

PREDICTIONS, written before running:

  P1  GATE, and the whole point. `Sigma_pred` reproduces the exact CME `Sigma` to within
      a few percent at every eps and Omega, with **no fitting**. If it does, the cost of
      restoration is a closed-form quantity and §38's G is a projection of it.
  P2  The ratio Sigma_pred/Sigma_exact is eps-INDEPENDENT to much better than §38's
      7-27%, because the traversal integral captures exactly what ln(theta/eps) missed --
      the drift's saturation near the attractor.
  P3  The residual, whatever it is, should shrink like 1/Omega: `Sigma_pred` is the
      deterministic path and `Sigma_exact` includes fluctuation contributions, which are
      subleading. A residual that does NOT shrink with Omega is a stochastic excess that
      the deterministic path cannot account for, and is itself a quantity worth naming.
  P4  With the closed form in hand the optimal drive is computable directly by minimising
      the integral over gamma, with no CME solve at all. It should land on §38's
      **gamma* ~ 0.20**. If it lands elsewhere, one of the two is wrong and the
      disagreement is the result.
  P5  If P1 FAILS -- if the deterministic path badly under- or over-predicts -- then the
      entropy of a restoring decision is dominated by fluctuations rather than by the
      mean path, which would be a strong and surprising statement, and would explain why
      no simple normalisation of Sigma has been preparation-free.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.networks.am_reversible import am_reversible, delta_star, reverse_pairing
from experiments.cost_of_reliability import cell
from experiments.slaving_axis import slaved

U = np.array([1.0, -1.0, 0.0])


def sigma_and_mu(net, delta, pairing):
    """Local entropy rate per molecule and drift of delta, on the slaved manifold."""
    st = slaved(net, delta)
    if st is None:
        return None
    f = net.fluxes(st)
    S = net.stoichiometry_matrix()
    sig = 0.0
    for r in range(net.n_reactions):
        rev = int(pairing[r])
        if rev < 0 or f[r] <= 0 or f[rev] <= 0:
            continue
        sig += float(f[r]) * float(np.log(f[r] / f[rev]))
    mu = float(U @ (S @ f))
    return sig, mu


def sigma_predicted(net, x0, xthr, pairing, n=4001):
    xs = np.linspace(x0, xthr, n)
    vals = []
    for x in xs:
        sm = sigma_and_mu(net, float(x), pairing)
        if sm is None or sm[1] <= 0:
            return float("nan")
        vals.append(sm[0] / sm[1])
    return float(np.trapezoid(vals, xs))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.07, 0.20, 0.30])
    ap.add_argument("--epss", type=float, nargs="+", default=[0.20, 0.30, 0.40, 0.50])
    ap.add_argument("--omegas", type=int, nargs="+", default=[150, 300, 450])
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/cost_absolute.json"))
    args = ap.parse_args()

    t0 = time.time()
    print("Sigma_pred = Omega * int sigma(d)/mu(d) dd   -- NO fitted parameter")
    print(f"{'gamma':>7}{'eps':>7}{'Omega':>7}{'Sigma exact':>14}{'Sigma pred':>13}"
          f"{'pred/exact':>12}")
    rows = []
    for g in args.gammas:
        net = am_reversible(g)
        pairing = reverse_pairing(net)
        ds = delta_star(g)
        for eps in args.epss:
            per_om = sigma_predicted(net, eps * ds, args.theta * ds, pairing)
            for om in args.omegas:
                try:
                    r = cell(g, om, eps, args.theta)
                except Exception:
                    continue
                if not np.isfinite(r["Sigma"]):
                    continue
                pred = om * per_om
                rows.append({"gamma": g, "eps": eps, "omega": om,
                             "exact": r["Sigma"], "pred": float(pred),
                             "ratio": float(pred / r["Sigma"])})
                print(f"{g:>7.2f}{eps:>7.2f}{om:>7}{r['Sigma']:>14.3f}"
                      f"{pred:>13.3f}{pred/r['Sigma']:>12.4f}")
        print()

    rr = np.array([r["ratio"] for r in rows])
    print(f"=== P1 gate: does the closed form reproduce the exact cost, unfitted?")
    print(f"  pred/exact over {len(rr)} cells: {rr.min():.4f}..{rr.max():.4f}"
          f"  mean {rr.mean():.4f}  sd {rr.std(ddof=1):.4f}")
    print(f"  -> P1 {'HOLDS' if abs(rr.mean()-1) < 0.05 else 'FAILS'}"
          f"   (mean within {100*abs(rr.mean()-1):.2f}% of 1)")

    print(f"\n=== P2: is the ratio eps-independent? (§38's G spread 7-27%)")
    print(f"{'gamma':>7}{'eps':>7}{'mean ratio':>13}")
    for g in args.gammas:
        per_eps = []
        for eps in args.epss:
            sel = [r["ratio"] for r in rows
                   if r["gamma"] == g and r["eps"] == eps]
            if sel:
                per_eps.append(np.mean(sel))
                print(f"{g:>7.2f}{eps:>7.2f}{np.mean(sel):>13.4f}")
        if len(per_eps) >= 2:
            pe = np.array(per_eps)
            print(f"{'':>7}{'spread':>7}{100*(pe.max()-pe.min())/pe.mean():>12.2f}%")

    print(f"\n=== P3: does the residual shrink with Omega?")
    print(f"{'Omega':>7}{'mean |ratio-1|':>17}")
    for om in args.omegas:
        sel = [abs(r["ratio"] - 1) for r in rows if r["omega"] == om]
        if sel:
            print(f"{om:>7}{np.mean(sel):>17.4f}")

    print(f"\n=== P4: the optimal drive from the closed form alone (no CME)")
    print(f"{'gamma':>8}{'int sigma/mu':>15}")
    gs = np.linspace(0.04, 0.44, 21)
    vals = []
    for g in gs:
        net = am_reversible(g)
        v = sigma_predicted(net, 0.35 * delta_star(g), args.theta * delta_star(g),
                            reverse_pairing(net), n=2001)
        vals.append(v)
    vals = np.array(vals)
    ok = np.isfinite(vals)
    for g, v in zip(gs[ok][::4], vals[ok][::4]):
        print(f"{g:>8.3f}{v:>15.5f}")
    i = int(np.argmin(vals[ok]))
    print(f"\n  minimum at gamma = {gs[ok][i]:.3f}, value {vals[ok][i]:.5f} k_B/molecule")
    print(f"  §38's CME-measured optimum was gamma* ~ 0.20  -> "
          f"{'AGREE' if abs(gs[ok][i]-0.20) < 0.05 else 'DISAGREE -- one of them is wrong'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"cells": rows,
                                    "gamma_scan": [(float(a), float(b))
                                                   for a, b in zip(gs, vals)]},
                                   indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
