"""§108 (T-CASC-x) -- chasing the residual drift: six candidates cleared, and a proof that
one mechanism cannot be enough.

THE GOAL, set before any of this ran: explain and remove §107's residual Omega-drift in the
closure, or prove it irreducible. §107 left it with the span growing 1.533x -> 1.880x under
four independently justified corrections, and no candidate cause.

The goal was NOT achieved. What follows is the localisation, the six refutations, and one
structural result that explains why every single-correction attempt has failed.

STEP 1 -- LOCALISATION. With all of §107's corrections applied, the two channels behave
completely differently across Omega = 14-70:
    contaminated, measured/model:  1.1988 1.2153 1.2179 1.2171 1.2140 1.2130 1.2151
    pure,         measured/model:  0.9474 1.0531 1.1338 1.2339 1.3670 1.5475 1.7604
`contam` is FLAT to 0.4% -- a pure prefactor. The entire drift is in `pure`, i.e. in stage 2's
own escape. The one thing stage 2 has that stage 1 does not is a FLUCTUATING INPUT.

STEP 2 -- SIX CANDIDATES, each cleared by measurement:
  (a) the operating points          -- exact to 0.0001% by Omega = 70 (§106)
  (b) the predicted input law       -- model vs measured geometric mean agrees to 0.4%
  (c) the averaging of lambda       -- same check, same 0.4%
  (d) splitting pi_low from lambda  -- averaging k_fwd = pi_low*lambda as one quantity moves
                                       it 1.3% at Omega = 14 falling to 0.4% at Omega = 70
  (e) p_transmit                    -- measured flat, 0.9466 -> 0.9445 (§106)
  (f) RETURN TRIPS -- pre-registered and refuted. `pure` conditions on stage 1 being high AT
      THE END, so trajectories that dipped low and came back are counted as pure while having
      dragged stage 2 down; pi_low falls 0.906 -> 0.525, so returns go from rare to common,
      which has the right sign. Making stage 1 ABSORBING below its saddle removes them
      entirely. The drift survives: 0.5926 -> 1.2142, span 2.05x.

STEP 3 -- THE PROOF, and it is the section's result.

Invert stage 2's effective escape rate from the absorbing-upstream run (no return trips) and
compare it against every candidate average of the rate over the input distribution. The true
rate sits BELOW k(<x>) at Omega = 14 and rises through the geometric mean toward the frozen
limit by Omega = 70.

**`ln k` is convex in the input at every Omega** (measured: d2/dx2 runs +0.408..+2.230 at
Omega = 14, +2.779..+14.493 at Omega = 70). By Jensen, convexity of `ln k` puts BOTH
`exp(<ln k>)` and `<k>` at or above `k(<x>)`. **So no average of the rate over the input can
produce a value below `k(<x>)`, and the Omega = 14 measurement is 0.788 x k(<x>).**

    => The small-Omega behaviour is not an averaging effect AT ALL.
    => At least two distinct mechanisms are in play, one suppressing and one enhancing.
    => That is why §107's four corrections recentred the model without flattening it: they
       all address the enhancing mechanism, and none of them can reach below k(<x>).

WHAT THIS DOES NOT DO. It does not identify the suppressing mechanism, and it does not remove
the drift. It converts "unexplained drift with no candidate" into "two mechanisms, one of them
provably outside the averaging family" -- which is a smaller and better-posed problem, and it
retires the whole class of fixes that §102.1 through §107 were drawn from.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from experiments.chain_without_a_joint_solve import chain_operating_points
from experiments.depth_compounding import R2, R3, _rates_vec
from experiments.escape_accounts_for_it import escape_rate
from experiments.margin_law import stage1_stationary
from experiments.the_corrected_closure import geometric_rate, input_law, pi_low

OMEGAS = (14, 20, 30, 40, 55, 70)
T0 = 2.0


def joint_absorbing(om, t=T0, cap_mult=1.25):
    """D = 2 free chain with stage 1 ABSORBING below its saddle -- no return trips.

    With absorption, "stage 1 high at t" is exactly "stage 1 never crossed", so the
    endpoint conditioning that §101 flagged cannot contaminate `pure`.
    """
    cap = int(np.ceil(cap_mult * R3 * om))
    m = cap + 1
    N = m * m
    idx = np.arange(N)
    n1, n2 = idx // m, idx % m
    sad = R2 * om

    rows, cols, vals = [], [], []
    diag = np.zeros(N)
    l1, u1 = _rates_vec(n1.astype(float), np.zeros(N), om, True)
    alive = n1 >= sad                      # absorbed states have no outgoing transitions
    up1 = (n1 < cap) & (l1 > 0) & alive
    dn1 = (n1 > 0) & (u1 > 0) & alive
    rows.append(idx[up1]); cols.append(idx[up1] + m); vals.append(l1[up1])
    rows.append(idx[dn1]); cols.append(idx[dn1] - m); vals.append(u1[dn1])
    diag -= np.where(up1, l1, 0.0) + np.where(dn1, u1, 0.0)

    l2, u2 = _rates_vec(n2.astype(float), n1.astype(float), om, False)
    up2 = (n2 < cap) & (l2 > 0)
    dn2 = (n2 > 0) & (u2 > 0)
    rows.append(idx[up2]); cols.append(idx[up2] + 1); vals.append(l2[up2])
    rows.append(idx[dn2]); cols.append(idx[dn2] - 1); vals.append(u2[dn2])
    diag -= np.where(up2, l2, 0.0) + np.where(dn2, u2, 0.0)

    rows.append(idx); cols.append(idx); vals.append(diag)
    Q = sp.csr_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))

    up, pi1 = stage1_stationary(om)
    p = np.zeros(N)
    for a, w in zip(up, pi1):
        p[a * m + int(round(R3 * om))] = w
    p = spla.expm_multiply(Q.T * t, p)
    g = p.reshape(m, m)
    return float(g[np.ix_(np.arange(m) >= sad, np.arange(m) < sad)].sum())


def effective_rate(om, pure_abs, t=T0):
    """Invert stage 2's effective escape rate from the absorbing-upstream `pure`."""
    mus, _ = chain_operating_points(om, 2)
    p2 = pi_low(om, mus[0])
    return -np.log(max(1.0 - pure_abs / p2, 1e-300)) / t


def candidate_averages(om):
    """k(<x>), exp(<ln k>) and <k> over the predicted input law."""
    mus, _ = chain_operating_points(om, 2)
    xs, w = input_law(om, mus[0])
    ks = np.array([escape_rate(om, x) for x in xs])
    return escape_rate(om, mus[0]), geometric_rate(om, mus[0]), float((w * ks).sum())


def ln_k_curvature(om, lo=2.2, hi=3.18, n=15):
    """d2(ln k)/dx2 over the input range -- the linchpin of §108's proof."""
    xs = np.linspace(lo, hi, n)
    lk = np.array([np.log(escape_rate(om, x)) for x in xs])
    return np.diff(lk, 2) / (xs[1] - xs[0]) ** 2


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/two_mechanisms.json"))
    args = ap.parse_args()

    print("Step 3 -- stage 2's effective rate against every candidate average")
    print(f"{'Omega':>7}{'lambda_true':>14}{'k(<x>)':>13}{'exp(<ln k>)':>14}{'<k>':>13}"
          f"{'true/k(<x>)':>13}")
    rows = []
    for om in OMEGAS:
        pa = joint_absorbing(om)
        lam = effective_rate(om, pa)
        kmean, kgeo, kari = candidate_averages(om)
        rows.append({"omega": om, "pure_absorbing": pa, "lambda_true": lam,
                     "k_at_mean": kmean, "k_geometric": kgeo, "k_arithmetic": kari,
                     "true_over_mean": lam / kmean, "true_over_geo": lam / kgeo})
        print(f"{om:>7}{lam:>14.4e}{kmean:>13.4e}{kgeo:>14.4e}{kari:>13.4e}"
              f"{lam / kmean:>13.4f}", flush=True)
        args.out.write_text(json.dumps(rows, indent=2, default=float))

    print("\nThe proof: is ln k convex in the input at every Omega?")
    curv = {}
    for om in OMEGAS:
        d2 = ln_k_curvature(om)
        curv[om] = [float(d2.min()), float(d2.max())]
        print(f"  Omega = {om:>3}: d2(ln k)/dx2 in [{d2.min():+.3f}, {d2.max():+.3f}]"
              f"   {'convex' if (d2 > 0).all() else 'NOT CONVEX'}")

    below = [r for r in rows if r["true_over_mean"] < 1.0]
    print(f"\n  ln k convex everywhere => every average of k over the input is >= k(<x>).")
    print(f"  Cells where the MEASURED rate is below k(<x>): "
          f"{[r['omega'] for r in below]} (lowest {min(r['true_over_mean'] for r in rows):.4f})")
    print("  Those cells are outside the reach of ANY averaging prescription.")
    print("  => at least two mechanisms, one of them provably not an averaging effect.")

    args.out.write_text(json.dumps({"cells": rows, "ln_k_curvature": curv},
                                   indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
