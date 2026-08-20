"""§102 (T-CASC-o) -- is the contaminated channel just upstream escape?

§101 measured that the failure channel is the majority of a free chain's error by D = 3,
and named a SUSPECT for why (rule 17): the cheapest route to a downstream error is one
upstream escape followed by near-certain propagation, rather than the last stage escaping
its own barrier. §101 also named the kill test, and this is it. The instrument is 1-D
generators only -- no joint master equation on the prediction side.

WHAT WAS SCOUTED BEFORE THIS FILE EXISTED, and is therefore NOT a prediction (rule 2). I
checked that the free stage-1 spectral gap reproduces the joint chain's P(stage 1 low at t)
through 1 - exp(-k t): at Omega = 30 it gives 7.5765e-3 against a measured 7.7739e-3 at
t = 2.0 (2.5% low), and misses by -34% at t = 0.5 and +13% at t = 8.0. So the exponential
is known to work at t = 2 and known to fail at both ends. Everything below about DEPTH and
about the CONTAMINATED/PURE SPLIT was written before running.

PREDICTIONS.

  P1  WIRING, and it is exact by construction. The Hill map is normalised so h(r3) = 1, so
      a downstream stage with its upstream pinned AT THE RAIL has propensities identical to
      stage 1's. Therefore escape_rate(om, R3) must equal §100's free stage-1 spectral gap
      to solver tolerance: 7.29805e-02 (Om=14), 3.80268e-03 (Om=30). A mismatch means the
      pinned generator is not the same element and nothing below means anything.

  P2  THE RATE IS A STEEP FUNCTION OF THE OPERATING POINT. escape_rate(om, x_up) must rise
      monotonically as x_up falls from r3 toward the saddle, and steeply -- this is the
      whole reason a degraded stage is more fragile than a fresh one. Reported as a curve;
      no threshold.

  P3  THE SUSPECT, STATED SO IT CAN FAIL. Take each stage's measured operating point in the
      FREE chain, read its escape rate off the 1-D curve, and predict

        P(stage i has failed by t)  ~  1 - exp(-k_i t)
        contaminated  ~  P(any upstream stage failed) x p_transmit
        pure          ~  P(no upstream failed) x (1 - exp(-k_last t))

      I predict the contaminated/pure ratio comes out within a FACTOR OF TWO of §101's
      measured value at t = 2.0, at both Omega and both depths. Within a factor of two the
      suspect survives as an account; an order of magnitude out and it is dead. That is a
      loose gate ON PURPOSE -- the exponential is already known to carry tens of percent of
      error at this window, so a tight gate would be testing the approximation rather than
      the mechanism (rule 20).

  P4  THE DISCRIMINATING PART. The suspect says the contaminated channel is governed by
      UPSTREAM escape and the pure channel by LAST-STAGE escape at a degraded operating
      point. Those two rates differ by a computable factor. If instead the contaminated
      channel were governed by the last stage's own barrier -- i.e. if propagation were
      irrelevant and both channels were the same physics -- contaminated/pure would be
      independent of the ratio k_upstream/k_last. **The kill:** k_upstream/k_last is a
      strong function of Omega (the barrier depths differ), so if the measured
      contaminated/pure does NOT track it across Omega = 14 vs 30, the suspect is wrong.

WHAT THIS CANNOT SETTLE. The operating points are taken from §101's joint solve rather than
predicted from §96, so this is not yet a CME-free prediction -- it tests whether ESCAPE
RATES AT THE OPERATING POINTS account for the split, not whether the whole chain can be
computed from single-element quantities. Closing that gap needs §96's predicted operating
points and is the follow-up, not this.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import scipy.sparse as sp

import experiments.chemical_cascade as cc
from experiments.depth_compounding import C, R2, R3
from experiments.free_upstream_depth import channel_split, last_low, solve
from experiments.what_reflection_costs import spectral_gap

P_TRANSMIT = 0.9376          # §100, Omega=30, t0=2.0 -- an endpoint co-occurrence


def escape_rate(om, x_up, cap_mult=1.25):
    """Slowest non-zero rate of ONE stage with its upstream pinned at concentration x_up."""
    cap = int(np.ceil(cap_mult * R3 * om))
    m = cap + 1
    rows, cols, vals = [], [], []
    for n in range(m):
        tot = 0.0
        lam, mu = cc.rates_stage(float(n), x_up * om, om, C, R3, False, "hill")
        if n < cap and lam > 0:
            rows.append(n); cols.append(n + 1); vals.append(lam); tot += lam
        if n > 0 and mu > 0:
            rows.append(n); cols.append(n - 1); vals.append(mu); tot += mu
        rows.append(n); cols.append(n); vals.append(-tot)
    Q = sp.csr_matrix((vals, (rows, cols)), shape=(m, m)).toarray()
    ev = np.sort(-np.real(np.linalg.eigvals(Q)))
    return float(ev[1])


def stage_mean_free(p, om, dims, strides, k):
    """Mean of stage k in a FULLY FREE chain, conditioned on that stage not having failed.

    §94's `stage_stats` cannot be used here: it infers the stage's indexing convention from
    its POSITION (every stage but the last reflected), which is right for §94's chain and
    wrong for this one, where every stage is free and indexed by its raw count. That
    mismatch is §96.1a's bug in the opposite direction, and it raised an IndexError rather
    than a wrong number only because the reflected range is shorter than the free one.
    """
    idx = np.arange(len(p))
    n = (idx // strides[k]) % dims[k]
    hi = n > R2 * om
    w = p[hi] / p[hi].sum()
    return float((w * n[hi].astype(float)).sum()) / om


def operating_points(om, D, t):
    """Each stage's mean concentration in the FREE chain, and the measured split."""
    p, ref, dims, strides, walled = solve(om, D, 0, t)
    assert not any(walled), "operating points are defined here for the fully free chain"
    mus = [stage_mean_free(p, om, dims, strides, k) for k in range(D)]
    tot, pure, contam = channel_split(p, om, dims, strides, ref, walled)
    return mus, tot, pure, contam


def predict(om, mus, t, p_transmit=P_TRANSMIT, legacy=False):
    """Rates at the measured operating points -> the contaminated/pure split.

    DEFAULT = the model that survives §106-§109. `legacy=True` restores the ORIGINAL model,
    which is KNOWN WRONG IN THREE WAYS and is kept only so §102's and §103's published
    numbers stay reproducible (rule 7). The maintained model is `the_corrected_closure.closure`.
    Every caller that reproduces a published table passes legacy=True explicitly; new work
    should not.

      (1) INDEXING (§106.3). `escape_rate(om, x_up)` is the rate of a stage whose UPSTREAM
          sits at x_up. The line below keys each stage to ITS OWN operating point. Stage 1
          has no upstream: its rate is escape_rate(om, R3), which equals the free spectral
          gap exactly at every Omega. As coded this is 12.7% high at Omega = 14.
      (2) ONE-WAY OCCUPANCY (§106.2). A free stage has no absorbing boundary, so
          P(low at t) = pi_low * (1 - exp(-lambda t)). pi_low falls 0.9057 -> 0.5247 over
          Omega = 14-70; taking it as 1 spans 1.71x and crosses 1 near Omega = 35, which is
          where §102 validated it.
      (3) RATE AT THE MEAN (§102.1, §107). The escape rate must be averaged over the
          fluctuating input; the fast-limit average is the GEOMETRIC mean exp(<ln k>), which
          is ~21% larger than the k(<x>) used here.

    A fourth defect is upstream of this function and cannot be fixed here: the joint chains
    that supply `mus` seed stage 2 as a delta at its rail and stage 1 from the REFLECTED
    stationary law rather than the free QSD (§109).
    """
    if not legacy:
        from experiments.the_corrected_closure import closure
        ratio, ks, p_t = closure(om, t, indexing=True, two_state=True, geometric=True)
        surv_contam = ratio / (1.0 + ratio)
        return ks, surv_contam, 1.0 - surv_contam
    ks = [escape_rate(om, x) for x in mus]          # legacy: see the three defects above
    surv = float(np.prod([np.exp(-k * t) for k in ks[:-1]]))   # no upstream stage failed
    contam = (1.0 - surv) * p_transmit
    pure = surv * (1.0 - np.exp(-ks[-1] * t))
    return ks, contam, pure


def rate_limits(om, D, t):
    """§102.1 -- bracket the LAST stage's effective escape rate by the two averaging limits.

    P3's residual is one-signed: the measured `pure` always exceeds a model that evaluates
    the escape rate AT the mean operating point. escape_rate is steeply convex in x_up, so
    the rate at the average and the average of the rate are different numbers, and they
    bracket the truth -- §92/§93's frozen/fast pair, one level up and applied to a rate
    rather than to a penalty.

    Returns (k_at_mean, k_averaged, k_effective, position) with position in log space,
    0 = rate at the mean (fast upstream), 1 = mean of the rate (frozen upstream).
    """
    p, ref, dims, strides, walled = solve(om, D, 0, t)
    idx = np.arange(len(p))
    up = (idx // strides[D - 2]) % dims[D - 2]          # the last stage's own input
    hi = up > R2 * om
    w = p[hi] / p[hi].sum()
    xs = up[hi].astype(float) / om

    k_mean = escape_rate(om, float((w * xs).sum()))
    tbl = {int(u): escape_rate(om, u / om) for u in np.unique(up[hi])}
    k_avg = float(sum(w[up[hi] == u].sum() * tbl[int(u)] for u in np.unique(up[hi])))

    # k_effective: invert the measured `pure` through the same survival model
    mus = [stage_mean_free(p, om, dims, strides, k) for k in range(D)]
    ks = [escape_rate(om, x) for x in mus]
    surv = float(np.prod([np.exp(-k * t) for k in ks[:-1]]))
    _, pure, _ = channel_split(p, om, dims, strides, ref, walled)
    frac = pure / surv
    k_eff = -np.log(max(1.0 - frac, 1e-300)) / t
    pos = np.log(k_eff / k_mean) / np.log(k_avg / k_mean)
    return k_mean, k_avg, float(k_eff), float(pos)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--deep", action="store_true", help="include Omega=30 D=3 (~25 min)")
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/escape_accounts_for_it.json"))
    args = ap.parse_args()

    print("P1 -- pinned at the rail, a downstream stage IS stage 1 (h(r3) = 1)")
    for om, ref_gap in ((14, 7.29805e-02), (30, 3.80268e-03)):
        k = escape_rate(om, R3)
        g, _ = spectral_gap(om, False)
        print(f"  Om={om:3d}  escape_rate(r3) = {k:.6e}   §100 free gap = {g:.6e}"
              f"   rel {abs(k/g-1):.2e}")

    print("\nP2 -- escape rate vs operating point")
    print(f"{'x_up':>8}" + "".join(f"{'Om='+str(o):>14}" for o in (14, 30)))
    curve = {}
    for x in (R3, 3.0, 2.8, 2.6, 2.4, 2.2, 2.0):
        row = [escape_rate(om, x) for om in (14, 30)]
        curve[f"{x:.4f}"] = row
        print(f"{x:>8.4f}" + "".join(f"{r:>14.4e}" for r in row))

    print("\nP3/P4 -- the suspect against §101's measured split, t0 = 2.0")
    cells = [(14, 2), (14, 3), (30, 2)] + ([(30, 3)] if args.deep else [])
    out = []
    print(f"{'Om':>4}{'D':>3}   {'operating points':>26}{'contam/pure pred':>18}"
          f"{'measured':>11}{'pred/meas':>11}")
    for om, D in cells:
        mus, tot, pure, contam = operating_points(om, D, 2.0)
        # legacy=True: this main() reproduces §102's PUBLISHED table (rule 7).
        ks, cp, pp = predict(om, mus, 2.0, legacy=True)
        r_pred, r_meas = cp / pp, contam / pure
        out.append({"omega": om, "D": D, "mus": mus, "ks": ks,
                    "pred_contam": cp, "pred_pure": pp, "pred_ratio": r_pred,
                    "meas_contam": contam, "meas_pure": pure, "meas_ratio": r_meas})
        print(f"{om:>4}{D:>3}   " + ", ".join(f"{m:.4f}" for m in mus).rjust(26)
              + f"{r_pred:>18.4f}{r_meas:>11.4f}{r_pred/r_meas:>11.4f}")

    print("\nP4 -- does the measured split track k_upstream/k_last across Omega?")
    for r in out:
        if r["D"] < 2:
            continue
        ratio_k = r["ks"][0] / r["ks"][-1]
        print(f"  Om={r['omega']:3d} D={r['D']}  k_up/k_last = {ratio_k:.4f}"
              f"   measured contam/pure = {r['meas_ratio']:.4f}")

    args.out.write_text(json.dumps({"curve": curve, "cells": out}, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
