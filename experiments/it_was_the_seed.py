"""§109 -- §108's central claim was a seeding artifact, and both seed mismatches are now fixed.

RULE 14, on a claim one commit old. §108 argued: stage 2's measured effective escape rate is
0.7880 x k(<x>) at Omega = 14; ln k is convex in the input at every Omega; therefore by Jensen
no average of the rate over the input can land that low; therefore at least two mechanisms are
in play, one of them provably outside the averaging family. The algebra is right. **The
premise is not a stationary escape rate at all.**

THE INSTRUMENT ERROR. §108 measured stage 2 in a chain where stage 2 is seeded as a DELTA at
its rail, while the averages it was compared against are built from stationary laws. In a
SINGLE stage with the upstream PINNED -- no fluctuating input, no averaging of any kind
available -- the delta seed carries a 25x transient:

    P(low at t)/t, Omega = 30, upstream pinned at 3.0319:
      t          0.5        1.0        2.0        4.0        8.0
      delta   1.25e-4    8.92e-4    2.02e-3    2.78e-3    3.15e-3
      QSD     3.11e-3    3.27e-3    3.40e-3    3.47e-3    3.49e-3

At t = 2 the delta seed alone suppresses by 1.68x. §108 compared that against a stationary
average and read the gap as physics.

THE SECOND SEED MISMATCH, which explains §108's other leftover. Stage 1 is seeded from
`stage1_stationary` -- the REFLECTED stage's stationary law -- while the model assumes the FREE
stage's QSD. The reflected law is not depleted near the saddle, so it escapes faster:
1.2460 / 1.2770 / 1.2768 / 1.2493 / 1.2402 across Omega = 14-70, flat. That is §108's
unexplained "flat 21.5% on the contam side", and §108 explicitly said a transient would have
the wrong sign -- true of stage 2's delta, false of stage 1's reflected law, which biases the
other way.

PREDICTIONS, WRITTEN BEFORE THE VERIFICATION RUNS BELOW (the two tables above were the
diagnostics that prompted this file and are disclosed as already computed, rule 2).

  P1  With QSD seeding and a pinned upstream, P(low at t)/t must converge to pi_low * lambda
      -- the two-state forward rate -- and be flat in t. If it still drifts in t the seed is
      not the whole transient.

  P2  With matched seeding, the two-state model must predict stage 1's escape to a FLAT
      few-percent at every Omega, with no drift. That is the claim that closes T-CASC-z.

  P3  §108's premise must FAIL under matched seeding: the measured stage-2 rate must sit AT OR
      ABOVE k(<x>) at every Omega, restoring the averaging family as a live candidate.

  P4  THE DRIFT MUST SURVIVE. Fixing a seed cannot remove a trend that runs across five
      barrier depths; if it does, then §106-§108 were measuring the seed all along and far more
      than §108 has to be withdrawn. I predict the span stays near 2x.

WHAT THIS DOES NOT DO. It does not explain the drift. It withdraws one wrong explanation of it,
identifies two real instrument defects, and closes T-CASC-z.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import experiments.chemical_cascade as cc
from experiments.chain_without_a_joint_solve import chain_operating_points
from experiments.depth_compounding import C, R2, R3, _rates_vec
from experiments.escape_accounts_for_it import escape_rate
from experiments.margin_law import stage1_stationary
from experiments.the_corrected_closure import pi_low
from experiments.two_mechanisms import candidate_averages
from experiments.what_reflection_costs import spectral_gap

OMEGAS = (14, 20, 30, 55, 70)


def _one_stage_generator(om, x_up, first=False, cap_mult=1.25):
    cap = int(np.ceil(cap_mult * R3 * om))
    m = cap + 1
    rows, cols, vals = [], [], []
    for s in range(m):
        tot = 0.0
        l, u = cc.rates_stage(float(s), x_up * om, om, C, R3, first, "hill")
        if s < cap and l > 0:
            rows.append(s); cols.append(s + 1); vals.append(l); tot += l
        if s > 0 and u > 0:
            rows.append(s); cols.append(s - 1); vals.append(u); tot += u
        rows.append(s); cols.append(s); vals.append(-tot)
    return sp.csr_matrix((vals, (rows, cols)), shape=(m, m)), m


def qsd(Q, m, om):
    """Quasi-stationary law above the saddle."""
    n = np.arange(m)
    keep = np.where(n > R2 * om)[0]
    Qs = Q[keep][:, keep].T.tocsc()
    _, v = spla.eigs(Qs, k=1, which="LR")
    p = np.abs(np.real(v[:, 0]))
    return keep, p / p.sum()


def one_stage_low(om, x_up, t, seed, first=False):
    """P(low at t) for a single stage with a pinned upstream, under a chosen seed."""
    Q, m = _one_stage_generator(om, x_up, first)
    n = np.arange(m)
    p = np.zeros(m)
    if seed == "qsd":
        keep, w = qsd(Q, m, om)
        p[keep] = w
    elif seed == "delta":
        p[int(round(R3 * om))] = 1.0
    elif seed == "reflected":
        up, pi = stage1_stationary(om)
        for a, wa in zip(up, pi):
            p[a] = wa
    else:
        raise ValueError(seed)
    p = spla.expm_multiply(Q.T * t, p)
    return float(p[n < R2 * om].sum())


def joint_absorbing_seeded(om, t=2.0, qsd_seed=True, cap_mult=1.25):
    """§108's absorbing-upstream chain, with stage 2 seeded either as a delta or from its QSD."""
    cap = int(np.ceil(cap_mult * R3 * om))
    m = cap + 1
    N = m * m
    idx = np.arange(N)
    n1, n2 = idx // m, idx % m
    sad = R2 * om
    rows, cols, vals = [], [], []
    diag = np.zeros(N)
    l1, u1 = _rates_vec(n1.astype(float), np.zeros(N), om, True)
    alive = n1 >= sad
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
    mus, _ = chain_operating_points(om, 2)
    p = np.zeros(N)
    if qsd_seed:
        Q2, m2 = _one_stage_generator(om, mus[0])
        keep, w2 = qsd(Q2, m2, om)
        for a, wa in zip(up, pi1):
            p[a * m + keep] = wa * w2
    else:
        for a, wa in zip(up, pi1):
            p[a * m + int(round(R3 * om))] = wa
    p = spla.expm_multiply(Q.T * t, p)
    g = p.reshape(m, m)
    return float(g[np.ix_(np.arange(m) >= sad, np.arange(m) < sad)].sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/it_was_the_seed.json"))
    args = ap.parse_args()
    out = {}

    print("P1 -- QSD seeding, pinned upstream: does P(low)/t converge to pi_low * lambda?")
    om = 30
    mus, _ = chain_operating_points(om, 2)
    target = pi_low(om, mus[0]) * escape_rate(om, mus[0])
    print(f"{'t':>7}{'delta':>13}{'QSD':>13}{'pi_low*lambda':>16}")
    p1 = []
    for t in (0.5, 1.0, 2.0, 4.0, 8.0):
        d = one_stage_low(om, mus[0], t, "delta") / t
        q = one_stage_low(om, mus[0], t, "qsd") / t
        p1.append({"t": t, "delta": d, "qsd": q})
        print(f"{t:>7.1f}{d:>13.4e}{q:>13.4e}{target:>16.4e}", flush=True)
    out["p1"] = {"rows": p1, "target": target}
    spread_q = max(r["qsd"] for r in p1) / min(r["qsd"] for r in p1)
    spread_d = max(r["delta"] for r in p1) / min(r["delta"] for r in p1)
    print(f"  spread over a 16x window: QSD {spread_q:.2f}x   delta {spread_d:.1f}x")

    print("\nP2 -- stage 1 with matched seeding, against the two-state model")
    print(f"{'Omega':>7}{'reflected':>13}{'QSD':>13}{'model':>13}{'QSD/model':>11}")
    p2 = []
    for o in OMEGAS:
        r = one_stage_low(o, 0.0, 2.0, "reflected", first=True)
        q = one_stage_low(o, 0.0, 2.0, "qsd", first=True)
        lam, _ = spectral_gap(o, False)
        mod = pi_low(o, R3) * (1 - np.exp(-lam * 2.0))
        p2.append({"omega": o, "reflected": r, "qsd": q, "model": mod,
                   "qsd_over_model": q / mod, "refl_over_qsd": r / q})
        print(f"{o:>7}{r:>13.4e}{q:>13.4e}{mod:>13.4e}{q / mod:>11.4f}", flush=True)
    out["p2"] = p2
    rr = [x["qsd_over_model"] for x in p2]
    print(f"  QSD/model spans {max(rr)/min(rr):.4f}x -- flat closes T-CASC-z")

    print("\nP3/P4 -- §108's premise under matched seeding, and does the drift survive?")
    print(f"{'Omega':>7}{'true/k(<x>) delta':>20}{'QSD':>10}{'bracket top':>14}")
    p3 = []
    for o in OMEGAS:
        pd = joint_absorbing_seeded(o, qsd_seed=False)
        pq = joint_absorbing_seeded(o, qsd_seed=True)
        p2l = pi_low(o, chain_operating_points(o, 2)[0][0])
        rate = lambda v: -np.log(max(1 - v / p2l, 1e-300)) / 2.0
        km, kg, ka = candidate_averages(o)
        p3.append({"omega": o, "delta": rate(pd) / km, "qsd": rate(pq) / km,
                   "bracket_top": ka / km})
        print(f"{o:>7}{rate(pd)/km:>20.4f}{rate(pq)/km:>10.4f}{ka/km:>14.4f}", flush=True)
    out["p3"] = p3
    below = [r["omega"] for r in p3 if r["qsd"] < 0.99]
    qs = [r["qsd"] for r in p3]
    print(f"  cells still below k(<x>) with matched seeding: {below or 'NONE'}"
          f"   -> §108's premise {'FAILS' if not below else 'survives'}")
    print(f"  drift span: delta {max(r['delta'] for r in p3)/min(r['delta'] for r in p3):.3f}x"
          f"   QSD {max(qs)/min(qs):.3f}x   -> P4 "
          f"{'HOLDS' if max(qs)/min(qs) > 1.7 else 'FAILS'}")

    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
