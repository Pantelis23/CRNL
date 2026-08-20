"""§101 (T-CASC-n) -- the depth ceiling when upstream stages are allowed to fail.

§100 found the channel §91-§98 never priced. Everything from §92 on reflects every stage
but the last at its saddle, so no upstream stage can fail; the penalty measured there is
the transfer of the upstream's RAIL FLUCTUATIONS. §100's refuted prediction is why this
matters: the Hill coupling's saturation attenuates a noisy upstream and does NOT attenuate
a failed one -- the downstream follows a fallen upstream at 0.73 rising to 0.98 across the
windows measured. A chain with free upstream stages therefore carries a second error
channel that transmits at near unity, and `D_max = c*/(penalty x eps)` was built without it.

At D = 2 that channel is worth 12.2% at t0 = 2 and 96% at t0 = 0.5 (§100.2). The question
here is what it does with DEPTH, which is the only thing that matters for a ceiling.

PREDICTIONS, WRITTEN BEFORE RUNNING.

  P1  WIRING, TWO WAYS, both exact. The generalised builder below must reproduce
      (a) §94's `depth_compounding.build_chain(om, D)` at n_reflected = D-1 and at
          n_reflected = D (its all_reflected=True), and
      (b) §100's independently written `what_reflection_costs.build_free` at D = 2,
          n_reflected = 0.
      Max |dQ| = 0 in every case. Two implementations written for different sections have
      to agree before anything below is worth reading.

  P2  THE MEASUREMENT. free/reflected error ratio at D = 2 and D = 3, swept over the
      window (rule 18 -- §100.2 was caught summarising one cell of exactly this trend).
      I predict the ratio GROWS with depth at every window, because D = 3 has two upstream
      stages that can fail where D = 2 has one. If the failure channels were independent
      and each contributed as at D = 2, then ratio(D=3) - 1 = 2 x (ratio(D=2) - 1), i.e.
      1.244 at t0 = 2. I predict LESS than that, because the wall's +49% inflation of the
      surviving branch (§100 P2) also compounds with depth and pushes the other way. Both
      numbers are reported; no threshold is attached to either (rule 20).

  P3  WHICH CHANNEL DOMINATES. Decompose free D = 3's P(stage 3 low) into
        - P(s3 lo, s1 hi, s2 hi)  -- pure fluctuation transfer, the thing §92-§98 measured
        - the remainder            -- failure-contaminated
      I predict fluctuation transfer still DOMINATES at t0 = 2, because §94's reflected
      penalty at D = 3 is 6.35 eps_iso while accumulated upstream failure is of order
      2 eps_iso -- and that the failure share RISES as the window shortens, following
      §100.2's ratio trend. Shares are reported at every window; either outcome prints.

  P4  THE DEPTH TREND, and what it does NOT support. Penalty sequences for both
      constructions at D = 2, 3 (Omega = 30, the published operating point) and D = 2, 3, 4
      (Omega = 14, where D = 4 fits). **No D_max is extrapolated from two or three points.**
      §90 withdrew a quoted precision for exactly this reason and rule 15 demands every
      candidate extrapolation or none. What is reported is the per-added-stage growth factor
      in each construction and whether the two differ.

SCOPE, stated up front. Omega = 14 is A*Omega = 2.66 -- the shallowest cell in §98's table
and the one whose penalty position read 1.3089, outside the two-limit bracket (§99(a)). Its
depth trend is a trend at a shallow barrier and is NOT transported to Omega = 30.

RESOURCE NOTE, stated rather than silent. A free D = 3 chain at Omega = 30 is 1,771,561
joint states and one window costs ~10 minutes, so Omega = 30 is run at THREE windows
(0.5, 2.0, 8.0, spanning the same 16x as §100.2) rather than five. Omega = 14 is cheap and
gets all five. No cell is dropped after being seen; the cut is by cost, decided here.

WHAT THIS CANNOT SETTLE. The seed is the same in both arms -- stage 1 from its
quasi-stationary law above the saddle, later stages at their rails -- so both arms start
CORRECT and this measures degradation from a correct start. A chain seeded from the true
joint stationary law would start with some stages already failed; that is a different
question and it is not asked here.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import experiments.chemical_cascade as cc
from experiments.depth_compounding import C, R2, R3, _rates_vec
from experiments.margin_law import stage1_stationary
from experiments.timescale_ratio import pinned_reference

WINDOWS = (0.5, 1.0, 2.0, 4.0, 8.0)


def build_gen(om, D, n_reflected, cap_mult=1.25):
    """D stages; the FIRST `n_reflected` of them are walled at their saddle, the rest free.

    n_reflected = D-1 is §94's construction, D is its all_reflected=True, 0 is a fully free
    chain. P1 checks all three against the implementations already in the tree.
    """
    if not 0 <= n_reflected <= D:
        raise ValueError(f"n_reflected must be in [0, {D}], got {n_reflected}")
    cap = int(np.ceil(cap_mult * R3 * om))
    nsad = int(np.ceil(R2 * om))
    ref = np.arange(nsad, cap + 1)
    walled = [i < n_reflected for i in range(D)]
    dims = [len(ref) if w else cap + 1 for w in walled]
    strides = [int(np.prod(dims[i + 1:])) for i in range(D)]
    N = int(np.prod(dims))
    idx = np.arange(N)
    counts = []
    for i in range(D):
        ni = (idx // strides[i]) % dims[i]
        counts.append(ref[ni].astype(float) if walled[i] else ni.astype(float))
    rows, cols, vals = [], [], []
    diag = np.zeros(N)
    for i in range(D):
        up = counts[i - 1] if i else np.zeros(N)
        lam, mu = _rates_vec(counts[i], up, om, i == 0)
        lo = float(ref[0]) if walled[i] else 0.0
        up_ok = (counts[i] < float(cap)) & (lam > 0)
        dn_ok = (counts[i] > lo) & (mu > 0)
        rows.append(idx[up_ok]); cols.append(idx[up_ok] + strides[i]); vals.append(lam[up_ok])
        rows.append(idx[dn_ok]); cols.append(idx[dn_ok] - strides[i]); vals.append(mu[dn_ok])
        diag -= np.where(up_ok, lam, 0.0) + np.where(dn_ok, mu, 0.0)
    rows.append(idx); cols.append(idx); vals.append(diag)
    Q = sp.csr_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))
    return Q, ref, dims, strides, cap, walled


def seed_gen(om, ref, dims, strides, walled, pi1):
    """Stage 1 from its quasi-stationary law above the saddle; later stages at their rails.

    §96.1a's bug was indexing the rail by the wrong convention for the stage type; here the
    convention is read off `walled` per stage rather than inferred from position.
    """
    p = np.zeros(int(np.prod(dims)))
    n_hi = int(round(R3 * om))
    pos = list(ref).index(n_hi)
    rest = sum(strides[i] * (pos if walled[i] else n_hi) for i in range(1, len(dims)))
    for a, w in enumerate(pi1):
        # pi1 is indexed by position above the saddle; a walled stage uses that position
        # directly, a free stage is indexed by its raw count ref[a].
        p[(a if walled[0] else int(ref[a])) * strides[0] + rest] = w
    return p


def last_low(p, om, dims, strides, walled, ref):
    idx = np.arange(len(p))
    n = (idx // strides[-1]) % dims[-1]
    counts = ref[n] if walled[-1] else n
    return float(p[counts < R2 * om].sum())


def solve(om, D, n_reflected, t):
    Q, ref, dims, strides, cap, walled = build_gen(om, D, n_reflected)
    _, pi1 = stage1_stationary(om)
    p = seed_gen(om, ref, dims, strides, walled, pi1)
    p = spla.expm_multiply(Q.T * t, p)
    return p, ref, dims, strides, walled


def channel_split(p, om, dims, strides, ref, walled):
    """Free chain only: split P(last low) by whether every upstream stage ended high."""
    idx = np.arange(len(p))
    D = len(dims)
    last = (idx // strides[-1]) % dims[-1]
    lo_last = (ref[last] if walled[-1] else last) < R2 * om
    up_hi = np.ones(len(p), bool)
    for k in range(D - 1):
        nk = (idx // strides[k]) % dims[k]
        ck = ref[nk] if walled[k] else nk
        up_hi &= ck >= R2 * om
    total = float(p[lo_last].sum())
    pure = float(p[lo_last & up_hi].sum())
    return total, pure, total - pure


def run(om, depths, windows=WINDOWS, out=None, sink=None):
    out = [] if out is None else out
    for D in depths:
        for t in windows:
            pf, ref, dims, strides, walled = solve(om, D, 0, t)
            free = last_low(pf, om, dims, strides, walled, ref)
            tot, pure, contam = channel_split(pf, om, dims, strides, ref, walled)
            pr, refr, dimr, strr, walr = solve(om, D, D - 1, t)
            refl = last_low(pr, om, dimr, strr, walr, refr)
            eps = pinned_reference(om, 1.0, t)
            out.append({"omega": om, "D": D, "t": t, "free": free, "refl": refl,
                        "ratio": free / refl, "pen_free": free / eps,
                        "pen_refl": refl / eps, "eps_iso": eps,
                        "pure": pure, "contam": contam,
                        "contam_share": contam / tot if tot > 0 else float("nan")})
            r = out[-1]
            print(f"  Om={om:3d} D={D} t={t:>4.1f}  free {free:.4e}  refl {refl:.4e}"
                  f"  ratio {r['ratio']:.4f}   pen {r['pen_free']:7.3f}/{r['pen_refl']:7.3f}"
                  f"   contam {100*r['contam_share']:5.1f}%", flush=True)
            if sink is not None:                       # incremental: a kill loses one cell
                sink.write_text(json.dumps(out, indent=2, default=float))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/free_upstream_depth.json"))
    ap.add_argument("--deep", action="store_true", help="add Omega=14 D=4 (10.6M states)")
    args = ap.parse_args()

    rows = []
    print("P4 -- Omega = 14, shallow barrier (A*Omega = 2.66), all five windows")
    run(14, (2, 3, 4) if args.deep else (2, 3), WINDOWS, rows, args.out)
    print("\nP2/P3 -- Omega = 30, the published operating point, three windows (see"
          " RESOURCE NOTE)")
    run(30, (2, 3), (0.5, 2.0, 8.0), rows, args.out)

    print("\nP2 summary: does the free/reflected ratio grow with depth?")
    print(f"{'Omega':>6}{'t':>6}{'D=2':>9}{'D=3':>9}{'D=4':>9}{'2x(D2-1)+1':>12}")
    for om in sorted({r["omega"] for r in rows}):
        for t in WINDOWS:
            cell = {r["D"]: r["ratio"] for r in rows if r["omega"] == om and r["t"] == t}
            if 2 not in cell:
                continue
            f = lambda d: f"{cell[d]:.4f}" if d in cell else "--"
            print(f"{om:>6}{t:>6.1f}{f(2):>9}{f(3):>9}{f(4):>9}"
                  f"{2*(cell[2]-1)+1:>12.4f}")

    print("\nP4 summary: penalty growth per added stage")
    for om in sorted({r["omega"] for r in rows}):
        for t in WINDOWS:
            cs = sorted([r for r in rows if r["omega"] == om and r["t"] == t],
                        key=lambda r: r["D"])
            if len(cs) < 2:
                continue
            gf = [cs[i + 1]["pen_free"] / cs[i]["pen_free"] for i in range(len(cs) - 1)]
            gr = [cs[i + 1]["pen_refl"] / cs[i]["pen_refl"] for i in range(len(cs) - 1)]
            print(f"  Om={om:3d} t={t:>4.1f}  free x" + ", x".join(f"{g:.3f}" for g in gf)
                  + "   reflected x" + ", x".join(f"{g:.3f}" for g in gr))
    print("\n  (no D_max is extrapolated from these -- see P4 in the docstring)")

    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
