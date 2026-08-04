"""T14-f: is §36's residual the 2-D minimum action, or is it instrument? — the definitive test

§36 measured §15's closed form against the exact collapse with the start ON the slaved
manifold and got pred/meas = 0.9993 / 0.9907 / 0.9822 / 1.0047 across gamma. It claimed
no mechanism for the 0.5-1.8% and named T14-f -- the 2-D minimum action -- as the only
candidate, while noting the residual is NOT one-signed where a path-minimisation
correction must be (a minimum over all paths cannot exceed the value along the slaved
one, so pred/meas >= 1 is required).

**TWO CANDIDATES HAVE SINCE BEEN ELIMINATED.**
  * The realised pool fraction `b` wobbles on the integer lattice exactly as eps does,
    and §36 showed the rate is far MORE sensitive to b than to eps -- so an uncontrolled
    b was the obvious suspect. Measured: the wobble is 0.00-0.53% and adding a b
    regressor moves the rates by under 0.03%. Dead.
  * **The residual is smaller than the fit's own window sensitivity.** Over six Omega
    windows spanning 200-1200 with 7 points each, pred/meas ranges by **1.30-2.29%** per
    gamma -- larger than the 0.5-1.8% being explained. At short windows the residual is
    entirely instrument.

**BUT THAT IS FIXABLE AND THE TEST IS THEREFORE NOT OVER.** §35.2 measured the rate `c`
stable to **0.12-0.19%** across windows when fitted over Omega = 150..2000 with 14
points -- an order of magnitude better than the short windows used above and in §36. So
the definitive test is §36's comparison redone at §35 grade: long windows, many points,
on-manifold start, eps- AND b-controlled.

PREDICTIONS, written before running:

  P1  GATE. At §35-grade windows the rate is determined to ~0.2%, an order better than
      the 1.30-2.29% short-window sensitivity. Verified by splitting each gamma's grid
      into halves and comparing. If it is not, the residual is unresolvable and T14-f is
      closed as unanswerable at any accessible precision.
  P2  THE TEST. The residual then either
        (a) becomes ONE-SIGNED and >= 1, confirming the 2-D minimum action and MEASURING
            it for the first time -- a correction of known sign whose size would be the
            last unexplained piece of the reliability half; or
        (b) stays two-signed at the sub-percent level, meaning §15's closed form is exact
            to better than 1% like-for-like and any 2-D correction is smaller still.
  P3  **(b) is the more likely outcome and is the stronger statement.** §36 already found
      the residual straddling 1, and halving the noise should shrink it rather than
      reveal a sign. Predicting the less interesting outcome in advance is the point.
  P4  If the residual is one-signed and >= 1, it must GROW as the timescale separation
      falls, since a lower separation means the true escape path can deviate further from
      the slaved manifold. sep = 3(1+2g)/(1-2g) runs 3.98 -> 17.0 over the gammas here,
      so the correction should be largest at gamma = 0.20 and smallest at 0.35. A
      one-signed residual with the WRONG ordering would not be the 2-D action.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
import scipy.sparse.linalg as spla

from crnl.cme import enumerate_states, generator
from crnl.networks.am_reversible import am_reversible, delta_star
from experiments.slaving_axis import rate_wkb, slaved

THETA = 0.80


def point(net, ds, x0, bt, omega):
    nb = int(round(bt * omega))
    rest = omega - nb
    d0 = max(1, int(round(x0 * omega)))
    if (rest - d0) % 2:
        d0 -= 1
    n0 = np.array([(rest + d0) // 2, (rest - d0) // 2, nb], dtype=np.int64)
    if min(n0) < 1:
        return None
    thr = max(2, int(round(THETA * ds * omega)))
    states, index = enumerate_states(3, int(omega))
    ab = np.array([abs(int(s[0]) - int(s[1])) >= thr for s in states])
    fav = np.array([int(s[1]) > int(s[0]) for s in states])[ab].astype(float)
    Q = generator(net, int(omega), float(omega))
    tr = np.where(~ab)[0]
    tm = {int(i): r for r, i in enumerate(tr)}
    A = Q[tr][:, tr].tocsr()
    b = -(Q[tr][:, np.where(ab)[0]].tocsr() @ fav)
    p = spla.spsolve(A, b)[tm[index[tuple(n0)]]]
    if not (0.0 < p < 1.0):
        return None
    return {"omega": float(omega), "lp": float(np.log(p)),
            "eps_r": float(d0 / (ds * omega)), "b_r": float(nb / omega)}


def fit(rows):
    om = np.array([r["omega"] for r in rows])
    lp = np.array([r["lp"] for r in rows])
    e = np.array([r["eps_r"] for r in rows])
    b = np.array([r["b_r"] for r in rows])
    cols = [om, np.log(om), om * (e - e.mean())]
    if b.std() > 0:
        cols.append(om * (b - b.mean()))
    cols.append(np.ones_like(om))
    A = np.vstack(cols).T
    c, *_ = np.linalg.lstsq(A, lp, rcond=None)
    res = lp - A @ c
    return float(c[0]), float(1 - res.var() / lp.var())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.20, 0.25, 0.30, 0.35])
    ap.add_argument("--eps-frac", type=float, default=0.35)
    ap.add_argument("--omegas", type=int, nargs="+",
                    default=[150, 300, 450, 600, 750, 900, 1050, 1200, 1400, 1600, 1800])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/manifold_residual.json"))
    args = ap.parse_args()

    t0 = time.time()
    print("§36 redone at §35 grade: long windows, on-manifold start, eps+b controlled")
    print(f"{'gamma':>7}{'sep':>7}{'cells':>7}{'decades':>9}{'measured':>12}"
          f"{'WKB pred':>12}{'pred/meas':>11}{'half-split':>12}{'R^2':>10}")
    out = []
    for g in args.gammas:
        net = am_reversible(g)
        ds = delta_star(g)
        x0 = args.eps_frac * ds
        bt = slaved(net, x0)[2]
        rows = [r for r in (point(net, ds, x0, bt, o) for o in args.omegas) if r]
        if len(rows) < 6:
            print(f"{g:>7.2f}   SKIPPED ({len(rows)} cells)")
            continue
        c, r2 = fit(rows)
        h = len(rows) // 2
        c1, _ = fit(rows[:h + 1])
        c2, _ = fit(rows[h:])
        split = abs(c1 - c2) / abs(c) * 100
        w = rate_wkb(net, x0)
        lp = np.array([r["lp"] for r in rows])
        dec = float((lp.max() - lp.min()) / np.log(10))
        sep = 3 * (1 + 2 * g) / (1 - 2 * g)
        out.append({"gamma": g, "sep": sep, "measured": c, "pred": w,
                    "ratio": w / c, "half_split_pct": split, "r2": r2,
                    "decades": dec, "n": len(rows)})
        print(f"{g:>7.2f}{sep:>7.2f}{len(rows):>7}{dec:>9.1f}{c:>12.6f}"
              f"{w:>12.6f}{w/c:>11.4f}{split:>11.2f}%{r2:>10.6f}")

    print(f"\n=== P1 gate: is the rate determined an order better than the short windows?")
    hs = np.array([r["half_split_pct"] for r in out])
    print(f"  half-split disagreement: {hs.min():.2f}%..{hs.max():.2f}%"
          f"   (short windows gave 1.30-2.29%)")
    print(f"  -> P1 {'HOLDS' if hs.max() < 0.6 else 'FAILS -- residual unresolvable'}")

    print(f"\n=== P2/P3: is the residual one-signed and >= 1?")
    rr = np.array([r["ratio"] for r in out])
    above = int((rr >= 1.0).sum())
    print(f"  ratios: " + "  ".join(f"{x:.4f}" for x in rr))
    print(f"  above 1: {above}/{len(rr)}   mean {rr.mean():.4f}  "
          f"spread {100*(rr.max()-rr.min()):.2f} points")
    if above == len(rr):
        print("  -> P2(a): ONE-SIGNED and >= 1. Consistent with the 2-D minimum action.")
    elif above == 0:
        print("  -> one-signed BELOW 1, which the 2-D minimum action FORBIDS "
              "(a minimum over paths cannot exceed the slaved value).")
    else:
        print("  -> P2(b)/P3: two-signed. No resolvable 2-D correction; the closed form "
              "is exact to this precision.")

    print(f"\n=== P4: if one-signed, does the correction grow as separation falls?")
    for r in sorted(out, key=lambda x: x["sep"]):
        print(f"  gamma={r['gamma']:.2f}  sep={r['sep']:>5.2f}  "
              f"residual {100*(r['ratio']-1):+.2f}%")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
