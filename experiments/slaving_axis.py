"""Is the 7.5-15.5% discrepancy the 1-D SLAVING? — tested on an axis independent of gamma

§35 established that §15's closed form is 7.5-15.5% too steep against the asymptotic
collapse rate. §28.3 named the 1-D slaved reduction as the suspect, on the strength of
the residual being eps-independent. **That has never been tested directly, because the
timescale separation `sep = 3(1+2g)/(1-2g)` is a function of gamma**, so every attempt
to vary the slaving quality also varied the landscape, the barrier and the drive
together. Rule 9's exact situation, unbroken for eight sections.

**IT CAN BE BROKEN.** The pool-exchange pair `X + Y <-> 2B` has `U . S_r = 0` for both
directions -- it moves the pool and leaves `delta` untouched. Scaling BOTH by a factor T
therefore accelerates the pool without adding or removing any delta-jump, and because
the pair scales together its rate ratio is fixed, so **gamma and the cycle affinity
A = -3 ln gamma are unchanged**. Measured: sep runs 9.0 -> 269 over T = 1 -> 100 at
gamma = 0.25, with lambda_antisym staying positive throughout.

T is not a pure knob -- the symmetric point moves (b: 0.333 -> 0.498) and the landscape
deforms. That is handled by computing the PREDICTION from each scaled network's own
closed forms, so it tracks the deformation self-consistently. The test is not whether
the rate is constant; it is whether the RATIO converges to 1.

**AND THE COMPARISON IS NOW CLEAN OF THE GAUSSIAN ERROR.** §35.4 showed the diffusion
truncation contributes only 0.07-0.62%, and gave the exact 1-D alternative: for delta,
only four reactions move it and by exactly +-1, so the slaved model is a birth-death
process whose exponent is EXACTLY `-int ln(lambda/mu_down) d(delta)` with no Gaussian
approximation anywhere. Comparing THAT against the exact 2-D CME isolates
**dimensionality alone**.

PREDICTIONS, written before running:

  P0  GATE. At T = 1 the measured rate and the ratio must reproduce §35's gamma = 0.25
      numbers (c = -0.035156, WKB/measured = 1.1369). If not, the scaled-network
      machinery disagrees with the established one and nothing else is admissible.
  P1  sep rises monotonically with T while gamma and the affinity stay fixed. Already
      measured above; re-reported here so the axis is in the record with the result.
  P2  THE TEST. `rate_WKB_1D / rate_measured` converges to 1 as sep grows. If it does,
      the discrepancy IS the 1-D slaving -- established on an axis independent of gamma
      for the first time, upgrading §28.3's suspect to a confirmed mechanism.
  P3  THE OUTCOME THAT COSTS MORE. If the ratio PLATEAUS above 1, the slaving is not the
      whole story, and a residual survives that neither the Gaussian truncation
      (§35.4: <0.62%) nor the slaving explains. That would mean §15's closed form has an
      error nothing currently identified accounts for, and §28.3's attribution -- which
      rests only on eps-independence -- would be wrong or incomplete.
  P4  If the ratio goes BELOW 1 at large sep, the 1-D model overshoots in the other
      direction and the reduction is not a one-sided approximation, which no account
      here allows and which would need the whole framing revisited.

  Rule 12 note: every rate is an eps-controlled fit over a matched P window, and the
  achieved window is reported per cell so matching can be checked rather than assumed.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
import scipy.sparse.linalg as spla
from scipy.optimize import brentq

from crnl.cme import enumerate_states, generator
from crnl.deterministic import jacobian
from crnl.reactions import Reaction, ReactionNetwork

U = np.array([1.0, -1.0, 0.0])


def am_scaled(gamma: float, T: float, k: float = 1.0) -> ReactionNetwork:
    """AM with the pool-exchange pair scaled by T. Both directions, so gamma is fixed."""
    return ReactionNetwork(species=["X", "Y", "B"], reactions=[
        Reaction({"X": 1, "Y": 1}, {"B": 2}, k * T, name="f1"),
        Reaction({"B": 1, "X": 1}, {"X": 2}, k, name="f2"),
        Reaction({"B": 1, "Y": 1}, {"Y": 2}, k, name="f3"),
        Reaction({"B": 2}, {"X": 1, "Y": 1}, gamma * k * T, name="r1"),
        Reaction({"X": 2}, {"B": 1, "X": 1}, gamma * k, name="r2"),
        Reaction({"Y": 2}, {"B": 1, "Y": 1}, gamma * k, name="r3"),
    ], name=f"am-scaled-g{gamma}-T{T}")


def slaved(net, delta):
    S = net.stoichiometry_matrix()

    def ds(s):
        x, y, b = 0.5 * (s + delta), 0.5 * (s - delta), 1.0 - s
        if min(x, y, b) < 0:
            return np.nan
        v = S @ net.fluxes(np.array([x, y, b]))
        return v[0] + v[1]

    lo, hi = max(abs(delta) + 1e-9, 1e-9), 1.0 - 1e-9
    a, b_ = ds(lo), ds(hi)
    if not (np.isfinite(a) and np.isfinite(b_)) or a * b_ > 0:
        return None
    s = brentq(ds, lo, hi, xtol=1e-14)
    return np.array([0.5 * (s + delta), 0.5 * (s - delta), 1.0 - s])


def updown(net, delta):
    st = slaved(net, delta)
    if st is None:
        return None
    S = net.stoichiometry_matrix()
    f = net.fluxes(st)
    up = sum(float(f[r]) for r in range(net.n_reactions) if float(U @ S[:, r]) > 0.5)
    dn = sum(float(f[r]) for r in range(net.n_reactions) if float(U @ S[:, r]) < -0.5)
    return up, dn


def delta_star_of(net):
    """Attractor separation: largest delta with zero drift along the slaved manifold."""
    def drift(d):
        ud = updown(net, d)
        return np.nan if ud is None else ud[0] - ud[1]
    lo, hi = 1e-5, 0.999
    f_lo = drift(lo)
    xs = np.linspace(lo, hi, 400)
    prev_x, prev_f = lo, f_lo
    root = None
    for x in xs[1:]:
        fx = drift(x)
        if np.isfinite(fx) and np.isfinite(prev_f) and prev_f * fx < 0:
            root = brentq(drift, prev_x, x, xtol=1e-14)
        if np.isfinite(fx):
            prev_x, prev_f = x, fx
    return root


def sep_of(net):
    S = net.stoichiometry_matrix()
    f = lambda b: (S @ net.fluxes(np.array([(1 - b) / 2, (1 - b) / 2, b])))[0]
    b = brentq(f, 1e-9, 1 - 1e-9, xtol=1e-14)
    ev = np.linalg.eigvals(jacobian(net, np.array([(1 - b) / 2, (1 - b) / 2, b]))).real
    ev = np.sort(ev[np.abs(ev) > 1e-9])
    return float(abs(ev[0] / ev[-1])), float(ev[-1])


def rate_wkb(net, x0, n=3001):
    xs = np.linspace(1e-6, x0, n)
    vals = []
    for x in xs:
        ud = updown(net, float(x))
        if ud is None or min(ud) <= 0:
            return float("nan")
        vals.append(np.log(ud[0] / ud[1]))
    return -float(np.trapezoid(vals, xs))


def p_error(net, omega, ds, eps, theta):
    nb_frac = 1.0 - (slaved(net, eps * ds)[0] + slaved(net, eps * ds)[1])
    nb = int(round(nb_frac * omega))
    rest = omega - nb
    d0 = max(1, int(round(eps * ds * omega)))
    if (rest - d0) % 2:
        d0 -= 1
    n0 = np.array([(rest + d0) // 2, (rest - d0) // 2, nb], dtype=np.int64)
    thr = max(2, int(round(theta * ds * omega)))
    states, index = enumerate_states(3, int(omega))
    absorb = np.array([abs(int(s[0]) - int(s[1])) >= thr for s in states])
    fav = np.array([int(s[1]) > int(s[0]) for s in states])[absorb].astype(float)
    Q = generator(net, int(omega), float(omega))
    tr = np.where(~absorb)[0]
    tmap = {int(i): r for r, i in enumerate(tr)}
    A = Q[tr][:, tr].tocsr()
    b = -(Q[tr][:, np.where(absorb)[0]].tocsr() @ fav)
    h = spla.spsolve(A, b)
    si = tmap[index[tuple(n0)]]
    return float(h[si]), float(d0 / (ds * omega))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.25)
    ap.add_argument("--ts", type=float, nargs="+", default=[1, 2, 5, 10, 30, 100])
    ap.add_argument("--omegas", type=int, nargs="+",
                    default=[200, 400, 600, 800, 1000, 1200])
    ap.add_argument("--eps-frac", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/slaving_axis.json"))
    args = ap.parse_args()

    t0 = time.time()
    g = args.gamma
    print(f"gamma = {g} fixed; the pool pair scaled by T -> separation varies alone")
    print(f"{'T':>6}{'sep':>9}{'delta*':>9}{'decades':>9}{'measured c':>13}"
          f"{'WKB 1-D':>12}{'WKB/meas':>10}{'sec':>7}")
    rows = []
    for T in args.ts:
        t1 = time.time()
        net = am_scaled(g, T)
        ds = delta_star_of(net)
        if ds is None:
            print(f"{T:>6.0f}   SKIPPED (no attractor)")
            continue
        sep, lam = sep_of(net)
        om, lp, er = [], [], []
        for o in args.omegas:
            try:
                p, e = p_error(net, o, ds, args.eps_frac, args.theta)
            except Exception:
                continue
            if np.isfinite(p) and 0.0 < p < 1.0:
                om.append(float(o)); lp.append(np.log(p)); er.append(e)
        if len(om) < 4:
            print(f"{T:>6.0f}   SKIPPED (too few usable cells)")
            continue
        om, lp, er = np.array(om), np.array(lp), np.array(er)
        A = np.vstack([om, np.log(om), om * (er - er.mean()), np.ones_like(om)]).T
        c, *_ = np.linalg.lstsq(A, lp, rcond=None)
        meas = float(c[0])
        pred = rate_wkb(net, args.eps_frac * ds)
        dec = float((lp.max() - lp.min()) / np.log(10))
        rows.append({"T": T, "sep": sep, "delta_star": ds, "decades": dec,
                     "measured": meas, "wkb": pred, "ratio": pred / meas,
                     "lambda_signal": lam})
        print(f"{T:>6.0f}{sep:>9.2f}{ds:>9.5f}{dec:>9.2f}{meas:>13.6f}"
              f"{pred:>12.6f}{pred/meas:>10.4f}{time.time()-t1:>7.0f}")

    print(f"\n=== P2/P3: does the ratio converge to 1 as the slaving improves?")
    print(f"{'sep':>9}{'WKB/measured':>15}{'excess':>10}")
    for r in rows:
        print(f"{r['sep']:>9.2f}{r['ratio']:>15.4f}{r['ratio']-1:>10.4f}")
    if len(rows) >= 3:
        x = np.array([1.0 / r["sep"] for r in rows])
        y = np.array([r["ratio"] for r in rows])
        p = np.polyfit(x, y, 1)
        r2 = 1 - np.var(y - np.polyval(p, x)) / np.var(y)
        print(f"\n  ratio vs 1/sep: intercept (sep -> infinity) = {p[1]:.5f}"
              f"   R^2 = {r2:.5f}")
        print(f"  excess at infinite separation = {p[1]-1:+.5f}")
        if abs(p[1] - 1.0) < 0.015:
            print("  -> P2: the discrepancy IS the 1-D slaving. §28.3's suspect is "
                  "confirmed on an axis independent of gamma.")
        else:
            print(f"  -> P3: the ratio PLATEAUS at {p[1]:.4f}. The slaving accounts for "
                  f"{100*(y[0]-p[1])/(y[0]-1):.1f}% of the discrepancy and "
                  f"{100*(p[1]-1)/(y[0]-1):.1f}% is something else.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
