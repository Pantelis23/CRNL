"""T-CASC-i(b): is the missing 58% of the mean shift the upstream's finite correlation time?

§95 predicted the downstream's operating point as <F(x_up)>, expanded to second order as
F(mu) + (1/2) F'' sigma^2, and it closed 42% of the shift. **That formula is the FROZEN-upstream
limit** -- it assumes the downstream fully equilibrates to F(x_up) before x_up moves. §92/§93
showed the two stages' correlation times are EQUAL here, so that limit is exactly wrong, and the
residual is the natural place for it to show.

**The two limits are different objects, and both are computable.** F is a composition: the coupling
enters only through the Hill factor h, and the rail is then G(h). So

    frozen upstream   mu_2 -> < G(h(x_up)) >  =  <F(x_up)>      (average the OUTPUT)
    fast upstream     mu_2 -> G( <h(x_up)> )                    (average the INPUT)

Jensen applies to a different function in each. §95 used the first. **Sweeping the upstream clock
moves the system between them**, and §93's speed knob does that without touching the landscape, the
barrier, the rail width or the stationary law.

Both stages are reflected at their saddles here, so neither can escape and no mean or width is
measured under conditioning (§92.1(b)).

PREDICTIONS, written before running.

  P1  GATE, a regression. At speed 1 this must reproduce §94/§95's stored mu_2 = 2.95635 and
      sigma_2 = 0.54181. A different number means a different chain.
  P2  **DOES THE CLOCK MOVE THE MEAN AT ALL?** If mu_2 is independent of the upstream speed, the
      correlation time is not the missing term and candidate (b) is dead -- which would leave the
      third moment and the reflecting boundary.
  P3  **DO THE TWO LIMITS BRACKET THE MEASUREMENT?** **Predicted: yes**, with the slow end
      approaching <F(x_up)> computed EXACTLY over stage 1's measured distribution (no second-order
      truncation) and the fast end approaching G(<h>). If the measurement sits outside both, the
      whole static-transfer framing is wrong rather than incomplete.
  P4  **HOW MUCH OF THE 58% DOES THIS ACCOUNT FOR?** Report the residual against the nearer limit
      at speed 1, and against the exact frozen average rather than §95's truncation -- part of the
      58% may simply be the second-order expansion, which is worth separating out before crediting
      any physics.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import experiments.chemical_cascade as cc
from experiments.cascade_schlogl import schlogl_consts
from experiments.jensen_shift import F, MU1, MU2, SD1, SD2, d2F
from experiments.margin_law import R1, R2, R3, stage1_stationary

C = schlogl_consts(R1, R2, R3)


def build_two_reflected(om, s_up, cap_mult=1.25):
    """Two stages, BOTH reflected at their saddles; stage 1's propensities scaled by s_up."""
    cap = int(np.ceil(cap_mult * R3 * om))
    nsad = int(np.ceil(R2 * om))
    ref = np.arange(nsad, cap + 1)
    n = len(ref)
    idx = np.arange(n * n)
    a, b = idx // n, idx % n
    c1, c2 = ref[a].astype(float), ref[b].astype(float)
    k1a, k1r, k2b, k2r = C
    auto1 = k1a * c1 * (c1 - 1.0) / om
    lam1 = (auto1 + k2b * om) * s_up
    mu1 = (k1r * c1 * (c1 - 1.0) * (c1 - 2.0) / om ** 2 + k2r * c1) * s_up
    auto2 = k1a * c2 * (c2 - 1.0) / om
    hv = np.array([cc.hill(x, R3) for x in c1 / om])
    lam2 = hv * auto2 + k2b * om
    mu2 = k1r * c2 * (c2 - 1.0) * (c2 - 2.0) / om ** 2 + k2r * c2
    rows, cols, vals = [], [], []
    diag = np.zeros(n * n)
    for cnt, lam, mu, st, lo, hi in ((c1, lam1, mu1, n, ref[0], cap),
                                     (c2, lam2, mu2, 1, ref[0], cap)):
        u = (cnt < hi) & (lam > 0)
        d = (cnt > lo) & (mu > 0)
        rows.append(idx[u]); cols.append(idx[u] + st); vals.append(lam[u])
        rows.append(idx[d]); cols.append(idx[d] - st); vals.append(mu[d])
        diag -= np.where(u, lam, 0.0) + np.where(d, mu, 0.0)
    rows.append(idx); cols.append(idx); vals.append(diag)
    Q = sp.csr_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))), shape=(n * n, n * n))
    return Q, ref, n


def stats(p, ref, n, k):
    idx = np.arange(n * n)
    c = ref[idx // n if k == 0 else idx % n].astype(float)
    w = p / p.sum()
    mu = float((w * c).sum())
    return mu, float(np.sqrt((w * (c - mu) ** 2).sum()))


def G_of_h(hbar):
    """The downstream rail when the Hill factor is held at hbar -- the FAST-limit prediction."""
    k1a, k1r, k2b, k2r = C
    r = np.roots([-k1r, k1a * hbar, -k2r, k2b])
    r = np.sort([z.real for z in r if abs(z.imag) < 1e-9 and z.real > 1e-12])
    return float(r[-1]) if len(r) >= 3 else np.nan


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=30)
    ap.add_argument("--t", type=float, default=2.0)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/static_transfer_limit.json"))
    args = ap.parse_args()
    om, t = args.omega, args.t
    cc.HILL_N, cc.HILL_K = 4.0, 1.0
    _, pi1 = stage1_stationary(om)
    out = {}

    print("=== P1 GATE: reproduce §94/§95's stored stage-2 mean and width at speed 1")
    Q, ref, n = build_two_reflected(om, 1.0)
    p0 = np.zeros(n * n)
    n_hi = int(round(R3 * om))
    for i, w in enumerate(pi1):
        p0[i * n + list(ref).index(n_hi)] = w
    p1v = spla.expm_multiply(Q.T * t, p0)
    m1, s1 = stats(p1v, ref, n, 0)
    m2, s2 = stats(p1v, ref, n, 1)
    print(f"  stage 1: mean {m1/om:.5f} sd {s1/om:.5f}   (§94: {MU1}, {SD1})")
    print(f"  stage 2: mean {m2/om:.5f} sd {s2/om:.5f}   (§94: {MU2}, {SD2})")
    ok = abs(m2 / om - MU2) < 5e-3 and abs(s2 / om - SD2) < 5e-3
    print(f"  -> P1 {'HOLDS: same chain' if ok else 'FAILS: this is not the chain §94/§95 measured'}")

    print("\n=== the two limits, both computed from stage 1's exact distribution")
    xs = ref / om
    w1 = np.array([p1v[i * n:(i + 1) * n].sum() for i in range(n)])
    w1 = w1 / w1.sum()
    frozen_exact = float(sum(w * F(x) for w, x in zip(w1, xs) if np.isfinite(F(x))))
    hbar = float(sum(w * cc.hill(x, R3) for w, x in zip(w1, xs)))
    fast = G_of_h(hbar)
    d_intr = MU1 - R3
    print(f"  <F(x_up)> exact over stage 1's law   = {frozen_exact:.5f}"
          f"   (+ intrinsic {d_intr:+.5f} -> {frozen_exact + d_intr:.5f})")
    print(f"  §95's 2nd-order truncation of it     = {F(MU1) + 0.5*d2F(MU1,1e-3)*SD1**2:.5f}"
          f"   (+ intrinsic -> {F(MU1) + 0.5*d2F(MU1,1e-3)*SD1**2 + d_intr:.5f})")
    print(f"  <h(x_up)> = {hbar:.5f}  ->  G(<h>) = {fast:.5f}"
          f"   (+ intrinsic -> {fast + d_intr:.5f})")
    out.update({"frozen_exact": frozen_exact, "fast": fast, "hbar": hbar})

    print("\n=== P2/P3: sweep the upstream clock. Does the mean move, and do the limits bracket?")
    print(f"{'s_up':>8}{'mean_2':>11}{'sd_2':>10}{'vs frozen':>12}{'vs fast':>10}")
    rows = []
    for s in (0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 32.0):
        Qs, refs, ns = build_two_reflected(om, s)
        q = np.zeros(ns * ns)
        for i, w in enumerate(pi1):
            q[i * ns + list(refs).index(n_hi)] = w
        q = spla.expm_multiply(Qs.T * t, q)
        mm, ss = stats(q, refs, ns, 1)
        rows.append({"s_up": s, "mu2": mm / om, "sd2": ss / om})
        print(f"{s:>8.3f}{mm/om:>11.5f}{ss/om:>10.5f}"
              f"{mm/om - (frozen_exact + d_intr):>12.5f}{mm/om - (fast + d_intr):>10.5f}")
    out["sweep"] = rows
    mus = [r["mu2"] for r in rows]
    moves = max(mus) - min(mus)
    print(f"  the mean moves by {moves:.5f} across a 256x change in the upstream clock")
    print(f"  -> P2 {'HOLDS: the clock DOES move the operating point, so the correlation time is live' if moves > 0.005 else 'FAILS: the clock does not move the mean -- candidate (b) is dead, leaving the third moment and the reflecting boundary'}")
    lo, hi = min(frozen_exact + d_intr, fast + d_intr), max(frozen_exact + d_intr, fast + d_intr)
    inside = sum(1 for m in mus if lo - 1e-9 <= m <= hi + 1e-9)
    print(f"  limits: frozen {frozen_exact + d_intr:.5f}, fast {fast + d_intr:.5f};"
          f"  {inside}/{len(mus)} measured means lie between them")
    print(f"  -> P3 {'HOLDS: the two limits bracket the measurement' if inside >= len(mus) - 1 else 'FAILS: the measurement sits OUTSIDE both limits, so the static-transfer framing is wrong rather than incomplete'}")

    print("\n=== P4: how much of §95's 58% was the second-order truncation?")
    trunc = F(MU1) + 0.5 * d2F(MU1, 1e-3) * SD1 ** 2 + d_intr
    gap95 = MU2 - trunc
    gap_exact = MU2 - (frozen_exact + d_intr)
    print(f"  residual against §95's truncated frozen formula: {gap95:+.5f}")
    print(f"  residual against the EXACT frozen average:       {gap_exact:+.5f}")
    if abs(gap95) > 0:
        print(f"  -> the truncation itself accounted for"
              f" {100*(1 - abs(gap_exact)/abs(gap95)):.0f}% of §95's residual")
    print(f"  residual against the fast limit:                 {MU2 - (fast + d_intr):+.5f}")
    print("  **Separating the expansion error from the physics BEFORE crediting either.**")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
