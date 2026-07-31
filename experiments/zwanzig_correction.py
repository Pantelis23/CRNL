"""§24.1b: is the cost of naive discarding exactly the term a proper projection supplies?

§24's arms zero a noise component and leave the drift alone. That is NOT a valid
reduced model: the rigorous construction (Zwanzig projection, arXiv:2512.03706)
returns a reduced system with a GENERALIZED FORCE and a STATE-DEPENDENT DIFFUSION --
eliminated degrees of freedom come back as modified drift and memory, not as
deletion. So §24 and §25 measure what NAIVE discarding costs, which is what a
practitioner does, and not what a correct reduction achieves.

This tests the difference, and the target is §24.1a's residual: `delta-only` -- which
deletes the pool's noise and keeps the pool as a variable -- overshoots P(error) by
+13.2% at gamma = 0.05 falling to +0.4% at gamma = 0.45, monotonically with the
timescale separation. That is the signature of a subleading term, and a proper
projection should name it.

WHICH TERM, and the reasoning matters because it rules out the obvious guess. The
`delta-only` projection preserves `a - b` exactly, so **delta's own diffusion is
UNCHANGED** -- the missing piece cannot be delta-diffusion. And the residual has a
sign: deleting the pool's noise makes errors MORE likely, so the pool's fluctuations
must be *strengthening* the restoring drift. Since `delta`'s drift depends on `n_B`
and is nonlinear in it, fluctuations of `s` shift the mean restoring force by the
curvature term

    generalized force:   1/2 * d2(b_delta)/ds2 * Var(s),      Var(s) = D_ss / (2 lam_s)

which is exactly the kind of object a Zwanzig projection produces. A second candidate
is the effective diffusion `s`-fluctuations induce on `delta` through the drift,

    induced diffusion:   D_extra = (d(b_delta)/ds)^2 * D_ss / lam_s^2

which acts in the OPPOSITE direction (more noise, more errors) and so cannot by
itself explain an overshoot. Both are computed per state by finite differences along
the pool direction `u = (0.5, 0.5, -1)` -- which raises `s` by one and leaves `delta`
untouched, conserving the total -- and both are run as separate arms, because I do
not know which dominates and reporting only the flattering one is what rule 15 exists
to prevent.

    delta-only          the §24.1 baseline: pool noise deleted, drift untouched
    +force              plus the curvature term in the drift
    +diffusion          plus the induced diffusion
    +both               the assembled reduction
    full                harness control

PREDICTIONS, written before running:

  P1  `+force` shrinks the residual, and shrinks it MOST where it is largest -- at
      small gamma, where the timescale separation is weakest. If the residual really
      is the leading neglected term, `+force` should take +13.2% at gamma = 0.05 down
      to a few percent.
  P2  `+diffusion` makes the residual WORSE, because it adds noise where the arm is
      already over-predicting errors. It is run anyway: it is the naive guess at
      "what did I delete", and showing it has the wrong sign is worth more than
      omitting it.
  P3  `+both` sits between them and is not the best arm. If it were, the two terms
      would be cancelling and neither would be identified.
  P4  WHAT WOULD REFUTE THE READING. If `+force` does not shrink the residual, or
      shrinks it uniformly across gamma rather than preferentially at small gamma,
      then §24.1a's identification of the residual as subleading transverse coupling
      is wrong, and the naive-vs-correct distinction is not what separates these
      models. The honest fallback would be that a proper projection needs the full
      memory kernel, not a Markovian curvature correction -- which this experiment
      cannot supply.
  P5  This does NOT test whether a proper projection rescues the CATEGORICAL failure
      (`s-only` reporting exactly 0). It cannot: `s-only` deletes the noise in the
      coordinate the observable is defined on, and no generalized force on a
      deterministic coordinate produces barrier crossings. Stated here so the scope
      of any positive result is not overread.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.approximations import propensities_batch
from crnl.networks.am_reversible import am_reversible
from crnl.vectorized import compile_network
from experiments.approximation_hierarchy import _setup, p_cme

ARMS = ("full", "delta-only", "+force", "+diffusion", "+both")
POOL_DIR = np.array([0.5, 0.5, -1.0])      # raises s by 1, leaves delta alone
H = 1.0


def _drift(comp, n, S):
    return propensities_batch(comp, n) @ S.T


def pool_terms(comp, n, S):
    """Curvature force and induced diffusion from eliminating the pool, per state.

    Finite differences along POOL_DIR. Everything is evaluated at the current state
    rather than at a nullcline, so the terms are state-dependent as a Zwanzig
    reduction requires -- a single constant would be a different (and weaker) model.
    """
    up = np.clip(n + H * POOL_DIR, 0.0, None)
    dn = np.clip(n - H * POOL_DIR, 0.0, None)
    b0, bu, bd = _drift(comp, n, S), _drift(comp, up, S), _drift(comp, dn, S)

    bd_delta = lambda b: b[:, 0] - b[:, 1]
    bd_s = lambda b: b[:, 0] + b[:, 1]

    d1 = (bd_delta(bu) - bd_delta(bd)) / (2 * H)
    d2 = (bd_delta(bu) - 2 * bd_delta(b0) + bd_delta(bd)) / (H ** 2)
    lam_s = -(bd_s(bu) - bd_s(bd)) / (2 * H)

    a = propensities_batch(comp, n)
    s_stoich = S[0, :] + S[1, :]                       # s-component per reaction
    D_ss = a @ (s_stoich ** 2)

    lam_s = np.where(lam_s > 1e-9, lam_s, np.nan)
    var_s = D_ss / (2.0 * lam_s)
    force = 0.5 * d2 * var_s                            # generalized force on delta
    D_extra = (d1 ** 2) * D_ss / (lam_s ** 2)           # induced diffusion on delta
    return (np.nan_to_num(force, nan=0.0),
            np.nan_to_num(np.clip(D_extra, 0.0, None), nan=0.0))


def run_arm(comp, n0, rng, *, dt, thr, trials, t_max, arm, max_steps=200_000) -> dict:
    S = comp.S.astype(float)
    n = np.tile(np.asarray(n0, dtype=float), (trials, 1))
    live = np.ones(trials, bool)
    t = np.zeros(trials)
    for _ in range(max_steps):
        if not live.any():
            break
        idx = np.where(live)[0]
        cur = n[idx]
        a = propensities_batch(comp, cur)
        mean = a * dt
        xi = np.sqrt(np.clip(mean, 0.0, None)) * rng.standard_normal(mean.shape)
        nz = xi @ S.T
        if arm != "full":
            h = 0.5 * (nz[:, 0] - nz[:, 1])             # keep delta's own noise only
            nz = np.stack([h, -h, np.zeros_like(h)], axis=1)
        drift = mean @ S.T
        if arm in ("+force", "+both"):
            force, _ = pool_terms(comp, cur, S)
            drift = drift + 0.5 * dt * np.stack(
                [force, -force, np.zeros_like(force)], axis=1)
        if arm in ("+diffusion", "+both"):
            _, dex = pool_terms(comp, cur, S)
            g = np.sqrt(np.clip(dex * dt, 0.0, None)) * rng.standard_normal(len(idx))
            nz = nz + 0.5 * np.stack([g, -g, np.zeros_like(g)], axis=1)
        cand = cur + drift + nz
        ok = (cand >= 0.0).all(axis=1)
        upd = idx[ok]
        n[upd] = cand[ok]
        t[upd] += dt
        done = (np.abs(n[idx, 0] - n[idx, 1]) >= thr) | (t[idx] >= t_max)
        live[idx[done]] = False
    fin = np.abs(n[:, 0] - n[:, 1]) >= thr
    nok = int(fin.sum())
    return {"p": float((n[fin, 0] <= n[fin, 1]).mean()) if nok else float("nan"),
            "n_ok": nok, "unfinished": int(live.sum())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.05, 0.15, 0.30, 0.45])
    ap.add_argument("--omegas", type=int, nargs="+", default=[40, 60, 80])
    ap.add_argument("--eps-frac", type=float, default=0.20)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--trials", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=20260731)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/zwanzig_correction.json"))
    args = ap.parse_args()

    t0 = time.time()
    print(f"eps={args.eps_frac} theta={args.theta} dt={args.dt} trials={args.trials}")
    print("  target: delta-only's residual, +13.2% at g=0.05 falling to +0.4% "
          "at g=0.45 (§24.1a)\n")
    print(f"{'gamma':>6} {'Om':>4} {'CME':>10} "
          + " ".join(f"{a:>14}" for a in ARMS))
    rows = []
    for g in args.gammas:
        net = am_reversible(g)
        for om in args.omegas:
            n0, thr, _ = _setup(g, om, args.eps_frac, args.theta)
            comp = compile_network(net, float(om))
            exact = p_cme(g, om, args.eps_frac, args.theta)
            got = {}
            for arm in ARMS:
                rng = np.random.default_rng(args.seed + 17 * om + int(1000 * g))
                got[arm] = run_arm(comp, n0, rng, dt=args.dt, thr=thr,
                                   trials=args.trials, t_max=4000.0, arm=arm)
            print(f"{g:>6.2f} {om:>4} {exact:>10.5f} "
                  + " ".join(f"{got[a]['p']:>7.5f}({got[a]['p']/exact-1:>+5.1%})"
                             for a in ARMS))
            rows.append({"gamma": g, "omega": om, "p_cme": exact,
                         **{f"p_{a}": got[a]["p"] for a in ARMS},
                         **{f"unfin_{a}": got[a]["unfinished"] for a in ARMS}})

    print(f"\n=== residual by arm and gamma (mean relative error vs exact)")
    print(f"{'gamma':>6} " + " ".join(f"{a:>14}" for a in ARMS))
    for g in args.gammas:
        rs = [r for r in rows if abs(r["gamma"] - g) < 1e-12]
        print(f"{g:>6.2f} " + " ".join(
            f"{np.mean([r[f'p_{a}']/r['p_cme']-1 for r in rs]):>+13.1%} "
            for a in ARMS))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
