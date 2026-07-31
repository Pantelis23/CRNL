"""T13-a: learn the generalized force §24.1b could not derive — and score it on a number it never saw

§24.1b established that the cost of naively deleting the pool's noise is a
GENERALIZED FORCE, measured directly (+0.0483 at gamma = 0.05, ~0 at gamma = 0.45,
tracking §24.1a's residual), and that both closed forms for it fail: the curvature
term is identically zero because `b_delta` is exactly linear in `s`, and the
cross-correlation term is wrong in sign at gamma = 0.05 and 11x-140x too large
elsewhere. That is the situation where deriving stops and learning starts -- and it is
what the coarse-graining literature does (arXiv:2512.03706 learns the reduced drift
and diffusion with random features and shallow networks).

**What CRNL can do that most of that literature cannot: validate against an EXACT
reference.** Learned coarse-grained models are normally checked against more
simulation. Here the learned force is scored against the exact CME.

THE DESIGN, and it is built so the network cannot be tuned into looking right.

    TRAIN on   the local missing force f(delta, s), measured by paired full-vs-
               projected runs under COMMON RANDOM NUMBERS. Dense, cheap, local.
    SCORE on   exact P(error) -- a tail probability, exponentially sensitive, with an
               exact CME reference, and **never in the loss**.

A network with enough capacity fits its training target by construction, so training
accuracy is not evidence of anything. The only informative number is the held-out
tail observable, which no part of the fitting procedure sees. This is rule 16 with a
universal approximator attached: "let it become what we need" is the failure mode,
"let it become something, then score it on an exact number it never saw" is the
experiment.

TWO PARTS.

  A. SATURATION -- T13-a's named kill test, RUN AND FOUND CONFOUNDED. The raw force
     estimate keeps growing with the window (0.021 at 0.5, 0.057 at 2.0, 0.202 at 8.0
     at gamma = 0.05), but that does NOT establish non-Markovian behaviour: two
     different dynamics started from the same state diverge in phase space, and once
     they do, local force differences compound. Growth of `d(delta)/tau` at long
     windows is what ANY difference between two dynamics produces, memory or not. So
     the saturation test as posed cannot answer T13-a and is reported as
     uninformative rather than as evidence. The usable window is bounded below by the
     pool correlation time (the effect must develop) and above by trajectory
     divergence.

  A'. CONSISTENCY -- the replacement test, and it is not confounded. If a Markovian
     closure exists, the fitted `f(delta, s)` must be THE SAME whichever window it was
     estimated at, within the valid range. Fit at several tau and compare the
     resulting P(error). Agreement across tau means a well-defined Markovian closure;
     systematic drift with tau means the closure is an artifact of the window and no
     state function reproduces the effect.
  B. CLOSURE -- fit f(delta, s) at the saturating window and run the projected arm with
     it added to the drift.

CONTROLS, fixed in advance:

  * A LINEAR model in (delta, s) is fitted alongside the network. If linear does as
     well, the network adds nothing and the force is simple -- reporting only the MLP
     would be dressing up a straight line.
  * The uncorrected `delta-only` arm and the `full` arm bracket the result; `full`
     is the harness control.
  * Training-set force MSE is reported for every model, so a model that fits the
     force well and misses P(error) is visible as such rather than as a failure of
     fitting.

PREDICTIONS, written before running:

  P1  The force SATURATES by roughly 3-4 pool-correlation times. §24.1b only sampled
      two windows, both short, so the observed growth is probably approach-to-plateau
      rather than genuine non-Markovian growth.
  P2  The learned force closes most of the +13.2% residual at gamma = 0.05, landing
      inside the few-percent band the full CLE achieves.
  P3  THE OUTCOME I MOST EXPECT TO HAVE TO REPORT. The learned force fits its training
      target well (low MSE) and still misses P(error) by more than the full CLE does,
      because a tail probability is exponentially sensitive to a drift the loss only
      constrained on average. That is a rule-16 result and it is worth more than a
      success: it would say a locally-accurate closure is not a tail-accurate closure.
  P4  Linear does nearly as well as the MLP. The measured force varied smoothly and
      almost monotonically with gamma, and `b_delta` is bilinear, so there is little
      reason to expect strong nonlinearity in (delta, s).
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


def _project_delta(nz):
    h = 0.5 * (nz[:, 0] - nz[:, 1])
    return np.stack([h, -h, np.zeros_like(h)], axis=1)


def paired_force(comp, states, dt, nsteps, reps, seed):
    """Missing force per state, via COMMON RANDOM NUMBERS on full vs projected.

    Both arms are driven by the SAME standard normal draws, each scaled by its own
    propensities. Without that coupling the difference of two means of order 3.6 with
    per-run spread ~4.6 is invisible at any affordable sample size; with it the
    difference is estimated directly.
    """
    S = comp.S.astype(float)
    K = len(states)
    nf = np.repeat(np.asarray(states, float), reps, axis=0)
    npj = nf.copy()
    rng = np.random.default_rng(seed)
    for _ in range(nsteps):
        af, ap = propensities_batch(comp, nf), propensities_batch(comp, npj)
        z = rng.standard_normal(af.shape)
        cf = nf + (af * dt) @ S.T + (np.sqrt(np.clip(af * dt, 0, None)) * z) @ S.T
        cp = npj + (ap * dt) @ S.T + _project_delta(
            (np.sqrt(np.clip(ap * dt, 0, None)) * z) @ S.T)
        nf = np.where((cf >= 0).all(axis=1)[:, None], cf, nf)
        npj = np.where((cp >= 0).all(axis=1)[:, None], cp, npj)
    d = ((nf[:, 0] - nf[:, 1]) - (npj[:, 0] - npj[:, 1])).reshape(K, reps)
    tau = nsteps * dt
    return d.mean(axis=1) / tau, d.std(axis=1, ddof=1) / np.sqrt(reps) / tau


def sample_states(comp, n0, thr, dt, rng, n_states, max_steps=20000):
    """States visited by the FULL dynamics before absorption -- the distribution the
    closure will actually be evaluated on, rather than a grid over unvisited space."""
    S = comp.S.astype(float)
    n = np.tile(np.asarray(n0, float), (n_states, 1))
    live = np.ones(n_states, bool)
    keep = n.copy()
    stop_at = rng.integers(1, 260, size=n_states)
    for step in range(max_steps):
        if not live.any():
            break
        idx = np.where(live)[0]
        a = propensities_batch(comp, n[idx])
        z = rng.standard_normal(a.shape)
        cand = n[idx] + (a * dt) @ S.T + (np.sqrt(np.clip(a * dt, 0, None)) * z) @ S.T
        ok = (cand >= 0).all(axis=1)
        n[idx[ok]] = cand[ok]
        keep[idx] = n[idx]
        live[idx[(np.abs(n[idx, 0] - n[idx, 1]) >= thr) | (step >= stop_at[idx])]] = False
    return keep


def run_with_closure(comp, n0, rng, *, dt, thr, trials, t_max, mode, force_fn=None,
                     max_steps=200_000):
    S = comp.S.astype(float)
    n = np.tile(np.asarray(n0, float), (trials, 1))
    live = np.ones(trials, bool)
    t = np.zeros(trials)
    for _ in range(max_steps):
        if not live.any():
            break
        idx = np.where(live)[0]
        cur = n[idx]
        a = propensities_batch(comp, cur)
        nz = (np.sqrt(np.clip(a * dt, 0, None))
              * rng.standard_normal(a.shape)) @ S.T
        if mode != "full":
            nz = _project_delta(nz)
        drift = (a * dt) @ S.T
        if force_fn is not None:
            f = force_fn(cur)
            drift = drift + 0.5 * dt * np.stack([f, -f, np.zeros_like(f)], axis=1)
        cand = cur + drift + nz
        ok = (cand >= 0).all(axis=1)
        n[idx[ok]] = cand[ok]
        t[idx[ok]] += dt
        live[idx[(np.abs(n[idx, 0] - n[idx, 1]) >= thr) | (t[idx] >= t_max)]] = False
    fin = np.abs(n[:, 0] - n[:, 1]) >= thr
    k = int(fin.sum())
    return float((n[fin, 0] <= n[fin, 1]).mean()) if k else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.05)
    ap.add_argument("--omegas", type=int, nargs="+", default=[40, 60, 80])
    ap.add_argument("--eps-frac", type=float, default=0.20)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--windows", type=float, nargs="+", default=[0.5, 1, 2, 4, 8])
    ap.add_argument("--taus", type=float, nargs="+", default=[1.0, 2.0, 4.0],
                    help="A': a Markovian closure must give the same answer at every tau")
    ap.add_argument("--n-states", type=int, default=400)
    ap.add_argument("--reps", type=int, default=3000)
    ap.add_argument("--trials", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=20260731)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/learned_closure.json"))
    args = ap.parse_args()

    g = args.gamma
    net = am_reversible(g)
    t0 = time.time()
    lam_s = 1.0 + 2.0 * g
    print(f"gamma={g}  pool relaxation lam_s ~ {lam_s:.2f}, "
          f"correlation time ~ {1/lam_s:.2f}")

    print(f"\n=== PART A: does the missing force SATURATE with the window? "
          f"(T13-a's kill test)")
    om = args.omegas[0]
    n0, thr, _ = _setup(g, om, args.eps_frac, args.theta)
    comp = compile_network(net, float(om))
    print(f"{'window':>8} {'/tau_s':>8} {'force':>10} {'SE':>9}")
    sat = []
    for w in args.windows:
        f, se = paired_force(comp, np.array([n0], float), args.dt,
                             int(round(w / args.dt)), 20000, args.seed + 3)
        sat.append({"window": w, "force": float(f[0]), "se": float(se[0])})
        print(f"{w:>8.2f} {w*lam_s:>8.2f} {f[0]:>10.5f} {se[0]:>9.5f}")

    print(f"\n=== PART B/A': learn f(delta, s) at each tau in {args.taus}; a Markovian\n"
          f"    closure must give the SAME P(error) at every tau. Scored on exact CME.")
    rows = []
    for om, tau in [(o, t) for o in args.omegas for t in args.taus]:
        n0, thr, _ = _setup(g, om, args.eps_frac, args.theta)
        comp = compile_network(net, float(om))
        rng = np.random.default_rng(args.seed + om)
        states = sample_states(comp, n0, thr, args.dt, rng, args.n_states)
        fmeas, fse = paired_force(comp, states, args.dt,
                                  int(round(tau / args.dt)), args.reps,
                                  args.seed + 11 * om)
        good = np.isfinite(fmeas)
        X = np.stack([states[good, 0] - states[good, 1],
                      states[good, 0] + states[good, 1]], axis=1) / om
        y = fmeas[good]

        from sklearn.linear_model import LinearRegression
        from sklearn.neural_network import MLPRegressor
        from sklearn.preprocessing import StandardScaler
        from sklearn.pipeline import make_pipeline

        models = {
            "linear": make_pipeline(StandardScaler(), LinearRegression()),
            "mlp": make_pipeline(StandardScaler(),
                                 MLPRegressor(hidden_layer_sizes=(64, 64),
                                              max_iter=20000, random_state=0,
                                              tol=1e-7)),
        }
        exact = p_cme(g, om, args.eps_frac, args.theta)
        rec = {"gamma": g, "omega": om, "tau": tau, "p_cme": exact,
               "n_train": int(good.sum())}
        print(f"\n  Omega={om} tau={tau}  exact P(error)={exact:.6f}  "
              f"train n={int(good.sum())}  force range "
              f"[{y.min():+.4f}, {y.max():+.4f}]")
        for name, mdl in models.items():
            mdl.fit(X, y)
            mse = float(np.mean((mdl.predict(X) - y) ** 2))
            var = float(np.var(y))
            def ff(cur, mdl=mdl, om=om):
                z = np.stack([cur[:, 0] - cur[:, 1], cur[:, 0] + cur[:, 1]],
                             axis=1) / om
                return mdl.predict(z)
            rr = np.random.default_rng(args.seed + 5 * om)
            p = run_with_closure(comp, n0, rr, dt=args.dt, thr=thr,
                                 trials=args.trials, t_max=4000.0,
                                 mode="proj", force_fn=ff)
            rec[f"p_{name}"] = p
            rec[f"mse_{name}"] = mse
            rec[f"r2_{name}"] = 1.0 - mse / var if var > 0 else float("nan")
            print(f"    {name:>7}: train R^2={rec[f'r2_{name}']:>6.3f}   "
                  f"P(error)={p:.6f}  ({p/exact-1:+.1%} vs exact)")
        for name, mode in (("delta-only", "proj"), ("full", "full")):
            rr = np.random.default_rng(args.seed + 5 * om)
            p = run_with_closure(comp, n0, rr, dt=args.dt, thr=thr,
                                 trials=args.trials, t_max=4000.0, mode=mode)
            rec[f"p_{name}"] = p
            print(f"    {name:>7}: {'':>16}   P(error)={p:.6f}  "
                  f"({p/exact-1:+.1%} vs exact)")
        rows.append(rec)

    print(f"\n=== summary: relative error on P(error), which was NEVER in the loss")
    for k in ("full", "delta-only", "linear", "mlp"):
        e = [r[f"p_{k}"] / r["p_cme"] - 1 for r in rows if np.isfinite(r.get(f"p_{k}", np.nan))]
        if e:
            print(f"  {k:>11}: {np.mean(e):>+7.1%} mean  "
                  f"[{min(e):+.1%}, {max(e):+.1%}]")
    print("\n=== A': is the closure the SAME at every tau? (Markovian iff yes)")
    for om in args.omegas:
        rs = [r for r in rows if r["omega"] == om]
        for k in ("linear", "mlp"):
            v = [(r["tau"], r[f"p_{k}"] / r["p_cme"] - 1) for r in rs]
            print(f"  Omega={om} {k:>7}: " + "  ".join(
                f"tau={t}:{e:+.1%}" for t, e in v))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"saturation": sat, "cells": rows},
                                   indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
