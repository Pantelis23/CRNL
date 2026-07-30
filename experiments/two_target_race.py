"""T12's last candidate: a race the POOL can win

§25 showed the required noise subspace varies by observable, but only DOWNWARD --
P(error) needs `span(delta)`, MFPT needs nothing. §25.1 then tested the sharpest
pure-noise candidate, `Var(T)`, and it too came back a `span(delta)` observable:
`delta-only` recovers it alone, `s-only` captures 16%. Three observables, one
subspace. The strong form of T12 -- an observable of this system that requires the
POOL's noise -- remains unfound.

This is the last candidate on the list: a RACE between two absorbing targets, one in
each coordinate.

    decision target   |n_X - n_Y| >= thr        (the usual restoration threshold)
    pool target       n_B <= m                  (the blank pool runs down)

and the observable is the splitting probability between them. `m` is chosen per cell
as the value minimising |P(pool first) - 1/2| against the exact CME, so the race is
genuinely competitive rather than a formality; at Omega = 40 that is m = 5 against a
starting pool of 9, i.e. the pool must fluctuate DOWN BY 4 to win. `cme.splitting_
probability` takes both predicates directly, so the reference is exact and needs no
new machinery. States satisfying both conditions are assigned to the decision target,
and the ambiguity is reported.

WHICH CELL IS INFORMATIVE, fixed before running because this design sits closer to
the definitional line than the review that proposed it acknowledged:

  * **`delta-only`'s failure would prove very little.** With the pool's own noise
    removed, `n_B` moves only by a drift that delta drives, so whether it ever
    reaches `m` is close to a deterministic question. If that arm collapses to ~0 or
    ~1 it is the mirror of §24.1's `delta-only` anchor and §24.2's `decision-only`
    trap -- an arm answering by construction. Predicted, and not counted as evidence.
  * **`s-only`'s ACCURACY is the whole test.** It has the pool's noise and none of
    delta's. If it reproduces the exact splitting probability, then this is an
    observable the pool's noise is SUFFICIENT for and the signal's is not -- the
    reversal T12 has been looking for, and the first such case in five sections.

PRE-COMMITTED CRITERION, fixed before any output is seen. The observable is a
probability near 1/2, so a ratio band is the wrong instrument (§25.1's 3x band would
be unreachable) and absolute deviation is used instead:

    recovers   |p_arm - p_exact| <= 0.05
    partial    0.05 < |p_arm - p_exact| <= 0.20
    fails      |p_arm - p_exact| > 0.20

PREDICTIONS, written before running:

  P1  `delta-only` fails, probably near-categorically. Stated so it cannot later be
      presented as a finding: it is semi-definitional and carries no weight.
  P2  THE INFORMATIVE ONE. `s-only` recovers, or at worst lands "partial". The pool
      target is a pool-fluctuation event and `s-only` is the arm that has pool
      fluctuations. If it recovers while being categorically wrong about P(error) on
      the same network, T12's strong form is established.
  P3  Harness control: `full` recovers in every cell. If it does not, `m` is in a
      regime where the CLE itself is unreliable and nothing else is readable.
  P4  THE OUTCOME THAT WOULD BE CLEANEST OF ALL, and which I would take over P2:
      NEITHER single-subspace arm recovers, so the race requires
      `span(delta) + span(s)`. That is not definitional in either direction -- both
      arms have a live mechanism -- and it would establish observable-dependence
      without leaning on an arm that answers by construction.
  P5  IF `s-only` ALSO FAILS while `delta-only` answers by construction, then this
      instrument has produced nothing, and the honest conclusion is about the
      INSTRUMENT rather than the chemistry: every observable reachable by
      first-passage on this simplex may be a question about delta, and finding a
      pool-requiring observable would need machinery this project does not have.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.approximations import propensities_batch
from crnl.cme import splitting_probability
from crnl.networks.am_reversible import am_reversible
from crnl.vectorized import compile_network
from experiments.approximation_hierarchy import _setup
from experiments.noise_placement import MODES, project

BAND_RECOVER = 0.05
BAND_FAIL = 0.20


def verdict(dp: float) -> str:
    if not np.isfinite(dp):
        return "fails"
    if abs(dp) <= BAND_RECOVER:
        return "recovers"
    if abs(dp) > BAND_FAIL:
        return "fails"
    return "partial"


def exact_race(net, om, n0, thr, m) -> dict:
    ab = lambda s: (abs(int(s[0]) - int(s[1])) >= thr) or (int(s[2]) <= m)
    fav = lambda s: (int(s[2]) <= m) and (abs(int(s[0]) - int(s[1])) < thr)
    return splitting_probability(net, om, float(om), n0, ab, fav)


def pick_m(net, om, n0, thr) -> tuple[int, float]:
    """m minimising |P(pool first) - 1/2|: the race must be competitive."""
    best, best_p = None, None
    for m in range(0, int(n0[2])):
        r = exact_race(net, om, n0, thr, m)
        if not r["valid"]:
            continue
        if best is None or abs(r["split"] - 0.5) < abs(best_p - 0.5):
            best, best_p = m, r["split"]
    return best, best_p


def run_arm(comp, n0, rng, *, dt, thr, m, trials, t_max, mode,
            max_steps=400_000) -> dict:
    S = comp.S.astype(float)
    n = np.tile(np.asarray(n0, dtype=float), (trials, 1))
    live = np.ones(trials, bool)
    t = np.zeros(trials)
    pool_win = np.zeros(trials, bool)
    done_any = np.zeros(trials, bool)
    both = 0
    for _ in range(max_steps):
        if not live.any():
            break
        idx = np.where(live)[0]
        a = propensities_batch(comp, n[idx])
        mean = a * dt
        xi = np.sqrt(np.clip(mean, 0.0, None)) * rng.standard_normal(mean.shape)
        nz = project(xi @ S.T, mode)
        cand = n[idx] + (mean @ S.T) + nz
        ok = (cand >= 0.0).all(axis=1)
        upd = idx[ok]
        n[upd] = cand[ok]
        t[upd] += dt
        dec = np.abs(n[idx, 0] - n[idx, 1]) >= thr
        pool = n[idx, 2] <= m
        both += int((dec & pool).sum())
        hit = dec | pool
        # states satisfying both are assigned to the DECISION target
        pool_win[idx[hit]] = pool[hit] & ~dec[hit]
        done_any[idx[hit]] = True
        live[idx[hit | (t[idx] >= t_max)]] = False
    nok = int(done_any.sum())
    return {"p_pool": float(pool_win[done_any].mean()) if nok else float("nan"),
            "n_ok": nok, "unfinished": int(live.sum()), "both_at_once": both}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.30)
    ap.add_argument("--omegas", type=int, nargs="+", default=[40, 50, 60, 70])
    ap.add_argument("--eps-fracs", type=float, nargs="+", default=[0.20])
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--trials", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=20260731)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/two_target_race.json"))
    args = ap.parse_args()

    net = am_reversible(args.gamma)
    t0 = time.time()
    print(f"gamma={args.gamma} theta={args.theta} dt={args.dt} trials={args.trials}")
    print(f"  criterion fixed in advance: recovers |dp|<={BAND_RECOVER}, "
          f"fails |dp|>{BAND_FAIL}\n")
    print(f"{'eps':>5} {'Om':>4} {'m':>3} {'nB0':>4} {'P(pool) exact':>14} "
          + " ".join(f"{m:>17}" for m in MODES))
    rows = []
    for eps in args.eps_fracs:
        for om in args.omegas:
            n0, thr, _ = _setup(args.gamma, om, eps, args.theta)
            m, pe = pick_m(net, om, n0, thr)
            if m is None:
                print(f"{eps:>5} {om:>4}  no valid m")
                continue
            comp = compile_network(net, float(om))
            got = {}
            for mode in MODES:
                rng = np.random.default_rng(args.seed + 17 * om + int(1000 * eps))
                got[mode] = run_arm(comp, n0, rng, dt=args.dt, thr=thr, m=m,
                                    trials=args.trials, t_max=4000.0, mode=mode)
            print(f"{eps:>5.2f} {om:>4} {m:>3} {int(n0[2]):>4} {pe:>14.4f} "
                  + " ".join(f"{got[k]['p_pool']:>8.4f}({verdict(got[k]['p_pool']-pe)[:4]:>4})"
                             for k in MODES))
            rows.append({"eps_frac": eps, "omega": om, "m": m, "nB0": int(n0[2]),
                         "thr": thr, "p_exact": pe,
                         **{f"{k}_{mm}": got[mm][k] for mm in MODES
                            for k in ("p_pool", "n_ok", "unfinished", "both_at_once")}})

    print(f"\n=== verdicts (deviation from exact), and P(error) behaviour for context")
    for mm in MODES:
        dp = [r[f"p_pool_{mm}"] - r["p_exact"] for r in rows
              if np.isfinite(r[f"p_pool_{mm}"])]
        vs = [verdict(x) for x in dp]
        print(f"  {mm:>14}: mean dev {np.mean(dp):+.4f}  "
              f"[{min(dp):+.4f}, {max(dp):+.4f}]   "
              f"{vs.count('recovers')} recovers, {vs.count('partial')} partial, "
              f"{vs.count('fails')} fails")
    tot_both = sum(r[f"both_at_once_{mm}"] for r in rows for mm in MODES)
    print(f"\n  states hitting BOTH targets in one step (assigned to decision): "
          f"{tot_both} across all arms and cells")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
