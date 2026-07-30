"""T12: is there an observable of this system that needs the POOL's noise?

§25 showed the required noise subspace varies by observable, but only DOWNWARD:
P(error) needs `span(delta)`, mean decision time needs nothing. An independent
review made the sharp objection -- an observable requiring the EMPTY subspace is
evidence that means are means, not that requirements differ by question. **The strong
form needs an observable whose requirement points at `s`,** the blank pool.

`Var(T)` at first passage is the review's best candidate and it is a good one:

  * It is a PURE-NOISE observable. The deterministic limit gives exactly zero, so
    there is no drift-dominated leading order to dilute the comparison against --
    the flaw that made MFPT a weak test.
  * TWO CREDIBLE MECHANISMS, and neither is guessable from the definition. Timing
    jitter can come from (a) diffusion along the CROSSING direction, which is
    delta's own noise, or from (b) fluctuations in the RATE OF ADVANCE: every
    recruitment propensity carries `n_B`, so a fluctuating pool modulates how fast
    delta is driven to the threshold. (a) points at `span(delta)`, (b) at `span(s)`.
  * It has an EXACT reference. `cme.first_passage_moments` solves `Qtt m2 = -2T`
    beside the existing mean, and `test_first_passage_variance_matches_the_ssa`
    pins it against sampled jump trajectories -- necessary because a wrong factor
    in that recursion still yields a positive, plausible variance.

At Omega = 40-100 the exact CV = sd(T)/E[T] runs 0.52-0.63, so the jitter is a large
effect, not a correction. And the pool is small -- `n_B ~ 9` at Omega = 40 -- so its
relative fluctuation is ~1/3, which is why mechanism (b) is credible rather than
notional.

NOT A DEFINITIONAL TRAP, checked the way §24.2 taught: under `s-only` delta reaches
the threshold by the drift every arm retains, and the TIME it takes varies because
`n_B` varies. So `s-only` returns a substantial, non-trivial Var(T) whose correctness
is an open question. Under `delta-only` the pool is deterministic but delta diffuses.
Both arms produce jitter by different routes; neither answer follows from the setup.

PRE-COMMITTED CRITERION, fixed before looking at any output because §24.1's exact
zeros will not recur and "categorical" must not be decided by eye:

    recovers  within 25% of the exact variance
    partial   between 25% and 3x off
    fails     more than 3x off, in either direction

PREDICTIONS, written before running:

  P1  BOTH arms produce substantial jitter and NEITHER recovers alone: `delta-only`
      captures mechanism (a) and misses (b), `s-only` the reverse, and only `full`
      lands inside the band. That is the STRONG form -- an observable requiring
      `span(delta) + span(s)` -- and it is what I expect.
  P2  THE STRONGEST POSSIBLE OUTCOME, and I do not expect it: `s-only` recovers
      Var(T) while remaining categorically wrong on P(error) from the SAME
      trajectories. That would be an observable of this system needing the pool's
      noise and not the signal's -- observable-dependence in both directions rather
      than merely downward.
  P3  WHAT WOULD SINK T12. If `s-only` returns a Var(T) near zero, the pool
      contributes no timing jitter, Var(T) is a third delta-observable, and the
      strong form is unfound. The honest reading would then be that all three
      observables tested so far are saddle-dominated -- which per the review does
      NOT establish "one stiff direction", but does mean I have been asking the same
      question three ways.
  P4  THE CONFOUND, again rule 12 and inherited from §25. The unconditional variance
      mixes correct and error paths, and error paths are slower, so part of it is
      bimodality rather than jitter. `s-only` has no error paths at all. Both are
      therefore reported: unconditional against the exact reference, and
      Var(T | correct) arm-to-arm, for which no exact reference exists here.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.approximations import propensities_batch
from crnl.cme import first_passage_moments
from crnl.networks.am_reversible import am_reversible
from crnl.vectorized import compile_network
from experiments.approximation_hierarchy import _setup
from experiments.noise_placement import MODES, project

BAND_RECOVER = 0.25
BAND_FAIL = 3.0


def verdict(ratio: float) -> str:
    if not np.isfinite(ratio) or ratio <= 0:
        return "fails"
    if abs(ratio - 1.0) <= BAND_RECOVER:
        return "recovers"
    if ratio > BAND_FAIL or ratio < 1.0 / BAND_FAIL:
        return "fails"
    return "partial"


def run_arm(comp, n0, rng, *, dt, thr, trials, t_max, mode,
            max_steps=400_000) -> dict:
    S = comp.S.astype(float)
    n = np.tile(np.asarray(n0, dtype=float), (trials, 1))
    live = np.ones(trials, bool)
    t = np.zeros(trials)
    t_abs = np.full(trials, np.nan)
    rejected = 0
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
        rejected += int((~ok).sum())
        upd = idx[ok]
        n[upd] = cand[ok]
        t[upd] += dt
        hit = np.abs(n[idx, 0] - n[idx, 1]) >= thr
        t_abs[idx[hit]] = t[idx[hit]]
        live[idx[hit | (t[idx] >= t_max)]] = False
    absorbed = np.isfinite(t_abs)
    if not absorbed.any():
        return {"var": float("nan"), "var_correct": float("nan"),
                "mean": float("nan"), "p_error": float("nan"),
                "n_ok": 0, "rejected": rejected}
    ta = t_abs[absorbed]
    wrong = n[absorbed, 0] <= n[absorbed, 1]
    right = ~wrong
    return {
        "var": float(ta.var(ddof=1)),
        "var_correct": float(ta[right].var(ddof=1)) if right.sum() > 1 else float("nan"),
        "mean": float(ta.mean()),
        "p_error": float(wrong.mean()),
        "n_ok": int(absorbed.sum()), "rejected": rejected,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.30)
    ap.add_argument("--omegas", type=int, nargs="+", default=[40, 60, 80, 100])
    ap.add_argument("--eps-fracs", type=float, nargs="+", default=[0.20, 0.35])
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--trials", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=20260731)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/timing_jitter.json"))
    args = ap.parse_args()

    net = am_reversible(args.gamma)
    t0 = time.time()
    print(f"gamma={args.gamma} theta={args.theta} dt={args.dt} trials={args.trials}")
    print(f"  verdict bands fixed in advance: recovers <= {BAND_RECOVER:.0%} off, "
          f"fails > {BAND_FAIL:.0f}x off\n")
    print(f"{'eps':>5} {'Om':>4} {'Var(T) exact':>13} "
          + " ".join(f"{m:>15}" for m in MODES))
    rows = []
    for eps in args.eps_fracs:
        for om in args.omegas:
            n0, thr, _ = _setup(args.gamma, om, eps, args.theta)
            comp = compile_network(net, float(om))
            ref = first_passage_moments(
                net, om, float(om), n0,
                lambda s, t=thr: abs(int(s[0]) - int(s[1])) >= t)
            got = {}
            for mode in MODES:
                rng = np.random.default_rng(args.seed + 17 * om + int(1000 * eps))
                got[mode] = run_arm(comp, n0, rng, dt=args.dt, thr=thr,
                                    trials=args.trials, t_max=4000.0, mode=mode)
            ve = ref["var_time"]
            print(f"{eps:>5.2f} {om:>4} {ve:>13.4f} "
                  + " ".join(f"{got[m]['var']:>7.3f}({verdict(got[m]['var']/ve)[:4]:>4})"
                             for m in MODES))
            rows.append({"eps_frac": eps, "omega": om, "var_cme": ve,
                         "mean_cme": ref["mean_time"], "cme_valid": ref["valid"],
                         **{f"{k}_{m}": got[m][k] for m in MODES
                            for k in ("var", "var_correct", "mean", "p_error",
                                      "n_ok", "rejected")}})

    print(f"\n=== verdicts on Var(T), and P(error) on the SAME trajectories")
    for m in MODES:
        rat = [r[f"var_{m}"] / r["var_cme"] for r in rows if r["var_cme"] > 0]
        vs = [verdict(x) for x in rat]
        zeros = sum(1 for r in rows if r[f"p_error_{m}"] == 0.0)
        print(f"  {m:>14}: Var(T) ratio {np.mean(rat):>6.3f} mean "
              f"[{min(rat):.3f}-{max(rat):.3f}]  "
              f"{vs.count('recovers')}/{len(vs)} recovers, "
              f"{vs.count('partial')} partial, {vs.count('fails')} fails"
              f"   |   P(error) exactly 0 in {zeros}/{len(rows)}")

    print(f"\n=== P4 control: Var(T | correct), arm-to-arm (no exact reference)")
    for m in MODES:
        rat = [r[f"var_correct_{m}"] / r["var_correct_full"] for r in rows
               if np.isfinite(r[f"var_correct_{m}"]) and r["var_correct_full"] > 0]
        print(f"  {m:>14}: {np.mean(rat):>6.3f} mean of full "
              f"[{min(rat):.3f}-{max(rat):.3f}]")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
