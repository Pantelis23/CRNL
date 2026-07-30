"""Is the required noise subspace a property of the SYSTEM or of the (system, OBSERVABLE) pair?

§24.1 established that for one observable -- P(error) at first passage -- the noise
in the blank pool is 80-88% of the variance and entirely discardable, while the noise
in the signal subspace is irreplaceable: removing it gives exactly 0 in all eight
cells, the ODE's own failure. §24's whole framing then ASSERTS that this is a
property of the (system, observable) pair, i.e. that a different question of the same
chemistry would nominate a different subspace.

**That has never been measured.** Every observable this project has tested, in
twenty-four sections, is the restoration error probability. If the required subspace
is the same for every observable then §24 collapses to "AM has one stiff direction" --
a fact about this network, not about simulation.

THE TEST. Identical machinery to §24.1 -- same projected-noise CLE arms, same drift,
same network, same start states -- but TWO observables read off the SAME trajectories:

    P(error)   probability the initially-favoured species loses.  A TAIL event: it
               requires delta to cross against its drift.
    MFPT       mean time to first passage at |n_X - n_Y| >= thr.  A BULK quantity,
               dominated by the drift carrying delta to the threshold.

Both have an EXACT reference and, better, they come out of the SAME `first_passage`
solve, so they are paired by construction. Projecting the noise does not rescale
time -- same `dt`, same drift -- so the arms share a clock and rule 11 is satisfied.

WHY THIS IS NOT A DEFINITIONAL TRAP, which is the failure mode this thread has hit
twice (§24.1's `delta-only`, §24.2's `decision-only`). Under `bookkeeping-only` the
signal still reaches the threshold: delta is carried there by the DRIFT, which every
arm keeps intact. So MFPT is finite and well-defined in every arm, and nothing about
its value is forced by the construction. The arm that reports EXACTLY ZERO for
P(error) will report a perfectly ordinary number for MFPT, and whether that number is
right is an open question the run answers.

PREDICTIONS, written before running:

  P1  Harness control: the full CLE recovers both observables to a few percent, as it
      did in §24.1 for P(error).
  P2  THE DECISIVE ONE. `bookkeeping-only` -- categorically wrong on P(error), 0 in
      all eight cells -- recovers MFPT to within ~15%. If so, the SAME system, the
      SAME arms and the SAME trajectories give a required subspace that DIFFERS BY
      OBSERVABLE, and §24's framing is measured rather than inferred.
  P3  `signal-only` also recovers MFPT. The asymmetry should be ONE-SIDED -- both
      arms fine for the bulk quantity, only one fine for the tail -- rather than a
      swap. A swap would be a stronger result and I do not expect it, because the
      drift that carries MFPT is present in every arm by construction.
  P4  WHAT KILLS THE CLAIM. If `bookkeeping-only` is also badly wrong on MFPT (say
      worse than 30%), then both observables need the same subspace, §24's
      observable-dependence is unsupported, and the honest retreat is that this
      network simply has one stiff direction. That outcome is entirely possible and
      would be worth as much as the other.
  P5  A CONFOUND I HAVE TO REPORT AROUND, and it is rule 12 again. The unconditional
      MFPT averages over trajectories that decide BOTH ways, and the wrong-way ones
      take longer. `bookkeeping-only` produces no wrong-way trajectories at all, so
      part of any MFPT agreement or disagreement is inherited from its P(error)
      failure rather than being an independent fact about timing. Both are therefore
      reported: the UNCONDITIONAL MFPT against the exact reference, and the MFPT
      CONDITIONED on the favoured outcome across arms. There is no exact reference
      for the conditional one -- it needs a Doob h-transform that is not implemented
      here -- so it is compared arm-to-arm only, and labelled as such.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.approximations import propensities_batch
from crnl.cme import first_passage
from crnl.networks.am_reversible import am_reversible, reverse_pairing
from crnl.vectorized import compile_network
from experiments.approximation_hierarchy import _setup
from experiments.noise_placement import MODES, project


def run_two_observables(comp, n0, rng, *, dt, thr, trials, t_max, mode,
                        max_steps=400_000) -> dict:
    """P(error) AND first-passage time from the same trajectories."""
    S = comp.S.astype(float)
    n = np.tile(np.asarray(n0, dtype=float), (trials, 1))
    live = np.ones(trials, bool)
    t = np.zeros(trials)
    t_abs = np.full(trials, np.nan)
    var_full = var_kept = 0.0
    for _ in range(max_steps):
        if not live.any():
            break
        idx = np.where(live)[0]
        a = propensities_batch(comp, n[idx])
        mean = a * dt
        xi = np.sqrt(np.clip(mean, 0.0, None)) * rng.standard_normal(mean.shape)
        nz_full = xi @ S.T
        nz = project(nz_full, mode)
        var_full += float((nz_full ** 2).sum())
        var_kept += float((nz ** 2).sum())
        cand = n[idx] + (mean @ S.T) + nz
        ok = (cand >= 0.0).all(axis=1)
        upd = idx[ok]
        n[upd] = cand[ok]
        t[upd] += dt
        hit = np.abs(n[idx, 0] - n[idx, 1]) >= thr
        t_abs[idx[hit]] = t[idx[hit]]
        live[idx[hit | (t[idx] >= t_max)]] = False
    absorbed = np.isfinite(t_abs)
    nok = int(absorbed.sum())
    if not nok:
        return {"p_error": float("nan"), "mfpt": float("nan"),
                "mfpt_correct": float("nan"), "n_ok": 0,
                "unfinished": int(live.sum()),
                "variance_kept_frac": var_kept / var_full if var_full else float("nan")}
    wrong = n[absorbed, 0] <= n[absorbed, 1]
    right = ~wrong
    return {
        "p_error": float(wrong.mean()),
        "mfpt": float(t_abs[absorbed].mean()),
        "mfpt_correct": float(t_abs[absorbed][right].mean()) if right.any() else float("nan"),
        "mfpt_wrong": float(t_abs[absorbed][wrong].mean()) if wrong.any() else float("nan"),
        "n_ok": nok, "unfinished": int(live.sum()),
        "variance_kept_frac": var_kept / var_full if var_full else float("nan"),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.30)
    ap.add_argument("--omegas", type=int, nargs="+", default=[40, 60, 80, 100])
    ap.add_argument("--eps-fracs", type=float, nargs="+", default=[0.20, 0.35])
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--trials", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=20260730)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/observable_dependence.json"))
    args = ap.parse_args()

    net = am_reversible(args.gamma)
    pair = reverse_pairing(net)
    t0 = time.time()
    print(f"gamma={args.gamma} theta={args.theta} dt={args.dt} trials={args.trials}")
    print("  two observables off the SAME trajectories, both exact-referenced "
          "from ONE first_passage solve\n")
    rows = []
    for eps in args.eps_fracs:
        print(f"=== eps/delta* = {eps}")
        print(f"{'Om':>4} {'obs':>10} {'exact':>10} "
              + " ".join(f"{m:>13}" for m in MODES))
        for om in args.omegas:
            n0, thr, _ = _setup(args.gamma, om, eps, args.theta)
            comp = compile_network(net, float(om))
            ref = first_passage(net, om, float(om), n0,
                                lambda s, t=thr: abs(int(s[0]) - int(s[1])) >= t,
                                pair)
            got = {}
            for mode in MODES:
                rng = np.random.default_rng(args.seed + 17 * om + int(1000 * eps))
                got[mode] = run_two_observables(comp, n0, rng, dt=args.dt, thr=thr,
                                                trials=args.trials, t_max=4000.0,
                                                mode=mode)
            pe = 1.0 - float(ref["split"])
            mt = float(ref["mean_time"])
            print(f"{om:>4} {'P(error)':>10} {pe:>10.5f} "
                  + " ".join(f"{got[m]['p_error']:>13.5f}" for m in MODES))
            print(f"{'':>4} {'MFPT':>10} {mt:>10.4f} "
                  + " ".join(f"{got[m]['mfpt']:>13.4f}" for m in MODES))
            print(f"{'':>4} {'MFPT|right':>10} {'--':>10} "
                  + " ".join(f"{got[m]['mfpt_correct']:>13.4f}" for m in MODES))
            rows.append({"eps_frac": eps, "omega": om, "thr": thr,
                         "p_cme": pe, "mfpt_cme": mt, "cme_valid": ref["valid"],
                         **{f"{k}_{m}": got[m][k] for m in MODES
                            for k in ("p_error", "mfpt", "mfpt_correct",
                                      "variance_kept_frac", "unfinished", "n_ok")}})
        print()

    print("=== the contrast, per arm: relative error on each observable")
    for m in MODES:
        pe_err, mt_err = [], []
        for r in rows:
            if r["p_cme"] > 0 and np.isfinite(r[f"p_error_{m}"]):
                pe_err.append(r[f"p_error_{m}"] / r["p_cme"] - 1.0)
            if np.isfinite(r[f"mfpt_{m}"]):
                mt_err.append(r[f"mfpt_{m}"] / r["mfpt_cme"] - 1.0)
        zeros = sum(1 for r in rows if r[f"p_error_{m}"] == 0.0)
        pe_txt = (f"{np.mean(pe_err):+.1%} mean" if pe_err else "no positive cells")
        print(f"  {m:>14}: P(error) {pe_txt}"
              + (f"  [{zeros}/{len(rows)} cells exactly 0]" if zeros else "")
              + f"   |   MFPT {np.mean(mt_err):+.1%} mean, "
                f"worst {max(mt_err, key=abs):+.1%}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
