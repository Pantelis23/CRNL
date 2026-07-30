"""T11-REFINED-a: how FEW noise directions does a restoration observable need?

§24.1 split `am_reversible`'s noise into the signal `delta = n_X - n_Y` and the
bookkeeping `s = n_X + n_Y`, and found the observable indifferent to the second
(88% of the variance, removable at a cost of 2-18%) and categorically dependent on
the first. But n = 2 has exactly one signal direction and one bookkeeping direction,
so the split had nowhere to hide -- "the bit" and "the pool" are the only two things
there are. THEORIES T11-REFINED-a names the test: go where the simplex has SEVERAL
signal directions and no single "the bit" coordinate.

At n = 3 the noise vector `(a1, a2, a3, ab)` sums to zero and decomposes into three
ORTHOGONAL free directions:

    d1 = (1, -1/2, -1/2)   CHAMPION vs RIVALS -- the direction the observable's
                           decision actually lives in (does X1 lose?)
    d2 = (0, 1, -1)        RIVAL vs RIVAL -- which of X2, X3 wins, a question the
                           observable never asks
    total (1, 1, 1)/3 + ab COMMITTED vs BLANK -- the pool, §24.1's `s`

`d1 . d2 = 0`, and both are orthogonal to `(1,1,1)`, so `signal-only` is exactly
`decision-only` + `rivals-only`. Every arm keeps the FULL DRIFT; only the noise is
projected. Reference is the exact CME via `cme.splitting_probability` (not
`first_passage`, whose favoured set is hardcoded `n[0] > n[1]` and is silently wrong
above n = 2 -- §21.4's own note).

WHAT IS INFORMATIVE HERE, stated before running. That killing `d1`'s noise kills the
observable is close to definitional, as it was in §24.1, and proves little on its own.
**The informative question is how much can be thrown away**: `decision-only` retains
ONE of three free directions. If it still recovers the exponent, then an n = 3 race
needs noise in a single direction out of three, which is a statement about simulation
cost rather than a tautology. And `rivals-only` tests whether the split is by SPECIES
ROLE (X2, X3 are both "signal" species) or by DIRECTION RELATIVE TO THE OBSERVABLE
(X2 - X3 is bookkeeping *for this question*).

PREDICTIONS, written before running:

  P1  `bookkeeping-only` fails categorically, reproducing §24.1's `s-only` at n = 3.
  P2  `signal-only` recovers the exponent to within roughly §24.1's 2-18%.
  P3  THE SHARP ONE. `decision-only` -- one direction of three -- recovers nearly as
      well as `signal-only`, and `rivals-only` fails. That would mean the useful split
      is not signal-species vs pool-species but **the direction the observable's
      decision lives in vs everything else**, and that a simulation may discard noise
      in signal-carrying species as long as it keeps the right combination of them.
  P4  IF `decision-only` FAILS WHILE `signal-only` WORKS, the decision needs the whole
      difference subspace and the split does NOT decompose direction-by-direction. The
      honest conclusion would then be that §24.1's result is a two-species convenience
      -- real at n = 2, not a general principle -- and T11-REFINED's wording would
      have to retreat from "the signal coordinate" to "the signal subspace".
  P5  IF `rivals-only` RECOVERS, the rival-rival direction feeds the champion's fate
      through the nonlinearity and the observable is not as directional as §24.1
      suggested. That would weaken the refinement rather than kill it.

OUTCOME (§24.2): P1 and P2 confirmed -- `bookkeeping-only` is categorically 0 in all
eight cells with 80% of the variance, `signal-only` carries the exponent to 6.6-8.5%.
P3 and P4 UNRESOLVED, because `decision-only` was invalid: with rival-vs-rival noise
zeroed, `X2 - X3` evolves deterministically from its start, and `setup` leaves the
rivals differing by 0 or 1 according to an Omega parity, so the arm returned 0 when
they started tied and 0.054 when they did not.

T11-REFINED-b, PREDICTIONS for the `--rival-skew` fix, written before running:

  P6  With a fixed skew the arm becomes smooth in Omega -- no more alternation
      between 0 and a finite value -- and `full` and `signal-only` barely move,
      since they always carried rival-vs-rival noise and never depended on the
      parity. If THEY move materially, the skew is not a neutral repair and the
      whole comparison has to be redone at matched start states.
  P7  THE SUBSTANTIVE ONE, and I expect `decision-only` to UNDER-estimate rather
      than recover. The champion loses to the BEST rival, not to the rival mean, and
      the best of two NOISY rivals sits higher than the best of two rivals whose
      separation is fixed by drift alone -- an order-statistic effect that
      `decision-only` throws away along with d2. So removing rival-vs-rival noise
      should LOWER the champion's failure probability, by a factor that grows with
      the number of rivals. If instead `decision-only` matches `signal-only`, then
      one direction of n suffices, and that is a real statement about how cheap a
      restoration simulation can be.
  P8  Sweeping the skew (2 vs 4) is the convergence check WITHIN the repair (rule
      13): `decision-only` should depend on it only weakly. A strong dependence
      means the arm is still being driven by the initial condition rather than by
      the dynamics, and the answer stays unresolved rather than becoming P7's.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.approximations import propensities_batch
from crnl.networks.n_winner_reversible import n_winner_reversible
from crnl.vectorized import compile_network
from crnl.networks.n_winner_reversible import gamma_critical
from experiments.approximation_hierarchy_nwinner import (
    _absorbing_batch, landscape_width, p_cme, setup,
)

MODES = ("full", "signal-only", "bookkeeping-only", "decision-only", "rivals-only")


def setup_skewed(n, gamma, omega, eps_frac, theta, width, skew):
    """Start state with the champion's margin AND a fixed rival asymmetry.

    T11-REFINED-b. `approximation_hierarchy_nwinner.setup` pins the rival maximum
    and takes the remainder off the other rivals, so at n = 3 the two rivals differ
    by 0 or 1 depending on an integer parity in Omega. That is harmless for arms
    that carry rival-vs-rival noise and FATAL for `decision-only`, which zeroes it:
    with X2 - X3 evolving deterministically from its initial value, rivals that
    start tied stay tied and the champion wins by construction (§24.2).

    Forcing `max(rivals) - min(rivals) = skew` exactly, independent of Omega,
    removes the parity. Same discipline as `setup`: pin the rival maximum, absorb
    the remainder by nudging the margin, and assert both invariants.
    """
    from crnl.networks.n_winner_reversible import symmetric_state
    _, b = symmetric_state(n, gamma)
    nb = int(round(b * omega))
    rest = omega - nb
    m = max(1, int(round(eps_frac * width * omega)))
    for dm in (0, 1, -1, 2, -2, 3, -3):
        mm = m + dm
        if mm < 1:
            continue
        num = rest - mm + (n - 2) * skew
        if num % n:
            continue
        R = num // n
        rivals = [R] + [R - skew] * (n - 2)
        c = R + mm
        n0 = np.array([c] + rivals + [nb], dtype=np.int64)
        if (n0 >= 0).all() and int(n0.sum()) == omega \
                and c - max(rivals) == mm \
                and (n == 2 or max(rivals) - min(rivals) == skew):
            thr = max(2, int(round(theta * width * omega)))
            return n0, thr, float(mm) / omega
    raise ValueError(f"no exact-margin start at Omega={omega}, n={n}, skew={skew}")


def _dirs(n: int) -> tuple[np.ndarray, np.ndarray]:
    d1 = np.zeros(n)
    d1[0] = 1.0
    d1[1:] = -1.0 / (n - 1)
    d2 = np.zeros(n)
    if n >= 3:
        d2[1], d2[2] = 1.0, -1.0
    return d1, d2


def project(nz: np.ndarray, mode: str, n: int) -> np.ndarray:
    """Project the species-space noise increment. Conservation is exact in all arms."""
    if mode == "full":
        return nz
    u = nz[:, :n]
    out = np.zeros_like(nz)
    tot = u.sum(axis=1, keepdims=True)
    if mode == "bookkeeping-only":
        out[:, :n] = tot / n
        out[:, n] = nz[:, n]
        return out
    diff = u - tot / n                       # committed differences, sums to zero
    if mode == "signal-only":
        out[:, :n] = diff
        return out
    d1, d2 = _dirs(n)
    d = d1 if mode == "decision-only" else d2
    if not np.any(d):
        raise ValueError(f"{mode} undefined for n={n}")
    out[:, :n] = ((diff @ d) / (d @ d))[:, None] * d[None, :]
    return out


def run_projected(comp, n0, rng, *, dt, stop, trials, t_max, mode, n,
                  max_steps=300_000) -> dict:
    S = comp.S.astype(float)
    m = np.tile(np.asarray(n0, dtype=float), (trials, 1))
    live = np.ones(trials, bool)
    t = np.zeros(trials)
    var_full = var_kept = 0.0
    rejected = 0
    for _ in range(max_steps):
        if not live.any():
            break
        idx = np.where(live)[0]
        a = propensities_batch(comp, m[idx])
        mean = a * dt
        xi = np.sqrt(np.clip(mean, 0.0, None)) * rng.standard_normal(mean.shape)
        nz_full = xi @ S.T
        nz = project(nz_full, mode, n)
        var_full += float((nz_full ** 2).sum())
        var_kept += float((nz ** 2).sum())
        cand = m[idx] + (mean @ S.T) + nz
        ok = (cand >= 0.0).all(axis=1)
        rejected += int((~ok).sum())
        upd = idx[ok]
        m[upd] = cand[ok]
        t[upd] += dt
        done = stop(m[idx]) | (t[idx] >= t_max)
        live[idx[done]] = False
    fin = stop(m)
    nok = int(fin.sum())
    wrong = int((np.argmax(m[fin][:, :n], axis=1) != 0).sum()) if nok else 0
    return {"p": wrong / nok if nok else float("nan"), "n_ok": nok,
            "unfinished": int(live.sum()), "rejected": rejected,
            "variance_kept_frac": var_kept / var_full if var_full > 0 else float("nan")}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--gamma-frac", type=float, default=0.60,
                    help="gamma as a fraction of gamma_c(n), matching §21.4. A raw "
                         "gamma copied from §21's AM run is WRONG here: gamma_c "
                         "falls with n, gamma_c(3) = 0.2023, and the first pass at "
                         "gamma = 0.30 measured a network with NO LANDSCAPE "
                         "(landscape_width = 0, threshold collapsed to its floor).")
    ap.add_argument("--omegas", type=int, nargs="+", default=[30, 45, 60, 80])
    ap.add_argument("--eps-fracs", type=float, nargs="+", default=[0.25, 0.40])
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--trials", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=20260730)
    ap.add_argument("--rival-skew", type=int, default=0,
                    help="fixed max(rivals)-min(rivals) in the start state. 0 uses "
                         "the original setup, whose 0-or-1 parity invalidated "
                         "decision-only (§24.2). >=2 is T11-REFINED-b's fix.")
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/noise_placement_nwinner.json"))
    args = ap.parse_args()

    n = args.n
    g = args.gamma_frac * gamma_critical(n)
    width = landscape_width(n, g)
    if not width > 1e-6:
        raise SystemExit(
            f"landscape_width = {width:.6g} at n={n}, gamma={g:.5f}: the landscape "
            f"is dead, so there is no restoration observable to measure. "
            f"gamma_c({n}) = {gamma_critical(n):.5f}.")
    t0 = time.time()
    print(f"n={n} gamma={g:.5f} (= {args.gamma_frac} x gamma_c={gamma_critical(n):.5f}) "
          f"theta={args.theta} dt={args.dt} trials={args.trials}")
    print(f"  landscape width (champion - best rival at the attractor) = {width:.5f}")
    print(f"  rival skew = {args.rival_skew}"
          + ("  [T11-REFINED-b: fixed, so decision-only is not decided by Omega parity]"
             if args.rival_skew else "  [original setup; decision-only INVALID here]"))
    print(f"\n{'eps':>5} {'Omega':>6} {'CME':>10} {'full':>10} {'signal':>10} "
          f"{'bookkeep':>10} {'decision':>10} {'rivals':>10}  variance kept")
    rows = []
    for eps in args.eps_fracs:
        for om in args.omegas:
            n0, thr, realised = (
                setup_skewed(n, g, om, eps, args.theta, width, args.rival_skew)
                if args.rival_skew else setup(n, g, om, eps, args.theta, width))
            comp = compile_network(n_winner_reversible(n, g), float(om))
            stop = _absorbing_batch(n, thr)
            exact = p_cme(n, g, om, n0, thr)
            got = {}
            for mode in MODES:
                rng = np.random.default_rng(args.seed + 17 * om + int(1000 * eps))
                got[mode] = run_projected(comp, n0, rng, dt=args.dt, stop=stop,
                                          trials=args.trials, t_max=6000.0,
                                          mode=mode, n=n)
            print(f"{eps:>5.2f} {om:>6} {exact:>10.5f} "
                  + " ".join(f"{got[m]['p']:>10.5f}" for m in MODES)
                  + "  " + "/".join(f"{got[m]['variance_kept_frac']:.3f}"
                                    for m in MODES[1:]))
            rows.append({"eps_frac": eps, "omega": om, "realised_eps": realised,
                         "p_cme": exact,
                         **{f"p_{m}": got[m]["p"] for m in MODES},
                         **{f"var_{m}": got[m]["variance_kept_frac"] for m in MODES},
                         **{f"unfin_{m}": got[m]["unfinished"] for m in MODES}})

    print(f"\n=== exponent in Omega of P(error)")
    summary = {}
    for eps in args.eps_fracs:
        rs = [r for r in rows if abs(r["eps_frac"] - eps) < 1e-12]
        om = np.array([r["omega"] for r in rs], float)
        line = {}
        for key in ["p_cme"] + [f"p_{m}" for m in MODES]:
            y = np.array([r[key] for r in rs])
            good = np.isfinite(y) & (y > 0)
            if good.sum() < 2:
                line[key] = None
                print(f"  eps={eps}: {key:>18} -> categorical failure "
                      f"({good.sum()}/{len(rs)} cells positive)")
                continue
            p = float(np.polyfit(np.log(om[good]), np.log(y[good]), 1)[0])
            line[key] = p
            print(f"  eps={eps}: {key:>18} -> exponent {p:>8.4f}"
                  f"   ({good.sum()}/{len(rs)} positive)")
        summary[str(eps)] = line

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"n": n, "width": width, "rows": rows,
                                    "exponents": summary}, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
