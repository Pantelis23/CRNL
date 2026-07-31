"""T15-a: does §29's multiplicative identity hold at n >= 3, and which arms does it cover?

§29 found that §24.1's `s-only` zero is a THEOREM at n = 2: the drift of the signal
coordinate carries no additive term,

    b_delta = delta * [ c_het*B - c_hom*(s-1) ]

so under any noise projection that leaves `delta` noise-free, `d(delta)/dt = delta*g(t)`
and sign(delta) is conserved exactly. No number of trajectories could have found a
crossing. That bounded what §24 may claim: the zero is a statement about AM's drift
structure, not about coarse-graining in general.

T15-a asks whether the same structure explains §24.2's zeros at n = 3. Worked out by
hand first, from the count-level propensities of `n_winner_reversible`, the pairwise
difference of ANY two committed species obeys

    d(n_i - n_j)/dt = (n_i - n_j) * (k/Omega) * [ n_B - sum_{l != i,j} n_l
                                                  - gamma*(n_i + n_j - 1) ]

for every n and every gamma. Three cancellations do it: the disagreement reaction of
the pair itself consumes both equally; its reverse `2B -> X_i + X_j` produces X_i and
X_j through the same (n-1) pairs each, so the whole gamma*n_B^2 term drops; and the
remaining terms are all bilinear in (n_i - n_j). At n = 2 the middle sum is empty and
this reduces to §29's identity exactly.

**IF THAT IS RIGHT IT IS STRONGER THAN §29 AND IT SORTS §24.2's ARMS.** The conserved
quantity is not "the signal coordinate" but sign(n_i - n_j) for EVERY pair
independently, so a projection is covered by the theorem exactly when it leaves some
pairwise difference direction noise-free:

  * `bookkeeping-only` (n = 3) zeroes all C(n,2) difference directions -> the ENTIRE
    ORDERING is frozen -> P(error) = 0 is a theorem, as `s-only` is at n = 2.
  * `decision-only` zeroes d2 = (0,1,-1) -> sign(n_2 - n_3) is frozen -> rivals that
    start tied stay tied forever. **That is §24.2's Omega-parity trap**, which was
    diagnosed empirically and blamed on integer rounding; the rounding only chose the
    initial condition, the conservation law is what made it fatal.
  * `rivals-only` keeps ONLY d2, whose projection has a NONZERO component along every
    one of the three pairwise differences (delta_12 gets -h, delta_13 gets +h,
    delta_23 gets 2h). **No pairwise sign is conserved, so the theorem does not cover
    it** -- and §24.2's stated reason for its zero, "X1 versus the rival mean is
    deterministic and X1 never loses", is a non-sequitur: the rival MEAN being
    noise-free does not stop one rival from breaking away from the other and passing
    X1, since delta_23 is exactly the direction that is being driven.

PREDICTIONS, written before running:

  P1  The identity holds to machine precision (relative residual < 1e-13) for
      n = 2..6 at several gamma, measured against the network's own stoichiometry and
      propensities -- NOT a hand-rolled copy of them, since the point is whether the
      shipped conventions satisfy it. And the bracket is independent of (n_i - n_j):
      varying the split at fixed n_i + n_j must leave b_ij/(n_i - n_j) constant.
  P2  Dynamical consequence at n = 3: under `bookkeeping-only` NO trajectory ever
      flips any pairwise sign -- not merely "none finishes wrong", which is what §24.2
      measured. This is the stronger statement and the one the theorem predicts.
  P3  THE DISCRIMINATING PREDICTION. Under `rivals-only` the theorem does NOT apply,
      so pairwise signs must be flippable. I expect a nonzero fraction of trajectories
      to drive delta_12 or delta_13 to zero, and P(error) to be nonzero once the trial
      count is raised far enough above §24.2's 40,000. If instead `rivals-only` stays
      exactly 0 with every pairwise sign intact at 10^6 trajectories, there is a
      conservation law here I have not found, my reading of d2 is wrong, and P3 fails
      in the way that costs the most.
  P4  `full` reproduces the exact CME at the same start state, or nothing above is
      admissible. Same gate as §24.1 and §29.

  Note on rule 10: the CLE step rejects any candidate with a negative count, which is
  a harness action the chemistry does not take. Rejections are counted and reported
  per arm, because a rejection could in principle mask a crossing.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.approximations import propensities_batch
from crnl.networks.n_winner_reversible import gamma_critical, n_winner_reversible
from crnl.vectorized import compile_network
from experiments.approximation_hierarchy_nwinner import (
    _absorbing_batch, landscape_width, p_cme,
)
from experiments.noise_placement_nwinner import project, setup_skewed

ARMS = ("full", "bookkeeping-only", "rivals-only", "decision-only")


def bracket(n: int, gamma: float, omega: float, counts: np.ndarray,
            i: int, j: int, k: float = 1.0) -> float:
    """The claimed multiplicative factor for pair (i, j), in count units."""
    n_b = counts[n]
    others = counts[:n].sum() - counts[i] - counts[j]
    return (k / omega) * (n_b - others - gamma * (counts[i] + counts[j] - 1.0))


def drift(comp, counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """dn/dt from the network's own stoichiometry and propensities.

    Returns the drift and, per species pair, the TOTAL ABSOLUTE TRAFFIC
    `sum_r |S_ir - S_jr| * a_r` that the identity requires to cancel down to
    something proportional to (n_i - n_j). That traffic is the only honest
    denominator for the residual: `max(|lhs|, |rhs|)` is degenerate exactly at
    n_i = n_j, where both sides are zero and any float noise reads as 100% error.
    """
    a = propensities_batch(comp, counts[None, :].astype(float))[0]
    S = comp.S
    return S @ a, a


def check_identity(n: int, gamma: float, omega: int, rng, n_states: int) -> dict:
    net = n_winner_reversible(n, gamma)
    comp = compile_network(net, float(omega))
    worst_abs, worst_rel, worst_where = 0.0, 0.0, None
    for _ in range(n_states):
        cuts = np.sort(rng.integers(0, omega + 1, size=n))
        counts = np.diff(np.concatenate([[0], cuts, [omega]])).astype(np.int64)
        assert counts.sum() == omega and len(counts) == n + 1
        b, a = drift(comp, counts)
        S = comp.S
        for i in range(n):
            for j in range(i + 1, n):
                d = float(counts[i] - counts[j])
                lhs = float(b[i] - b[j])
                rhs = d * bracket(n, gamma, omega, counts, i, j)
                err = abs(lhs - rhs)
                traffic = float(np.abs(S[i] - S[j]) @ a)
                scale = max(traffic, 1e-300)
                if err > worst_abs:
                    worst_abs = err
                if err / scale > worst_rel:
                    worst_rel, worst_where = err / scale, (counts.tolist(), i, j)
    return {"n": n, "gamma": gamma, "omega": omega, "states": n_states,
            "worst_abs": worst_abs, "worst_rel": worst_rel, "at": worst_where}


def check_independence(n: int, gamma: float, omega: int, rng, n_states: int) -> dict:
    """b_ij/(n_i - n_j) must not move when the split changes at fixed n_i + n_j."""
    net = n_winner_reversible(n, gamma)
    comp = compile_network(net, float(omega))
    worst = 0.0
    for _ in range(n_states):
        cuts = np.sort(rng.integers(0, omega + 1, size=n))
        base = np.diff(np.concatenate([[0], cuts, [omega]])).astype(np.int64)
        i, j = 0, 1
        pair_total = int(base[i] + base[j])
        if pair_total < 4:
            continue
        ratios = []
        for split in range(1, pair_total):
            if split * 2 == pair_total:
                continue                       # b_ij = 0 there; ratio is 0/0
            c = base.copy()
            c[i], c[j] = split, pair_total - split
            b, _ = drift(comp, c)
            ratios.append(float(b[i] - b[j]) / float(c[i] - c[j]))
        if len(ratios) >= 2:
            r = np.array(ratios)
            worst = max(worst, float(np.ptp(r) / max(abs(r.mean()), 1e-300)))
    return {"n": n, "gamma": gamma, "omega": omega, "worst_rel_spread": worst}


def run_tracked(comp, n0, rng, *, dt, stop, trials, t_max, mode, n,
                max_steps=300_000) -> dict:
    """Projected-noise CLE that also watches every pairwise sign, not just the winner."""
    S = comp.S.astype(float)
    m = np.tile(np.asarray(n0, dtype=float), (trials, 1))
    live = np.ones(trials, bool)
    t = np.zeros(trials)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    sgn0 = np.array([np.sign(float(n0[i] - n0[j])) for i, j in pairs])
    if (sgn0 == 0).any():
        raise ValueError(f"start state has tied species; pair signs undefined: {n0}")
    # smallest value reached by sgn0 * (n_i - n_j), per trajectory per pair
    worst = np.tile(np.abs(np.array([float(n0[i] - n0[j]) for i, j in pairs])),
                    (trials, 1))
    flipped = np.zeros((trials, len(pairs)), bool)
    var_full = var_kept = 0.0
    rejected = steps_taken = 0
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
        steps_taken += len(idx)
        upd = idx[ok]
        m[upd] = cand[ok]
        t[upd] += dt
        cur = np.stack([m[upd][:, i] - m[upd][:, j] for i, j in pairs], axis=1) * sgn0
        worst[upd] = np.minimum(worst[upd], cur)
        flipped[upd] |= (cur <= 0.0)
        done = stop(m[idx]) | (t[idx] >= t_max)
        live[idx[done]] = False
    fin = stop(m)
    nok = int(fin.sum())
    wrong = int((np.argmax(m[fin][:, :n], axis=1) != 0).sum()) if nok else 0
    gap0 = np.abs(np.array([float(n0[i] - n0[j]) for i, j in pairs]))
    frac = worst / gap0
    champ = [k for k, (i, j) in enumerate(pairs) if i == 0]   # champion vs each rival
    rival = [k for k, (i, j) in enumerate(pairs) if i != 0]   # rival vs rival
    return {"p": wrong / nok if nok else float("nan"), "n_ok": nok, "wrong": wrong,
            "unfinished": int(live.sum()), "rejected": rejected,
            "reject_frac": rejected / steps_taken if steps_taken else float("nan"),
            "pairs": [f"{i+1}-{j+1}" for i, j in pairs],
            "flipped_any": int(flipped.any(axis=1).sum()),
            "flipped_champion": int(flipped[:, champ].any(axis=1).sum()),
            "flipped_rival": int(flipped[:, rival].any(axis=1).sum()) if rival else 0,
            "flipped_per_pair": [int(c) for c in flipped.sum(axis=0)],
            "closest_champion": float(frac[:, champ].min()),
            "closest_rival": float(frac[:, rival].min()) if rival else float("nan"),
            "median_closest_champion": float(np.median(frac[:, champ].min(axis=1))),
            "variance_kept_frac": var_kept / var_full if var_full > 0 else float("nan")}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ns", type=int, nargs="+", default=[2, 3, 4, 5, 6])
    ap.add_argument("--gamma-fracs", type=float, nargs="+",
                    default=[0.10, 0.40, 0.60, 0.90])
    ap.add_argument("--id-omega", type=int, default=97)
    ap.add_argument("--id-states", type=int, default=40)
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--gamma-frac", type=float, default=0.60)
    ap.add_argument("--omega", type=int, default=30)
    ap.add_argument("--eps-frac", type=float, default=0.25)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--rival-skew", type=int, default=2)
    ap.add_argument("--trials", type=int, default=40000)
    ap.add_argument("--deep-trials", type=int, default=0,
                    help="extra trials for rivals-only alone (P3). 0 skips.")
    ap.add_argument("--seed", type=int, default=20260801)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/pairwise_identity.json"))
    args = ap.parse_args()

    t0 = time.time()
    rng = np.random.default_rng(args.seed)

    print("=== P1: d(n_i - n_j)/dt = (n_i - n_j) * (1/Omega) * "
          "[n_B - sum_{l!=i,j} n_l - gamma*(n_i+n_j-1)]")
    print(f"{'n':>3}{'gamma':>10}{'pairs':>8}{'worst abs':>13}{'worst rel':>13}")
    ident = []
    for n in args.ns:
        gc = gamma_critical(n)
        for gf in args.gamma_fracs:
            g = gf * gc
            r = check_identity(n, g, args.id_omega, rng, args.id_states)
            r["gamma_frac"] = gf
            ident.append(r)
            npair = args.id_states * n * (n - 1) // 2
            print(f"{n:>3}{g:>10.6f}{npair:>8}{r['worst_abs']:>13.3e}"
                  f"{r['worst_rel']:>13.3e}")
    worst_rel = max(r["worst_rel"] for r in ident)
    print(f"\n  worst relative residual over all {len(ident)} (n, gamma) cells: "
          f"{worst_rel:.3e}   -> P1 {'HOLDS' if worst_rel < 1e-13 else 'FAILS'}")

    print("\n=== P1b: is the bracket independent of the difference itself?")
    print(f"{'n':>3}{'gamma':>10}{'worst spread of b_ij/(n_i-n_j)':>34}")
    indep = []
    for n in args.ns:
        g = 0.60 * gamma_critical(n)
        r = check_independence(n, g, args.id_omega, rng, max(4, args.id_states // 4))
        indep.append(r)
        print(f"{n:>3}{g:>10.6f}{r['worst_rel_spread']:>34.3e}")
    worst_sp = max(r["worst_rel_spread"] for r in indep)
    print(f"\n  worst: {worst_sp:.3e}  -> the factor does not depend on the coordinate "
          f"it multiplies")

    n = args.n
    g = args.gamma_frac * gamma_critical(n)
    width = landscape_width(n, g)
    if not width > 1e-6:
        raise SystemExit(f"landscape_width = {width:.6g}: no landscape, nothing to do.")
    n0, thr, realised = setup_skewed(n, g, args.omega, args.eps_frac, args.theta,
                                     width, args.rival_skew)
    comp = compile_network(n_winner_reversible(n, g), float(args.omega))
    stop = _absorbing_batch(n, thr)
    exact = p_cme(n, g, args.omega, n0, thr)
    print(f"\n=== P2/P3/P4: n={n} gamma={g:.5f} Omega={args.omega} "
          f"eps={args.eps_frac} skew={args.rival_skew}")
    print(f"  start {n0.tolist()}  threshold {thr}  realised eps {realised:.5f}")
    print(f"  exact CME P(error) = {exact:.6f}")
    print(f"\n  flips are reported SEPARATELY for champion-vs-rival pairs (1-2, 1-3)"
          f" and the\n  rival-vs-rival pair (2-3): only the first can make the "
          f"observable wrong.")
    print(f"\n{'arm':>18}{'var':>7}{'P(err)':>10}{'wrong':>7}{'flip champ':>12}"
          f"{'flip rival':>12}{'min champ/d0':>14}{'rej/step':>10}{'unfin':>7}")
    dyn = {}
    for mode in ARMS:
        r = np.random.default_rng(args.seed + 991)
        res = run_tracked(comp, n0, r, dt=args.dt, stop=stop, trials=args.trials,
                          t_max=6000.0, mode=mode, n=n)
        dyn[mode] = res
        print(f"{mode:>18}{res['variance_kept_frac']:>7.3f}"
              f"{res['p']:>10.6f}{res['wrong']:>7}{res['flipped_champion']:>12}"
              f"{res['flipped_rival']:>12}{res['closest_champion']:>14.5f}"
              f"{res['reject_frac']:>10.4f}{res['unfinished']:>7}")

    bk = dyn["bookkeeping-only"]
    print(f"\n  P2: bookkeeping-only flipped {bk['flipped_any']} signs of any kind in "
          f"{args.trials} trajectories.\n      Closest any pair came to crossing: "
          f"{min(bk['closest_champion'], bk['closest_rival']):.2e} of its initial gap "
          f"-- it approaches zero\n      but cannot reach it, which is the theorem and "
          f"not mere distance.   -> {'HOLDS' if bk['flipped_any'] == 0 else 'FAILS'}")
    do = dyn["decision-only"]
    print(f"  P2b: decision-only froze the RIVAL pair only -- rival flips "
          f"{do['flipped_rival']}, champion flips {do['flipped_champion']}. That is "
          f"§24.2's parity trap\n       restated as a conservation law.")
    ro = dyn["rivals-only"]
    print(f"  P3: rivals-only champion flips {ro['flipped_champion']} of "
          f"{args.trials}, closest {ro['closest_champion']:.4f}, "
          f"median trajectory {ro['median_closest_champion']:.4f}")
    print(f"  P4 gate: full {dyn['full']['p']:.6f} against exact {exact:.6f}  "
          f"(ratio {dyn['full']['p']/exact:.3f})")

    deep = None
    if args.deep_trials:
        print(f"\n=== P3 deep: rivals-only at {args.deep_trials} trajectories")
        r = np.random.default_rng(args.seed + 7777)
        deep = run_tracked(comp, n0, r, dt=args.dt, stop=stop,
                           trials=args.deep_trials, t_max=6000.0,
                           mode="rivals-only", n=n)
        print(f"  P(error) = {deep['p']:.3e}  ({deep['wrong']} wrong of "
              f"{deep['n_ok']} absorbed)   champion flips {deep['flipped_champion']}")
        print(f"  closest champion approach {deep['closest_champion']:.4f}, "
              f"rejections/step {deep['reject_frac']:.4f}, "
              f"unfinished {deep['unfinished']}")
        if deep["wrong"] == 0 and deep["flipped_champion"] == 0:
            print("  -> P3 FAILS: no crossing at all. The theorem does not cover this "
                  "arm, so something else conserves the ordering here and the reading "
                  "of d2 above is wrong.")
        else:
            print("  -> P3 holds: rivals-only is a SAMPLING zero, not a theorem zero, "
                  "unlike bookkeeping-only.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(
        {"identity": ident, "independence": indep,
         "cell": {"n": n, "gamma": g, "omega": args.omega, "eps_frac": args.eps_frac,
                  "skew": args.rival_skew, "start": n0.tolist(), "thr": int(thr),
                  "realised_eps": realised, "p_cme": exact},
         "dynamics": dyn, "deep": deep}, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
