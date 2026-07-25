"""Is the radix penalty convention-dependent? (FINDINGS open question 2)

§3 measured a radix penalty under a FIXED PAIRWISE MARGIN — the champion leads
each rival by δ = 0.10 — and flagged its own scope: "A different convention
(fixed champion *share*, or symmetric plurality) could well be unbounded; that is
untested." This tests it.

    margin: f_champ = 1/n + δ(n−1)/n,  f_other = 1/n − δ/n   → champ share → δ
    share : f_champ = s,               f_other = (1−s)/(n−1) → champ share = s

Both anchored at n=2 to the same 0.55 / 0.45 start, so δ = 0.10 and s = 0.55.

THE ANSWER IS YES, AND THE MECHANISM IS MUNDANE. Under fixed share the penalty
does not merely weaken, it vanishes: P(champion wins) = 1.0000 at every n ≥ 3.
But fixed share is not an alternative way of holding the contest difficulty
constant — it hands the champion a pairwise lead that GROWS with n:

    n            2      3      4      8     16     24
    margin conv  0.100  0.100  0.100  0.100  0.100  0.100
    share conv   0.100  0.325  0.400  0.486  0.520  0.530

a 5.3× larger lead by n=24. So the two conventions do not disagree about physics;
they ask different questions, and only the fixed-margin one isolates the effect of
the alphabet size. §3's choice was the right one and its finding stands.

A NOTE ON A PREDICTION THAT FAILED. `THEORIES.md` T3 predicted the opposite —
that fixed share would give an UNBOUNDED penalty, "because the champion's
per-rival lead keeps shrinking". The lead is `s − (1−s)/(n−1)`, which grows
toward `s`. The reasoning was wrong before any code ran, which is why the
prediction is now written down before the sweep rather than after.

Nor is share itself the governing variable: at a FIXED share of 0.50, P(win) is
0.606 at n=2 and 0.997 at n=3, because the same share means a pairwise lead of
~0 at n=2 and 0.25 at n=3.

    python -m experiments.radix_convention
    python -m experiments.radix_convention --quick
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

from crnl.networks import n_winner
from crnl.stochastic import seed_for
from crnl.vectorized import compile_network, gillespie_fast


def _quantise(exact: np.ndarray, omega: int) -> np.ndarray:
    """Largest-remainder rounding, then guarantee a strict champion lead.

    Same guard as radix_wall.champion_counts: floor is monotonic, so rivals can
    at most tie the champion; if one does, move a single molecule. Without it a
    tie at t=0 makes 'champion wins' a coin flip that looks like physics.
    """
    counts = np.floor(exact).astype(np.int64)
    short = int(omega - counts.sum())
    for i in np.argsort(-(exact - counts))[:short]:
        counts[i] += 1
    if counts[1:].size and counts[1:].max() >= counts[0]:
        donor = 1 + int(np.argmax(counts[1:]))
        counts[donor] -= 1
        counts[0] += 1
    return np.concatenate([counts, [0]])


def counts_margin(n: int, omega: int, delta: float = 0.10) -> np.ndarray:
    """Champion leads EACH rival by delta. Champion's share -> delta as n grows."""
    exact = np.array([1 / n + delta * (n - 1) / n]
                     + [1 / n - delta / n] * (n - 1)) * omega
    return _quantise(exact, omega)


def counts_share(n: int, omega: int, share: float = 0.55) -> np.ndarray:
    """Champion holds a fixed share. Its pairwise lead GROWS toward `share`."""
    exact = np.array([share] + [(1 - share) / (n - 1)] * (n - 1)) * omega
    return _quantise(exact, omega)


def pairwise_lead(counts: np.ndarray, n: int, omega: int) -> float:
    return float((counts[0] - counts[1:n].max()) / omega)


def p_champion_wins(n: int, omega: int, counts: np.ndarray, trials: int,
                    seed: int = 0) -> dict:
    """P(champion wins | a decision was reached). Conditional, because the
    all-blank bin is a genuine outcome of irreversible AM (FINDINGS §1)."""
    net = n_winner(n)
    comp = compile_network(net, float(omega))
    names = list(net.species)
    wins = decided = blank = 0
    for t in range(trials):
        r = gillespie_fast(comp, counts, seed_for(omega, t, base=seed),
                           species=names)
        final = r.n_final[:n]
        if final.sum() == 0:
            blank += 1
            continue
        decided += 1
        wins += int(np.argmax(final) == 0)
    return {"p_win": wins / decided if decided else float("nan"),
            "decided_frac": decided / trials, "blank": blank, "trials": trials}


def main():
    here = os.path.dirname(__file__)
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ns", type=int, nargs="+", default=[2, 3, 4, 6, 8, 12, 16, 24])
    p.add_argument("--omega", type=int, default=120)
    p.add_argument("--delta", type=float, default=0.10)
    p.add_argument("--share", type=float, default=0.55)
    p.add_argument("--trials", type=int, default=3000)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--data", default=os.path.join(
        here, os.pardir, "results", "radix_convention.json"))
    args = p.parse_args()
    if args.quick:
        args.ns, args.trials = [2, 4, 8], 800

    print(f"Radix penalty vs convention  (Ω={args.omega}, {args.trials} trials, "
          f"δ={args.delta}, s={args.share}; both anchored to 0.55/0.45 at n=2)")
    print(f"\n{'n':>4} {'P(win) margin':>14} {'lead':>6}  {'P(win) share':>13} "
          f"{'lead':>6}")
    rows = []
    for n in args.ns:
        cm = counts_margin(n, args.omega, args.delta)
        cs = counts_share(n, args.omega, args.share)
        rm = p_champion_wins(n, args.omega, cm, args.trials, args.seed)
        rs = p_champion_wins(n, args.omega, cs, args.trials, args.seed)
        lm, ls = pairwise_lead(cm, n, args.omega), pairwise_lead(cs, n, args.omega)
        rows.append({"n": n, "margin": {**rm, "lead": lm,
                                        "champ_share": float(cm[0] / args.omega)},
                     "share": {**rs, "lead": ls,
                               "champ_share": float(cs[0] / args.omega)}})
        print(f"{n:>4} {rm['p_win']:>14.4f} {lm:>6.3f}  {rs['p_win']:>13.4f} "
              f"{ls:>6.3f}")

    m0, m1 = rows[0]["margin"]["p_win"], rows[-1]["margin"]["p_win"]
    s0, s1 = rows[0]["share"]["p_win"], rows[-1]["share"]["p_win"]
    print(f"\nfixed margin: P(win) {m0:.4f} -> {m1:.4f} over n={args.ns[0]}..{args.ns[-1]}")
    print(f"fixed share : P(win) {s0:.4f} -> {s1:.4f}")
    print("\nThe conventions disagree because only ONE holds the contest fixed. The")
    print("share convention's pairwise lead grows "
          f"{rows[-1]['share']['lead'] / rows[0]['share']['lead']:.1f}x across "
          "this range, so it\nis asking an easier question at every n, not "
          "measuring the same one differently.")

    os.makedirs(os.path.dirname(os.path.abspath(args.data)), exist_ok=True)
    with open(args.data, "w") as fh:
        json.dump({"omega": args.omega, "delta": args.delta, "share": args.share,
                   "trials": args.trials, "rows": rows}, fh, indent=2)
    print(f"\nwrote data -> {args.data}")


if __name__ == "__main__":
    main()
