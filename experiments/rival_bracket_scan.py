"""T15-a-ii: is `rivals-only` gated by the rival bracket, or just by the barrier?

§30.1 proposed that `rivals-only` can only fail when the rival-vs-rival bracket

    G_23 = (1/Omega) * [ n_B - n_1 - gamma*(n_2 + n_3 - 1) ]

is positive, since by §30's identity delta_23 grows exactly when G_23 > 0, and by the
u-identity the champion's noise-free mean margin is eroded at a rate proportional to
delta_23^2. The eps sweep ordered `rivals-only` over four decades in exactly that
sequence -- +5.51 -> 0.935, +3.75 -> 0.426, +2.00 -> 0.109, +0.24 -> 0.012,
-1.52 -> 0.0002, -5.03 -> 0.

**THAT SWEEP CANNOT SUPPORT THE CLAIM AND I SHOULD NOT HAVE READ IT AS IF IT COULD.**
eps sets the champion's margin, which sets BOTH the barrier and G_23, so the two move
together by construction and the ordering is equally consistent with "`rivals-only`
fails when the barrier is low, like every other arm". This is rule 9 -- constancy
along the axis you happened to sweep is not constancy -- and the project has been
wrong this way three times.

**The one accidental matched-barrier pair in §30.1 already argues against the claim:**
Omega=60/eps=0.03 and Omega=90/eps=0.05 have barriers matched to 0.42% (CME 0.48003
against 0.47804) and G_23 differing by 3.1x (+1.21 against +3.75), yet the paired
ratio `rivals-only / full` is 0.8908 against 0.8940 -- a difference of 0.36%. One pair
is not a test, and it also varies Omega, so this experiment separates the axes on
purpose.

THE SEPARATION. G_23 depends on n_B, which the champion's margin does not. Holding the
margin m = n_1 - max(rivals) FIXED at 7 and the skew fixed at 2, the start state is
determined by n_B alone:

    3*R = Omega - n_B - m + skew

so raising n_B lowers R, and n_B - R rises about twice as fast as n_B. Over
n_B = 13..58 at Omega = 90 that drives G_23*Omega from about -23 to about +40 -- a
swing four times wider than the whole eps sweep produced -- **with the champion's
margin never changing**. The barrier does still move (a start far from the attractor
is a different problem), so the exact CME is measured per cell and every comparison is
the paired within-cell ratio `rivals-only / full`, never a raw probability across
cells (rule 18).

PREDICTIONS, written before running:

  P1  IF §30.1's mechanism is right, the paired ratio rises monotonically with G_23
      and collapses toward 0 at the negative end, even between cells whose CME
      barriers are comparable. G_23, not the barrier, would then be the gate.
  P2  IF the ratio instead tracks the CME barrier and is flat or unordered in G_23,
      §30.1's mechanism is DEAD: `rivals-only` is then simply a 6.6%-variance arm that
      fails when failure is easy, and its zero in §24.2 and §30 was barrier height all
      along. **This is what the matched-barrier pair predicts, so it is the outcome I
      now expect**, and saying so before the run is the point.
  P3  The control, again: `bookkeeping-only` returns exactly 0 with zero pairwise
      flips in EVERY cell, including any cell where the champion is barely ahead. §30
      showed this at a CME of 0.597; a conservation law does not care about n_B either.
  P4  If the ratio exceeds 1 anywhere -- `rivals-only` failing MORE than full noise --
      that is not a bug. §30.1 measured 1.5604 at Omega=90, eps=0.03, where the
      champion led by a single count: removing d1's noise removes restoring
      fluctuations as well as destroying ones. It is reported, not clipped.

  Rule 10: the CLE's negative-count rejection is reported per cell. At high n_B the
  committed counts get small and rejections should rise; if the ratio turns over
  exactly where rejections spike, the harness is the suspect and the cell is not
  admissible.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.networks.n_winner_reversible import gamma_critical, n_winner_reversible
from crnl.vectorized import compile_network
from experiments.approximation_hierarchy_nwinner import (
    _absorbing_batch, landscape_width, p_cme,
)
from experiments.pairwise_identity import bracket, run_tracked

ARMS = ("full", "bookkeeping-only", "rivals-only")


def state_for_nb(omega: int, n_b: int, margin: int, skew: int) -> np.ndarray:
    """Start state with n_B chosen and the champion's margin held fixed."""
    num = omega - n_b - margin + skew
    if num <= 0 or num % 3:
        raise ValueError("no integer state")
    r = num // 3
    n0 = np.array([r + margin, r, r - skew, n_b], dtype=np.int64)
    if (n0 <= 0).any():
        raise ValueError("non-positive count")
    assert int(n0.sum()) == omega
    assert int(n0[0] - max(n0[1], n0[2])) == margin
    assert int(n0[1] - n0[2]) == skew
    return n0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--gamma-frac", type=float, default=0.60)
    ap.add_argument("--omega", type=int, default=90)
    ap.add_argument("--margin", type=int, default=7)
    ap.add_argument("--skew", type=int, default=2)
    ap.add_argument("--nbs", type=int, nargs="+",
                    default=[13, 22, 31, 40, 49, 58])
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--trials", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=20260802)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/rival_bracket_scan.json"))
    args = ap.parse_args()

    t0 = time.time()
    n, om = args.n, args.omega
    g = args.gamma_frac * gamma_critical(n)
    width = landscape_width(n, g)
    thr = max(2, int(round(args.theta * width * om)))
    comp = compile_network(n_winner_reversible(n, g), float(om))
    stop = _absorbing_batch(n, thr)

    print(f"n={n} gamma={g:.5f} Omega={om} threshold={thr} "
          f"margin FIXED at {args.margin}, skew {args.skew}, trials {args.trials}")
    print(f"\n{'n_B':>5}{'start':>18}{'G23*Om':>9}{'CME':>9}{'full':>9}"
          f"{'bookkeep':>10}{'rivals':>9}{'rivals/full':>12}{'bk flip':>9}"
          f"{'rej ro':>8}")
    rows = []
    for n_b in args.nbs:
        try:
            n0 = state_for_nb(om, n_b, args.margin, args.skew)
        except ValueError as e:
            print(f"{n_b:>5}   SKIPPED ({e})")
            continue
        g23 = bracket(n, g, om, n0, 1, 2) * om
        exact = p_cme(n, g, om, n0, thr)
        got = {}
        for mode in ARMS:
            r = np.random.default_rng(args.seed + 13 * n_b)
            got[mode] = run_tracked(comp, n0, r, dt=args.dt, stop=stop,
                                    trials=args.trials, t_max=6000.0, mode=mode, n=n)
        f, ro, bk = got["full"]["p"], got["rivals-only"]["p"], got["bookkeeping-only"]
        ratio = ro / f if f > 0 else float("nan")
        print(f"{n_b:>5}{str(n0.tolist()):>18}{g23:>9.2f}{exact:>9.5f}{f:>9.5f}"
              f"{bk['p']:>10.5f}{ro:>9.5f}{ratio:>12.4f}"
              f"{bk['flipped_champion']:>9}{got['rivals-only']['reject_frac']:>8.4f}")
        rows.append({"n_b": n_b, "start": n0.tolist(), "G23_times_omega": float(g23),
                     "p_cme": exact, "ratio": ratio, **{m: got[m] for m in ARMS}})

    print(f"\n=== P3 control: bookkeeping-only")
    bad = [r for r in rows if r["bookkeeping-only"]["flipped_champion"] > 0
           or r["bookkeeping-only"]["wrong"] > 0]
    print(f"  cells with any flip or error: {len(bad)} of {len(rows)}"
          f"   -> P3 {'HOLDS' if not bad else 'FAILS'}")

    print(f"\n=== P1 vs P2: what does the paired ratio follow?")
    gg = np.array([r["G23_times_omega"] for r in rows])
    cme = np.array([r["p_cme"] for r in rows])
    rat = np.array([r["ratio"] for r in rows])
    ok = np.isfinite(rat)
    print(f"  G_23*Omega swing: {gg.min():+.2f} to {gg.max():+.2f}  "
          f"(the eps sweep managed {-5.03:+.2f} to {5.51:+.2f})")
    print(f"  CME barrier range over the same cells: {cme.min():.5f} to {cme.max():.5f}")
    if ok.sum() >= 3:
        c_g = float(np.corrcoef(gg[ok], rat[ok])[0, 1])
        c_c = float(np.corrcoef(cme[ok], rat[ok])[0, 1])
        print(f"  corr(ratio, G_23)  = {c_g:+.4f}")
        print(f"  corr(ratio, CME)   = {c_c:+.4f}")
        mono_g = bool(np.all(np.diff(rat[ok][np.argsort(gg[ok])]) > 0))
        print(f"  ratio monotone increasing in G_23? {'yes' if mono_g else 'NO'}")
        print(f"\n  Cells with comparable CME but different G_23 are the decisive ones:")
        order = np.argsort(cme)
        for a in range(len(order)):
            for b in range(a + 1, len(order)):
                ia, ib = order[a], order[b]
                if abs(cme[ia] - cme[ib]) / max(cme[ib], 1e-12) < 0.15:
                    print(f"    n_B {rows[ia]['n_b']:>3} vs {rows[ib]['n_b']:>3}:  "
                          f"CME {cme[ia]:.5f}/{cme[ib]:.5f}  "
                          f"G_23 {gg[ia]:+.2f}/{gg[ib]:+.2f}  "
                          f"ratio {rat[ia]:.4f}/{rat[ib]:.4f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"gamma": g, "omega": om, "thr": int(thr),
                                    "margin": args.margin, "rows": rows},
                                   indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
