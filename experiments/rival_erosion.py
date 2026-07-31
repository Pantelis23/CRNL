"""T15-a-i: why `rivals-only` never fails, and the cell where it should.

§30 (`pairwise_identity.py`) confirmed the pairwise identity at n = 2..6 to 4.4e-16
and with it P2: under `bookkeeping-only` no trajectory flipped any pairwise sign in
40,000 attempts, approaching to 9.1e-06 of the initial gap without ever crossing.

**P3 FAILED and this experiment is the follow-up.** `rivals-only` was predicted to be
a mere sampling zero, since d2 = (0,1,-1) has a nonzero component along every pairwise
difference. It is not: in 440,000 trajectories the champion's margin never flipped and
never came closer than 0.51 of its initial gap. That is a barrier, not rarity, so
there is a mechanism the sign-conservation argument missed.

THE MECHANISM, derived from the same identity rather than guessed. Under `rivals-only`
the noise is (0, +h, -h, 0), so n_1, n_2 + n_3 and n_B are all noise-free and the
champion's MEAN margin u = n_1 - (n_2 + n_3)/2 carries no noise at all. Its drift
follows from differencing two brackets:

    G_13 - G_12 = (1/Omega)[n_3 - n_2 - gamma*(n_3 - n_2)] = -(1 - gamma)*delta_23/Omega

and therefore, with Gbar = (G_12 + G_13)/2 and delta_23 = n_2 - n_3,

    du/dt = u * Gbar  -  (1 - gamma) * delta_23^2 / (4 Omega)

**u has an ADDITIVE term, so u is NOT sign-conserved and the champion is not protected
by a conservation law.** The term is negative for every gamma < 1: rival spread erodes
the champion's mean margin at a rate QUADRATIC in that spread. This is §24.2's P7
order-statistic intuition -- "the champion loses to the best rival, not to the rival
mean" -- as an exact identity rather than a hand-wave.

So `rivals-only` fails only if delta_23 can GROW, and by the identity delta_23 grows
exactly when its own bracket is positive:

    G_23 = (1/Omega) * [ n_B - n_1 - gamma*(n_2 + n_3 - 1) ]

which is negative whenever the champion is well ahead of the blank pool -- the case in
every cell §24.2 and §30 happened to run. **The cells were chosen by the barrier, and
the sign of G_23 came along for the ride** (rule 9: constancy along the axis you
happened to sweep is not constancy). At small eps the champion is closer to the pool
and G_23 turns POSITIVE: at Omega = 90 it runs +5.51, +3.75, +2.00, +0.24, -1.52 across
eps = 0.03 .. 0.18, changing sign inside a range this project routinely sweeps.

PREDICTIONS, written before running:

  P1  The u-identity above is exact -- relative residual < 1e-13 against the network's
      own drift, at n = 3 over random states and several gamma. If it is not, the
      derivation is wrong and nothing below means anything.
  P2  THE CONTROL, and the strongest form of §30's theorem. `bookkeeping-only` returns
      EXACTLY 0 with zero pairwise flips in EVERY cell, including the shallowest one
      where the champion leads by a single count and full noise fails ~40% of the time.
      A conservation law does not care how low the barrier is. If bookkeeping-only ever
      returns nonzero, the theorem is false and §29/§30 both fall.
  P3  THE TEST. `rivals-only` becomes NONZERO in cells with G_23 > 0 and stays 0 in
      cells with G_23 < 0 -- the SIGN OF THE BRACKET predicts it, not the retained
      variance (which barely moves across these cells) and not the barrier height.
  P4  If `rivals-only` is nonzero at G_23 < 0 too, the sign criterion is wrong and its
      zero was only ever barrier height; §24.2's arm would then be uninformative rather
      than structural. If it stays 0 even at G_23 = +5.5, there is a further constraint
      on the d2-only projection that neither the pairwise identity nor the u-identity
      captures, and I do not have a candidate for one.
  P5  Rule 18, paired: every arm and the CME reference run on the SAME start state per
      cell, and ratios are quoted within a cell. Nothing is compared across cells.

  Rule 10 note: as in §30 the CLE step rejects negative candidates. At small eps the
  counts are smaller and rejections should RISE, so the rejection fraction is reported
  per arm per cell -- if `rivals-only` turns nonzero exactly where rejections spike,
  the harness is a suspect for the effect and the reading is not admissible.
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
from experiments.noise_placement_nwinner import setup_skewed
from experiments.pairwise_identity import bracket, drift, run_tracked

ARMS = ("full", "bookkeeping-only", "rivals-only")


def u_residual(gamma: float, omega: int, rng, n_states: int) -> float:
    """Worst relative residual of du/dt = u*Gbar - (1-gamma)*delta_23^2/(4 Omega)."""
    n = 3
    comp = compile_network(n_winner_reversible(n, gamma), float(omega))
    worst = 0.0
    for _ in range(n_states):
        cuts = np.sort(rng.integers(0, omega + 1, size=n))
        c = np.diff(np.concatenate([[0], cuts, [omega]])).astype(np.int64)
        b, a = drift(comp, c)
        lhs = float(b[0] - 0.5 * (b[1] + b[2]))
        u = float(c[0]) - 0.5 * float(c[1] + c[2])
        d23 = float(c[1] - c[2])
        gbar = 0.5 * (bracket(n, gamma, omega, c, 0, 1)
                      + bracket(n, gamma, omega, c, 0, 2))
        rhs = u * gbar - (1.0 - gamma) * d23 ** 2 / (4.0 * omega)
        scale = max(abs(b[0]), abs(b[1]), abs(b[2]), 1e-300)
        worst = max(worst, abs(lhs - rhs) / scale)
    return worst


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--gamma-frac", type=float, default=0.60)
    ap.add_argument("--cells", type=str, nargs="+",
                    default=["90:0.03", "90:0.05", "90:0.08", "90:0.12",
                             "90:0.18", "90:0.25", "60:0.03", "60:0.25"],
                    help="Omega:eps_frac pairs, chosen to straddle G_23 = 0")
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--rival-skew", type=int, default=2)
    ap.add_argument("--trials", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=20260801)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/rival_erosion.json"))
    args = ap.parse_args()

    t0 = time.time()
    n = args.n
    g = args.gamma_frac * gamma_critical(n)
    width = landscape_width(n, g)
    rng = np.random.default_rng(args.seed)

    print("=== P1: du/dt = u*Gbar - (1-gamma)*delta_23^2/(4 Omega)")
    ures = {}
    for gf in (0.10, 0.40, 0.60, 0.90):
        gg = gf * gamma_critical(3)
        r = u_residual(gg, 97, rng, 60)
        ures[f"{gf}"] = {"gamma": gg, "worst_rel": r}
        print(f"  gamma = {gg:.6f} ({gf} x gamma_c)   worst relative residual {r:.3e}")
    wu = max(v["worst_rel"] for v in ures.values())
    print(f"  -> P1 {'HOLDS' if wu < 1e-13 else 'FAILS'} (worst {wu:.3e})\n")

    print(f"n={n} gamma={g:.5f} width={width:.5f} skew={args.rival_skew} "
          f"trials={args.trials}")
    print(f"\n{'Omega':>6}{'eps':>6}{'start':>18}{'G23*Om':>9}{'CME':>9}"
          f"{'full':>9}{'bookkeep':>10}{'rivals':>9}"
          f"{'bk flips':>10}{'ro flips':>10}{'ro min':>9}{'rej ro':>8}")
    rows = []
    for spec in args.cells:
        om_s, eps_s = spec.split(":")
        om, eps = int(om_s), float(eps_s)
        n0, thr, realised = setup_skewed(n, g, om, eps, args.theta, width,
                                         args.rival_skew)
        c = n0.astype(np.int64)
        g23 = bracket(n, g, om, c, 1, 2) * om
        comp = compile_network(n_winner_reversible(n, g), float(om))
        stop = _absorbing_batch(n, thr)
        exact = p_cme(n, g, om, n0, thr)
        got = {}
        for mode in ARMS:
            r = np.random.default_rng(args.seed + 31 * om + int(1000 * eps))
            got[mode] = run_tracked(comp, n0, r, dt=args.dt, stop=stop,
                                    trials=args.trials, t_max=6000.0, mode=mode, n=n)
        bk, ro = got["bookkeeping-only"], got["rivals-only"]
        print(f"{om:>6}{eps:>6.2f}{str(n0.tolist()):>18}{g23:>9.2f}{exact:>9.5f}"
              f"{got['full']['p']:>9.5f}{bk['p']:>10.5f}{ro['p']:>9.5f}"
              f"{bk['flipped_champion']:>10}{ro['flipped_champion']:>10}"
              f"{ro['closest_champion']:>9.4f}{ro['reject_frac']:>8.4f}")
        rows.append({"omega": om, "eps_frac": eps, "start": n0.tolist(),
                     "thr": int(thr), "realised_eps": realised,
                     "G23_times_omega": float(g23), "p_cme": exact,
                     **{m: got[m] for m in ARMS}})

    print(f"\n=== P2 (control): bookkeeping-only must be exactly 0 everywhere")
    bad = [r for r in rows if r["bookkeeping-only"]["flipped_champion"] > 0
           or r["bookkeeping-only"]["wrong"] > 0]
    shallow = min(rows, key=lambda r: r["p_cme"] if r["p_cme"] > 0 else 1.0)
    deepest = max(rows, key=lambda r: r["p_cme"])
    print(f"  cells with any bookkeeping-only flip or error: {len(bad)}  "
          f"-> P2 {'HOLDS' if not bad else 'FAILS'}")
    print(f"  hardest case for it: CME = {deepest['p_cme']:.4f} at "
          f"Omega={deepest['omega']} eps={deepest['eps_frac']}, full noise fails "
          f"{deepest['full']['p']:.4f} of the time and bookkeeping-only still returns "
          f"{deepest['bookkeeping-only']['p']:.5f}")

    print(f"\n=== P3: does sign(G_23) predict whether rivals-only can fail?")
    print(f"{'Omega':>6}{'eps':>6}{'G23*Om':>9}{'sign':>7}{'rivals P':>10}"
          f"{'champ flips':>13}{'closest':>9}")
    for r in rows:
        s = "+" if r["G23_times_omega"] > 0 else "-"
        print(f"{r['omega']:>6}{r['eps_frac']:>6.2f}{r['G23_times_omega']:>9.2f}"
              f"{s:>7}{r['rivals-only']['p']:>10.5f}"
              f"{r['rivals-only']['flipped_champion']:>13}"
              f"{r['rivals-only']['closest_champion']:>9.4f}")
    pos = [r for r in rows if r["G23_times_omega"] > 0]
    neg = [r for r in rows if r["G23_times_omega"] < 0]
    pos_live = [r for r in pos if r["rivals-only"]["flipped_champion"] > 0]
    neg_live = [r for r in neg if r["rivals-only"]["flipped_champion"] > 0]
    print(f"\n  G_23 > 0 cells with champion flips: {len(pos_live)}/{len(pos)}")
    print(f"  G_23 < 0 cells with champion flips: {len(neg_live)}/{len(neg)}")
    if pos and neg:
        if len(pos_live) == len(pos) and not neg_live:
            print("  -> P3 HOLDS: the sign of the rival bracket, not the variance and "
                  "not the barrier, decides whether rivals-only can fail at all.")
        elif neg_live:
            print("  -> P3 FAILS as stated: rivals-only also fails where G_23 < 0, so "
                  "its zero was barrier height rather than structure (P4).")
        else:
            print("  -> P3 FAILS the other way: rivals-only stays frozen even where "
                  "delta_23 must grow. Neither identity accounts for that (P4).")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"gamma": g, "width": width,
                                    "u_identity": ures, "rows": rows},
                                   indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
