"""A relic ABUNDANCE, which the fixed-drive model could not have.

FINDINGS 5.1 showed the expanding SSA is exactly ordinary SSA stopped at internal
time 1/H. That is forced: every AM reaction is 2 -> 2, so dilution scales every
propensity identically and the landscape is invariant under expansion. Freeze-out
there is a deadline, and what survives is only a SIGN (which species won) plus a
binary decided/undecided -- there is no abundance, because the equilibrium never
moves.

`crnl/cooling.py` lets the medium cool: gamma(s) = gamma0 ** ((1-s)^-w) on the
rescaled internal clock s = H*tau in [0,1). The forward reactions are untouched
and only the reverses are suppressed, so the equilibrium minority abundance falls
toward zero while the reaction rate also falls. Freeze-out is then the classic
race, and it leaves a relic ABOVE the equilibrium value -- which is what a relic
abundance means.

TWO ARMS, differing only in whether the medium cools:

  cooling   gamma0 = 0.55 (just ABOVE gamma_c = 1/2, so there is no landscape at
            the start) with w > 0. Cooling carries gamma down through gamma_c,
            pitchfork happens at a definite time, and the system must choose. The
            symmetry breaking is driven by the expansion itself.
  fixed     gamma = 0.05 held constant (w = 0), the FINDINGS 5.1 model. The
            landscape is there from the start and never changes.

PREDICTIONS, written before running:

  P1  In the cooling arm the relic minority fraction RISES with H: a faster
      expansion freezes the annihilation earlier, leaving more of the minority
      species behind. This is the standard cosmological direction (faster
      expansion -> larger relic).
  P2  The relic EXCEEDS the equilibrium minority fraction at the gamma where it
      froze. If it merely tracked equilibrium there would be no freeze-out, just
      adiabatic following.
  P3  In the fixed arm the relic is set by the chemistry (gamma) and is
      essentially H-INDEPENDENT once the system decides at all. That contrast is
      the whole point: abundance set by expansion versus abundance set by
      chemistry.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.cooling import CoolingSchedule, gillespie_cooling, reverse_mask
from crnl.networks.am_reversible import (
    am_reversible, delta_star, reverse_pairing,
)
from crnl.vectorized import compile_network


def equilibrium_minority(gamma: float) -> float:
    """y*/(x*+y*) at the attractor: the abundance the system would relax to.

    x*+y* = 1/(1+gamma) and x*-y* = delta*(gamma), so the minority share of the
    committed population is (1 - (1+gamma) delta*)/2. It is 1/2 at gamma_c (no
    landscape, equal populations) and falls to 0 as gamma -> 0.
    """
    if gamma >= 0.5:
        return 0.5
    return float((1.0 - (1.0 + gamma) * delta_star(gamma)) / 2.0)


def run_cell(gamma0: float, w: float, omega: int, hubble: float,
             trials: int, seed: int) -> dict:
    net = am_reversible(gamma0)
    mask = reverse_mask(net, reverse_pairing(net))
    sched = CoolingSchedule.build(gamma0, w)
    comp = compile_network(net, float(omega))
    third = omega // 3
    n0 = np.array([third, third, omega - 2 * third], dtype=np.int64)
    rng = np.random.default_rng(seed)

    relics, gammas, undecided, blank = [], [], 0, 0
    for _ in range(trials):
        r = gillespie_cooling(comp, n0, rng, hubble, sched, mask,
                              species=["X", "Y", "B"])
        nx, ny = int(r.n_final[0]), int(r.n_final[1])
        if nx + ny == 0:
            blank += 1
            continue
        minority = min(nx, ny) / (nx + ny)
        # "undecided" keeps FINDINGS 5.1's meaning: never resolved at all. Those
        # runs must NOT enter the abundance -- an undecided run sits near 0.5 and
        # a mean over both populations reports "large relic" for a system that
        # simply never annihilated. The first run of this experiment did exactly
        # that and reported a relic of 0.42 at 65% undecided.
        if minority > 0.40:
            undecided += 1
            continue
        relics.append(minority)
        gammas.append(r.gamma_final)
    if not relics:
        return {}
    gm = float(np.mean(gammas))
    return {
        "gamma0": gamma0, "w": w, "omega": omega, "hubble": hubble,
        "trials": trials,
        "relic_mean": float(np.mean(relics)),
        "n_decided": len(relics),
        "relic_median": float(np.median(relics)),
        "relic_sem": float(np.std(relics) / np.sqrt(len(relics))),
        "gamma_freeze": gm,
        "equilibrium_at_freeze": equilibrium_minority(gm),
        "p_undecided": undecided / trials,
        "p_blank": blank / trials,
    }


def report(rows, title):
    print(f"\n=== {title}")
    print(f"{'H':>7} {'relic':>10} {'+-sem':>9} {'median':>9} {'g_freeze':>10} "
          f"{'eq(g_frz)':>10} {'relic/eq':>9} {'undec':>7}")
    for r in rows:
        eq = r["equilibrium_at_freeze"]
        ratio = r["relic_mean"] / eq if eq > 0 else float("inf")
        print(f"{r['hubble']:>7.3f} {r['relic_mean']:>10.5f} {r['relic_sem']:>9.5f} "
              f"{r['relic_median']:>9.5f} {r['gamma_freeze']:>10.2e} "
              f"{eq:>10.2e} {ratio:>9.1f} {r['p_undecided']:>7.3f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=300)
    ap.add_argument("--w", type=float, default=1.0 / 3.0)
    ap.add_argument("--gamma0-cooling", type=float, default=0.55)
    ap.add_argument("--gamma-fixed", type=float, default=0.05)
    ap.add_argument("--hubbles", type=float, nargs="+",
                    default=[0.01, 0.02, 0.05, 0.10, 0.20])
    ap.add_argument("--trials", type=int, default=400)
    ap.add_argument("--seed", type=int, default=20260728)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/cooling_relic.json"))
    args = ap.parse_args()

    t0 = time.time()
    sch = CoolingSchedule.build(args.gamma0_cooling, args.w)
    print(f"Omega={args.omega}  w={args.w:.4f}  cooling arm gamma0="
          f"{args.gamma0_cooling} (above gamma_c=0.5: no landscape at s=0)")
    print(f"  the pitchfork gamma_c=0.5 is crossed at s={sch.s_of_gamma(0.5):.4f} "
          f"of the available internal time, for every H")

    cool = [run_cell(args.gamma0_cooling, args.w, args.omega, h,
                     args.trials, args.seed + i)
            for i, h in enumerate(args.hubbles)]
    cool = [c for c in cool if c]
    report(cool, "COOLING arm: does the relic rise with H?  (P1, P2)")

    fixed = [run_cell(args.gamma_fixed, 0.0, args.omega, h,
                      args.trials, args.seed + 100 + i)
             for i, h in enumerate(args.hubbles)]
    fixed = [f for f in fixed if f]
    report(fixed, "FIXED arm (w=0, FINDINGS 5.1's model): H-independent?  (P3)")

    # Summarise ONLY the cells that actually decided. Comparing the arms over the
    # full H range makes them look identical (1162x vs 1318x) because at large H
    # both are dominated by marginal decisions -- conditioning on "decided" does
    # not rescue a cell that was 62% undecided. That summary line supported the
    # opposite of the correct conclusion and is the reason for this filter.
    UNDEC_MAX = 0.10
    for rows, name in ((cool, "cooling"), (fixed, "fixed  ")):
        clean = [r for r in rows if r["p_undecided"] <= UNDEC_MAX]
        if len(clean) < 2:
            print(f"\n{name} arm: fewer than 2 clean cells "
                  f"(undecided <= {UNDEC_MAX:.0%}); nothing quotable")
            continue
        hs = [c["hubble"] for c in clean]
        rs = [c["relic_mean"] for c in clean]
        print(f"\n{name} arm, clean cells only (undecided <= {UNDEC_MAX:.0%}): "
              f"H {min(hs):g}-{max(hs):g} ({max(hs)/min(hs):.0f}x) -> relic "
              f"{rs[0]:.5f} to {rs[-1]:.5f}  ({rs[-1]/max(rs[0], 1e-12):.0f}x)")
        eqs = [c["relic_mean"] / c["equilibrium_at_freeze"] for c in clean
               if c["equilibrium_at_freeze"] > 0]
        if eqs:
            print(f"{' ' * len(name)}      relic/equilibrium over those cells: "
                  f"{min(eqs):.3g} to {max(eqs):.3g}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"cooling": cool, "fixed": fixed}, indent=2))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
