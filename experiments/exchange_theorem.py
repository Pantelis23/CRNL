"""T15-d: exchange symmetry gives DIVISIBILITY; §30's constant ratio is a degree accident

T15-c (§42, run first) found that §30's identity survives two and four conservation laws
and dies the moment exchange symmetry is broken by a single rate constant. Reading that
back exposes a distinction §30 never had to make, because in the AM family the two
properties coincide.

**AN ANTISYMMETRIC POLYNOMIAL IN TWO VARIABLES IS DIVISIBLE BY THEIR DIFFERENCE.** If a
mass-action network is invariant under swapping species i and j, then

    b_i(sigma n) = b_j(n)     =>     (b_i - b_j) is antisymmetric under n_i <-> n_j

and propensities are polynomials in the counts (falling factorials are polynomials), so

    b_i - b_j = (n_i - n_j) * P(n),      P symmetric under i <-> j.                   (*)

That is a THEOREM, not a measurement, and it holds at ANY reaction order and for ANY number
of conservation laws. Its consequence is the whole restoration claim: `delta = 0` is an
invariant manifold of the deterministic flow, so **the sign of a lead is a deterministic
invariant and every reversal is a fluctuation**.

§30's identity claims something STRICTLY STRONGER: that P does not depend on the split at
fixed `s = n_i + n_j`. P is symmetric, so P = P(s, delta^2); independence of delta needs P
to have degree 0 in delta^2, i.e. the pair to enter every reaction at total degree <= 2.
AM satisfies that -- its highest pair term is `2X -> B + X`, giving P ~ (s - 1) -- so §30
could not tell the two apart. **A cubic pair term separates them.**

Take the symmetric pair `X + 2Y -> 3B` and `2X + Y -> 3B`. In the difference,

    b_X - b_Y  ~  a(X+2Y) - a(2X+Y)  ~  n_X n_Y (n_Y - 1) - n_X (n_X - 1) n_Y
               =  -n_X n_Y * delta

so P ~ n_X n_Y = (s^2 - delta^2)/4: **still divisible, no longer constant in delta.**

PREDICTIONS, written before running:

  P1  GATE. `am_reversible` passes BOTH tests; `am_asymmetric` fails BOTH. If a network
      that violates §30's ratio still passed divisibility here for the wrong reason, the
      probe would have no power and nothing below is admissible.
  P2  THE SEPARATION. The cubic network **passes divisibility and FAILS the constant
      ratio**. The two properties are distinct, §30 measured the stronger one, and only
      the weaker one is what the restoration claim actually needs.
  P3  The cubic network's ratio is EXACTLY affine in delta^2 -- residual at machine
      precision, not merely a good fit -- because P is symmetric of degree 2. A fit that
      is merely good would mean the algebra above is wrong somewhere.
  P4  THE THEOREM. On randomly generated networks that are symmetrised by construction --
      random species counts, random orders up to 4, random rate constants -- divisibility
      holds in **every** case, to machine precision. One counterexample refutes (*).
  P5  On the same random networks with symmetrisation REMOVED, divisibility generically
      FAILS. Without this the P4 result would be untestable: a probe that passes
      everything measures nothing.
  P6  Divisibility is indifferent to the number of conservation laws, as T15-c found for
      the stronger property. Reported per law-count so the claim is not read off a
      single value.
"""
from __future__ import annotations

import argparse
import json
import pathlib
from collections import Counter

import numpy as np

from crnl.networks.am_asymmetric import am_asymmetric
from crnl.networks.am_reversible import am_reversible
from crnl.reactions import Reaction, ReactionNetwork


def am_cubic(gamma: float, k: float = 1.0) -> ReactionNetwork:
    """AM whose disagreement channel is CUBIC in the pair, symmetrised.

    X + 2Y -> 3B and 2X + Y -> 3B with equal constants: exchange-symmetric, but the pair
    enters at total degree 3, so P acquires a delta^2.
    """
    rx = [
        Reaction({"X": 1, "Y": 2}, {"B": 3}, k, name="dis-xyy"),
        Reaction({"B": 3}, {"X": 1, "Y": 2}, gamma * k, name="rev-dis-xyy"),
        Reaction({"X": 2, "Y": 1}, {"B": 3}, k, name="dis-xxy"),
        Reaction({"B": 3}, {"X": 2, "Y": 1}, gamma * k, name="rev-dis-xxy"),
        Reaction({"B": 1, "X": 1}, {"X": 2}, k, name="rec-X"),
        Reaction({"X": 2}, {"B": 1, "X": 1}, gamma * k, name="rev-rec-X"),
        Reaction({"B": 1, "Y": 1}, {"Y": 2}, k, name="rec-Y"),
        Reaction({"Y": 2}, {"B": 1, "Y": 1}, gamma * k, name="rev-rec-Y"),
    ]
    return ReactionNetwork(species=["X", "Y", "B"], reactions=rx, name=f"am-cubic-g{gamma}")


def swap_counts(d: dict, a: str, b: str) -> dict:
    return {(b if s == a else a if s == b else s): v for s, v in d.items()}


def random_network(rng, n_extra=2, n_rx=5, max_order=4, symmetrise=True):
    """Random mass-action network on X, Y and n_extra spectators.

    With `symmetrise`, every reaction is accompanied by its X<->Y image at the SAME rate
    constant, which is exactly the hypothesis of (*). Reactions that are their own image
    are added once.
    """
    species = ["X", "Y"] + [f"S{i}" for i in range(n_extra)]
    rx, seen = [], set()
    for r in range(n_rx):
        k = float(rng.uniform(0.3, 2.0))
        left, right = {}, {}
        for side in (left, right):
            for _ in range(int(rng.integers(1, max_order + 1))):
                s = species[int(rng.integers(0, len(species)))]
                side[s] = side.get(s, 0) + 1
        pairs = [(left, right)]
        if symmetrise:
            pairs.append((swap_counts(left, "X", "Y"), swap_counts(right, "X", "Y")))
        for li, ri in pairs:
            if li == ri:
                continue
            key = (tuple(sorted(li.items())), tuple(sorted(ri.items())))
            if key in seen:
                continue
            seen.add(key)
            rx.append(Reaction(dict(li), dict(ri), k, name=f"r{r}-{len(rx)}"))
    if not rx:
        return None
    return ReactionNetwork(species=species, reactions=rx, name="random")


def n_conservation_laws(net) -> int:
    S = net.stoichiometry_matrix().astype(float)
    return int(S.shape[0] - np.linalg.matrix_rank(S))


def probe(net, counts, omega, i=0, j=1):
    """Return (divisibility residual, ratio spread, deltas, ratios) at fixed n_i+n_j.

    Divisibility is checked where n_i == n_j: b_i - b_j must be EXACTLY zero there,
    normalised by the total absolute traffic through the pair so the test has meaning
    when the drift itself is tiny (§30's normalisation; a bare 0/0 read 'fail' there).
    """
    cs = net.stochastic_constants(float(omega))
    S = net.stoichiometry_matrix().astype(np.int64)
    total = int(counts[i] + counts[j])
    if total % 2:
        total -= 1
    scale = np.abs(S[i, :] - S[j, :]).astype(float)

    def at(ni, nj):
        c = np.array(counts, dtype=np.int64)
        c[i], c[j] = ni, nj
        a = net.propensities(c, cs)
        return float((S[i, :] - S[j, :]) @ a), float(scale @ np.abs(a))

    val, traffic = at(total // 2, total // 2)
    div = abs(val) / max(traffic, 1e-300)

    deltas, ratios = [], []
    for s in range(2, total - 1):
        if 2 * s == total:
            continue
        v, _ = at(s, total - s)
        d = float(2 * s - total)
        deltas.append(d)
        ratios.append(v / d)
    if len(ratios) < 5:
        return div, None, None, None
    r = np.array(ratios)
    spread = float(np.ptp(r) / max(abs(r.mean()), 1e-300))
    return div, spread, np.array(deltas), r


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.25)
    ap.add_argument("--beta", type=float, default=0.20)
    ap.add_argument("--omega", type=int, default=90)
    ap.add_argument("--states", type=int, default=25)
    ap.add_argument("--random-nets", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20260809)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/exchange_theorem.json"))
    args = ap.parse_args()

    g, om, be = args.gamma, args.omega, args.beta
    rng = np.random.default_rng(args.seed)

    named = [
        ("am_reversible", am_reversible(g), "symmetric, pair degree 2"),
        ("am_cubic", am_cubic(g), "symmetric, pair degree 3"),
        ("am_asymmetric", am_asymmetric(g, be), "asymmetric"),
    ]

    print(f"gamma={g} beta={be} Omega={om}")
    print(f"{'network':>16}{'laws':>6}{'divisibility':>15}{'div?':>7}"
          f"{'ratio spread':>15}{'§30 ratio?':>12}   note")
    rows, keep = [], {}
    for name, net, note in named:
        worst_div, worst_sp, sample = 0.0, 0.0, None
        for _ in range(args.states):
            cuts = np.sort(rng.integers(6, om - 6, size=2))
            c = [int(cuts[0]), int(cuts[1] - cuts[0]), int(om - cuts[1])]
            if min(c) < 6:
                continue
            c += [int(rng.integers(2, om)) for _ in range(len(net.species) - 3)]
            d, sp, dl, rt = probe(net, c, om)
            worst_div = max(worst_div, d)
            if sp is not None:
                worst_sp = max(worst_sp, sp)
                if sample is None:
                    sample = (dl, rt)
        dv, rr = worst_div < 1e-12, worst_sp < 1e-12
        keep[name] = sample
        rows.append({"network": name, "cons_laws": n_conservation_laws(net),
                     "divisibility": worst_div, "divisible": bool(dv),
                     "ratio_spread": worst_sp, "constant_ratio": bool(rr), "note": note})
        print(f"{name:>16}{n_conservation_laws(net):>6}{worst_div:>15.3e}"
              f"{'YES' if dv else 'NO':>7}{worst_sp:>15.3e}"
              f"{'YES' if rr else 'NO':>12}   {note}")

    by = {r["network"]: r for r in rows}
    gate = (by["am_reversible"]["divisible"] and by["am_reversible"]["constant_ratio"]
            and not by["am_asymmetric"]["divisible"]
            and not by["am_asymmetric"]["constant_ratio"])
    print(f"\n=== P1 gate: {'OK' if gate else 'FAILED -- nothing below is admissible'}")

    cu = by["am_cubic"]
    print("\n=== P2: do the two properties separate?")
    print(f"  am_cubic (symmetric, degree 3): divisible {'YES' if cu['divisible'] else 'NO'}"
          f", §30 ratio {'YES' if cu['constant_ratio'] else 'NO'}")
    if cu["divisible"] and not cu["constant_ratio"]:
        print("  -> SEPARATED. Divisibility is the weaker property, it is what the")
        print("     restoration claim needs, and §30's constant ratio is the special")
        print("     case where the pair enters at total degree <= 2.")

    print("\n=== P3: is the cubic ratio EXACTLY affine in delta^2?")
    dl, rt = keep["am_cubic"]
    A = np.vstack([np.ones_like(dl), dl ** 2]).T
    coef, *_ = np.linalg.lstsq(A, rt, rcond=None)
    resid = float(np.max(np.abs(A @ coef - rt)) / max(np.abs(rt).max(), 1e-300))
    print(f"  ratio = {coef[0]:.6e} + {coef[1]:.6e} * delta^2"
          f"   max rel residual {resid:.3e}"
          f"   -> {'EXACT' if resid < 1e-12 else 'NOT exact'}")

    print(f"\n=== P4/P5: the theorem on {args.random_nets} random networks")
    stats = {True: [], False: []}
    laws, laws_ok = Counter(), Counter()
    for sym in (True, False):
        for t in range(args.random_nets):
            net = random_network(rng, n_extra=int(rng.integers(1, 4)),
                                 n_rx=int(rng.integers(3, 8)),
                                 max_order=int(rng.integers(2, 5)), symmetrise=sym)
            if net is None:
                continue
            c = [int(rng.integers(6, om)) for _ in net.species]
            try:
                d, _, _, _ = probe(net, c, om)
            except Exception:
                continue
            stats[sym].append(d)
            if sym:
                nl = n_conservation_laws(net)
                laws[nl] += 1
                laws_ok[nl] += int(d < 1e-12)
        arr = np.array(stats[sym])
        n_div = int((arr < 1e-12).sum())
        tag = "SYMMETRISED" if sym else "not symmetrised"
        print(f"  {tag:>16}: {n_div}/{len(arr)} divisible"
              f"   worst residual {arr.max():.3e}   median {np.median(arr):.3e}")

    sym_arr, asym_arr = np.array(stats[True]), np.array(stats[False])
    p4 = bool((sym_arr < 1e-12).all())
    p5 = bool((asym_arr < 1e-12).mean() < 0.5)
    print(f"\n  P4 {'HOLDS -- no counterexample to (*)' if p4 else 'REFUTED'};"
          f"  P5 {'HOLDS -- the probe has power' if p5 else 'FAILS -- probe is vacuous'}")

    print("\n=== P6: divisibility vs number of conservation laws (symmetrised only)")
    for nl in sorted(laws):
        print(f"  {nl:>2} laws: {laws_ok[nl]:>4}/{laws[nl]:<4} divisible")

    out = {"named": rows,
           "cubic_fit": {"const": float(coef[0]), "delta2": float(coef[1]),
                         "max_rel_resid": resid},
           "random": {"symmetrised_worst": float(sym_arr.max()),
                      "symmetrised_n": len(sym_arr),
                      "symmetrised_divisible": int((sym_arr < 1e-12).sum()),
                      "plain_n": len(asym_arr),
                      "plain_divisible": int((asym_arr < 1e-12).sum())},
           "by_laws": {str(k): [laws_ok[k], laws[k]] for k in sorted(laws)}}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
