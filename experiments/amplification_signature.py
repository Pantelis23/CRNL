"""T15-f: sign(P) is a c-weighted sum with NON-NEGATIVE weights — so topology decides the
unanimous cases and rate constants decide only the mixed ones

§53 closed T15-e negatively: sign(P) is not determined by the stoichiometry (17 of 188
random topologies flip it under rate changes alone), and amplification is rare -- 21/200
against divisibility's 200/200. That leaves the obvious question unasked: **what do the
10.5% have in common?**

**P DECOMPOSES EXACTLY, AND THE DECOMPOSITION SEPARATES TOPOLOGY FROM RATES.** Group the
reactions of an exchange-symmetric network into mirror pairs {r, rbar}. Let reaction r have
X-power p and Y-power q in its reactants with p > q, all other species entering through a
symmetric factor O_r. Then a_r = c_r O_r x^p y^q and a_rbar = c_r O_r x^q y^p, while
symmetry forces S_X(rbar) = S_Y(r), so with

    **d_r = S_X(r) - S_Y(r)**

the pair contributes d_r c_r O_r (x^p y^q - x^q y^p) = d_r c_r O_r (xy)^q (x - y) * SUM,
where SUM = sum_{m} x^m y^{p-q-1-m} >= 0. A self-mirror reaction (p = q) forces S_X = S_Y
and contributes nothing. Hence

    **P = SUM_pairs d_r * c_r * [ O_r (xy)^q SUM_r ]**,  every bracket NON-NEGATIVE on x,y >= 0

**So sign(P) is a rate-weighted sum of the integers d_r with non-negative weights.** Three
regimes follow immediately, and two of them are purely combinatorial:

  * **all d_r <= 0  =>  P <= 0 everywhere.** Restoration impossible, from stoichiometry
    alone, whatever the rate constants.
  * **all d_r >= 0, some > 0  =>  P > 0 everywhere** (on x, y > 0). Restoration guaranteed
    from stoichiometry alone.
  * **mixed signs  =>  the rate constants decide**, and only here.

**d_r > 0 IS NOT THE SAME AS AUTOCATALYSIS**, which is the obvious guess and is wrong.
d_r = S_X - S_Y > 0 holds for AM's `B + X -> 2X` (X makes more X), but equally for
`2X + Y -> 2X + B`, where S_X = 0 and S_Y = -1: **X catalysing Y's DESTRUCTION amplifies a
lead exactly as well as X catalysing its own production.** The right notion is positive
feedback on the *difference*, not on either species.

PREDICTIONS, written before running:

  P1  GATE, exact. The decomposition reproduces `P_at` to machine precision on AM, on
      `am_cubic`, and on random symmetrised networks. If it does not, the algebra above is
      wrong and nothing else here is admissible.
  P2  **all d_r <= 0 => P <= 0 at every sampled state**, over random symmetrised networks
      and over many random rate-constant draws each. One counterexample refutes.
  P3  **all d_r >= 0 with some d_r > 0 => P > 0 at every sampled state with x, y > 0.**
      Same standard: one counterexample refutes.
  P4  **THE CROSS-CHECK AGAINST §53, and it is the sharp one.** §53 found 17/188 topologies
      flipping sign(P) under rate changes. **Those must be exactly the MIXED ones**, and no
      unanimous topology may flip. This links a combinatorial classification computed here
      to a measurement made in a different experiment for a different purpose.
  P5  The 10.5% amplifying fraction of §53 should be roughly the fraction of topologies that
      are all-non-negative or mixed-and-favourably-weighted -- reported as a decomposition
      of that number, not as a new one.
  P6  If P2 or P3 fails, the non-negativity of the bracket is wrong somewhere -- most likely
      where a species appears on BOTH sides so that O_r is not what I claim, or where the
      CME's falling factorials differ from the concentration powers used here. That
      distinction is checked explicitly rather than assumed.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from crnl.networks.am_reversible import am_reversible
from experiments.amplification_sign import P_at
from experiments.exchange_theorem import am_cubic, random_network, swap_counts


def mirror_pairs(net, i="X", j="Y"):
    """[(d_r, r_index)] over mirror pairs, taking the X-heavy member of each."""
    specs = list(net.species)
    si, sj = specs.index(i), specs.index(j)
    S = net.stoichiometry_matrix()
    key = lambda r: (tuple(sorted(r.reactants.items())), tuple(sorted(r.products.items())))
    byk = {key(r): idx for idx, r in enumerate(net.reactions)}
    seen, out = set(), []
    for idx, r in enumerate(net.reactions):
        if idx in seen:
            continue
        mk = (tuple(sorted(swap_counts(r.reactants, i, j).items())),
              tuple(sorted(swap_counts(r.products, i, j).items())))
        jdx = byk.get(mk)
        if jdx is None:
            return None                      # not symmetric; caller must skip
        if jdx == idx:
            seen.add(idx)
            continue                          # self-mirror: contributes nothing
        seen.update({idx, jdx})
        p, q = r.reactants.get(i, 0), r.reactants.get(j, 0)
        heavy = idx if p > q else jdx
        out.append((int(S[si, heavy] - S[sj, heavy]), heavy))
    return out


def P_decomposed(net, x, i="X", j="Y"):
    """P rebuilt as sum_pairs d_r * c_r * (non-negative bracket)."""
    specs = list(net.species)
    xi, xj = x[specs.index(i)], x[specs.index(j)]
    pairs = mirror_pairs(net, i, j)
    if pairs is None:
        return None, None
    total, brackets = 0.0, []
    for d, idx in pairs:
        r = net.reactions[idx]
        p, q = r.reactants.get(i, 0), r.reactants.get(j, 0)
        O = 1.0
        for sp, n in r.reactants.items():
            if sp not in (i, j):
                O *= x[specs.index(sp)] ** n
        ssum = sum(xi ** m * xj ** (p - q - 1 - m) for m in range(p - q))
        bracket = O * (xi * xj) ** q * ssum
        brackets.append(bracket)
        total += d * r.k * bracket
    return float(total), brackets


def classify(net, i="X", j="Y"):
    """Classify by the d_r of pairs that actually CONTRIBUTE.

    A p == q pair has an empty bracket sum and contributes identically zero, but
    `mirror_pairs` records its d as an arbitrary +-1 because "the X-heavy member" is
    undefined there (§55). Counting those misclassifies genuinely all-d<=0 networks as
    "mixed" -- 148 of 232 in the draw that caught it.
    """
    pairs = mirror_pairs(net, i, j)
    if pairs is None:
        return None
    ds = [d for d, idx in pairs
          if net.reactions[idx].reactants.get(i, 0)
          != net.reactions[idx].reactants.get(j, 0)]
    if not ds or all(d == 0 for d in ds):
        return "trivial"
    if all(d <= 0 for d in ds):
        return "all<=0"
    if all(d >= 0 for d in ds):
        return "all>=0"
    return "mixed"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--nets", type=int, default=300)
    ap.add_argument("--seed", type=int, default=20260811)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/amplification_signature.json"))
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    print("=== P1 GATE: does the decomposition reproduce P exactly?")
    worst = 0.0
    for gamma in (0.05, 0.20, 0.35, 0.55):
        for net in (am_reversible(gamma), am_cubic(gamma)):
            for s in (0.3, 0.6, 0.9):
                for d in (0.05, 0.2):
                    x = [(s + d) / 2, (s - d) / 2, 1 - s]
                    a = P_at(net, x)
                    b, _ = P_decomposed(net, x)
                    worst = max(worst, abs(a - b) / max(abs(a), 1e-30))
    nrand = 0
    for _ in range(120):
        net = random_network(rng, n_extra=int(rng.integers(1, 4)),
                             n_rx=int(rng.integers(3, 7)), max_order=4, symmetrise=True)
        if net is None:
            continue
        x = rng.uniform(0.1, 1.5, len(net.species))
        a = P_at(net, x)
        b, _ = P_decomposed(net, x)
        if b is None:
            continue
        nrand += 1
        worst = max(worst, abs(a - b) / max(abs(a), 1e-12))
    print(f"  worst relative error over AM/am_cubic + {nrand} random networks:"
          f" {worst:.3e}   -> P1 {'HOLDS' if worst < 1e-9 else 'FAILS'}")

    print("\n=== P2/P3: do the unanimous classes behave as the algebra demands?")
    tally = {"all<=0": [0, 0], "all>=0": [0, 0], "mixed": [0, 0], "trivial": [0, 0]}
    viol = {"all<=0": 0, "all>=0": 0}
    for _ in range(args.nets):
        net = random_network(rng, n_extra=int(rng.integers(1, 4)),
                             n_rx=int(rng.integers(3, 8)), max_order=4, symmetrise=True)
        if net is None:
            continue
        cls = classify(net)
        if cls is None:
            continue
        pos = neg = False
        for _ in range(15):
            x = rng.uniform(0.05, 1.5, len(net.species))
            p = P_at(net, x)
            if p > 1e-9:
                pos = True
            if p < -1e-9:
                neg = True
        tally[cls][0] += 1
        tally[cls][1] += int(pos)
        if cls == "all<=0" and pos:
            viol["all<=0"] += 1
        if cls == "all>=0" and neg:
            viol["all>=0"] += 1
    print(f"{'class':>10}{'networks':>10}{'P>0 somewhere':>16}{'violations':>12}")
    for k in ("all<=0", "all>=0", "mixed", "trivial"):
        n, a = tally[k]
        v = viol.get(k, 0)
        print(f"{k:>10}{n:>10}{a:>16}{v:>12}")
    print(f"  -> P2 {'HOLDS: all<=0 never amplifies' if viol['all<=0'] == 0 else 'REFUTED'}")
    print(f"  -> P3 {'HOLDS: all>=0 never contracts' if viol['all>=0'] == 0 else 'REFUTED'}")

    print("\n=== P4 CROSS-CHECK: do only MIXED topologies flip under rate changes?")
    flips = {"all<=0": 0, "all>=0": 0, "mixed": 0, "trivial": 0}
    counts = {k: 0 for k in flips}
    for _ in range(args.nets):
        base = random_network(rng, n_extra=int(rng.integers(1, 3)),
                              n_rx=int(rng.integers(3, 6)), max_order=3, symmetrise=True)
        if base is None:
            continue
        cls = classify(base)
        if cls is None:
            continue
        x = rng.uniform(0.2, 1.0, len(base.species))
        signs = set()
        for _ in range(12):
            scale, rx = {}, []
            for r in base.reactions:
                k = r.name.split("-")[0]
                scale.setdefault(k, float(rng.uniform(0.05, 20.0)))
                rx.append(type(r)(dict(r.reactants), dict(r.products),
                                  r.k * scale[k], name=r.name))
            p = P_at(type(base)(species=list(base.species), reactions=rx, name="p"), x)
            if abs(p) > 1e-9:
                signs.add(p > 0)
        counts[cls] += 1
        if len(signs) == 2:
            flips[cls] += 1
    print(f"{'class':>10}{'networks':>10}{'flipped':>10}")
    for k in ("all<=0", "all>=0", "mixed", "trivial"):
        print(f"{k:>10}{counts[k]:>10}{flips[k]:>10}")
    unan = flips["all<=0"] + flips["all>=0"] + flips["trivial"]
    print(f"  unanimous topologies that flipped: {unan}"
          f"  -> P4 {'HOLDS -- only mixed topologies flip' if unan == 0 else 'REFUTED'}")
    if counts["mixed"]:
        print(f"  mixed topologies that flipped: {flips['mixed']}/{counts['mixed']}"
              f" = {100*flips['mixed']/counts['mixed']:.0f}%")

    print("\n=== P5: decomposing §53's 10.5%")
    tot = sum(tally[k][0] for k in tally)
    amp = sum(tally[k][1] for k in tally)
    print(f"  amplifying overall: {amp}/{tot} = {100*amp/max(tot,1):.1f}%"
          f"   (§53 measured 10.5% on a different draw)")
    for k in ("all>=0", "mixed"):
        n, a = tally[k]
        if n:
            print(f"    {k:>7}: {a}/{n} = {100*a/n:.0f}% amplify, contributing"
                  f" {100*a/max(tot,1):.1f} points")

    out = {"p1_worst": worst, "tally": tally, "violations": viol,
           "flips": flips, "counts": counts}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
