"""T15-h: the non-restoring rate constants form a CONVEX CONE — capability is combinatorial

§53 concluded "restoration is tuned, not topological". §54 gave the decomposition
P(x) = sum_r d_r c_r B_r(x) with every bracket B_r >= 0, and measured that mixed networks
amplify only 18% of the time at sampled rates. **Both of those readings are incomplete, and
the decomposition already contains the correction.**

**P IS LINEAR IN THE RATE CONSTANTS.** Write v_r(x) = d_r B_r(x), so P(x) = <c, v(x)>. Then

    NON-RESTORING  <=>  <c, v(x)> <= 0 for every accessible x

and that condition is closed under addition and positive scaling:

    <c1 + c2, v(x)> = <c1, v(x)> + <c2, v(x)> <= 0,     <t c, v(x)> = t <c, v(x)> <= 0

> **THE NON-RESTORING RATE VECTORS FORM A CONVEX CONE** -- the polar of the cone generated
> by {v(x) : x accessible}.

Two corollaries follow that §54 did not state, and the second contradicts the natural
reading of §54's own 18%:

  **(a)** all d_r <= 0 makes v(x) <= 0 componentwise, so the WHOLE positive orthant of rate
  constants lies in the cone. §54's `all<=0` class is a special case of the theorem.

  **(b)** some d_r > 0 means one can load c onto that reaction and choose x where its
  bracket dominates, so **the network restores for SOME c. CAPABILITY IS COMBINATORIAL**,
  even though realisation is not. §54's 18% would then be a fact about the RATES SAMPLED,
  not about the networks.

**THE SCOPE LIMIT IS THE INTERESTING PART.** The domination argument in (b) needs room to
scale x until one monomial beats the others. On the open positive orthant there is always
room. **On a COMPACT accessible set -- a simplex, which is what a conservation law imposes
-- there may not be**, and then capability can fail even with d_r > 0. That is the
difference between "restoration is possible in principle" and "possible at accessible
concentrations", and it is exactly where a conservation law bites.

PREDICTIONS, written before running:

  P1  GATE. P is exactly linear in c: P(a*c1 + b*c2) = a*P(c1) + b*P(c2) to machine
      precision, at fixed x. If not, the decomposition is not what §54 says it is.
  P2  **THE CONE THEOREM.** If c1 and c2 are both non-restoring then so is c1 + c2. Tested
      on mixed networks over many random pairs. **One counterexample refutes the theorem**,
      not merely the experiment.
  P3  **CAPABILITY (b), constructively.** For every network with some d_r > 0, loading c on
      the positive-d reactions and starving the negative-d ones restores. On the open
      orthant this must work in **every** case.
  P4  §54's `all<=0` networks are non-restoring for EVERY c, not merely for the ones §54
      sampled. Corollary (a), tested over wide random c.
  P5  **THE REFINEMENT OF MY OWN §54.** Its mixed class amplified 18% of the time. Under
      (b), **100% of them are capable**. If that holds, §54's 18% is reported for what it
      is -- a property of the rate distribution sampled -- and the network-level statement
      is that mixed means capable.
  P6  **THE SCOPE LIMIT.** On networks with a conservation law, where the accessible set is
      a compact simplex, capability may FAIL despite d_r > 0. Predicted: a strictly smaller
      fraction than on the orthant. If instead it is also 100%, compactness does not bite
      and the theorem is stronger than the proof sketch requires -- reported either way.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from crnl.networks.am_reversible import am_reversible
from experiments.amplification_signature import classify, mirror_pairs
from experiments.exchange_theorem import random_network


def v_of(net, x, i="X", j="Y"):
    """v_r(x) = d_r * B_r(x) over mirror pairs, and the matching rate vector c."""
    specs = list(net.species)
    xi, xj = x[specs.index(i)], x[specs.index(j)]
    pairs = mirror_pairs(net, i, j)
    if pairs is None:
        return None, None
    v, c = [], []
    for d, idx in pairs:
        r = net.reactions[idx]
        p, q = r.reactants.get(i, 0), r.reactants.get(j, 0)
        O = 1.0
        for sp, m in r.reactants.items():
            if sp not in (i, j):
                O *= x[specs.index(sp)] ** m
        ssum = sum(xi ** a * xj ** (p - q - 1 - a) for a in range(p - q))
        v.append(d * O * (xi * xj) ** q * ssum)
        c.append(r.k)
    return np.array(v, dtype=float), np.array(c, dtype=float)


def restores(net, c, states, i="X", j="Y", tol=1e-12):
    """Is <c, v(x)> > 0 for any sampled accessible x?"""
    for x in states:
        v, _ = v_of(net, x, i, j)
        if v is None:
            return None
        if float(np.dot(c, v)) > tol:
            return True
    return False


def orthant_states(rng, m, n=60):
    """Wide log-spread samples of the open positive orthant."""
    return [np.exp(rng.uniform(-6, 6, m)) for _ in range(n)]


def simplex_states(rng, m, n=60):
    """Samples of the compact simplex sum(x) = 1 -- what a conservation law imposes."""
    out = []
    for _ in range(n):
        z = rng.dirichlet(np.ones(m))
        out.append(z)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--nets", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20260812)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/restoration_cone.json"))
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    print("=== P1 GATE: is P exactly linear in the rate constants?")
    worst = 0.0
    for _ in range(60):
        net = random_network(rng, n_extra=2, n_rx=5, max_order=3, symmetrise=True)
        if net is None:
            continue
        x = np.exp(rng.uniform(-2, 2, len(net.species)))
        v, _ = v_of(net, x)
        if v is None or v.size == 0:
            continue
        c1, c2 = rng.uniform(0.1, 5, v.size), rng.uniform(0.1, 5, v.size)
        a, b = 0.37, 2.1
        lhs = float(np.dot(a * c1 + b * c2, v))
        rhs = a * float(np.dot(c1, v)) + b * float(np.dot(c2, v))
        worst = max(worst, abs(lhs - rhs) / max(abs(lhs), 1e-30))
    print(f"  worst relative deviation: {worst:.3e}   -> P1 {'HOLDS' if worst < 1e-12 else 'FAILS'}")

    print("\n=== P2 THE CONE THEOREM: c1, c2 non-restoring => c1 + c2 non-restoring")
    tested = viol = 0
    for _ in range(args.nets):
        net = random_network(rng, n_extra=int(rng.integers(1, 3)),
                             n_rx=int(rng.integers(3, 6)), max_order=3, symmetrise=True)
        if net is None or classify(net) != "mixed":
            continue
        states = orthant_states(rng, len(net.species), 40)
        v0, _ = v_of(net, states[0])
        if v0 is None or v0.size == 0:
            continue
        bad = []
        for _ in range(40):
            c = np.exp(rng.uniform(-4, 4, v0.size))
            if restores(net, c, states) is False:
                bad.append(c)
            if len(bad) >= 2:
                break
        if len(bad) < 2:
            continue
        tested += 1
        if restores(net, bad[0] + bad[1], states):
            viol += 1
    print(f"  {tested} networks with two non-restoring rate vectors found;"
          f" {viol} violations")
    print(f"  -> P2 {'HOLDS -- the non-restoring set is closed under addition' if viol == 0 else 'REFUTED, the theorem is wrong'}")

    print("\n=== P3/P5: is CAPABILITY combinatorial? (load c on the positive-d reactions)")
    stats = {"mixed": [0, 0], "all<=0": [0, 0], "all>=0": [0, 0], "trivial": [0, 0]}
    for _ in range(args.nets):
        net = random_network(rng, n_extra=int(rng.integers(1, 4)),
                             n_rx=int(rng.integers(3, 8)), max_order=4, symmetrise=True)
        if net is None:
            continue
        cls = classify(net)
        if cls is None:
            continue
        states = orthant_states(rng, len(net.species), 80)
        v0, _ = v_of(net, states[0])
        if v0 is None or v0.size == 0:
            continue
        pairs = mirror_pairs(net)
        ds = np.array([d for d, _ in pairs], dtype=float)
        c = np.where(ds > 0, 1e4, np.where(ds < 0, 1e-4, 1.0))
        ok = restores(net, c, states)
        stats[cls][0] += 1
        stats[cls][1] += int(bool(ok))
    print(f"{'class':>10}{'networks':>10}{'capable':>10}{'fraction':>11}"
          f"   (§54 measured at SAMPLED rates)")
    for k in ("all<=0", "all>=0", "mixed", "trivial"):
        n, a = stats[k]
        if n:
            print(f"{k:>10}{n:>10}{a:>10}{100*a/n:>10.1f}%"
                  + ("   <- §54 found 18% amplifying here" if k == "mixed" else "")
                  + ("   <- §54 found 0/113" if k == "all<=0" else ""))
    m_n, m_a = stats["mixed"]
    print(f"  -> P3/P5 {'HOLD: mixed means CAPABLE' if m_n and m_a == m_n else 'FAIL'}")
    z_n, z_a = stats["all<=0"]
    print(f"  -> P4 {'HOLDS: all<=0 is non-restoring for every c' if z_a == 0 else 'REFUTED'}")

    print("\n=== P6 SCOPE LIMIT: does capability survive on a COMPACT simplex?")
    sstats = {"mixed": [0, 0]}
    for _ in range(args.nets):
        net = random_network(rng, n_extra=int(rng.integers(1, 4)),
                             n_rx=int(rng.integers(3, 8)), max_order=4, symmetrise=True)
        if net is None or classify(net) != "mixed":
            continue
        states = simplex_states(rng, len(net.species), 80)
        v0, _ = v_of(net, states[0])
        if v0 is None or v0.size == 0:
            continue
        ds = np.array([d for d, _ in mirror_pairs(net)], dtype=float)
        c = np.where(ds > 0, 1e4, np.where(ds < 0, 1e-4, 1.0))
        sstats["mixed"][0] += 1
        sstats["mixed"][1] += int(bool(restores(net, c, states)))
    n, a = sstats["mixed"]
    if n:
        print(f"  mixed networks on the simplex: {a}/{n} = {100*a/n:.1f}% capable"
              f"   (orthant: {100*m_a/max(m_n,1):.1f}%)")
        print(f"  -> {'compactness does NOT bite -- capability survives' if a == n else 'compactness BITES: a conservation law can forbid restoration a network is otherwise capable of'}")

    print("\n=== AM as a worked case")
    for gamma in (0.2, 0.5, 0.8):
        net = am_reversible(gamma)
        ds = [d for d, _ in mirror_pairs(net)]
        st = simplex_states(rng, 3, 200)
        r = restores(net, np.array([net.reactions[i].k for _, i in mirror_pairs(net)]), st)
        print(f"  gamma={gamma}: d = {ds}, class = {classify(net)}, restores = {r}")
    print("  AM is mixed, hence CAPABLE; gamma decides whether it realises it.")

    out = {"p1_worst": worst, "cone_tested": tested, "cone_violations": viol,
           "orthant": stats, "simplex": sstats}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
