"""T15-e: is sign(P) structural, or is only the cannot-reverse half a theorem?

§43 proved that for ANY exchange-symmetric mass-action network

    b_i - b_j = (n_i - n_j) * P(n),   P symmetric

so delta = 0 is an invariant manifold and **the sign of a lead is a deterministic
invariant** -- 200/200 random symmetrised networks, 157 of them conserving nothing at all.
§43 was explicit that this is only half of restoration: divisibility gives NO-REVERSAL,
while AMPLIFICATION additionally needs P > 0, and §43's own tables show P going negative
(the bracket ranges start at -0.16). **"Restoration is divisibility plus sign(P) > 0, and
only the first half is a theorem."** This asks whether the second half can be made one.

**P HAS AN IDENTITY THAT §43 DID NOT NAME.** Since b_i - b_j = delta * P, differentiating
at delta = 0 gives

    **P(symmetric state) = d(b_i - b_j)/d(delta)**,

which at a symmetric FIXED POINT is exactly the symmetry-breaking eigenvalue -- the
quantity THEORIES T7 and §14 built the n-winner barrier on, where lambda(n) = 1/(2n-1) at
gamma = 0. So §43's P and T7's lambda are the same object evaluated at the same place, and
the question "does this network amplify" is the question "is the symmetric state unstable".

For AM, P = k(1 - (1+gamma)s) explicitly (§51), so P = 0 at s = 1/(1+gamma), i.e. at
b = gamma/(1+gamma) -- **the off-symmetric fixed points lie on P's zero set**, since
d(delta)/dt = delta*P vanishes there with delta != 0. P's zero set is where the attractors
are, not where the separatrix is: §43's theorem already makes the separatrix delta = 0.

**THE SHARP QUESTION IS WHETHER sign(P) IS COMBINATORIAL.** Divisibility is: it follows
from the stoichiometry and the exchange symmetry alone, with no reference to rate constants.
If sign(P) were too, restoration would be structural end to end.

PREDICTIONS, written before running:

  P1  GATE, exact. At any state with n_i = n_j, `P` equals the antisymmetric directional
      derivative of (b_i - b_j), to machine precision, on AM and on random symmetrised
      networks. And at AM's symmetric fixed point it equals the Jacobian's eigenvalue along
      (1,-1,0), the symmetry-breaking eigenvalue. If these disagree the identity is wrong.
  P2  For AM, P = 0 exactly at the off-symmetric fixed points -- checked against
      `delta_star(gamma)` on the slaved manifold, which is an independent route to the same
      points.
  P3  **THE TEST, and it is a kill test in the strict sense.** Hold the TOPOLOGY of a random
      symmetric network fixed and vary ONLY the rate constants. **If sign(P) flips, sign(P)
      is not combinatorial** and §43's theorem covers the cannot-reverse half only, for
      good. AM itself does this -- P at the symmetric point is positive for gamma < 1/2 and
      negative above it, which is gamma_c -- so the prediction is that it flips, and the
      point of running it on random networks is to show AM is not special.
  P4  THE PARTIAL STRUCTURAL STATEMENT, if any. Over random symmetrised networks, what
      FRACTION amplify at all (P > 0 somewhere)? A number near 0 or 1 would suggest the
      topology at least biases the answer; a number near 1/2 would say it carries no
      information. Reported either way.
  P5  Does P keep one sign between the symmetric state and the attractor, or can it flip in
      between -- a network that amplifies a small lead, contracts a medium one, then
      amplifies again? Nothing forbids it, and if it happens the phrase "restoring" needs a
      region attached to it.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from crnl.deterministic import jacobian
from crnl.networks.am_reversible import GAMMA_C, am_reversible, delta_star
from experiments.exchange_theorem import random_network
from experiments.slaving_axis import slaved


def P_at(net, x, i=0, j=1, h=1e-6):
    """P = (b_i - b_j)/(x_i - x_j), or its limit when x_i == x_j."""
    S = net.stoichiometry_matrix()
    if abs(x[i] - x[j]) > 1e-12:
        b = S @ net.fluxes(np.asarray(x, dtype=float))
        return float((b[i] - b[j]) / (x[i] - x[j]))
    xp, xm = np.array(x, float), np.array(x, float)
    xp[i] += h; xp[j] -= h
    xm[i] -= h; xm[j] += h
    bp = S @ net.fluxes(xp)
    bm = S @ net.fluxes(xm)
    return float(((bp[i] - bp[j]) - (bm[i] - bm[j])) / (4 * h))


def antisym_eig(net, x, i=0, j=1):
    """Jacobian eigenvalue along the antisymmetric direction e_i - e_j."""
    J = jacobian(net, np.asarray(x, dtype=float))
    v = np.zeros(len(x)); v[i], v[j] = 1.0, -1.0
    Jv = J @ v
    return float(Jv[i] - Jv[j]) / 2.0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--nets", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20260810)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/amplification_sign.json"))
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    print("=== P1 GATE: P at a symmetric state == the antisymmetric derivative")
    worst = 0.0
    for gamma in (0.05, 0.20, 0.35, 0.49):
        net = am_reversible(gamma)
        for s in (0.4, 0.6, 0.8):
            x = [s / 2, s / 2, 1 - s]
            a, b = P_at(net, x), antisym_eig(net, x)
            worst = max(worst, abs(a - b) / max(abs(a), 1e-30))
    for _ in range(40):
        net = random_network(rng, n_extra=2, n_rx=5, max_order=3, symmetrise=True)
        x = rng.uniform(0.1, 1.0, len(net.species))
        x[1] = x[0]
        a, b = P_at(net, x), antisym_eig(net, x)
        if abs(a) > 1e-9:
            worst = max(worst, abs(a - b) / abs(a))
    print(f"  worst relative disagreement over 52 states: {worst:.3e}"
          f"   -> P1 {'HOLDS' if worst < 1e-6 else 'FAILS'}")

    print("\n=== P1b: at AM's symmetric fixed point, P == the symmetry-breaking eigenvalue")
    print(f"{'gamma':>7}{'P':>14}{'lambda':>14}{'sign':>8}")
    for gamma in (0.05, 0.20, 0.35, 0.45, 0.49, 0.55, 0.70):
        net = am_reversible(gamma)
        st = slaved(net, 0.0)
        if st is None:
            continue
        p, lam = P_at(net, st), antisym_eig(net, st)
        print(f"{gamma:>7.2f}{p:>14.8f}{lam:>14.8f}{'+' if p > 0 else '-':>8}")

    print(f"\n  gamma_c = {GAMMA_C}: P changes sign there, which is the pitchfork.")

    print("\n=== P2: are AM's attractors on P's zero set?")
    print(f"{'gamma':>7}{'delta*':>11}{'P at attractor':>17}")
    for gamma in (0.05, 0.20, 0.35, 0.45):
        net = am_reversible(gamma)
        ds = delta_star(gamma)
        st = slaved(net, ds)
        if st is None:
            continue
        print(f"{gamma:>7.2f}{ds:>11.6f}{P_at(net, st):>17.3e}")

    print("\n=== P3 KILL TEST: fix the TOPOLOGY, vary only the rate constants")
    flips, tested = 0, 0
    examples = []
    for t in range(args.nets):
        base = random_network(rng, n_extra=int(rng.integers(1, 3)),
                              n_rx=int(rng.integers(3, 6)),
                              max_order=3, symmetrise=True)
        if base is None:
            continue
        x = rng.uniform(0.2, 1.0, len(base.species))
        x[1] = x[0]
        signs = set()
        vals = []
        for _ in range(12):
            # same reactions, same stoichiometry, NEW rate constants
            rx = []
            scale = {}
            for r in base.reactions:
                key = r.name.split("-")[0]
                if key not in scale:
                    scale[key] = float(rng.uniform(0.05, 20.0))
                rx.append(type(r)(dict(r.reactants), dict(r.products),
                                  r.k * scale[key], name=r.name))
            net = type(base)(species=list(base.species), reactions=rx, name="perturbed")
            p = P_at(net, x)
            if abs(p) > 1e-9:
                signs.add(p > 0)
                vals.append(p)
        if len(vals) >= 6:
            tested += 1
            if len(signs) == 2:
                flips += 1
                if len(examples) < 3:
                    examples.append((min(vals), max(vals)))
    print(f"  {flips}/{tested} topologies show sign(P) FLIPPING under rate changes alone")
    for lo, hi in examples:
        print(f"    example: P ranges {lo:+.4f} .. {hi:+.4f} at one fixed state")
    p3 = flips > 0.5 * tested
    print(f"  -> P3 {'sign(P) is NOT combinatorial' if p3 else 'sign(P) looks topology-determined'}")

    print("\n=== P4: what fraction of random symmetric networks amplify at all?")
    amp, tot = 0, 0
    for _ in range(args.nets):
        net = random_network(rng, n_extra=int(rng.integers(1, 4)),
                             n_rx=int(rng.integers(3, 8)), max_order=4, symmetrise=True)
        if net is None:
            continue
        pos = False
        for _ in range(12):
            x = rng.uniform(0.05, 1.5, len(net.species))
            x[1] = x[0]
            try:
                if P_at(net, x) > 1e-9:
                    pos = True
                    break
            except Exception:
                pass
        tot += 1
        amp += int(pos)
    print(f"  {amp}/{tot} = {100*amp/max(tot,1):.1f}% have P > 0 somewhere")

    print("\n=== P5: does P keep one sign from the symmetric state out to the attractor?")
    print(f"{'gamma':>7}{'P(0)':>12}{'P(0.5d*)':>12}{'P(d*)':>12}{'monotone?':>11}")
    for gamma in (0.05, 0.20, 0.35, 0.45):
        net = am_reversible(gamma)
        ds = delta_star(gamma)
        ps = []
        for f in (0.0, 0.5, 1.0):
            st = slaved(net, f * ds)
            ps.append(P_at(net, st) if st is not None else float("nan"))
        mono = all(ps[i] >= ps[i + 1] - 1e-12 for i in range(len(ps) - 1))
        print(f"{gamma:>7.2f}{ps[0]:>12.6f}{ps[1]:>12.6f}{ps[2]:>12.3e}"
              f"{'yes' if mono else 'NO':>11}")

    out = {"p1_worst": worst, "p3_flips": flips, "p3_tested": tested,
           "p4_amplifying": amp, "p4_total": tot}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
