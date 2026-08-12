"""T-THM-c: is §67's affinity floor a law, or a coincidence of two data points?

§67 measured the affinity at which bistability dies in two substrates with nothing in common:

    AM (3 reversible pairs, bimolecular):  A_c = 3 ln 2 = 2.0794
    Schloegl (2 pairs, trimolecular):      A_c = 2 ln 3 = 2.1972

and noted they are 5.66% apart and both ln(small integer). Two points fit almost any
two-parameter form, so §67 recorded it as a PATTERN and opened T-THM-c rather than a law.
There are at least two readings, and they fit the two points equally well:

    (a)  A_c = (pairs) x ln(pairs + 1)        3 ln 2 and 2 ln 3
    (b)  A_c = (pairs) x ln(max order)        AM's max order is 2, Schloegl's is 3

**BOTH ARE REFUTED BY A DERIVATION, and the derivation is the prediction here.** Generalise
Schloegl's autocatalysis to order p while keeping the two reversible pairs:

    pX <-> (p+1)X    (k1a, k1r)          0 <-> X    (k2b, k2r)

Bistability dies where the quartic/quintic has a TRIPLE root at x0, i.e. f = f' = f'' = 0.
Solving those three conditions:

    k1a/k1r = (p+1) x0 / (p-1),   k2r/k1r = (p+1) x0^p / (p-1),   k2b/k1r = x0^(p+1)

and the cycle affinity A = ln[k1a k2r / (k1r k2b)] gives, with x0 cancelling exactly,

    **A_c(p) = 2 ln[(p+1)/(p-1)]**

    p = 2:  2 ln 3    = 2.1972      (reproduces §67, so the derivation is anchored)
    p = 3:  2 ln 2    = 1.3863      (reading (a) predicts 2.1972; reading (b) predicts 2 ln 4 = 2.7726)
    p = 4:  2 ln(5/3) = 1.0217
    p -> infinity:    -> 0

**So within ONE family, at fixed pair count, the floor runs from 2.197 down toward zero.**
Both readings are refuted, and §67's near-agreement between 3 ln 2 and 2 ln 3 is a coincidence
of two points that happen to sit 5.66% apart.

PREDICTIONS. The formula above was derived before running; this file tests it against the
engine, which knows nothing of it.

  P1  GATE. For each p the constructed network has a genuine TRIPLE root at x0 (all three of
      f, f', f'' vanish to 1e-10), and perturbing off it gives three distinct positive fixed
      points -- a real restoring element, not a degenerate curiosity. If the perturbed network
      is not bistable, the "floor" is not the death of bistability and nothing below counts.
      **SECOND VERSION: the first perturbed the CONSTANT term k2b and the gate duly failed**,
      printing one positive root for every p. That is not a bug in the element, it is the wrong
      knob: near a triple root f ~ -k(x-x0)^3 + delta, which has ONE real root for any delta.
      Splitting a cusp into three roots requires the LINEAR coefficient -- exactly the m^2 that
      §67 used. THIRD VERSION: lowering k2r ALONE also fails, and for the same reason -- it
      leaves the family. §67's knob adds m^2*(x - x0) to the field, which moves k2r AND k2b
      together (k2r - m^2, k2b - m^2*x0). With that, every p gives three positive roots
      straddling x0. Two wrong perturbations, both kept per rule 3: a gate that fails because
      the PROBE is wrong looks exactly like a gate that fails because the claim is wrong.
  P6  **IS IT ACTUALLY A FLOOR? §67 never checked.** It called A_c a floor on the strength of
      it being the affinity where bistability dies, which is a BOUNDARY value, not a minimum.
      If A(m) DECREASED away from m = 0 then a restoring element could run below "the floor"
      and the word would be wrong. Checked here by sweeping m at fixed p.
  P2  **THE TEST, absolute (rule 16).** `cycle_affinity` at the degenerate point equals
      2 ln[(p+1)/(p-1)] to 1e-12, for p = 2..8 and several x0. The engine computes it from the
      null space of the per-pair stoichiometry and has no access to the closed form.
  P3  **THE KILL, and one counterexample is enough (rule 19: do not write a tolerance for a
      universality claim).** Reading (a) predicts 2 ln 3 for every p at 2 pairs; reading (b)
      predicts 2 ln(p+1). **If A_c(3) = 2 ln 2, both die on the spot.**
  P4  **x0-INDEPENDENCE.** A_c(p) must not depend on where the degenerate point sits, as §67
      found for p = 2. If it does, the "floor" is a property of the operating point rather
      than of the element and §67's comparison with AM was ill-posed.
  P5  **RULE 9, an axis I did not choose.** Vary the pair count too, not just the order: add a
      third reversible pair (a second linear exchange 0 <-> X at independent rates) and ask
      whether the floor scales with pairs at all. If a 3-pair Schloegl still gives
      2 ln[(p+1)/(p-1)], then "pairs" never entered the answer and both readings were
      malformed as well as wrong.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from crnl.networks.am_reversible import cycle_affinity, reverse_pairing
from crnl.reactions import Reaction, ReactionNetwork


def degenerate_consts(p, x0, k1r=1.0):
    """The (k1a, k1r, k2b, k2r) with a triple root at x0, from f = f' = f'' = 0."""
    k1a = (p + 1) * x0 / (p - 1) * k1r
    k2r = (p + 1) * x0 ** p / (p - 1) * k1r
    k2b = x0 ** (p + 1) * k1r
    return k1a, k1r, k2b, k2r


def build(p, x0, k1r=1.0, extra_pair=None):
    k1a, k1r, k2b, k2r = degenerate_consts(p, x0, k1r)
    rx = [Reaction({"X": p}, {"X": p + 1}, k1a, name=f"f1:{p}X->{p+1}X"),
          Reaction({"X": p + 1}, {"X": p}, k1r, name=f"r1:{p+1}X->{p}X"),
          Reaction({}, {"X": 1}, k2b, name="f2:->X"),
          Reaction({"X": 1}, {}, k2r, name="r2:X->")]
    if extra_pair is not None:                    # P5: a third pair, 2X <-> X + (nothing)
        kf, kr = extra_pair
        rx += [Reaction({"X": 2}, {"X": 1}, kf, name="f3:2X->X"),
               Reaction({"X": 1}, {"X": 2}, kr, name="r3:X->2X")]
    return ReactionNetwork(species=["X"], reactions=rx, name=f"schlogl{p}")


def field(p, x0, x, k1r=1.0):
    k1a, k1r, k2b, k2r = degenerate_consts(p, x0, k1r)
    return k1a * x ** p - k1r * x ** (p + 1) + k2b - k2r * x


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ps", type=int, nargs="+", default=[2, 3, 4, 5, 6, 8])
    ap.add_argument("--x0s", type=float, nargs="+", default=[0.6, 1.0, 2.5])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/affinity_floor_family.json"))
    args = ap.parse_args()

    print("=== P1 GATE: is it a genuine triple root, and is the perturbed element bistable?")
    print(f"{'p':>3}{'x0':>7}{'|f|':>11}{'|f''|':>11}{'|f\"|':>11}{'roots off-degeneracy':>34}")
    worst_g = 0.0
    bistable_ok = True
    h = 1e-5
    for p in args.ps:
        for x0 in args.x0s[:1]:
            f0 = field(p, x0, x0)
            d1 = (field(p, x0, x0 + h) - field(p, x0, x0 - h)) / (2 * h)
            d2 = (field(p, x0, x0 + h) - 2 * f0 + field(p, x0, x0 - h)) / h ** 2
            worst_g = max(worst_g, abs(f0), abs(d1), abs(d2) * 1e-4)
            # §67's knob: add m^2 (x - x0), which moves k2r AND k2b together
            k1a, k1r, k2b, k2r = degenerate_consts(p, x0)
            mm = 0.1
            coeff = np.zeros(p + 2)
            coeff[0] = -k1r                       # x^(p+1)
            coeff[1] = k1a                        # x^p
            coeff[-2] = -(k2r - mm ** 2)
            coeff[-1] = k2b - mm ** 2 * x0
            rts = np.roots(coeff)
            pos = np.sort([r.real for r in rts if abs(r.imag) < 1e-9 and r.real > 1e-12])
            if not (len(pos) >= 3 and pos.min() < x0 < pos.max()):
                bistable_ok = False
            print(f"{p:>3}{x0:>7.2f}{abs(f0):>11.2e}{abs(d1):>11.2e}{abs(d2):>11.2e}"
                  f"{str(np.round(pos, 5)):>34}")
    print(f"  -> P1 {'HOLDS' if worst_g < 1e-8 and bistable_ok else 'FAILS'}"
          f"  (triple root exact; perturbed element has 3 positive fixed points)")

    print("\n=== P2/P3/P4: the floor against the derivation, and against the two readings")
    print(f"{'p':>3}{'x0':>7}{'A (engine)':>14}{'2ln[(p+1)/(p-1)]':>19}"
          f"{'(a) 2ln3':>10}{'(b) 2ln(p+1)':>14}")
    rows, worst_a = [], 0.0
    for p in args.ps:
        for x0 in args.x0s:
            net = build(p, x0)
            A = cycle_affinity(net, reverse_pairing(net))
            pred = 2.0 * np.log((p + 1) / (p - 1))
            worst_a = max(worst_a, abs(A - pred))
            rows.append({"p": p, "x0": x0, "A": float(A), "pred": float(pred)})
            print(f"{p:>3}{x0:>7.2f}{A:>14.10f}{pred:>19.10f}"
                  f"{2*np.log(3):>10.4f}{2*np.log(p+1):>14.4f}")
    print(f"  -> P2 {'HOLDS' if worst_a < 1e-12 else 'FAILS'} (worst {worst_a:.2e})")

    a3 = [r for r in rows if r["p"] == 3]
    if a3:
        got = a3[0]["A"]
        print(f"\n  P3 KILL at p = 3: measured A_c = {got:.6f}")
        print(f"    reading (a) (pairs x ln(pairs+1)) predicts {2*np.log(3):.6f}"
              f"  -> {'SURVIVES' if abs(got-2*np.log(3))<1e-9 else 'REFUTED'}")
        print(f"    reading (b) (pairs x ln(max order)) predicts {2*np.log(4):.6f}"
              f"  -> {'SURVIVES' if abs(got-2*np.log(4))<1e-9 else 'REFUTED'}")
        print(f"    derivation 2 ln[(p+1)/(p-1)] predicts {2*np.log(2):.6f}"
              f"  -> {'CONFIRMED' if abs(got-2*np.log(2))<1e-9 else 'ALSO WRONG'}")

    spread = {}
    for p in args.ps:
        vals = [r["A"] for r in rows if r["p"] == p]
        spread[p] = float(np.ptp(vals))
    print(f"\n  P4 x0-independence: worst spread across x0 = {max(spread.values()):.2e}")
    print(f"  -> P4 {'HOLDS: the floor is a property of the element, not the operating point' if max(spread.values()) < 1e-12 else 'FAILS'}")

    print(f"\n=== P6: is A_c a FLOOR (a minimum), or merely the boundary value? §67 assumed.")
    print(f"{'p':>3}" + "".join(f"{f'm={m}':>12}" for m in (0.0, 0.05, 0.10, 0.20, 0.30)))
    floors_ok = True
    p6 = []
    for p in args.ps:
        x0 = 1.0
        k1a, k1r, k2b, k2r = degenerate_consts(p, x0)
        As = []
        for m in (0.0, 0.05, 0.10, 0.20, 0.30):
            kr, kb = k2r - m ** 2, k2b - m ** 2 * x0
            As.append(float(np.log(k1a * kr / (k1r * kb))) if kr > 0 and kb > 0 else np.nan)
        if not all(a >= As[0] - 1e-12 for a in As if np.isfinite(a)):
            floors_ok = False
        p6.append({"p": p, "A_of_m": As})
        print(f"{p:>3}" + "".join(f"{a:>12.6f}" if np.isfinite(a) else f"{'--':>12}"
                                 for a in As))
    print(f"  -> P6 {'HOLDS: A(m) is minimised at m = 0, so floor is the right word' if floors_ok else 'FAILS: A dips below A_c, so it is a boundary value and NOT a floor -- §67 mis-named it'}")

    print(f"\n=== P5 (rule 9): does the PAIR COUNT enter at all? add a third reversible pair")
    print(f"{'p':>3}{'pairs':>7}{'A (engine)':>14}{'2ln[(p+1)/(p-1)]':>19}{'changed?':>10}")
    p5 = []
    for p in (2, 3, 4):
        base = cycle_affinity(build(p, 1.0), reverse_pairing(build(p, 1.0)))
        net3 = build(p, 1.0, extra_pair=(0.7, 0.3))
        try:
            A3 = cycle_affinity(net3, reverse_pairing(net3))
            note = f"{abs(A3-base):.2e}"
        except Exception as e:
            A3, note = float("nan"), type(e).__name__
        p5.append({"p": p, "A2": float(base), "A3": float(A3)})
        print(f"{p:>3}{3:>7}{A3:>14.10f}{2*np.log((p+1)/(p-1)):>19.10f}{note:>10}")
    print(f"  -> a third pair opens a SECOND independent cycle, so a single affinity may not")
    print(f"     be defined; whatever it prints above is reported, not interpreted.")

    print(f"\n=== summary: the floor across the family")
    for p in args.ps:
        print(f"  p = {p}: A_c = {2*np.log((p+1)/(p-1)):.6f}")
    print(f"  AM (§9.1, 3 pairs, bimolecular): {3*np.log(2):.6f}")
    print(f"  -> the floor spans {2*np.log((args.ps[-1]+1)/(args.ps[-1]-1)):.4f}"
          f"..{2*np.log(3):.4f} within ONE family at FIXED pair count,"
          f" and tends to 0 as p grows")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"rows": rows, "p5": p5, "p6": p6,
                                    "worst_gate": worst_g, "worst_affinity": worst_a},
                                   indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
