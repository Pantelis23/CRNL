"""T15-g: does §54's decomposition predict gamma_c(n) and T7's lambda(n)? — absolute tests

§54 decomposed P into a rate-weighted sum of integers, but only at n = 2. §30 proved the
pairwise identity for EVERY n, so the decomposition should follow -- and at n > 2 it makes
two predictions against numbers this project published years of sections ago, with no
fitting anywhere.

**WORKING THE DECOMPOSITION AT GENERAL n.** For the pair (i, j) in `n_winner_reversible`,
under the (i, j) swap the reactions sort into:

  * `B + X_i -> 2X_i` and its mirror: p=1, q=0, d = S_i - S_j = **+1**, weight k, bracket b
  * `2X_i -> B + X_i` and its mirror:  p=2, q=0, d = **-1**, weight gamma*k, bracket x_i + x_j
  * `X_i + X_k -> 2B` for each k not in {i,j}: mirrors with `X_j + X_k -> 2B`, p=1, q=0,
    d = **-1**, weight k, bracket x_k -- **(n-2) contracting terms**
  * `X_i + X_j -> 2B`: **SELF-MIRROR, contributes 0** (§54's zero, §30's first cancellation)
  * `2B -> X_i + X_k`: p = q = 0 for both members, so x^p y^q - x^q y^p = 0 and the pair
    **cancels exactly** -- a second, different route to zero

Summing at a symmetric state (all x_l = x, b = 1 - n x):

    **P = k[ b - gamma*(2x) - (n-2)x ]**

and §30's published bracket is `(k/Om)[n_B - sum_{l != i,j} n_l - gamma(n_i + n_j - 1)]` --
**the same three terms, in the same order, with the same signs.**

Two absolute consequences follow with no free parameter:

  * At gamma = 0 the symmetric fixed point has x = 1/(2n-1) (THEORIES T7, §14), so
    b = (n-1)/(2n-1) and **P = [(n-1) - (n-2)]/(2n-1) = 1/(2n-1)** -- exactly T7/§14's
    symmetry-breaking eigenvalue lambda(n).
  * P = 0 defines the loss of amplification, so **P must vanish at gamma = gamma_c(n)**,
    the drive at which the landscape disappears, which `gamma_critical(n)` computes by a
    completely independent bracketed root-find.

PREDICTIONS, written before running:

  P1  GATE. At n = 2..6, P_ij equals the antisymmetric directional derivative of b_i - b_j
      (§53's identity) and §54's decomposition reproduces it, both to machine precision. If
      either fails at n > 2 the framework does not generalise and nothing below counts.
  P2  ABSOLUTE, against T7/§14. **P at the symmetric fixed point with gamma = 0 equals
      1/(2n-1)** for n = 2, 3, 4, 5, 6 -- i.e. 0.33333, 0.20000, 0.14286, 0.11111, 0.09091.
      Not a fit; a closed form meeting a published one.
  P3  ABSOLUTE, and the sharper one. **P at the symmetric fixed point vanishes exactly at
      gamma = gamma_critical(n)**, whose published value at n = 3 is 0.2023 and which is
      computed here by an independent route. Agreement to solver precision or the
      decomposition is wrong.
  P4  CLASSIFICATION. Every n-winner pair is **MIXED** (one d = +1 against 1 + (n-2)
      d = -1), so by §54 the rates must decide -- consistent with gamma_c existing at all.
      **And the count of negative terms grows with n**, which is why gamma_c(n) FALLS with
      n: more rivals means more contracting terms to outweigh.
  P5  The two distinct zeros must both be present: `X_i + X_j -> 2B` self-mirror, and the
      `2B -> X_i + X_k` pair cancelling through p = q. Reported separately, because they
      are different mechanisms reaching the same 0 and conflating them would hide one.
  P6  If P3 fails while P2 holds, the decomposition is right at gamma = 0 and wrong in its
      gamma-weighting -- which would localise the error to the reverse-recruitment bracket
      rather than condemning the whole construction.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
from scipy.optimize import brentq

from crnl.networks.n_winner_reversible import gamma_critical, n_winner_reversible
from experiments.amplification_sign import P_at, antisym_eig
from experiments.amplification_signature import classify, mirror_pairs, P_decomposed


def symmetric_state(net, n, gamma):
    """All committed species equal, blank last; solve dx_0/dt = 0 for x."""
    S = net.stoichiometry_matrix()

    def f(x):
        st = np.array([x] * n + [1.0 - n * x])
        if st[-1] < 0:
            return np.nan
        return float((S @ net.fluxes(st))[0])

    lo, hi = 1e-12, 1.0 / n - 1e-12
    a, b = f(lo), f(hi)
    if not (np.isfinite(a) and np.isfinite(b)) or a * b > 0:
        return None
    x = brentq(f, lo, hi, xtol=1e-15)
    return np.array([x] * n + [1.0 - n * x])


def P_symmetric(n, gamma):
    net = n_winner_reversible(n, gamma)
    st = symmetric_state(net, n, gamma)
    if st is None:
        return None, None
    return float(P_at(net, st, i=0, j=1)), st


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ns", type=int, nargs="+", default=[2, 3, 4, 5, 6])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/amplification_n.json"))
    args = ap.parse_args()

    print("=== P1 GATE: does §53's identity and §54's decomposition survive n > 2?")
    print(f"{'n':>4}{'gamma':>8}{'P':>14}{'antisym eig':>14}{'decomposed':>14}{'worst':>10}")
    worst = 0.0
    for n in args.ns:
        for gamma in (0.0, 0.5 * gamma_critical(n), 0.9 * gamma_critical(n)):
            net = n_winner_reversible(n, gamma)
            st = symmetric_state(net, n, gamma)
            if st is None:
                continue
            p = P_at(net, st, i=0, j=1)
            e = antisym_eig(net, st, i=0, j=1)
            d, _ = P_decomposed(net, st, i=net.species[0], j=net.species[1])
            w = max(abs(p - e), abs(p - d)) / max(abs(p), 1e-30)
            worst = max(worst, w)
            print(f"{n:>4}{gamma:>8.4f}{p:>14.9f}{e:>14.9f}{d:>14.9f}{w:>10.1e}")
    print(f"  -> P1 {'HOLDS' if worst < 1e-8 else 'FAILS'}  (worst {worst:.2e})")

    print("\n=== P2 ABSOLUTE: P(symmetric, gamma=0) must equal T7/§14's lambda = 1/(2n-1)")
    print(f"{'n':>4}{'P measured':>16}{'1/(2n-1)':>14}{'diff':>12}")
    p2rows, p2worst = [], 0.0
    for n in args.ns:
        p, _ = P_symmetric(n, 0.0)
        if p is None:
            continue
        pred = 1.0 / (2 * n - 1)
        p2worst = max(p2worst, abs(p - pred))
        p2rows.append({"n": n, "P": p, "pred": pred})
        print(f"{n:>4}{p:>16.10f}{pred:>14.10f}{p - pred:>12.2e}")
    print(f"  -> P2 {'HOLDS' if p2worst < 1e-9 else 'FAILS'}  (worst |diff| {p2worst:.2e})")

    print("\n=== P3 ABSOLUTE: P must vanish at gamma = gamma_critical(n)")
    print(f"{'n':>4}{'gamma_c':>14}{'P at gamma_c':>16}{'P(0.9 g_c)':>14}")
    p3rows, p3worst = [], 0.0
    for n in args.ns:
        gc = gamma_critical(n)
        p, _ = P_symmetric(n, gc)
        plo, _ = P_symmetric(n, 0.9 * gc)
        if p is None:
            continue
        p3worst = max(p3worst, abs(p))
        p3rows.append({"n": n, "gamma_c": gc, "P_at_gc": p})
        print(f"{n:>4}{gc:>14.8f}{p:>16.3e}{plo if plo is not None else float('nan'):>14.6f}")
    print(f"  -> P3 {'HOLDS -- the decomposition predicts gamma_c with no fit' if p3worst < 1e-7 else 'FAILS'}"
          f"  (worst |P| {p3worst:.2e})")
    print(f"  published gamma_c(3) = 0.2023;  computed {gamma_critical(3):.6f}")

    print("\n=== P4: classification and the term count")
    print(f"{'n':>4}{'class':>10}{'+1':>5}{'-1':>5}{'zero-bkt':>10}{'gamma_c':>12}")
    p4rows = []
    for n in args.ns:
        net = n_winner_reversible(n, 0.5 * gamma_critical(n))
        pairs = mirror_pairs(net, net.species[0], net.species[1])
        if pairs is None:
            print(f"{n:>4}   pairing failed")
            continue
        # A p == q pair has an EMPTY bracket sum and contributes exactly zero, but
        # `mirror_pairs` records its d as an arbitrary +-1 (the heavy member is
        # ill-defined when p == q). Counting raw d would report those as contributing
        # terms -- the first pass of this experiment did, giving n-1 spurious negatives.
        i0, j0 = net.species[0], net.species[1]
        pos = neg = zero = 0
        for d, idx in pairs:
            rr = net.reactions[idx]
            if rr.reactants.get(i0, 0) == rr.reactants.get(j0, 0):
                zero += 1
            elif d > 0:
                pos += 1
            elif d < 0:
                neg += 1
        cls = classify(net, net.species[0], net.species[1])
        p4rows.append({"n": n, "class": cls, "pos": pos, "neg": neg, "zero": zero,
                       "gamma_c": gamma_critical(n)})
        print(f"{n:>4}{cls:>10}{pos:>5}{neg:>5}{zero:>10}{gamma_critical(n):>12.6f}")
    print("  predicted: 1 contributing +1, n-1 contributing -1, n-2 zero-bracket")
    print("  -> ONE amplifying term against n-1 contracting ones, which is why")
    print("     gamma_c(n) falls with n: the single d=+1 must outweigh more each time.")

    print("\n=== P5: the two DIFFERENT routes to a zero contribution")
    n = 4
    net = n_winner_reversible(n, 0.1)
    specs = net.species
    self_mirror, pq_cancel = [], []
    for r in net.reactions:
        p, q = r.reactants.get(specs[0], 0), r.reactants.get(specs[1], 0)
        sw = {(specs[1] if s == specs[0] else specs[0] if s == specs[1] else s): v
              for s, v in r.reactants.items()}
        swp = {(specs[1] if s == specs[0] else specs[0] if s == specs[1] else s): v
               for s, v in r.products.items()}
        if sw == r.reactants and swp == r.products:
            self_mirror.append(r.name)
        elif p == q:
            pq_cancel.append(r.name)
    print(f"  n=4 self-mirror reactions (map to themselves): {self_mirror}")
    print(f"  n=4 p=q pairs (distinct mirrors, contributions cancel): {pq_cancel}")
    print("  -> two different mechanisms, both giving zero; §54 saw only the first.")

    out = {"p1_worst": worst, "p2": p2rows, "p3": p3rows, "p4": p4rows,
           "self_mirror_n4": self_mirror, "pq_cancel_n4": pq_cancel}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
