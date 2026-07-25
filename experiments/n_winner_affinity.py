"""What does an n-symbol landscape cost in drive? (THEORIES Q3)

FINDINGS 9.1 measured the two-symbol answer: the landscape dies at gamma_c = 1/2,
so it costs a minimum cycle affinity A = 3 ln 2. This sweeps n.

The elementary cycle is three reactions for EVERY n (X_i+X_j -> 2B, B+X_i -> 2X_i,
B+X_j -> 2X_j returns every count to its start), so the affinity per cycle stays
A = -3 ln gamma and the whole question is where gamma_c(n) sits.

Two guesses that the sweep kills: A_c = n ln 2 (THEORIES' original) and
gamma_c = 1/n giving A_c = 3 ln n (predicted before running, because it
reproduces the exact n=2 value). Measured, gamma_c falls far faster than 1/n --
by n=32 it is ~1000x below it.

    python -m experiments.n_winner_affinity
    python -m experiments.n_winner_affinity --quick
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

from crnl.networks.n_winner_reversible import (
    affinity_critical, gamma_critical, lambda_breaking, symmetric_state,
)

LN2 = np.log(2.0)


def main():
    here = os.path.dirname(__file__)
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ns", type=int, nargs="+",
                   default=[2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256])
    p.add_argument("--tail-from", type=int, default=64,
                   help="fit the asymptotic exponent using n >= this")
    p.add_argument("--quick", action="store_true")
    p.add_argument("--data", default=os.path.join(
        here, os.pardir, "results", "n_winner_affinity.json"))
    args = p.parse_args()
    if args.quick:
        args.ns, args.tail_from = [2, 4, 8, 16, 32], 8

    print("Minimum drive for an n-symbol landscape (exact: analytic Jacobian at "
          "the symmetric point)")
    print(f"\n{'n':>5} {'gamma_c':>13} {'A_c':>9} {'local exp':>10} "
          f"{'A_c/ln n':>9} {'A_c/ln2 per bit':>16}")
    rows, prev = [], None
    for n in args.ns:
        gc = gamma_critical(n)
        ac = affinity_critical(n)
        x, b = symmetric_state(n, gc)
        exp_local = (float(np.log(gc / prev[1]) / np.log(n / prev[0]))
                     if prev else None)
        rows.append({"n": n, "gamma_c": gc, "A_c": ac,
                     "local_exponent": exp_local,
                     "A_c_over_ln_n": ac / np.log(n),
                     "per_bit": ac / np.log2(n),
                     "symmetric_x": x, "symmetric_b": b})
        le = "" if exp_local is None else f"{exp_local:>10.3f}"
        print(f"{n:>5} {gc:>13.6e} {ac:>9.4f} {le:>10} "
              f"{ac / np.log(n):>9.4f} {ac / np.log2(n):>16.4f}")
        prev = (n, gc)

    tail = [r for r in rows if r["n"] >= args.tail_from]
    fit = None
    if len(tail) >= 3:
        slope, icept = np.polyfit([np.log(r["n"]) for r in tail],
                                  [np.log(r["gamma_c"]) for r in tail], 1)
        fit = {"exponent": float(slope), "intercept": float(icept),
               "A_c_coefficient": float(-3 * slope), "from_n": args.tail_from}
        print(f"\nfit over n >= {args.tail_from}:  gamma_c ~ n^({slope:.4f})"
              f"   =>   A_c ~ {-3 * slope:.3f} ln n")
        print(f"  n^-3 would give exponent -3.000 and A_c = 9 ln n = "
              f"{9 * LN2:.3f} k_BT per bit of alphabet")

    print(f"\nn=2 is NOT on that asymptote: A_c(2) = {rows[0]['A_c']:.4f} "
          f"= 3 ln 2, where 9 ln 2 = {9 * LN2:.4f}.")
    print("The famous case is the special one; the law is approached from below.")

    os.makedirs(os.path.dirname(os.path.abspath(args.data)), exist_ok=True)
    with open(args.data, "w") as fh:
        json.dump({"rows": rows, "fit": fit}, fh, indent=2)
    print(f"\nwrote data -> {args.data}")


if __name__ == "__main__":
    main()
